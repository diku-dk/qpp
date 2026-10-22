-- | Rewrite passes on flat circuits: flattening to elementary gates,
--   compaction to support, absorbing permutations into wire labels, fusing
--   adjacent gates, and collecting global phases. Every pass preserves the
--   dense matrix of 'toQOp'; 'simplify' runs them all in a sensible order.
--   Input is assumed well-formed, as for 'fromQOp'.
--
--   Diagonal gates form a gate class of their own. A composite diagonal
--   operator is kept as one gate and put into a canonical form: a balanced
--   direct-sum tree of phases over ascending wires, with its global phase
--   split off into a gate on no wires. Equal diagonals are then equal terms,
--   and two diagonal gates on the same wires fuse by adding phases. Atoms
--   under controls, such as @C (R Z t)@ or @C Z@, are kept as written.
module QPP.QCirc.Rewrite
  ( flatten, elementary
  , compact
  , absorbPermutations
  , fuse, orderWires
  , collectPhases
  , simplify
    -- * Diagonal gates
  , diagonal, diagonalPhases, canonicalDiagonal, maxDiagonalWires
  ) where

import Data.Array (Array, listArray, (!))
import Data.Bits (bit, testBit, xor)
import Data.List (elemIndex, partition, sort, sortOn)
import Data.Maybe (fromMaybe, isJust, mapMaybe, maybeToList)
import Data.Ratio (denominator, numerator)
import qualified Data.Set as S
import QPP.QCirc.Syntax
import QPP.Rewrite (assocTensors, cleanOnes, doComp, fixpoint)
import QPP.Syntax
import QPP.Util (reduceAngle)

--------------------------------------------------------------------------------
-- Flattening to elementary gates
--------------------------------------------------------------------------------

-- | Split every gate into elementary gates: an atom (@X@, @Y@, @Z@, @H@,
--   @SX@, @R@, @Permute@, @Phase@) under zero or more controls, or a
--   canonical diagonal gate. Compositions, tensor products and adjoints are
--   unfolded, @C@ distributes over them, and a non-diagonal direct sum
--   @a ⊕ b@ becomes @X@-conjugated controls of @a@ followed by controls of
--   @b@ on the selector wire. A composite diagonal operator, controls
--   included, becomes one canonical diagonal gate instead, unless it is wider
--   than 'maxDiagonalWires'. Identity padding disappears, so idle wires
--   inside a control are dropped as well. Rotation axes are kept as written.
flatten :: QCirc -> QCirc
flatten (QCirc n gs) = QCirc n (concatMap (\(Gate g ws) -> unfold [] g ws) gs)

-- Controls accumulate outermost first; the target wires follow them.
unfold :: [Wire] -> QOp -> [Wire] -> [Gate]
unfold cs g ws = case g of
  Id _          -> []
  Adjoint a     -> unfold cs (dagger a) ws
  C a           -> unfold (cs ++ [c]) a targets
  _ | not (atom g), diagonal g, Just gates <- canonical -> gates
  Compose a b   -> unfold cs b ws ++ unfold cs a ws
  Tensor a b    -> let (wa, wb) = splitAt (op_qubits a) ws
                   in unfold cs a wa ++ unfold cs b wb
  DirectSum a b -> conjugated (unfold (cs ++ [c]) a targets) ++ unfold (cs ++ [c]) b targets
  _             -> maybeToList (compactGate controlled)  -- An atom, kept as written.
  where
    controlled = Gate (iterate C g !! length cs) (cs ++ ws)
    -- The support of the controlled gate bounds the width of its phase list.
    -- Compaction first, so that identity padding does not widen the list.
    canonical
      | length cs + S.size (op_support g) > maxDiagonalWires = Nothing
      | otherwise = case compactGate controlled of
          Nothing            -> Just []
          Just (Gate g' ws') -> canonicalDiagonal ws' <$> diagonalPhases g'
    (c, targets) = case ws of
      (w : more) -> (w, more)
      []         -> error "QPP.QCirc.Rewrite.flatten: control without a wire"
    conjugated [] = []
    conjugated hs = Gate X [c] : hs ++ [Gate X [c]]

-- | An atom under zero or more controls, or a tree of phases, as 'flatten'
--   produces.
elementary :: Gate -> Bool
elementary (Gate g _) = atom (snd (peelControls g)) || tree g
  where
    tree (Id _) = False
    tree a      = phaseTree a

-- | 'cleanOnes' to a fixpoint: identities under controls, zero rotations and
--   full turns collapse to a single 'Id', and angles are reduced.
clean :: QOp -> QOp
clean = fixpoint 100 cleanOnes

atom :: QOp -> Bool
atom op = case op of
  X         -> True
  Y         -> True
  Z         -> True
  H         -> True
  SX        -> True
  R _ _     -> True
  Permute _ -> True
  Phase _   -> True
  _         -> False

-- | A balanced direct-sum tree of phases, with identity subtrees as 'Id' and 'C'.
phaseTree :: QOp -> Bool
phaseTree op = case op of
  Phase _       -> True
  Id _          -> True
  C a           -> phaseTree a
  DirectSum a b -> phaseTree a && phaseTree b
  _             -> False

--------------------------------------------------------------------------------
-- Diagonal gates
--------------------------------------------------------------------------------

-- | Composite diagonal gates wider than this are not canonicalized, since a
--   phase list has 2^k entries; they fall back to structural unfolding, and
--   wider diagonal pairs are not fused by phases.
maxDiagonalWires :: Int
maxDiagonalWires = 12

-- | Structurally diagonal: identities, phases, @Z@, rotations about a
--   diagonal axis, and sums, controls, products, compositions and adjoints
--   of diagonals.
diagonal :: QOp -> Bool
diagonal op = case op of
  Id _          -> True
  Phase _       -> True
  Z             -> True
  R a _         -> diagonal a
  C a           -> diagonal a
  DirectSum a b -> diagonal a && diagonal b
  Tensor a b    -> diagonal a && diagonal b
  Compose a b   -> diagonal a && diagonal b
  Adjoint a     -> diagonal a
  _             -> False

-- | The exponents θ_x of a diagonal operator diag(e^{iπθ_x}), indexed by basis
--   state with local qubit 0 as the most significant bit. 'Nothing' for a
--   non-diagonal operator, or a rotation about a non-Hermitian diagonal axis.
diagonalPhases :: QOp -> Maybe [Rational]
diagonalPhases op = case op of
  Id n          -> Just (replicate (2 ^ n) 0)
  Phase t       -> Just [t]
  Z             -> Just [0, 1]
  R a 0         -> Just (replicate (2 ^ op_qubits a) 0)
  R a t         -> diagonalPhases a >>= mapM (rotation t)
  C a           -> (\ps -> map (const 0) ps ++ ps) <$> diagonalPhases a
  DirectSum a b -> (++) <$> diagonalPhases a <*> diagonalPhases b
  Tensor a b    -> (\pa pb -> [ p + q | p <- pa, q <- pb ]) <$> diagonalPhases a <*> diagonalPhases b
  Compose a b   -> zipWith (+) <$> diagonalPhases a <*> diagonalPhases b
  Adjoint a     -> map negate <$> diagonalPhases a
  _             -> Nothing
  where
    -- A Hermitian diagonal axis has entries ±1, and exp(-iπt/2 · (±1)) = e^{∓iπt/2}.
    rotation t p
      | denominator p == 1 = Just (if even (numerator p) then negate (t / 2) else t / 2)
      | otherwise          = Nothing

-- | Canonical gates for the diagonal with exponents @ps@ over wires @ws@, in
--   that bit order: a global 'Phase' gate on no wires when the exponent at
--   |0…0⟩ is nonzero, then one gate holding a balanced direct-sum tree of
--   phases over the wires in ascending order. Exponents are reduced into
--   (-1, 1], wires the diagonal does not depend on are dropped, and identity
--   subtrees collapse to 'Id' and 'C'. An identity yields no gates.
canonicalDiagonal :: [Wire] -> [Rational] -> [Gate]
canonicalDiagonal ws ps =
  [ Gate (Phase global) [] | global /= 0 ] ++ [ Gate tree kept | not (isId tree) ]
  where
    sorted = sort ws
    reordered = reorderPhases ws sorted ps
    base = case reordered of
      (p : _) -> p
      []      -> 0
    global = reduceAngle 2 base
    (kept, relative) = dropIdleWires sorted (map (\p -> reduceAngle 2 (p - base)) reordered)
    tree = build relative
    build [p] = if p == 0 then Id 0 else Phase p
    build qs  = let (l, r) = splitAt (length qs `div` 2) qs in node (build l) (build r)
    node (Id m) (Id k) | m == k = Id (m + 1)
    node (Id _) t = C t
    node l r      = DirectSum l r
    isId (Id _) = True
    isId _      = False

-- | Re-index a phase list from wire order @from@ to wire order @to@, the same
--   set of wires. The first wire is the most significant index bit.
reorderPhases :: [Wire] -> [Wire] -> [Rational] -> [Rational]
reorderPhases from to ps
  | from == to = ps
  | otherwise  = [ table ! source x | x <- [0 .. 2 ^ n - 1] ]
  where
    n = length to
    table :: Array Int Rational
    table = listArray (0, 2 ^ n - 1) ps
    position w = fromMaybe (error "QPP.QCirc.Rewrite.reorderPhases: wire missing") (elemIndex w from)
    source :: Int -> Int
    source x = sum [ bit (n - 1 - position w) | (j, w) <- zip [0 ..] to, testBit x (n - 1 - j) ]

-- | Drop every wire whose bit never changes the exponent, from the last wire
--   backwards so that earlier bit positions stay valid.
dropIdleWires :: [Wire] -> [Rational] -> ([Wire], [Rational])
dropIdleWires ws ps = foldr check (ws, ps) [0 .. length ws - 1]
  where
    check j (vs, qs)
      | all (\x -> table ! x == table ! (x `xor` mask)) [0 .. size - 1] =
          (take j vs ++ drop (j + 1) vs, [ q | (x, q) <- zip [0 :: Int ..] qs, not (testBit x b) ])
      | otherwise = (vs, qs)
      where
        n = length vs
        size = 2 ^ n :: Int
        b = n - 1 - j
        mask = bit b
        table = listArray (0, size - 1) qs

--------------------------------------------------------------------------------
-- Compaction to support
--------------------------------------------------------------------------------

-- | Restrict every gate to its structural support ('op_support'): identity
--   factors inside the operator are removed and the wire list keeps only the
--   wires the gate acts on, in order. Gates that reduce to an identity are
--   dropped. @C (Id 2 ⊗ X)@ on four wires becomes @C X@ on the control and
--   the target; @R (X ⊗ I ⊗ Z) t@ becomes @R (X ⊗ Z) t@ on two wires.
compact :: QCirc -> QCirc
compact (QCirc n gs) = QCirc n (mapMaybe compactGate gs)

-- | 'compact' for one gate; 'Nothing' when it reduces to an identity. The
--   operator is cleaned before its support is read, so that a full turn
--   vanishes and a half turn leaves only its phase.
compactGate :: Gate -> Maybe Gate
compactGate (Gate g ws) = case clean (dropIdle idle cleaned) of
  Id _ -> Nothing
  g'   -> Just (Gate g' (map (ws !!) support))
  where
    cleaned = clean g
    support = S.toAscList (op_support cleaned)
    idle = S.fromList [0 .. op_qubits cleaned - 1] `S.difference` S.fromList support

-- | Remove the given local qubits, which must lie outside 'op_support'.
dropIdle :: S.Set Nat -> QOp -> QOp
dropIdle idle g
  | S.null idle = g
  | otherwise = case g of
      Id m          -> Id (m - S.size idle)
      R _ 0         -> Id (op_qubits g - S.size idle)
      R a t         -> R (dropIdle idle a) t
      C a           -> C (dropIdle (below idle) a)
      DirectSum a b -> DirectSum (dropIdle (below idle) a) (dropIdle (below idle) b)
      Compose a b   -> Compose (dropIdle idle a) (dropIdle idle b)
      Adjoint a     -> Adjoint (dropIdle idle a)
      Tensor a b    -> let na = op_qubits a
                           (ia, ib) = S.partition (< na) idle
                       in tensor (dropIdle ia a) (dropIdle (S.map (subtract na) ib) b)
      Permute ks    -> Permute (restrict ks)
      _             -> g  -- Atoms have full support, so their idle set is empty.
  where
    below = S.map (subtract 1)  -- The selector wire of C and ⊕ is never idle.
    tensor (Id 0) b = b
    tensor a (Id 0) = a
    tensor a b      = Tensor a b
    -- Idle positions of a permutation are fixed points, so the kept positions
    -- map among themselves; renumber them in order.
    restrict ks =
      let keep = [ i | i <- [0 .. length ks - 1], i `S.notMember` idle ]
          index i = length (takeWhile (/= i) keep)
      in [ index (ks !! i) | i <- keep ]

--------------------------------------------------------------------------------
-- Permutations as wire relabeling
--------------------------------------------------------------------------------

-- | Remove 'Permute' gates by relabeling the wires of every later gate. The
--   accumulated relabeling is emitted as one final 'Permute' on all wires
--   when it is not the identity, so 'toQOp' is preserved; drop that last
--   gate to get the circuit up to a renaming of its output wires.
absorbPermutations :: QCirc -> QCirc
absorbPermutations (QCirc n gs) = QCirc n (go [0 .. n - 1] gs)
  where
    -- tau !! j is the wire of the output circuit that holds what the input
    -- circuit has on wire j at this point.
    go tau [] = [ Gate (Permute tau) [0 .. n - 1] | tau /= [0 .. n - 1] ]
    go tau (Gate (Permute ks) ws : rest) = go (relabel tau ws ks) rest
    go tau (Gate g ws : rest) = Gate g (map (tau !!) ws) : go tau rest

    -- After Permute ks on wires ws, wire (ws !! j) holds what wire
    -- (ws !! (ks !! j)) held before.
    relabel tau ws ks =
      [ maybe (tau !! w) (\j -> tau !! (ws !! (ks !! j))) (elemIndex w ws)
      | w <- [0 .. n - 1] ]

--------------------------------------------------------------------------------
-- Gate fusion
--------------------------------------------------------------------------------

-- | Fuse or cancel pairs of gates with no gate on their wires in between.
--   Wires are first put in the order of 'orderWires'. Gates on the same wires
--   then fuse structurally: adjoint pairs and self-inverse gates cancel,
--   rotations about one axis add their angles (a half turn becoming a global
--   or controlled phase), phases add, permutations compose, and controls fuse
--   when their targets do. Two diagonal gates on the same set of wires fuse
--   by adding their phases into a canonical diagonal gate, and a diagonal
--   gate may pass earlier diagonal gates to find its partner, since diagonals
--   commute. Repeats until nothing changes.
fuse :: QCirc -> QCirc
fuse (QCirc n gs) = QCirc n (converge (map orderWires gs))
  where
    converge xs =
      let ys = reverse (foldl step [] xs)
      in if ys == xs then xs else converge ys
    -- The accumulator holds the output so far, most recent gate first.
    step recent g = fromMaybe (g : recent) (fuseInto g recent)

-- | Order the wires of a gate where its operator allows it for free, so that
--   equal gates listed differently become equal terms: controls ascend, since
--   controls commute; the factors of a tensor product, or of a rotation about
--   a tensor product of one-qubit axes, follow their wires; a permutation is
--   re-expressed on its ascending wires. Any other operator keeps its order.
orderWires :: Gate -> Gate
orderWires (Gate g ws) = Gate (iterate C core' !! k) (sort cs ++ ts')
  where
    (k, core) = peelControls g
    (cs, ts) = splitAt k ws
    (core', ts') = case core of
      Permute p -> (Permute (onAscending ts p), sort ts)
      R axis t | Just (axis', us) <- factorsByWire axis ts -> (R axis' t, us)
      _ | Just (g', us) <- factorsByWire core ts -> (g', us)
        | otherwise -> (core, ts)

-- | A tensor product of one-qubit factors on shuffled wires, re-expressed
--   with its factors in wire order.
factorsByWire :: QOp -> [Wire] -> Maybe (QOp, [Wire])
factorsByWire op ws
  | length factors > 1, all ((== 1) . op_qubits) factors, ws /= sorted =
      Just (foldr1 Tensor (map snd ordered), sorted)
  | otherwise = Nothing
  where
    factors = assocTensors op
    ordered = sortOn fst (zip ws factors)
    sorted = map fst ordered

-- | The permutation that, on the ascending wires, equals @Permute p@ on wires
--   @ws@, where wire @ws !! i@ receives from wire @ws !! (p !! i)@.
onAscending :: [Wire] -> [Int] -> [Int]
onAscending ws p = [ position sorted (ws !! (p !! position ws w)) | w <- sorted ]
  where
    sorted = sort ws
    position vs v = fromMaybe (error "QPP.QCirc.Rewrite.orderWires: wire missing") (elemIndex v vs)

-- | Place a gate in the output (most recent first) by fusing it with the
--   nearest earlier gate it can fuse with, if any. A gate on no wires shares
--   no wire with anything and never fuses here; 'collectPhases' merges those.
fuseInto :: Gate -> [Gate] -> Maybe [Gate]
fuseInto (Gate _ []) _ = Nothing
fuseInto g@(Gate op ws) recent = go recent
  where
    go [] = Nothing
    go (p@(Gate h vs) : older)
      | not (any (`elem` ws) vs)   = (p :) <$> go older
      | Just fused <- fusePair p g = Just (reverse fused ++ older)  -- Most recent first.
      | diagonal op && diagonal h  = (p :) <$> go older
      | otherwise                  = Nothing

-- | The gates replacing @h@ followed by @g@, if the pair fuses.
fusePair :: Gate -> Gate -> Maybe [Gate]
fusePair (Gate h vs) (Gate g ws)
  | vs == ws, ([], c) <- doComp g h = Just (maybeToList (compactGate (Gate c ws)))
  | sort vs == sorted, length ws <= maxDiagonalWires, diagonal h, diagonal g
  , Just ph <- diagonalPhases h, Just pg <- diagonalPhases g =
      Just (canonicalDiagonal sorted
              (zipWith (+) (reorderPhases vs sorted ph) (reorderPhases ws sorted pg)))
  | otherwise = Nothing
  where
    sorted = sort ws

--------------------------------------------------------------------------------
-- Phase collection
--------------------------------------------------------------------------------

-- | Merge every gate on no wires into one global 'Phase' at the start of the
--   circuit, with its angle reduced to (-1, 1] modulo two; a zero total
--   leaves no gate. Phases inside other gates are relative and are kept.
collectPhases :: QCirc -> QCirc
collectPhases (QCirc n gs) = QCirc n ([ Gate (Phase t) [] | t /= 0 ] ++ rest)
  where
    (scalars, rest) = partition scalar gs
    scalar (Gate g ws) = null ws && isJust (scalarPhase g)
    t = reduceAngle 2 (sum [ p | Gate g _ <- scalars, Just p <- [scalarPhase g] ])

-- | The phase parameter of a zero-qubit operator, when it is syntactically a
--   product of phases.
scalarPhase :: QOp -> Maybe Rational
scalarPhase op = case op of
  Id 0        -> Just 0
  Phase t     -> Just t
  Tensor a b  -> (+) <$> scalarPhase a <*> scalarPhase b
  Compose a b -> (+) <$> scalarPhase a <*> scalarPhase b
  Adjoint a   -> negate <$> scalarPhase a
  R a t       -> do  -- exp(-i*pi*t/2 * e^{i*pi*q}) with e^{i*pi*q} = ±1 for a Hermitian scalar.
    q <- scalarPhase a
    if denominator q == 1
      then Just (if even (numerator q) then negate (t / 2) else t / 2)
      else Nothing
  _           -> Nothing

--------------------------------------------------------------------------------
-- Everything, in order
--------------------------------------------------------------------------------

-- | 'flatten', 'compact', 'absorbPermutations', 'fuse', 'collectPhases'.
simplify :: QCirc -> QCirc
simplify = collectPhases . fuse . absorbPermutations . compact . flatten
