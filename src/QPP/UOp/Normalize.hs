-- | Symbolic simplification of well-formed UOp terms; no matrices are built.
--
-- This is a structural normal form, not a decision procedure for equality.
-- Input must already satisfy the UOp typing conditions (for example, come
-- from 'QPP.UOp.FromQOp.fromQOp'). In particular, every nonzero rotation
-- needs a Hermitian axis BEFORE its angle is reduced modulo four.
module QPP.UOp.Normalize (normalize) where

import Data.List (transpose)
import QPP.UOp.Syntax
import QPP.Util (permApply)

-- | Flatten associations, remove units, collect scalar phases, and simplify
-- rotations and adjacent composition factors. Interchange lifts
-- compositions out of tensors and direct sums, producing a sequence of
-- parallel layers. Tensor products of sums stay factored; no basis states
-- or Cartesian products of branches are enumerated.
normalize :: UOp -> UOp
normalize a
  | b == a    = a
  | otherwise = normalize b  -- A local rewrite can expose another simplification.
  where b = pass a

pass :: UOp -> UOp
pass op = case op of
  R a t             -> rotation (pass a) t
  Tensor xs         -> tensorOp (map pass xs)
  DirectSum xs      -> sumOp (map pass xs)
  Compose d xs      -> composeOp d (map pass xs)
  PermuteTensor ds p -> permutation PermuteTensor (product ds) ds p
  PermuteSum ds p    -> permutation PermuteSum (sum ds) ds p

-- Keep small negative angles readable: phase parameters use (-1,1],
-- non-scalar rotations (-2,2]. 'phase' then encodes its parameter as R I (-2*t).
angle :: Rational -> Rational -> Rational
angle period t = t - period * fromInteger (ceiling ((t - period/2) / period))

isIdentity :: UOp -> Bool
isIdentity (Compose _ []) = True
isIdentity _              = False

finish :: ([UOp] -> UOp) -> Dim -> [UOp] -> UOp
finish con d xs
  | d == 0 || all isIdentity xs = identity d
  | [a] <- xs                  = a
  | otherwise                  = con xs

tensorOp :: [UOp] -> UOp
tensorOp xs = scale (sum phases) body
  where
    flat = concatMap unpack xs
    unpack (Tensor ys) = ys
    unpack a           = [a]
    parts = map splitPhase flat
    phases = map fst parts
    factors = filter (/= identity 1) (map snd parts)
    body = parallel Tensor (product (map dim factors))
             (mergeIdentities (*) factors)

sumOp :: [UOp] -> UOp
sumOp xs = parallel DirectSum (sum (map dim blocks)) (mergeIdentities (+) blocks)
  where
    flat = concatMap unpack xs
    unpack (DirectSum ys) = ys
    unpack a              = [a]
    blocks = filter (/= identity 0) flat

-- Interchange, oriented towards composition on the outside. For example,
-- (A B) tensor C = (A tensor C) (B tensor I), and likewise for direct sum.
-- Pad shorter columns with identities: C must occur once, not in every layer.
parallel :: ([UOp] -> UOp) -> Dim -> [UOp] -> UOp
parallel con d xs
  | d == 0 = identity 0
  | any isComposition xs = composeOp d (map (finish con d) layers)
  | otherwise = finish con d xs
  where
    isComposition (Compose _ (_:_)) = True
    isComposition _                 = False
    column (Compose _ ys) = ys
    column a              = [a]
    columns = map column xs
    height = maximum (0 : map length columns)
    padded a ys = ys ++ replicate (height - length ys) (identity (dim a))
    layers = transpose (zipWith padded xs columns)

-- Only adjacent identities can merge: nontrivial tensor factors and blocks
-- retain their order. I_2 tensor I_3 = I_6, whereas I_2 direct-sum I_3 = I_5.
mergeIdentities :: (Dim -> Dim -> Dim) -> [UOp] -> [UOp]
mergeIdentities combine = foldr join []
  where
    join (Compose d []) (Compose e [] : rest) = identity (combine d e) : rest
    join a rest = a : rest

-- Scalar factors commute through tensor and composition, but phases of
-- different direct-sum blocks remain relative and must not be added.
-- R I t = exp(-i*pi*t/2), so the scalar phase parameter is -t/2.
splitPhase :: UOp -> (Rational, UOp)
splitPhase (R (Compose 1 []) t) = (-t/2, identity 1)
splitPhase (Tensor (R (Compose 1 []) t : xs)) =
  (-t/2, finish Tensor (product (map dim xs)) xs)
splitPhase a = (0, a)

scale :: Rational -> UOp -> UOp
scale t a
  | dim b == 0 || r == 0 = b
  | b == identity 1     = phase r
  | DirectSum xs <- b   = sumOp (map (scale r) xs)  -- Apply the SAME phase to every block.
  | Compose d (x:xs) <- b = Compose d (scale r x : xs)  -- Apply the phase once, keeping composition outside.
  | Tensor xs <- b
  , (prefix, block@DirectSum{} : suffix) <- break isSum xs =
      Tensor (prefix ++ scale r block : suffix)  -- Absorb into one block factor, even past identity padding.
  | Tensor xs <- b      = Tensor (phase r : xs)
  | otherwise           = Tensor [phase r, b]
  where
    (s, b) = splitPhase a
    r = angle 2 (t + s)
    isSum DirectSum{} = True
    isSum _           = False

rotation :: UOp -> Rational -> UOp
rotation a t
  | d == 0 || r == 0 = identity d
  | r == 2          = scale 1 (identity d)  -- R(A,2) = -I for a Hermitian unitary axis.
  | (1, b) <- splitPhase a = rotation b (-r)  -- R(-B,t) = R(B,-t).
  | abs r == 1      = scale (-r/2) a        -- R(A,+-1) = -/+i*A.
  | isIdentity a    = scale (-r/2) a
  | DirectSum xs <- a = sumOp [rotation x r | x <- xs]
  | otherwise       = R a r
  where
    d = dim a
    r = angle 4 t

composeOp :: Dim -> [UOp] -> UOp
composeOp d xs = scale (sum (map fst parts)) body
  where
    flat = concatMap unpack xs
    unpack (Compose _ ys) = ys  -- Well-formed input has the same dimension d throughout.
    unpack a              = [a]
    parts = map splitPhase flat
    body = finish (Compose d) d (foldr reduce [] (map snd parts))
    reduce a rest
      | isIdentity a = rest
    reduce a (b : rest)
      | Just c <- fuse a b = reduce c rest  -- Revisit the new neighbour after cancellation.
    reduce a rest = a : rest

fuse :: UOp -> UOp -> Maybe UOp
fuse (R a s) (R b t)
  | a == b = Just (rotation a (s + t))
fuse (PermuteTensor left q) (PermuteTensor ds p)
  | left == permApply p ds =
      Just (permutation PermuteTensor (product ds) ds (permApply q p))
fuse (PermuteSum left q) (PermuteSum ds p)
  | left == permApply p ds =
      Just (permutation PermuteSum (sum ds) ds (permApply q p))
fuse a b
  | a == dagger b = Just (identity (dim a))  -- Syntactic adjoint matching is sufficient, not complete.
  | otherwise    = Nothing

permutation :: ([Dim] -> Permutation -> UOp) -> Dim -> [Dim] -> Permutation -> UOp
permutation con d ds p
  | d <= 1 || p == [0 .. length ds - 1] = identity d
  | otherwise = con ds p
