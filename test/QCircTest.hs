-- | Flat circuits (QPP.QCirc): extraction from QOp, embedding into QOp,
--   Operator instances, and validation.
--
--   Fixed cases pin gate order and wire windows. Properties over random QOps
--   check that the round trip preserves the dense matrix, that extracted
--   circuits are well-formed layer atoms on consecutive wires, that extraction
--   is idempotent, and that the circuit instances of <>, ⊗, ⊕ and adj agree
--   with QOp's. Hand-built circuits with shuffled, non-consecutive wires
--   exercise the permutation conjugation in toQOp, including against an
--   independent CNOT specification built from the bit rule. Layering is
--   checked for completeness, disjointness and depth bounds; the text and
--   quantikz diagrams are pinned on small circuits and checked for shape.
module Main where

import Data.Bits (bit, testBit, xor)
import Data.List (intercalate, isPrefixOf, nub, partition, sort)
import Data.Ratio ((%), denominator, numerator)
import qualified Data.Set as S
import Generators
import qualified Numeric.LinearAlgebra as H
import Programs.QFT (qft)
import QPP
import QPP.QCirc (Gate(..), QCirc(..))
import qualified QPP.QCirc as QC
import QPP.Semantics.Matrix (CMat)
import qualified QPP.Semantics.Matrix as MS
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck

main :: IO ()
main = defaultMain $ qcOpts $ testGroup "QCircTest"
  [ fromQOpTests, toQOpTests, extractionTests, instanceTests, circuitTests, rewriteTests
  , layerTests, diagramTests ]

qcOpts :: TestTree -> TestTree
qcOpts = localOption (QuickCheckTests 50) . localOption (QuickCheckMaxSize 4)

------------------------------------------------------------------------
-- Fixed cases
------------------------------------------------------------------------

fromQOpTests :: TestTree
fromQOpTests = testGroup "fromQOp fixed cases"
  [ testCase "parallel gates in wire order" $
      QC.fromQOp (H ⊗ X) @?= QCirc 2 [Gate H [0], Gate X [1]]
  , testCase "sequential gates in time order" $
      QC.fromQOp ((H ⊗ I) ∘ C X) @?= QCirc 2 [Gate (C X) [0,1], Gate H [0]]
  , testCase "a control with an idle wire inside is one gate on its window" $
      QC.fromQOp (C (Id 2 ⊗ X)) @?= QCirc 4 [Gate (C (Id 2 ⊗ X)) [0,1,2,3]]
  , testCase "adjoints are pushed into the gates" $
      QC.fromQOp (Adjoint (H ⊗ X)) @?= QCirc 2 [Gate H [0], Gate X [1]]
  , testCase "identity has no gates" $
      QC.fromQOp (Id 3) @?= QCirc 3 []
  , testCase "a global phase is a gate on no wires" $
      QC.fromQOp (Phase (1 % 4) ⊗ Id 2) @?= QCirc 2 [Gate (Phase (1 % 4)) []]
  , testCase "an identity-axis rotation keeps its phase" $
      QC.fromQOp (R (Id 1) (1 % 2)) @?= QCirc 1 [Gate (Phase ((-1) % 4)) []]
  , testCase "a buried composition is lifted into layers" $
      QC.fromQOp (X ⊗ (H ∘ Z)) @?= QCirc 2 [Gate Z [1], Gate H [1], Gate X [0]]
  , testCase "qft 3: three Hadamards, three controlled phases, one permutation" $ do
      let QCirc n gs = QC.fromQOp (qft 3)
      n @?= 3
      length gs @?= 7
      length [ () | Gate H _ <- gs ] @?= 3
      length [ () | Gate (C _) _ <- gs ] @?= 3
      last gs @?= Gate (Permute [2,1,0]) [0,1,2]
  ]

toQOpTests :: TestTree
toQOpTests = testGroup "toQOp fixed cases"
  [ testCase "consecutive wires use identity padding" $
      QC.toQOp (QCirc 4 [Gate H [2]]) @?= Id 2 ⊗ (H ⊗ Id 1)
  , testCase "no zero-width identities at the edges" $ do
      QC.toQOp (QCirc 3 [Gate H [0]]) @?= H ⊗ Id 2
      QC.toQOp (QCirc 3 [Gate H [2]]) @?= Id 2 ⊗ H
  , testCase "a phase gate on no wires" $
      QC.toQOp (QCirc 3 [Gate (Phase (1 % 4)) []]) @?= Phase (1 % 4) ⊗ Id 3
  , testCase "the empty circuit is the identity" $
      QC.toQOp (QCirc 2 []) @?= Id 2
  , testCase "time order becomes mathematical order" $
      QC.toQOp (QCirc 1 [Gate H [0], Gate X [0]]) @?= X ∘ H
  , testCase "CNOT across idle wires" $
      assertClose "CNOT 0->3"
        (MS.evalOp (QC.toQOp (QCirc 4 [Gate (C X) [0,3]])))
        (MS.evalOp (C (Id 2 ⊗ X)))
  , testCase "control below target, across idle wires" $
      assertClose "CNOT 3->0"
        (MS.evalOp (QC.toQOp (QCirc 4 [Gate (C X) [3,0]])))
        (MS.evalOp ((H ⊗ Id 2 ⊗ H) ∘ C (Id 2 ⊗ X) ∘ (H ⊗ Id 2 ⊗ H)))
  , testCase "control below target, adjacent" $
      assertClose "CNOT 1->0"
        (MS.evalOp (QC.toQOp (QCirc 2 [Gate (C X) [1,0]])))
        (MS.evalOp ((H ⊗ H) ∘ C X ∘ (H ⊗ H)))
  ]

------------------------------------------------------------------------
-- Extraction properties over random QOps
------------------------------------------------------------------------

prop_round_trip :: RandomQOp -> Property
prop_round_trip (RandomQOp _ op) =
  MS.evalOp (QC.toQOp (QC.fromQOp op)) ~~ MS.evalOp op

prop_width :: RandomQOp -> Property
prop_width (RandomQOp n op) = circ_qubits (QC.fromQOp op) === n

prop_well_formed :: RandomQOp -> Property
prop_well_formed (RandomQOp _ op) = QC.checkQCirc (QC.fromQOp op) === Right ()

-- | Every extracted gate is a layer atom on a consecutive, increasing window.
prop_layer_atoms :: RandomQOp -> Property
prop_layer_atoms (RandomQOp _ op) =
  conjoin [ counterexample (QC.showQCirc c) (atom g) | g <- circ_gates c ]
  where
    c = QC.fromQOp op
    atom (Gate g ws) = atomic g && ws == window ws
    atomic g = case g of
      Id _        -> False
      Compose _ _ -> False
      Tensor _ _  -> False
      _           -> True
    window []       = []
    window ws@(k:_) = [k .. k + length ws - 1]

prop_idempotent :: RandomQOp -> Property
prop_idempotent (RandomQOp _ op) =
  let c = QC.fromQOp op
  in QC.fromQOp (QC.toQOp c) === c

extractionTests :: TestTree
extractionTests = testGroup "Extraction properties"
  [ testProperty "toQOp . fromQOp preserves evalOp" prop_round_trip
  , testProperty "fromQOp keeps the qubit count"      prop_width
  , testProperty "fromQOp output passes checkQCirc"  prop_well_formed
  , testProperty "gates are layer atoms on windows"  prop_layer_atoms
  , testProperty "fromQOp . toQOp . fromQOp = fromQOp" prop_idempotent
  ]

------------------------------------------------------------------------
-- Operator instances agree with QOp's
------------------------------------------------------------------------

-- | Generate a same-arity companion op for a `RandomQOp`.
sameArity :: Int -> Gen QOp
sameArity n = sized $ \sz -> genQOpAt n (max 1 (min sz 3))

prop_compose :: RandomQOp -> Property
prop_compose (RandomQOp n a) =
  forAll (sameArity n) $ \b ->
    MS.evalOp (QC.toQOp (QC.fromQOp a <> QC.fromQOp b)) ~~ MS.evalOp (a ∘ b)

prop_tensor :: RandomQOp -> RandomQOp -> Property
prop_tensor (RandomQOp _ a) (RandomQOp _ b) =
  op_qubits a + op_qubits b <= 5 ==>
    MS.evalOp (QC.toQOp (QC.fromQOp a ⊗ QC.fromQOp b)) ~~ MS.evalOp (a ⊗ b)

prop_direct_sum :: RandomQOp -> Property
prop_direct_sum (RandomQOp n a) =
  forAll (sameArity n) $ \b ->
    MS.evalOp (QC.toQOp (QC.fromQOp a ⊕ QC.fromQOp b)) ~~ MS.evalOp (a ⊕ b)

prop_adjoint :: RandomQOp -> Property
prop_adjoint (RandomQOp _ a) =
  MS.evalOp (QC.toQOp (adj (QC.fromQOp a))) ~~ MS.evalOp (Adjoint a)

instanceTests :: TestTree
instanceTests = testGroup "Operator instances"
  [ testProperty "a <> b"  prop_compose
  , testProperty "a ⊗ b"   prop_tensor
  , testProperty "a ⊕ b"   prop_direct_sum
  , testProperty "adj a"   prop_adjoint
  ]

------------------------------------------------------------------------
-- Hand-built circuits: arbitrary wire lists through toQOp
------------------------------------------------------------------------

-- | CNOT on any two wires against the bit rule: flip the target bit of every
--   basis index whose control bit is set. Wire 0 is the most significant index
--   bit, as `ket` orders it. Exercises non-self-inverse wire permutations.
prop_cnot_any_wires :: Property
prop_cnot_any_wires =
  forAll (choose (2, 4)) $ \n ->
  forAll (genTwoWires n) $ \(c, t) ->
    let d = 2 ^ n :: Int
        wireBit w = n - 1 - w
        image j = if testBit j (wireBit c) then j `xor` bit (wireBit t) else j
        expected = H.assoc (d, d) 0 [ ((image j, j), 1) | j <- [0 .. d - 1] ] :: CMat
    in counterexample ("control " ++ show c ++ ", target " ++ show t) $
         MS.evalOp (QC.toQOp (QCirc n [Gate (C X) [c, t]])) ~~ expected

-- | Two distinct wires of an n-wire register, in either order.
genTwoWires :: Int -> Gen (Int, Int)
genTwoWires n = do
  c <- choose (0, n - 1)
  t <- choose (0, n - 1) `suchThat` (/= c)
  return (c, t)

-- | Random circuit on 1–4 wires: up to four gates, each on up to two distinct
--   wires in shuffled order, with a depth-1 random payload of matching arity.
newtype RandomQCirc = RandomQCirc QCirc
instance Show RandomQCirc where
  show (RandomQCirc c) = QC.showQCirc c

instance Arbitrary RandomQCirc where
  arbitrary = do
    n <- choose (1, 4)
    k <- choose (0, 4)
    gs <- vectorOf k (genGate n)
    return (RandomQCirc (QCirc n gs))
  shrink (RandomQCirc (QCirc n gs)) =
       [ RandomQCirc (QCirc n (front ++ back)) | (front, _, back) <- splits gs ]
    ++ [ RandomQCirc (QCirc n (front ++ Gate g' ws : back))
       | (front, Gate g ws, back) <- splits gs, g' <- shrinkQOp (length ws) g ]

genGate :: Int -> Gen Gate
genGate n = do
  m <- choose (0, min 2 n)
  ws <- take m <$> genPermutation n
  g <- genQOpAt m 1
  return (Gate g ws)

splits :: [a] -> [([a], a, [a])]
splits xs = [ (take i xs, x, drop (i + 1) xs) | (i, x) <- zip [0 ..] xs ]

-- Size measures for the shrinker check, as in AlgebraTest.
ratMeasure :: Rational -> Integer
ratMeasure r = 2 * (abs (numerator r) + denominator r) + (if r < 0 then 1 else 0)

qopMeasure :: QOp -> Integer
qopMeasure op = case op of
  Id _          -> 1
  X             -> 2
  Permute _     -> 2
  Y             -> 3
  Z             -> 3
  H             -> 3
  SX            -> 4
  Phase q       -> 2 + ratMeasure q
  R ax θ        -> 2 + qopMeasure ax + ratMeasure θ
  C a           -> 1 + qopMeasure a
  Adjoint a     -> 1 + qopMeasure a
  Tensor a b    -> 1 + qopMeasure a + qopMeasure b
  DirectSum a b -> 1 + qopMeasure a + qopMeasure b
  Compose a b   -> 1 + qopMeasure a + qopMeasure b

circMeasure :: QCirc -> Integer
circMeasure (QCirc _ gs) = sum [ 1 + qopMeasure g | Gate g _ <- gs ]

prop_random_well_formed :: RandomQCirc -> Property
prop_random_well_formed (RandomQCirc c) = QC.checkQCirc c === Right ()

prop_random_unitary :: RandomQCirc -> Property
prop_random_unitary (RandomQCirc c) =
  let u = MS.evalOp (QC.toQOp c)
  in (u H.<> H.tr u) ~~ identCMat (H.rows u)

prop_random_round_trip :: RandomQCirc -> Property
prop_random_round_trip (RandomQCirc c) =
  let q = QC.toQOp c
  in MS.evalOp (QC.toQOp (QC.fromQOp q)) ~~ MS.evalOp q

prop_shrink_decreases :: RandomQCirc -> Property
prop_shrink_decreases (RandomQCirc c) =
  conjoin [ counterexample (QC.showQCirc c') (circMeasure c' < circMeasure c)
          | RandomQCirc c' <- shrink (RandomQCirc c) ]

circuitTests :: TestTree
circuitTests = testGroup "Hand-built circuits"
  [ testProperty "CNOT on any two wires matches the bit rule" prop_cnot_any_wires
  , testProperty "random circuits pass checkQCirc"            prop_random_well_formed
  , testProperty "toQOp of a random circuit is unitary"       prop_random_unitary
  , testProperty "toQOp . fromQOp . toQOp preserves evalOp"   prop_random_round_trip
  , testProperty "shrink candidates are strictly smaller"     prop_shrink_decreases
  ]

------------------------------------------------------------------------
-- Rewrite passes
------------------------------------------------------------------------

-- | A pass preserves the dense matrix of the circuit.
preserves :: (QCirc -> QCirc) -> QCirc -> Property
preserves pass c = MS.evalOp (QC.toQOp (pass c)) ~~ MS.evalOp (QC.toQOp c)

passTests :: String -> (QCirc -> QCirc) -> [TestTree] -> TestTree
passTests name pass extra = testGroup name $
  [ testProperty "preserves evalOp on extracted circuits" $
      \(RandomQOp _ op) -> preserves pass (QC.fromQOp op)
  , testProperty "preserves evalOp on hand-built circuits" $
      \(RandomQCirc c) -> preserves pass c
  , testProperty "preserves evalOp on permutation circuits" $
      \(PermCirc c) -> preserves pass c
  , testProperty "preserves evalOp on diagonal circuits" $
      \(DiagCirc c) -> preserves pass c
  ] ++ extra

-- | Circuits mixing permutations on shuffled wire subsets with one-qubit
--   gates, so that absorbing permutations meets cycles and nested relabelings.
newtype PermCirc = PermCirc QCirc
instance Show PermCirc where
  show (PermCirc c) = QC.showQCirc c

instance Arbitrary PermCirc where
  arbitrary = do
    n <- choose (2, 4)
    k <- choose (0, 5)
    gs <- vectorOf k (oneof [permGate n, oneQubit n])
    return (PermCirc (QCirc n gs))
    where
      permGate n = do
        m <- choose (2, n)
        ws <- take m <$> genPermutation n
        ks <- genPermutation m
        return (Gate (Permute ks) ws)
      oneQubit n = do
        w <- choose (0, n - 1)
        g <- genQOpBase 1
        return (Gate g [w])
  shrink (PermCirc (QCirc n gs)) =
    [ PermCirc (QCirc n (front ++ back)) | (front, _, back) <- splits gs ]

-- | Circuits of composite diagonal gates on two to four shuffled wires, so
--   that phase reordering, idle-wire detection and fusion across wire orders
--   meet three and four wires, not only the two of 'RandomQCirc'.
newtype DiagCirc = DiagCirc QCirc
instance Show DiagCirc where
  show (DiagCirc c) = QC.showQCirc c

instance Arbitrary DiagCirc where
  arbitrary = do
    n <- choose (2, 4)
    k <- choose (1, 4)
    gs <- vectorOf k (diagonalGate n)
    return (DiagCirc (QCirc n gs))
    where
      diagonalGate n = do
        m <- choose (2, n)
        ws <- take m <$> genPermutation n
        g <- genDiagonal m 2
        return (Gate g ws)
  shrink (DiagCirc (QCirc n gs)) =
    [ DiagCirc (QCirc n (front ++ back)) | (front, _, back) <- splits gs ]

-- | A structurally diagonal operator on @m@ qubits of syntactic depth at most
--   @d@: identities, phases, @Z@, rotations about tensor products of @Z@ and
--   @Id 1@, and controls, sums, products, compositions and adjoints of those.
--   Every axis is Hermitian, so the phase list of the operator exists.
genDiagonal :: Int -> Int -> Gen QOp
genDiagonal m 0 = leaf
  where
    leaf = oneof $
         [ return (Id m) ]
      ++ [ Phase <$> genRat | m == 0 ]
      ++ [ return Z | m == 1 ]
      ++ [ R . foldr1 Tensor <$> vectorOf m (elements [Z, Id 1]) <*> genRat | m > 0 ]
genDiagonal m d = frequency $
     [ (1, genDiagonal m 0)
     , (2, Compose <$> smaller m <*> smaller m)
     , (1, Adjoint <$> smaller m)
     , (2, do k <- choose (0, m); Tensor <$> smaller k <*> smaller (m - k))
     ]
  ++ [ (3, C <$> smaller (m - 1)) | m > 0 ]
  ++ [ (2, DirectSum <$> smaller (m - 1) <*> smaller (m - 1)) | m > 0 ]
  where
    smaller k = genDiagonal k (d - 1)

rewriteTests :: TestTree
rewriteTests = testGroup "Rewrite passes"
  [ passTests "flatten" QC.flatten flattenCases
  , passTests "compact" QC.compact compactCases
  , passTests "absorbPermutations" QC.absorbPermutations absorbCases
  , passTests "fuse" QC.fuse fuseCases
  , passTests "collectPhases" QC.collectPhases phaseCases
  , passTests "simplify" QC.simplify simplifyCases
  ]

flattenCases :: [TestTree]
flattenCases =
  [ testCase "controls distribute over composition" $
      flat [Gate (C (H ∘ X)) [0,1]] @?= [Gate (C X) [0,1], Gate (C H) [0,1]]
  , testCase "idle wires inside a control disappear" $
      flat [Gate (C (Id 2 ⊗ X)) [0,1,2,3]] @?= [Gate (C X) [0,3]]
  , testCase "a direct sum becomes X-conjugated controls" $
      flat [Gate (DirectSum X Z) [0,1]]
        @?= [Gate X [0], Gate (C X) [0,1], Gate X [0], Gate (C Z) [0,1]]
  , testCase "an identity branch needs no conjugation" $
      flat [Gate (DirectSum (Id 1) X) [0,1]] @?= [Gate (C X) [0,1]]
  , testCase "nested controls distribute over tensor factors" $
      flat [Gate (C (C (X ⊗ Z))) [0,1,2,3]]
        @?= [Gate (C (C X)) [0,1,2], Gate (C (C Z)) [0,1,3]]
  , testCase "a controlled phase gate is one canonical controlled phase" $
      flat [Gate (C (Phase (1 % 4) ⊗ R Z (1 % 2))) [2,0]]
        @?= [Gate (C (C (Phase (1 % 2)))) [0,2]]
  , testCase "a sum of phases splits a global phase off a P gate" $
      flat [Gate (Phase (1 % 3) ⊕ Phase (1 % 2)) [0]]
        @?= [Gate (Phase (1 % 3)) [], Gate (C (Phase (1 % 6))) [0]]
  , testCase "equal branches make the wire idle" $
      flat [Gate (Phase (1 % 3) ⊕ Phase (1 % 3)) [0]] @?= [Gate (Phase (1 % 3)) []]
  , testCase "composed rotations canonicalize" $
      flat [Gate (R Z (1 % 3) ∘ R Z (1 % 6)) [0]]
        @?= [Gate (Phase ((-1) % 4)) [], Gate (C (Phase (1 % 2))) [0]]
  , testCase "a diagonal with a nontrivial first branch" $
      flat [Gate (DirectSum Z (Id 1)) [0,1]] @?= [Gate (DirectSum (C (Phase 1)) (Id 1)) [0,1]]
  , testCase "atoms under controls stay as written" $
      flat [Gate (C (R Z (1 % 3))) [3,0]] @?= [Gate (C (R Z (1 % 3))) [3,0]]
  , testCase "a controlled phase past an idle wire, on shuffled wires" $
      flat [Gate (C (Id 1 ⊗ C (Phase (1 % 3)))) [2,0,1]] @?= [Gate (C (C (Phase (1 % 3)))) [1,2]]
  , testCase "an asymmetric diagonal on reversed wires is re-indexed" $
      flat [Gate (C (Phase (1 % 3) ⊕ Phase (1 % 2))) [1,0]]
        @?= [Gate (DirectSum (C (Phase (1 % 3))) (C (Phase (1 % 2)))) [0,1]]
  , testCase "a wide composite diagonal falls back to structural unfolding" $ do
      let wide = QC.flatten (QCirc 13 [Gate (iterate C (Phase (1 % 4) ⊗ R Z (1 % 2)) !! 12) [0 .. 12]])
      length (circ_gates wide) @?= 2
      assertBool "elementary" (all QC.elementary (circ_gates wide))
  , testCase "adjoints are resolved" $
      flat [Gate (Adjoint (C (R Z (1 % 3)))) [0,1]] @?= [Gate (C (R Z ((-1) % 3))) [0,1]]
  , testProperty "every gate is elementary" $
      \(RandomQOp _ op) -> all QC.elementary (circ_gates (QC.flatten (QC.fromQOp op)))
  ]
  where flat gs = circ_gates (QC.flatten (QCirc 4 gs))

compactCases :: [TestTree]
compactCases =
  [ testCase "identity padding inside a control" $
      comp [Gate (C (Id 2 ⊗ X)) [0,1,2,3]] @?= [Gate (C X) [0,3]]
  , testCase "identity factors in a rotation axis" $
      comp [Gate (R (X ⊗ Id 1 ⊗ Z) (1 % 3)) [0,1,2]] @?= [Gate (R (X ⊗ Z) (1 % 3)) [0,2]]
  , testCase "fixed points of a permutation" $
      comp [Gate (Permute [0,2,1]) [0,1,2]] @?= [Gate (Permute [1,0]) [1,2]]
  , testCase "identity gates are dropped" $
      comp [Gate (Id 3) [0,1,2], Gate (R X 0) [3]] @?= []
  , testCase "wires keep their order" $
      comp [Gate (C (Id 1 ⊗ X)) [3,1,0]] @?= [Gate (C X) [3,0]]
  , testProperty "every gate has full support" $
      \(RandomQOp _ op) -> all fullSupport (circ_gates (QC.compact (QC.fromQOp op)))
  ]
  where
    comp gs = circ_gates (QC.compact (QCirc 4 gs))
    fullSupport (Gate g ws) = S.toAscList (op_support g) == [0 .. length ws - 1]

absorbCases :: [TestTree]
absorbCases =
  [ testCase "a swap relabels the later gate and trails the circuit" $
      QC.absorbPermutations (QCirc 2 [Gate (Permute [1,0]) [0,1], Gate X [0]])
        @?= QCirc 2 [Gate X [1], Gate (Permute [1,0]) [0,1]]
  , testCase "a swap and its inverse cancel around a gate" $
      QC.absorbPermutations
        (QCirc 3 [Gate (Permute [1,0]) [0,2], Gate H [2], Gate (Permute [1,0]) [0,2]])
        @?= QCirc 3 [Gate H [0]]
  , testCase "a trailing permutation on all wires is unchanged" $
      let c = QC.fromQOp (qft 3) in QC.absorbPermutations c @?= c
  , testProperty "no permutation remains except a last one on all wires" $
      \(PermCirc c) -> noPermutes (QC.absorbPermutations c)
  ]
  where
    noPermutes (QCirc n gs) = case reverse gs of
      Gate (Permute _) ws : earlier -> ws == [0 .. n - 1] && not (any isPermute earlier)
      earlier                       -> not (any isPermute earlier)
    isPermute (Gate (Permute _) _) = True
    isPermute _                    = False

fuseCases :: [TestTree]
fuseCases =
  [ testCase "self-inverse gates cancel" $
      fu 1 [Gate X [0], Gate X [0]] @?= []
  , testCase "cancellation across a gate on other wires" $
      fu 2 [Gate X [0], Gate H [1], Gate X [0]] @?= [Gate H [1]]
  , testCase "rotations about the same axis add" $
      fu 1 [Gate (R Z (1 % 3)) [0], Gate (R Z (1 % 3)) [0]] @?= [Gate (R Z (2 % 3)) [0]]
  , testCase "a half turn becomes a global phase" $
      fu 1 [Gate (R Z 1) [0], Gate (R Z 1) [0]] @?= [Gate (Phase 1) []]
  , testCase "a controlled half turn becomes a phase on the control" $
      fu 2 [Gate (C (R Z 1)) [0,1], Gate (C (R Z 1)) [0,1]] @?= [Gate (C (Phase 1)) [0]]
  , testCase "controlled gates cancel" $
      fu 2 [Gate (C X) [0,1], Gate (C X) [0,1]] @?= []
  , testCase "permutations compose" $
      fu 3 [Gate (Permute [1,2,0]) [0,1,2], Gate (Permute [1,2,0]) [0,1,2]]
        @?= [Gate (Permute [2,0,1]) [0,1,2]]
  , testCase "adjoint pairs cancel" $
      fu 1 [Gate SX [0], Gate (Adjoint SX) [0]] @?= []
  , testCase "square roots of X fuse" $
      fu 1 [Gate SX [0], Gate SX [0]] @?= [Gate X [0]]
  , testCase "controls fuse in either order" $
      fu 3 [Gate (C (C X)) [0,1,2], Gate (C (C X)) [1,0,2]] @?= []
  , testCase "swaps on reversed wires cancel" $
      fu 2 [Gate (Permute [1,0]) [0,1], Gate (Permute [1,0]) [1,0]] @?= []
  , testCase "a permutation is re-expressed on ascending wires" $
      QC.orderWires (Gate (Permute [1,2,0]) [2,0,1]) @?= Gate (Permute [1,2,0]) [0,1,2]
  , testCase "tensor factors follow their wires" $
      QC.orderWires (Gate (C (X ⊗ Z)) [2,1,0]) @?= Gate (C (Z ⊗ X)) [2,0,1]
  , testCase "rotation axes follow their wires" $
      fu 2 [Gate (R (X ⊗ Z) (1 % 3)) [0,1], Gate (R (Z ⊗ X) (1 % 3)) [1,0]]
        @?= [Gate (R (X ⊗ Z) (2 % 3)) [0,1]]
  , testCase "phase gates on the same wire add structurally" $
      fu 1 [Gate (C (Phase (1 % 4))) [0], Gate (C (Phase (1 % 4))) [0]]
        @?= [Gate (C (Phase (1 % 2))) [0]]
  , testCase "Z and its P-gate spelling cancel" $
      fu 1 [Gate Z [0], Gate (C (Phase 1)) [0]] @?= []
  , testCase "a rotation and a phase gate fuse into canonical form" $
      fu 1 [Gate (R Z (1 % 2)) [0], Gate (C (Phase (1 % 2))) [0]]
        @?= [Gate (Phase ((-1) % 4)) [], Gate (C (Phase 1)) [0]]
  , testCase "controlled phases fuse across wire order" $
      fu 4 [Gate (C (C (Phase (1 % 4)))) [0,3], Gate (C (C (Phase (1 % 4)))) [3,0]]
        @?= [Gate (C (C (Phase (1 % 2)))) [0,3]]
  , testCase "a diagonal gate passes earlier diagonal gates" $
      fu 2 [Gate (C Z) [0,1], Gate (C (Phase 1)) [0], Gate (C Z) [0,1]]
        @?= [Gate (C (Phase 1)) [0]]
  , testCase "a non-diagonal gate blocks" $
      fu 1 [Gate Z [0], Gate H [0], Gate Z [0]] @?= [Gate Z [0], Gate H [0], Gate Z [0]]
  , testCase "no rule, no change" $
      fu 1 [Gate H [0], Gate X [0], Gate H [0]] @?= [Gate H [0], Gate X [0], Gate H [0]]
  , testCase "control and target do not swap" $
      fu 2 [Gate (C X) [0,1], Gate (C X) [1,0]] @?= [Gate (C X) [0,1], Gate (C X) [1,0]]
  , testProperty "diagonal gates never share a wire set after fusing" $
      \(DiagCirc c) ->
        let sets = [ sort ws | Gate _ ws <- circ_gates (QC.fuse c), not (null ws) ]
        in sets === nub sets
  ]
  where fu n gs = circ_gates (QC.fuse (QCirc n gs))

phaseCases :: [TestTree]
phaseCases =
  [ testCase "phases add and move to the front" $
      ph [Gate (Phase (1 % 4)) [], Gate X [0], Gate (Phase (3 % 4)) []]
        @?= [Gate (Phase 1) [], Gate X [0]]
  , testCase "a full turn vanishes" $
      ph [Gate (Phase 1) [], Gate (Phase 1) []] @?= []
  , testCase "angles reduce into (-1, 1]" $
      ph [Gate (Phase (3 % 2)) []] @?= [Gate (Phase ((-1) % 2)) []]
  , testCase "controlled phases are relative and stay" $
      ph [Gate (C (Phase (1 % 2))) [0]] @?= [Gate (C (Phase (1 % 2))) [0]]
  , testProperty "at most one gate on no wires, and it comes first" $
      \(RandomQOp _ op) ->
        let gs = circ_gates (QC.collectPhases (QC.fromQOp op))
        in length (filter scalar gs) <= 1 && not (any scalar (drop 1 gs))
  ]
  where
    ph gs = circ_gates (QC.collectPhases (QCirc 1 gs))
    scalar (Gate _ ws) = null ws

simplifyCases :: [TestTree]
simplifyCases =
  [ testCase "qft 4: elementary gates and one trailing permutation" $ do
      let QCirc n gs = QC.simplify (QC.fromQOp (qft 4))
      n @?= 4
      length gs @?= 11
      assertBool "elementary" (all QC.elementary gs)
      last gs @?= Gate (Permute [3,2,1,0]) [0,1,2,3]
  , testProperty "idempotent on extracted circuits" $
      \(RandomQOp _ op) -> let c = QC.simplify (QC.fromQOp op) in QC.simplify c === c
  , testCase "an adjoint square root of X fuses with its inverse pair" $
      QC.simplify (QCirc 1 [Gate (Adjoint SX) [0], Gate SX [0], Gate SX [0]])
        @?= QCirc 1 [Gate SX [0]]
  , testProperty "idempotent on hand-built circuits" $
      \(RandomQCirc c) -> let c' = QC.simplify c in QC.simplify c' === c'
  , testProperty "idempotent on diagonal circuits" $
      \(DiagCirc c) -> let c' = QC.simplify c in QC.simplify c' === c'
  ]

------------------------------------------------------------------------
-- Layers
------------------------------------------------------------------------

-- | The gates on a wire, in the order given.
onWire :: Int -> [Gate] -> [Gate]
onWire w = filter ((w `elem`) . QC.gateWires)

prop_layers_complete :: RandomQCirc -> Property
prop_layers_complete (RandomQCirc c) =
  length flat === length gs .&&.
    conjoin [ onWire w flat === onWire w gs | w <- [0 .. circ_qubits c - 1] ]
  where
    gs = circ_gates c
    flat = concat (QC.layers c)

prop_layers_disjoint :: RandomQCirc -> Bool
prop_layers_disjoint (RandomQCirc c) =
  all (distinct . concatMap QC.gateWires) (QC.layers c)
    && all (distinct . concatMap QC.gateSpan) (QC.layersBy QC.gateSpan (circ_gates c))
  where
    distinct ws = nub ws == ws

prop_depth_bounds :: RandomQCirc -> Property
prop_depth_bounds (RandomQCirc c) =
  counterexample ("depth " ++ show d) (busiest <= d && d <= length (filter (not . null . QC.gateWires) gs))
  where
    gs = circ_gates c
    d = QC.depth c
    busiest = maximum (0 : [ length (onWire w gs) | w <- [0 .. circ_qubits c - 1] ])

layerTests :: TestTree
layerTests = testGroup "Layers"
  [ testCase "parallel gates share a layer" $
      QC.layers (QCirc 2 [Gate H [0], Gate H [1]]) @?= [[Gate H [0], Gate H [1]]]
  , testCase "a gate follows every earlier gate on its wires; a scalar joins the first layer" $
      QC.layers (QCirc 2 [Gate (Phase (1 % 2)) [], Gate H [1], Gate (C X) [0,1], Gate H [0]])
        @?= [[Gate (Phase (1 % 2)) [], Gate H [1]], [Gate (C X) [0,1]], [Gate H [0]]]
  , testCase "a scalar alone has depth 0" $
      QC.depth (QCirc 1 [Gate (Phase (1 % 2)) []]) @?= 0
  , testCase "qft 3 simplified has depth 6" $
      QC.depth (QC.simplify (QC.fromQOp (qft 3))) @?= 6
  , testCase "spans block the wires a gate crosses" $
      QC.layersBy QC.gateSpan [Gate (C X) [0,2], Gate H [1]] @?= [[Gate (C X) [0,2]], [Gate H [1]]]
  , testProperty "layers keep every gate and the order on each wire" prop_layers_complete
  , testProperty "gates of a layer are on disjoint wires, spans on disjoint rows" prop_layers_disjoint
  , testProperty "depth lies between the busiest wire and the gate count" prop_depth_bounds
  ]

------------------------------------------------------------------------
-- Diagrams
------------------------------------------------------------------------

block :: [String] -> String
block = intercalate "\n"

textCases :: [TestTree]
textCases =
  [ testCase "a Bell pair" $
      QC.showQCircDiagram (QCirc 2 [Gate H [0], Gate (C X) [0,1]]) @?= block
        [ "q0: ──H──●──"
        , "         │"
        , "q1: ─────⊕──" ]
  , testCase "a control below its target crosses the wire between" $
      QC.showQCircDiagram (QCirc 3 [Gate (C X) [2,0]]) @?= block
        [ "q0: ──⊕──"
        , "      │"
        , "q1: ──┼──"
        , "      │"
        , "q2: ──●──" ]
  , testCase "a swap, a controlled swap and a Toffoli with its target between its controls" $
      QC.showQCircDiagram (QCirc 3 [ Gate (Permute [1,0]) [0,1], Gate (C (Permute [1,0])) [0,1,2]
                                   , Gate (C (C X)) [0,2,1] ]) @?= block
        [ "q0: ──×──●──●──"
        , "      │  │  │"
        , "q1: ──×──×──⊕──"
        , "         │  │"
        , "q2: ─────×──●──" ]
  , testCase "phase gates and controlled phases in either wire order" $
      QC.showQCircDiagram (QCirc 2 [ Gate (C (C (Phase (1 % 4)))) [0,1], Gate (C (Phase (1 % 2))) [1]
                                   , Gate (C (C (Phase ((-3) % 4)))) [1,0] ]) @?= block
        [ "q0: ────●─────────────P(-3π/4)──"
        , "        │                │"
        , "q1: ──P(π/4)──P(π/2)─────●──────" ]
  , testCase "a box on non-consecutive wires numbers its other inputs" $
      QC.showQCircDiagram (QCirc 3 [Gate (R (X ⊗ Y) (1 % 3)) [2,0]]) @?= block
        [ "q0: ───────#1────────"
        , "            │"
        , "q1: ────────┼────────"
        , "            │"
        , "q2: ──R(X ⊗ Y, π/3)──" ]
  , testCase "a control inside its target's rows makes the whole gate a box" $
      QC.showQCircDiagram (QCirc 3 [Gate (C (R (X ⊗ Y) (1 % 3))) [1,0,2]]) @?= block
        [ "q0: ─────────#1─────────"
        , "             │"
        , "q1: ──C(R(X ⊗ Y, π/3))──"
        , "             │"
        , "q2: ─────────#2─────────" ]
  , testCase "scalars come first" $
      QC.showQCircDiagram (QCirc 1 [Gate (Phase (3 % 4)) [], Gate H [0], Gate (Phase 1 ∘ Phase (1 % 2)) []])
        @?= block [ "scalar: Phase(3π/4)", "scalar: Phase(π) ∘ Phase(π/2)", "q0: ──H──" ]
  , testCase "the empty circuit is bare wires" $
      QC.showQCircDiagram (QCirc 2 []) @?= block [ "q0: ──", "", "q1: ──" ]
  ]

verticalCases :: [TestTree]
verticalCases =
  [ testCase "a Bell pair" $
      QC.showQCircDiagramVertical (QCirc 2 [Gate H [0], Gate (C X) [0,1]]) @?= block
        [ "q0  q1"
        , "│   │"
        , "H   │"
        , "│   │"
        , "●───⊕"
        , "│   │" ]
  , testCase "a connector crosses the wire between; columns take their widest cell" $
      QC.showQCircDiagramVertical (QCirc 3 [Gate (C (C (Phase (1 % 4)))) [0,2], Gate H [1]]) @?= block
        [ "q0  q1    q2"
        , "│   │     │"
        , "●───┼───P(π/4)"
        , "│   │     │"
        , "│   H     │"
        , "│   │     │" ]
  ]

-- | After the scalar lines, a header and a bare-wire row, then a layer row
--   and a bare-wire row per span layer, all bare-wire rows alike.
prop_vertical_shape :: RandomQCirc -> Property
prop_vertical_shape (RandomQCirc c) =
  counterexample diagram $
    length rows === 2 + 2 * length (QC.layersBy QC.gateSpan wired)
      .&&. nub [ r | (i, r) <- zip [0 :: Int ..] rows, odd i ] === take 1 (drop 1 rows)
  where
    diagram = QC.showQCircDiagramVertical c
    (scalars, wired) = partition (null . QC.gateWires) (circ_gates c)
    rows = drop (length scalars) (lines diagram)

quantikzCases :: [TestTree]
quantikzCases =
  [ testCase "a Bell pair" $
      QC.showQuantikz (QCirc 2 [Gate H [0], Gate (C X) [0,1]]) @?= block
        [ "\\begin{quantikz}"
        , "\\lstick{$q_{0}$} & \\gate{H} & \\ctrl{1} & \\qw \\\\"
        , "\\lstick{$q_{1}$} & \\qw & \\targ{} & \\qw"
        , "\\end{quantikz}" ]
  , testCase "controls link towards the target from both sides" $
      QC.showQuantikz (QCirc 3 [ Gate (Permute [1,0]) [0,1], Gate (C (Permute [1,0])) [0,1,2]
                               , Gate (C (C X)) [0,2,1] ]) @?= block
        [ "\\begin{quantikz}"
        , "\\lstick{$q_{0}$} & \\swap{1} & \\ctrl{1} & \\ctrl{1} & \\qw \\\\"
        , "\\lstick{$q_{1}$} & \\targX{} & \\swap{1} & \\targ{} & \\qw \\\\"
        , "\\lstick{$q_{2}$} & \\qw & \\targX{} & \\ctrl{-1} & \\qw"
        , "\\end{quantikz}" ]
  , testCase "a box on non-consecutive wires spans them with a line through" $
      QC.showQuantikz (QCirc 3 [Gate (R (X ⊗ Y) (1 % 3)) [2,0]]) @?= block
        [ "\\begin{quantikz}"
        , "\\lstick{$q_{0}$} & \\gate[wires=3]{R_{X \\otimes Y}(\\pi/3)}\\gateinput{1} & \\qw \\\\"
        , "\\lstick{$q_{1}$} & \\linethrough & \\qw \\\\"
        , "\\lstick{$q_{2}$} & \\gateinput{0} & \\qw"
        , "\\end{quantikz}" ]
  , testCase "a phase gate and a scalar comment" $
      QC.showQuantikz (QCirc 1 [Gate (Phase (3 % 4)) [], Gate (C (Phase (1 % 4))) [0]]) @?= block
        [ "% scalar: e^{i\\,3\\pi/4}"
        , "\\begin{quantikz}"
        , "\\lstick{$q_{0}$} & \\phase{\\pi/4} & \\qw"
        , "\\end{quantikz}" ]
  ]

-- | A row per wire and a gap row between neighbours after the scalar lines;
--   wire rows start with their name and have one width, gap rows are shorter.
prop_text_shape :: RandomQCirc -> Property
prop_text_shape (RandomQCirc c) =
  counterexample diagram $
    length rows === 2 * n - 1
      .&&. and [ ("q" ++ show i ++ ":") `isPrefixOf` row | (i, row) <- zip [0 :: Int ..] wireRows ]
      .&&. all ((== width) . length) wireRows
      .&&. all ((< width) . length) gapRows
  where
    n = circ_qubits c
    diagram = QC.showQCircDiagram c
    scalars = filter (null . QC.gateWires) (circ_gates c)
    rows = drop (length scalars) (lines diagram)
    wireRows = [ r | (i, r) <- zip [0 :: Int ..] rows, even i ]
    gapRows = [ r | (i, r) <- zip [0 :: Int ..] rows, odd i ]
    width = maximum (map length wireRows)

-- | A comment per scalar, a row per wire, and in every row one cell per span
--   layer plus the stick and the trailing wire.
prop_quantikz_shape :: RandomQCirc -> Property
prop_quantikz_shape (RandomQCirc c) =
  counterexample fragment $
    length comments === length scalars
      .&&. length rows === circ_qubits c
      .&&. map cells rows === replicate (circ_qubits c) (length (QC.layersBy QC.gateSpan wired) + 2)
  where
    fragment = QC.showQuantikz c
    (scalars, wired) = partition (null . QC.gateWires) (circ_gates c)
    (comments, rest) = break (== "\\begin{quantikz}") (lines fragment)
    rows = takeWhile (/= "\\end{quantikz}") (drop 1 rest)
    cells = (+ 1) . length . filter (== '&')

diagramTests :: TestTree
diagramTests = testGroup "Diagrams"
  [ testGroup "text" textCases
  , testGroup "vertical" verticalCases
  , testGroup "quantikz" quantikzCases
  , testProperty "text diagrams have a row per wire and gap, of one width" prop_text_shape
  , testProperty "vertical diagrams have a row per layer and bare wires between" prop_vertical_shape
  , testProperty "quantikz rows have one cell per column"               prop_quantikz_shape
  ]
