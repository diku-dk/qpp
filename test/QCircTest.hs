-- | Flat circuits (QPP.QCirc): extraction from QOp, embedding into QOp,
--   Operator instances, and validation.
--
--   Fixed cases pin gate order and wire windows. Properties over random QOps
--   check that the round trip preserves the dense matrix, that extracted
--   circuits are well-formed layer atoms on consecutive wires, that extraction
--   is idempotent, and that the circuit instances of <>, ⊗, ⊕ and adj agree
--   with QOp's. Hand-built circuits with shuffled, non-consecutive wires
--   exercise the permutation conjugation in toQOp, including against an
--   independent CNOT specification built from the bit rule.
module Main where

import Data.Bits (bit, testBit, xor)
import Data.Ratio ((%), denominator, numerator)
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
  [ fromQOpTests, toQOpTests, extractionTests, instanceTests, circuitTests ]

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
