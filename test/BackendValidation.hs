{-# LANGUAGE ScopedTypeVariables #-}
module Main where

import QPP
import QPP.Util (permInvert)
import QPP.Semantics (Convertible(..))
import qualified QPP.Semantics.Matrix      as MS
import qualified QPP.Semantics.MPS         as MPS
import qualified QPP.Semantics.Statevector as SV
import QPP.Semantics.Matrix (CMat)
import Programs.QFT (qft)

import Numeric.LinearAlgebra (rows)
import qualified Numeric.LinearAlgebra as H

import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck

import Generators

------------------------------------------------------------------------
-- Properties
------------------------------------------------------------------------

-- | The MPS validator produces the same matrix as the Matrix backend's evalOp.
prop_mps_matches_ms :: RandomQOp -> Property
prop_mps_matches_ms (RandomQOp _ op) =
  let mref = MS.evalOp op
      mmps = to (MPS.evalOp op) :: CMat
  in mref ~~ mmps

-- | The Statevector validator produces the same matrix as the Matrix backend's evalOp.
prop_sv_matches_ms :: RandomQOp -> Property
prop_sv_matches_ms (RandomQOp _ op) =
  let mref = MS.evalOp op
      msv  = to (SV.evalOp op) :: CMat
  in mref ~~ msv

-- | MPS evalOp output is unitary: M · M† = I.
prop_unitary_mps :: RandomQOp -> Property
prop_unitary_mps (RandomQOp _ op) =
  let m = to (MPS.evalOp op) :: CMat
      n = rows m
  in (m H.<> H.tr m) ~~ identCMat n

-- | op ∘ adj op = I, evaluated under MPS.
prop_inverse_mps :: RandomQOp -> Property
prop_inverse_mps (RandomQOp n op) =
  let m = to (MPS.evalOp (Compose op (Adjoint op))) :: CMat
  in m ~~ identCMat (2 ^ n)

-- | adj (adj op) = op, evaluated under MPS.
prop_double_adjoint_mps :: RandomQOp -> Property
prop_double_adjoint_mps (RandomQOp _ op) =
  let m1 = to (MPS.evalOp op) :: CMat
      m2 = to (MPS.evalOp (Adjoint (Adjoint op))) :: CMat
  in m1 ~~ m2

-- | R ax 0 = Id (for any Pauli string ax).
prop_R_zero :: RotGen -> Property
prop_R_zero (RotGen n ax _) =
  let m = to (MPS.evalOp (R ax 0)) :: CMat
  in m ~~ identCMat (2 ^ n)

-- | R ax (2θ) = (R ax θ) ∘ (R ax θ).
prop_R_doubling :: RotGen -> Property
prop_R_doubling (RotGen _ ax θ) =
  let lhs = to (MPS.evalOp (R ax (2*θ))) :: CMat
      rhs = to (MPS.evalOp (Compose (R ax θ) (R ax θ))) :: CMat
  in lhs ~~ rhs

-- | R ax θ ∘ R ax (-θ) = Id.
prop_R_inverse :: RotGen -> Property
prop_R_inverse (RotGen n ax θ) =
  let m = to (MPS.evalOp (Compose (R ax θ) (R ax (-θ)))) :: CMat
  in m ~~ identCMat (2 ^ n)

-- | Permute π ∘ Permute (permInvert π) = Id.
prop_permute_inverse :: Perm -> Property
prop_permute_inverse (Perm p) =
  let n = length p
      m = to (MPS.evalOp (Compose (Permute p) (Permute (permInvert p)))) :: CMat
  in m ~~ identCMat (2 ^ n)

------------------------------------------------------------------------
-- Unit tests (fixed cases for regression coverage)
------------------------------------------------------------------------

-- | Compare op's matrix across all three backends.
caseAllBackends :: String -> QOp -> TestTree
caseAllBackends lbl op = testCase lbl $ do
  let mref = MS.evalOp op
      mmps = to (MPS.evalOp op) :: CMat
      msv  = to (SV.evalOp op)  :: CMat
  assertClose ("MPS vs MS: " ++ lbl) mref mmps
  assertClose ("SV  vs MS: " ++ lbl) mref msv

primTests :: TestTree
primTests = testGroup "Primitives"
  [ caseAllBackends "I"  (Id 1)
  , caseAllBackends "X"  X
  , caseAllBackends "Y"  Y
  , caseAllBackends "Z"  Z
  , caseAllBackends "H"  H
  , caseAllBackends "SX" SX
  ]

phaseAngles :: [Rational]
phaseAngles = [0, 1/4, 1/2, 1, 2/3, -1/4, 5/4]

rotAngles :: [Rational]
rotAngles = [0, 1/8, 1/4, 1/3, 1/2, 2/3, 3/4, 1, 5/4, 2, -1/4, -1/2]

phaseTests :: TestTree
phaseTests = testGroup "Phase"
  [ caseAllBackends ("Phase " ++ show q) (Phase q) | q <- phaseAngles ]

idTests :: TestTree
idTests = testGroup "Identity"
  [ caseAllBackends ("Id " ++ show n) (Id n) | n <- [0, 1, 2, 3] ]

rotTests :: TestTree
rotTests = testGroup "Single-qubit rotations"
  [ caseAllBackends ("R " ++ axn ++ " " ++ show q) (R ax q)
  | (axn, ax) <- [("X", X), ("Y", Y), ("Z", Z)]
  , q <- rotAngles
  ]

multiRotTests :: TestTree
multiRotTests = testGroup "Multi-qubit rotations" $
  [ caseAllBackends ("R (X⊗Z) "   ++ show q) (R (Tensor X Z) q)
    | q <- [1/4, 1/3, 1/2, -1/4 :: Rational] ]
  ++
  [ caseAllBackends ("R (X⊗Y⊗Z) " ++ show q) (R (Tensor X (Tensor Y Z)) q)
    | q <- [1/4, 1/3, 1/2, -1/4 :: Rational] ]

permTests :: TestTree
permTests = testGroup "Permutations"
  [ caseAllBackends ("Permute " ++ show p) (Permute p)
  | p <- [[1,0], [2,0,1], [0,2,1], [3,1,2,0]]
  ]

tensorCases :: [(String, QOp)]
tensorCases =
  [ ("X ⊗ H",            Tensor X H)
  , ("X ⊗ Y ⊗ Z",        Tensor X (Tensor Y Z))
  , ("H ⊗ Permute[1,0]", Tensor H (Permute [1,0]))
  , ("(R Z 1/4) ⊗ H",    Tensor (R Z (1/4)) H)
  ]

composeCases :: [(String, QOp)]
composeCases =
  [ ("H ∘ X",             Compose H X)
  , ("X ∘ H ∘ X",         Compose X (Compose H X))
  , ("(R Z 1/3) ∘ H",     Compose (R Z (1/3)) H)
  , ("(H ⊗ I) ∘ (I ⊗ H)", Compose (Tensor H (Id 1)) (Tensor (Id 1) H))
  ]

controlledCases :: [(String, QOp)]
controlledCases =
  [ ("C X",              C X)
  , ("C H",              C H)
  , ("C (X ⊗ Y)",        C (Tensor X Y))
  , ("C (R Z 1/4)",      C (R Z (1/4)))
  , ("C (Permute [1,0])", C (Permute [1,0]))
  ]

dsumCases :: [(String, QOp)]
dsumCases =
  [ ("H ⊕ X",           DirectSum H X)
  , ("(R Z 1/3) ⊕ H",   DirectSum (R Z (1/3)) H)
  ]

adjCases :: [(String, QOp)]
adjCases =
  [ ("Adjoint H",                 Adjoint H)
  , ("Adjoint SX",                Adjoint SX)
  , ("Adjoint (R Z 1/3)",         Adjoint (R Z (1/3)))
  , ("Adjoint (R X 1/4)",         Adjoint (R X (1/4)))
  , ("Adjoint (R (X⊗Z) 1/4)",     Adjoint (R (Tensor X Z) (1/4)))
  , ("Adjoint (C X)",             Adjoint (C X))
  , ("Adjoint (Permute [2,0,1])", Adjoint (Permute [2,0,1]))
  ]

qftCases :: [(String, QOp)]
qftCases = [("qft " ++ show n, qft n) | n <- [1, 2, 3, 4, 5]]

cases :: String -> [(String, QOp)] -> TestTree
cases groupName cs = testGroup groupName [ caseAllBackends lbl op | (lbl, op) <- cs ]

fixedTests :: TestTree
fixedTests = testGroup "Fixed cases"
  [ primTests
  , phaseTests
  , idTests
  , rotTests
  , multiRotTests
  , permTests
  , cases "Tensors"      tensorCases
  , cases "Compositions" composeCases
  , cases "Controlled"   controlledCases
  , cases "Direct sums"  dsumCases
  , cases "Adjoints"     adjCases
  , cases "QFT circuits" qftCases
  ]

------------------------------------------------------------------------
-- Property suite
------------------------------------------------------------------------

-- | Use a modest number of tests since each evaluates 3 backends + ≤ 16 columns.
qcOpts :: TestTree -> TestTree
qcOpts = localOption (QuickCheckTests 50) . localOption (QuickCheckMaxSize 5)

propTests :: TestTree
propTests = qcOpts $ testGroup "Properties"
  [ testProperty "MPS matches MS"           prop_mps_matches_ms
  , testProperty "SV  matches MS"           prop_sv_matches_ms
  , testProperty "MPS output is unitary"    prop_unitary_mps
  , testProperty "op ∘ adj op = I (MPS)"    prop_inverse_mps
  , testProperty "adj² (op) = op (MPS)"     prop_double_adjoint_mps
  , testProperty "R ax 0 = Id (MPS)"        prop_R_zero
  , testProperty "R ax 2θ = (R ax θ)²"      prop_R_doubling
  , testProperty "R ax θ ∘ R ax (-θ) = I"   prop_R_inverse
  , testProperty "Permute π ∘ π⁻¹ = Id"    prop_permute_inverse
  ]

------------------------------------------------------------------------

main :: IO ()
main = defaultMain $ testGroup "BackendValidation" [fixedTests, propTests]
