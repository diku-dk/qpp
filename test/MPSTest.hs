{-# LANGUAGE ScopedTypeVariables #-}

-- | MPS-internal invariants that the cross-backend suites can't see:
--
--   * Canonical-center validity: every site left of `center_site` must be a
--     left isometry, every site right of it a right isometry. All probability
--     computations (`measure1`, `sampleAll`) silently assume this; a stale
--     center gives *wrong probabilities* while the state vector itself stays
--     correct (so `BackendValidation` passes).
--   * Truncation that actually binds: `approxCfg` with maxBond smaller than
--     the exact bond dimension. (`SamplingTest`'s truncation case uses a
--     maxBond above the state's exact bond dim, so the truncating branch of
--     `svd_compact` was never exercised.)
--   * Profiling stats under binding truncation.
--   * Deterministic program-level oracle vs the Matrix backend: same RNG in,
--     same outcomes and same collapsed state out.
--   * `R` with an `Adjoint`/`Phase`-carrying Pauli-string axis (axisPaulis
--     corner cases the `genPauli`-based properties never generate).
module Main where

import QPP
import QPP.Semantics (Convertible(..))
import qualified QPP.Semantics.Matrix as MS
import qualified QPP.Semantics.MPS    as MPS
import QPP.Semantics.Matrix (CMat)

import qualified Data.Vector as V
import qualified Numeric.LinearAlgebra as H

import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck

import Generators

------------------------------------------------------------------------
-- Canonical-center validity
------------------------------------------------------------------------

isoTol :: Double
isoTol = 1e-8

-- | ‖A₀†A₀ + A₁†A₁ − I‖₁ — zero for a left isometry.
leftIsoDefect :: MPS.Site -> Double
leftIsoDefect (MPS.Site x0 x1) =
  let g = (H.tr x0 H.<> x0) + (H.tr x1 H.<> x1)
  in mDiff g (identCMat (H.cols x0))

-- | ‖A₀A₀† + A₁A₁† − I‖₁ — zero for a right isometry.
rightIsoDefect :: MPS.Site -> Double
rightIsoDefect (MPS.Site x0 x1) =
  let g = (x0 H.<> H.tr x0) + (x1 H.<> H.tr x1)
  in mDiff g (identCMat (H.rows x0))

-- | All isometry violations of the claimed canonical form, as (site, side,
--   defect) triples. Empty iff `center_site` is a valid orthogonality center.
centerDefects :: MPS.MPS -> [(Int, String, Double)]
centerDefects st =
  let c  = MPS.center_site st
      sV = MPS.sites st
      n  = V.length sV
      defL = [ (p, "left-iso",  d) | p <- [0 .. c-1],   let d = leftIsoDefect  (sV V.! p), d > isoTol ]
      defR = [ (p, "right-iso", d) | p <- [c+1 .. n-1], let d = rightIsoDefect (sV V.! p), d > isoTol ]
  in defL ++ defR

centerProp :: MPS.MPS -> Property
centerProp st =
  let defects = centerDefects st
  in counterexample
       ("center_site=" ++ show (MPS.center_site st) ++ " invalid; defects: " ++ show defects)
       (null defects)

assertCenterValid :: String -> MPS.MPS -> Assertion
assertCenterValid msg st =
  let defects = centerDefects st
  in assertBool (msg ++ ": center_site=" ++ show (MPS.center_site st)
                     ++ " invalid; defects: " ++ show defects)
                (null defects)

------------------------------------------------------------------------
-- Generators
------------------------------------------------------------------------

zerosKet :: Int -> MPS.StateT
zerosKet n = MPS.ket (replicate n 0)

applyOp :: QOp -> MPS.StateT -> MPS.StateT
applyOp op = MPS.apply (MPS.evalOp op)

-- | Two independent ops on the same register, applied in sequence. The first
--   op parks the center inside its own support; the second op's support may
--   be disjoint from it, which is exactly where stale-center bugs live.
data OpPair = OpPair Int QOp QOp
instance Show OpPair where
  show (OpPair n a b) =
    "OpPair{n=" ++ show n ++ "} " ++ showOp a ++ "  THEN  " ++ showOp b
instance Arbitrary OpPair where
  arbitrary = do
    n <- choose (2, 4)
    d <- choose (1, 3)
    OpPair n <$> genQOpAt n d <*> genQOpAt n d
  shrink (OpPair n a b) =
    [ OpPair n a' b | a' <- shrinkQOp n a ] ++
    [ OpPair n a b' | b' <- shrinkQOp n b ]

-- | A unitary/measure/unitary/measure-all program plus a fixed RNG prefix.
data RandProg = RandProg Int Program [Double]
instance Show RandProg where
  show (RandProg n prog rs) =
    "RandProg{n=" ++ show n ++ ", rng=" ++ show (take 4 rs) ++ "...}\n"
    ++ unlines (map show prog)
instance Arbitrary RandProg where
  arbitrary = do
    n   <- choose (2, 4)
    op1 <- genQOpAt n 3
    op2 <- genQOpAt n 3
    ks  <- sublistOf [0 .. n-1]
    -- Stay away from 0/1 so that FP noise in the computed probabilities
    -- cannot flip an outcome between backends.
    rs  <- vectorOf (4*n) (choose (0.02, 0.98))
    return (RandProg n [ Unitary op1, Measure ks, Unitary op2, Measure [0 .. n-1] ] rs)

------------------------------------------------------------------------
-- Center-validity properties
------------------------------------------------------------------------

prop_center_after_op :: RandomQOp -> Property
prop_center_after_op (RandomQOp n op) =
  centerProp (applyOp op (zerosKet n))

prop_center_after_op_pair :: OpPair -> Property
prop_center_after_op_pair (OpPair n a b) =
  centerProp (applyOp b (applyOp a (zerosKet n)))

-- | The state *vector* must be right even when the center bookkeeping is the
--   only thing under suspicion (mpsToDenseVec is center-agnostic).
prop_dense_after_op_pair :: OpPair -> Property
prop_dense_after_op_pair (OpPair n a b) =
  let dMPS = MPS.mpsToDenseVec (applyOp b (applyOp a (zerosKet n)))
      dMS  = MS.evalOp b H.<> MS.evalOp a H.<> MS.ket (replicate n 0)
  in dMPS ~~ dMS

prop_center_after_measure :: OpPair -> Property
prop_center_after_measure (OpPair n a b) =
  forAll (choose (0, n-1)) $ \k ->
  forAll (choose (0.02, 0.98)) $ \r ->
    let st          = applyOp b (applyOp a (zerosKet n))
        (st', _, _) = MPS.measure1 (st, [], [r]) k
    in centerProp st'

------------------------------------------------------------------------
-- Program-level oracle vs the Matrix backend (deterministic, same RNG)
------------------------------------------------------------------------

prop_program_oracle :: RandProg -> Property
prop_program_oracle (RandProg n prog rs) =
  let (stM, outsM, _) = MS.evalProg  prog (MS.ket (replicate n 0)) rs
      (stP, outsP, _) = MPS.evalProg prog (zerosKet n)             rs
  in counterexample ("outcomes: MS=" ++ show outsM ++ " MPS=" ++ show outsP)
       (outsM == outsP)
     .&&. (MPS.mpsToDenseVec stP ~~ stM)
     .&&. centerProp stP

------------------------------------------------------------------------
-- Binding truncation
------------------------------------------------------------------------

-- | 4-qubit circuit with exact bond dimension 4 at the middle cut: two layers
--   of overlapping two-qubit rotations with incommensurate angles.
entangler :: QOp
entangler = foldr1 Compose
  [      R (Tensor X X) (1/3) @> 2
  , 1 <@ R (Tensor Y Y) (1/4) @> 1
  , 2 <@ R (Tensor Z Z) (1/5)
  ,      R (Tensor Y X) (2/7) @> 2
  , 1 <@ R (Tensor X Z) (3/8) @> 1
  , 2 <@ R (Tensor Y Z) (1/7)
  ]

entangledKet :: MPS.EvalCfg -> MPS.StateT
entangledKet c = applyOp entangler ((zerosKet 4) { MPS.cfg = c })

truncTests :: TestTree
truncTests = testGroup "Binding truncation (approxCfg)"
  [ testCase "entangler needs bond > 2 exactly" $ do
      let st = entangledKet MPS.defaultCfg
      assertBool ("max bond = " ++ show (MPS.maxBondDimension st))
                 (MPS.maxBondDimension st > 2)

  , testCase "Exact matches dense reference" $
      assertClose "entangler" (MS.evalOp entangler H.<> MS.ket [0,0,0,0])
                              (MPS.mpsToDenseVec (entangledKet MPS.defaultCfg))

  , testCase "loose Truncate (maxBond 64) matches dense reference" $
      assertClose "entangler/loose"
        (MS.evalOp entangler H.<> MS.ket [0,0,0,0])
        (MPS.mpsToDenseVec (entangledKet (MPS.approxCfg 64 1e-12)))

  , testCase "binding Truncate (maxBond 2): runs, bond ≤ 2, norm ≤ 1" $ do
      let st   = entangledKet (MPS.approxCfg 2 1e-12)
          nrm  = norm st
      assertBool ("max bond = " ++ show (MPS.maxBondDimension st))
                 (MPS.maxBondDimension st <= 2)
      assertBool ("norm = " ++ show nrm) (nrm <= 1 + 1e-9 && nrm > 0.1)
      assertCenterValid "post-truncation" st

  , testCase "binding Truncate (maxBond 2): profile populated, error bounded" $ do
      let cfg0    = MPS.approxCfg 2 1e-12
          (st, p) = MPS.applyProfiled (MPS.evalOp entangler)
                                      ((zerosKet 4) { MPS.cfg = cfg0 })
      assertBool ("nSVDs = "      ++ show (MPS.nSVDs p))      (MPS.nSVDs p > 0)
      assertBool ("maxBondDim = " ++ show (MPS.maxBondDim p)) (MPS.maxBondDim p <= 2)
      assertBool ("discardedWeight = " ++ show (MPS.discardedWeight p))
                 (MPS.discardedWeight p > 0)
      -- Cauchy–Schwarz: (Σ√δ)² ≥ Σδ
      assertBool "normError² ≥ discardedWeight"
                 (MPS.normError p * MPS.normError p >= MPS.discardedWeight p - 1e-15)
      -- normError is a 2-norm bound on the truncation error of this
      -- (unitary-only) run, so it bounds the fidelity loss vs the exact state.
      let dExact  = MPS.mpsToDenseVec (entangledKet MPS.defaultCfg)
          dTrunc  = MPS.mpsToDenseVec st
          overlap = H.norm_2 (H.flatten (H.tr dExact H.<> dTrunc))
      assertBool ("normError = " ++ show (MPS.normError p) ++ " not usefully small")
                 (MPS.normError p < 0.9)
      assertBool ("overlap = " ++ show overlap ++ " < 1 - normError = "
                  ++ show (1 - MPS.normError p))
                 (overlap >= 1 - MPS.normError p - 1e-9)
  ]

-- | Random ops under a binding bond cap: must not crash, must respect the cap
--   both in the state and in the reported profile.
prop_trunc_bond_capped :: RandomQOp -> Property
prop_trunc_bond_capped (RandomQOp n op) =
  let (st, p) = MPS.applyProfiled (MPS.evalOp op)
                                  ((zerosKet n) { MPS.cfg = MPS.approxCfg 2 1e-12 })
  in counterexample ("max bond = " ++ show (MPS.maxBondDimension st)
                     ++ ", profile = " ++ show p)
       (MPS.maxBondDimension st <= 2 && MPS.maxBondDim p <= 2)

------------------------------------------------------------------------
-- Profiling: statistics must follow the evaluation, not one branch of it
------------------------------------------------------------------------

-- | Componentwise comparison with FP slack on the accumulated doubles.
profApprox :: MPS.Profile -> MPS.Profile -> Property
profApprox p q =
  counterexample (show p ++ " /≈ " ++ show q) $
        MPS.nSVDs p == MPS.nSVDs q
    &&  MPS.maxBondDim p == MPS.maxBondDim q
    &&  abs (MPS.discardedWeight p - MPS.discardedWeight q) < 1e-12
    &&  abs (MPS.normError p - MPS.normError q) < 1e-12

-- | The profile is a pure annotation: the profiled entry point computes the
--   very same state as the plain one.
prop_profiled_matches_plain :: RandomQOp -> Property
prop_profiled_matches_plain (RandomQOp n op) =
  let st = zerosKet n
  in fst (MPS.applyProfiled (MPS.evalOp op) st) === applyOp op st

-- | Under Exact evaluation nothing is ever discarded.
prop_exact_run_zero_error :: RandomQOp -> Property
prop_exact_run_zero_error (RandomQOp n op) =
  let (_, p) = MPS.applyProfiled (MPS.evalOp op) (zerosKet n)
  in counterexample (show p)
       (MPS.discardedWeight p == 0 && MPS.normError p == 0)

-- | Profiles accumulate along the execution path: evaluating `b ∘ a` reports
--   exactly the two segments' profiles combined.
prop_path_additive :: OpPair -> Property
prop_path_additive (OpPair n a b) =
  let st0       = zerosKet n
      (st1, p1) = MPS.applyProfiled (MPS.evalOp a) st0
      (st2, p2) = MPS.applyProfiled (MPS.evalOp b) st1
      (stC, pC) = MPS.applyProfiled (MPS.evalOp (Compose b a)) st0
  in stC === st2 .&&. profApprox pC (p1 <> p2)

-- | Evaluating `C a` includes evaluating `a` (on the projected true-branch)
--   plus the projection/merge work, so its SVD count must dominate the count
--   for `a` placed on the same register. This is the property the old
--   cfg-embedded counters violated (they reported the *fewer* SVDs of
--   whichever branch's state survived the merge).
profileTests :: TestTree
profileTests = testGroup "Profiling"
  [ testCase "C a counts the controlled branch's SVD work" $ do
      let inner = Compose (R (Tensor Y Y) (1/3)) (R (Tensor X Z) (1/4))
          nSVDsOf op = MPS.nSVDs (snd (MPS.applyProfiled (MPS.evalOp op) (zerosKet 3)))
          nC = nSVDsOf (C inner)
          nA = nSVDsOf (Tensor (Id 1) inner)
      assertBool ("nSVDs: C=" ++ show nC ++ " standalone=" ++ show nA)
                 (nC >= nA && nA > 0)
  ]

------------------------------------------------------------------------
-- axisPaulis corner cases: Phase/Adjoint inside a rotation axis
------------------------------------------------------------------------

axisCases :: TestTree
axisCases = testGroup "R with Phase/Adjoint in the axis"
  [ testCase lbl $ assertClose lbl (MS.evalOp op) (to (MPS.evalOp op) :: CMat)
  | (lbl, op) <-
      [ ("R (Phase (1/3) ⊗ X) (1/2)",
          R (Tensor (Phase (1/3)) X) (1/2))
      , ("R (Adjoint (Phase (1/3) ⊗ X)) (1/2)",
          R (Adjoint (Tensor (Phase (1/3)) X)) (1/2))
      , ("R (Adjoint (Phase (1/2) ⊗ X ⊗ Z)) (2/3)",
          R (Adjoint (Tensor (Phase (1/2)) (Tensor X Z))) (2/3))
      , ("R (Adjoint (Adjoint (Phase (1/3) ⊗ Y))) (1/2)",
          R (Adjoint (Adjoint (Tensor (Phase (1/3)) Y))) (1/2))
      ]
  ]

------------------------------------------------------------------------

qcOpts :: TestTree -> TestTree
qcOpts = localOption (QuickCheckTests 100) . localOption (QuickCheckMaxSize 5)

main :: IO ()
main = defaultMain $ testGroup "MPSTest"
  [ qcOpts $ testGroup "Canonical-center validity"
      [ testProperty "valid after a random op"        prop_center_after_op
      , testProperty "valid after two ops"            prop_center_after_op_pair
      , testProperty "dense state right after two ops" prop_dense_after_op_pair
      , testProperty "valid after measure1"           prop_center_after_measure
      ]
  , qcOpts $ testGroup "Program oracle vs the Matrix backend"
      [ testProperty "same RNG ⇒ same outcomes + same state" prop_program_oracle ]
  , truncTests
  , qcOpts $ testGroup "Truncation properties"
      [ testProperty "bond ≤ maxBond under binding cap" prop_trunc_bond_capped ]
  , profileTests
  , qcOpts $ testGroup "Profiling properties"
      [ testProperty "profiled state ≡ plain state"     prop_profiled_matches_plain
      , testProperty "Exact run has zero error"         prop_exact_run_zero_error
      , testProperty "profiles add along Compose paths" prop_path_additive
      ]
  , axisCases
  ]
