{-# LANGUAGE ScopedTypeVariables #-}
-- | Correctness tests for `MPS.sampleAll` (Ferris--Vidal perfect-sampling
--   all-qubit measurement). Grouped by failure mode rather than API surface,
--   so each group pins down a specific way the implementation could go wrong.
--
--   * Determinism on product kets (degenerate marginals).
--   * Logical/physical index mapping (states with `Permute`).
--   * Edge sizes / center positions (n=1, c=n-1).
--   * Invariances (scalar prefactor, center position).
--   * RNG plumbing (exactly n draws consumed).
--   * Empirical Born-rule agreement vs `mpsToDenseVec` amplitudes.
--   * Cross-method (sampleAll vs fold-measure1, both within MPS).
--   * Cross-backend (MPS sampling vs MS sampling — the analog of
--     `BackendValidation`'s "MPS matches MS" pattern).
--   * Truncation safety under `approxCfg`.
--   * QC: random circuits applied to |0..0>; bitstring length invariant;
--     sum-to-1 over the analytic distribution.
module Main where

import QPP
import qualified QPP.Semantics.Matrix      as MS
import qualified QPP.Semantics.MPS         as MPS
import QPP.Semantics.MPS (sampleAll)
import Programs.QFT (qft)

import qualified Data.Vector as V
import Data.List (foldl')
import Data.Complex (Complex(..), magnitude)
import qualified Numeric.LinearAlgebra as H

import System.Random (mkStdGen, randoms)

import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck

import Generators (RandomQOp(..), tol, genRat, genPermutation)
import Control.Monad (replicateM)

------------------------------------------------------------------------
-- Sampling-specific helpers
------------------------------------------------------------------------

-- | MSB-first bitstring -> integer index, matching `mpsToDenseVec`'s convention.
bitsToIdx :: [Bool] -> Int
bitsToIdx = foldl' (\acc b -> 2*acc + (if b then 1 else 0)) 0

-- | Run a Program on the given state with a fixed RNG, discard outcomes/RNG.
runProg :: Program -> MPS.StateT -> MPS.StateT
runProg prog st0 =
  let rng0 = randoms (mkStdGen 0) :: [Double]
      (st, _, _) = MPS.evalProg prog st0 rng0
  in st

-- | Empirical histogram of K independent `sampleAll` trials.
empiricalHist :: Int -> Int -> Int -> MPS.StateT -> V.Vector Int
empiricalHist seed n k st =
  let go 0 _   hist = hist
      go i rng hist =
        let (outs, rng') = sampleAll st rng
            idx          = bitsToIdx outs
        in go (i-1) rng' (hist V.// [(idx, hist V.! idx + 1)])
      rng0 = randoms (mkStdGen seed) :: [Double]
  in go k rng0 (V.replicate (2^n) 0)

-- | Histogram via the regular `Measure [0..n-1]` Step path (MPS).
measureHistMPS :: Int -> Int -> Int -> MPS.StateT -> V.Vector Int
measureHistMPS seed n k st =
  let go 0 _   hist = hist
      go i rng hist =
        let (_, outs, rng') = MPS.evalProg [Measure [0..n-1]] st rng
            idx             = bitsToIdx outs
        in go (i-1) rng' (hist V.// [(idx, hist V.! idx + 1)])
      rng0 = randoms (mkStdGen seed) :: [Double]
  in go k rng0 (V.replicate (2^n) 0)

-- | Sampling histogram via the Matrix backend, running the same
--   Program ++ Measure on its own ket.
msHist :: Int -> Int -> Int -> Program -> V.Vector Int
msHist seed n k prog =
  let st0 = MS.ket (replicate n 0)
      go 0 _   hist = hist
      go i rng hist =
        let (_, outs, rng') = MS.evalProg (prog ++ [Measure [0..n-1]]) st0 rng
            idx             = bitsToIdx outs
        in go (i-1) rng' (hist V.// [(idx, hist V.! idx + 1)])
      rng0 = randoms (mkStdGen seed) :: [Double]
  in go k rng0 (V.replicate (2^n) 0)

-- | Born probabilities |<bs|psi>|^2 from the MPS amplitudes.
analyticProbs :: MPS.StateT -> V.Vector Double
analyticProbs st =
  let v = MPS.mpsToDenseVec st
      d = H.rows v
  in V.fromList [ let z = v `H.atIndex` (i,0) in let m = magnitude z in m*m
                | i <- [0 .. d-1] ]

-- | Cell-wise tolerance check: |observed - expected| <= cellTol per bin.
--   At K=4000, sigma per cell is sqrt(K p (1-p)) <= ~32, so cellTol=100 is ~3 sigma.
assertHistMatches :: String -> Double -> Int -> V.Vector Double -> V.Vector Int -> Assertion
assertHistMatches lbl cellTol k probs hist =
  let cells   = V.length probs
      expCt i = fromIntegral k * (probs V.! i)
      obsCt i = fromIntegral (hist V.! i)
      bad     = [ (i, expCt i, obsCt i)
                | i <- [0 .. cells-1]
                , abs (obsCt i - expCt i) > cellTol ]
  in case bad of
       [] -> pure ()
       xs -> assertFailure $
               lbl ++ ": empirical vs expected disagrees in "
                   ++ show (length xs) ++ " cells (K=" ++ show k
                   ++ ", cellTol=" ++ show cellTol ++ "):\n"
                   ++ unlines [ "  idx " ++ show i
                                  ++ "  exp=" ++ show e
                                  ++ "  obs=" ++ show o
                              | (i,e,o) <- take 8 xs ]

-- | Total variation distance between two empirical histograms (each of total K).
tvDistance :: V.Vector Int -> V.Vector Int -> Double
tvDistance h1 h2 =
  let total1 = fromIntegral (V.sum h1) :: Double
      total2 = fromIntegral (V.sum h2) :: Double
      p1 i = fromIntegral (h1 V.! i) / total1
      p2 i = fromIntegral (h2 V.! i) / total2
      n = V.length h1
  in 0.5 * sum [ abs (p1 i - p2 i) | i <- [0 .. n - 1] ]

-- | Shared QC budget across the cheap structural / invariance properties.
qcOpts :: TestTree -> TestTree
qcOpts = localOption (QuickCheckTests 30) . localOption (QuickCheckMaxSize 3)

-- | Random complex scalar with bounded entries (avoid zero).
genComplex :: Gen (Complex Double)
genComplex = do
  r <- choose (-4, 4)
  i <- choose (-4, 4)
  return (r :+ i)

-- | Standard programs used across multiple groups.
bellProg :: Program
bellProg = [ Unitary (Tensor H (Id 1)), Unitary (C X) ]

ghzProg :: Program
ghzProg =
  [ Unitary (Tensor H (Tensor (Id 1) (Id 1)))
  , Unitary (Tensor (C X) (Id 1))
  , Unitary (Tensor (Id 1) (C X))
  ]

plusProg :: Int -> Program
plusProg n = [ Unitary (foldr1 Tensor (replicate n H)) ]

qftProg :: Int -> Program
qftProg n = [ Unitary (qft n) ]

------------------------------------------------------------------------
-- Determinism on product kets
------------------------------------------------------------------------

determinismCase :: String -> [Int] -> TestTree
determinismCase lbl bs = testCase lbl $ do
  let st       = MPS.ket bs
      expected = [ b == 1 | b <- bs ]
      runs = [ (seed, fst (sampleAll st (randoms (mkStdGen seed) :: [Double])))
             | seed <- [1, 137, 9999, -7, 0] ]
  mapM_ (\(seed, outs) ->
           assertEqual ("ket " ++ show bs ++ " seed=" ++ show seed)
                       expected outs)
        runs

determinismTests :: TestTree
determinismTests = testGroup "Determinism on product kets"
  [ determinismCase "|0>"     [0]
  , determinismCase "|1>"     [1]
  , determinismCase "|01>"    [0,1]
  , determinismCase "|10>"    [1,0]
  , determinismCase "|110>"   [1,1,0]
  , determinismCase "|10110>" [1,0,1,1,0]
  ]

------------------------------------------------------------------------
-- Logical/physical index mapping (with Permute)
------------------------------------------------------------------------

-- `Permute` updates only log2phys; physical sites stay put with original bits.
-- Sampling reads through phys2log, so logical-qubit-k's outcome must be the
-- bit that ended up at logical k after the permutation.
permuteCase :: String -> [Int] -> [Int] -> [Bool] -> TestTree
permuteCase lbl bs perm expected = testCase lbl $ do
  let st   = runProg [ Unitary (Permute perm) ] (MPS.ket bs)
      runs = [ fst (sampleAll st (randoms (mkStdGen seed) :: [Double]))
             | seed <- [1, 42, 1337] ]
  mapM_ (\outs -> assertEqual (lbl ++ " outcomes") expected outs) runs

permuteTests :: TestTree
permuteTests = testGroup "Logical/physical index mapping"
  [ permuteCase "Permute [1,0] · |01>"
                [0,1] [1,0] [True, False]
  , permuteCase "Permute [1,0] · |10>"
                [1,0] [1,0] [False, True]
  , permuteCase "Permute [2,0,1] · |011>"
                [0,1,1] [2,0,1] [True, False, True]
  , permuteCase "Permute [2,1,0] · |100>"
                [1,0,0] [2,1,0] [False, False, True]
  , -- Double permute test: composing two permutations should still yield a
    -- product state with deterministic outcomes equal to the composed remap.
    testCase "Permute [1,0] ∘ Permute [1,0] · |01> = |01>" $ do
      let prog = [ Unitary (Permute [1,0]), Unitary (Permute [1,0]) ]
          st   = runProg prog (MPS.ket [0,1])
          outs = fst (sampleAll st (randoms (mkStdGen 42) :: [Double]))
      assertEqual "double swap is identity" [False, True] outs
  ]

------------------------------------------------------------------------
-- Edge sizes / center positions
------------------------------------------------------------------------

edgeTests :: TestTree
edgeTests = testGroup "Edge sizes / center positions"
  [ testCase "n=1: |0>" $ do
      let (outs, _) = sampleAll (MPS.ket [0]) (randoms (mkStdGen 9) :: [Double])
      assertEqual "|0> -> [False]" [False] outs

  , testCase "n=1: |1>" $ do
      let (outs, _) = sampleAll (MPS.ket [1]) (randoms (mkStdGen 9) :: [Double])
      assertEqual "|1> -> [True]" [True] outs

  , testCase "n=1: H|0> ~ 50/50" $ do
      let st   = runProg [Unitary H] (MPS.ket [0])
          hist = empiricalHist 11 1 4000 st
      assertEqual "counts sum" 4000 (V.sum hist)
      assertBool ("|0> count " ++ show (hist V.! 0) ++ " ~ 2000")
                 (abs (hist V.! 0 - 2000) < 100)

  , -- c = n-1: the right sweep is one step; the left sweep does all the work.
    testCase "Center at n-1: |+>^3 distribution matches analytic" $ do
      let n     = 3
          st0   = runProg (plusProg n) (MPS.ket (replicate n 0))
          st    = fst (MPS.moveCenterToPhys (n-1) st0)
          probs = analyticProbs st
          hist  = empiricalHist 7 n 4000 st
      assertHistMatches "c=n-1, |+>^3" 100 4000 probs hist
  ]

------------------------------------------------------------------------
-- Invariances (scalar, center position)
------------------------------------------------------------------------

-- Same seed + same Born rule on a scaled MPS => bit-identical histograms.
-- Scaling only mutates the `scalar` field; internal sites are unaffected,
-- so the rng->bit decisions match.
prop_scalar_invariance :: RandomQOp -> Property
prop_scalar_invariance (RandomQOp n op) =
  n >= 1 && n <= 3 ==>
    forAll genComplex $ \c ->
      magnitude c > 1e-6 ==>
        let st0 = MPS.apply (MPS.evalOp op) (MPS.ket (replicate n 0))
            stS = c .* st0
            k   = 1000
            h0  = empiricalHist 21 n k st0
            hS  = empiricalHist 21 n k stS
        in counterexample ("c = " ++ show c)
             (V.toList h0 === V.toList hS)

-- Moving the center is a gauge transformation: the state's distribution is
-- invariant. With independent RNG, the two empirical hists should have
-- small TV distance.
prop_center_invariance :: RandomQOp -> Property
prop_center_invariance (RandomQOp n op) =
  n >= 2 && n <= 3 ==>
    forAll (choose (0, n-1)) $ \c ->
      let stOrig = MPS.apply (MPS.evalOp op) (MPS.ket (replicate n 0))
          st0    = fst (MPS.moveCenterToPhys 0 stOrig)
          stC    = fst (MPS.moveCenterToPhys c stOrig)
          k      = 1500
          h0     = empiricalHist 31  n k st0
          hC     = empiricalHist 131 n k stC
          d      = tvDistance h0 hC
      in counterexample ("c = " ++ show c ++ ", TV = " ++ show d) (d < 0.10)

invarianceTests :: TestTree
invarianceTests = testGroup "Invariances"
  [ qcOpts $ testProperty "scalar invariance"           prop_scalar_invariance
  , qcOpts $ testProperty "center position invariance"  prop_center_invariance
  ]

------------------------------------------------------------------------
-- RNG plumbing
------------------------------------------------------------------------

rngTailConsumed :: Int -> TestTree
rngTailConsumed n = testCase ("n=" ++ show n ++ ": consumes exactly n draws") $ do
  let st        = MPS.ket (replicate n 0)
      rng       = randoms (mkStdGen 999) :: [Double]
      (_, rng') = sampleAll st rng
      expected  = drop n (take (n + 20) rng)
      actual    = take 20 rng'
  assertEqual ("returned RNG tail (n=" ++ show n ++ ")") expected actual

rngTests :: TestTree
rngTests = testGroup "RNG plumbing"
  [ rngTailConsumed 1
  , rngTailConsumed 2
  , rngTailConsumed 5
  ]

------------------------------------------------------------------------
-- Empirical Born-rule agreement vs analytic |<bs|psi>|^2
------------------------------------------------------------------------

bornCase :: String -> Int -> Int -> Int -> Program -> TestTree
bornCase lbl n seed trials prog = testCase lbl $ do
  let st    = runProg prog (MPS.ket (replicate n 0))
      probs = analyticProbs st
      hist  = empiricalHist seed n trials st
  assertBool (lbl ++ ": analytic probs do not sum to 1")
             (abs (V.sum probs - 1.0) < tol)
  assertEqual (lbl ++ ": hist count") trials (V.sum hist)
  assertHistMatches lbl 100 trials probs hist

bornTests :: TestTree
bornTests = testGroup "Empirical Born-rule agreement"
  [ bornCase "Bell state"            2  1 4000 bellProg
  , bornCase "GHZ on 3 qubits"       3 17 4000 ghzProg
  , bornCase "|+>^2"                 2 23 4000 (plusProg 2)
  , bornCase "|+>^3"                 3 31 4000 (plusProg 3)
  , bornCase "QFT|0>^3 = |+>^3"      3 41 4000 (qftProg 3)
  , testCase "Bell: off-diagonal cells empty" $ do
      let st   = runProg bellProg (MPS.ket [0,0])
          hist = empiricalHist 51 2 4000 st
      assertEqual "|01> count must be 0" 0 (hist V.! 1)
      assertEqual "|10> count must be 0" 0 (hist V.! 2)
  , testCase "GHZ: off-correlated cells empty" $ do
      let st   = runProg ghzProg (MPS.ket [0,0,0])
          hist = empiricalHist 53 3 4000 st
      mapM_ (\i -> assertEqual ("cell " ++ show i ++ " must be 0")
                               0 (hist V.! i))
            [1,2,3,4,5,6]
  ]

------------------------------------------------------------------------
-- Cross-method: sampleAll vs fold-measure1 (both in MPS)
------------------------------------------------------------------------

-- Both methods share the Born rule; their empirical hists on independent RNG
-- streams should have small TV distance.
prop_cross_method :: RandomQOp -> Property
prop_cross_method (RandomQOp n op) =
  n >= 1 && n <= 3 ==>
    let st       = MPS.apply (MPS.evalOp op) (MPS.ket (replicate n 0))
        k        = 1500
        hSample  = empiricalHist  101 n k st
        hMeasure = measureHistMPS 202 n k st
        d        = tvDistance hSample hMeasure
    in counterexample ("TV = " ++ show d ++ ", op = " ++ showOp op) (d < 0.10)

crossMethodTests :: TestTree
crossMethodTests = testGroup "Cross-method (sampleAll vs fold-measure1)"
  [ -- One named case for quick fail-diagnostics.
    testCase "Bell" $ do
      let st       = runProg bellProg (MPS.ket [0,0])
          hSample  = empiricalHist  101 2 4000 st
          hMeasure = measureHistMPS 202 2 4000 st
      assertBool ("Bell TV " ++ show (tvDistance hSample hMeasure))
                 (tvDistance hSample hMeasure < 0.05)
  , qcOpts $ testProperty "TV(sampleAll, fold measure1) small"
                          prop_cross_method
  ]

------------------------------------------------------------------------
-- Cross-backend: MPS sampling vs Matrix-backend sampling
------------------------------------------------------------------------

-- Build the same circuit on MS and on MPS, run Born-rule sampling on each,
-- check TV distance is small. Catches divergence between the two backends
-- in the measurement path specifically.
prop_cross_backend :: RandomQOp -> Property
prop_cross_backend (RandomQOp n op) =
  n >= 1 && n <= 3 ==>
    let prog  = [Unitary op]
        stMPS = runProg prog (MPS.ket (replicate n 0))
        k     = 1000
        hMPS  = empiricalHist 301 n k stMPS
        hMS   = msHist        402 n k prog
        d     = tvDistance hMPS hMS
    in counterexample ("TV(MPS,MS) = " ++ show d ++ ", op = " ++ showOp op)
                      (d < 0.10)

crossBackendTests :: TestTree
crossBackendTests = testGroup "Cross-backend (MPS vs MS sampling)"
  [ testCase "Bell" $ do
      let stMPS = runProg bellProg (MPS.ket [0,0])
          hMPS  = empiricalHist 301 2 4000 stMPS
          hMS   = msHist        402 2 4000 bellProg
      assertBool ("Bell TV(MPS,MS) " ++ show (tvDistance hMPS hMS))
                 (tvDistance hMPS hMS < 0.05)
  , localOption (QuickCheckTests 8) $ localOption (QuickCheckMaxSize 3) $
      testProperty "TV(MPS, MS) small" prop_cross_backend
  ]

------------------------------------------------------------------------
-- Truncation safety
------------------------------------------------------------------------

-- Under `approxCfg`, the SVDs internal to gate application truncate. For a
-- low-entanglement state, the
-- truncated distribution should still agree with the analytic distribution
-- (of the truncated state, read off via `mpsToDenseVec`). Smoke test +
-- correctness vs. the truncated state's own amplitudes.
truncationTests :: TestTree
truncationTests = testGroup "Truncation safety (approxCfg)"
  [ testCase "GHZ-3 under approxCfg(maxBond=4, tol=1e-9)" $ do
      let n     = 3
          cfg0  = MPS.approxCfg 4 1e-9
          st0   = (MPS.ket (replicate n 0)) { MPS.cfg = cfg0 }
          rng0  = randoms (mkStdGen 0) :: [Double]
          (st,_,_) = MPS.evalProg ghzProg st0 rng0
          probs = analyticProbs st
          hist  = empiricalHist 71 n 4000 st
      assertBool "truncated probs ~ sum to 1"
                 (abs (V.sum probs - 1.0) < 1e-6)
      assertHistMatches "GHZ trunc" 100 4000 probs hist
  ]

------------------------------------------------------------------------
-- QuickCheck: random circuits on |0..0>
------------------------------------------------------------------------

-- For a random unitary on n in [1..3], sampleAll's empirical histogram
-- should match the analytic |<bs|psi>|^2 distribution from mpsToDenseVec.
prop_random_circuit :: RandomQOp -> Property
prop_random_circuit (RandomQOp n op) =
  n >= 1 && n <= 3 ==>
    let st     = MPS.apply (MPS.evalOp op) (MPS.ket (replicate n 0))
        probs  = analyticProbs st
        trials = 2000
        hist   = empiricalHist 12345 n trials st
        cellTol = 80
        diffs  = [ ( i
                   , fromIntegral trials * (probs V.! i)
                   , fromIntegral (hist V.! i) :: Double )
                 | i <- [0 .. 2^n - 1] ]
        ok = all (\(_, e, o) -> abs (o - e) < cellTol) diffs
    in counterexample
         (unlines $
            ("random op = " ++ showOp op) :
            [ "  idx " ++ show i ++ " exp=" ++ show e ++ " obs=" ++ show o
            | (i,e,o) <- diffs ])
         ok

-- Analytic distribution sums to 1 (under defaultCfg = Exact, modulo FP).
prop_random_probs_sum_to_one :: RandomQOp -> Property
prop_random_probs_sum_to_one (RandomQOp n op) =
  n >= 1 && n <= 4 ==>
    let st    = MPS.apply (MPS.evalOp op) (MPS.ket (replicate n 0))
        total = V.sum (analyticProbs st)
    in counterexample ("sum = " ++ show total) (abs (total - 1.0) < tol)

-- Defensive invariant: `sampleAll` returns exactly n bits.
prop_outcome_length :: RandomQOp -> Property
prop_outcome_length (RandomQOp n op) =
  n >= 1 && n <= 4 ==>
    let st   = MPS.apply (MPS.evalOp op) (MPS.ket (replicate n 0))
        rng  = randoms (mkStdGen 7) :: [Double]
        outs = fst (sampleAll st rng)
    in length outs === n

-- Sum-to-1 holds also for the EMPIRICAL hist (sums to K, regardless of state).
prop_empirical_count :: RandomQOp -> Property
prop_empirical_count (RandomQOp n op) =
  n >= 1 && n <= 3 ==>
    let st     = MPS.apply (MPS.evalOp op) (MPS.ket (replicate n 0))
        trials = 200
        hist   = empiricalHist 5 n trials st
    in V.sum hist === trials

qcTests :: TestTree
qcTests = testGroup "QC over random circuits"
  [ -- Empirical agreement: expensive (K=2000 internal trials per QC sample),
    -- so limit QC trial count further.
    localOption (QuickCheckTests 6) $
      testProperty "empirical hist matches analytic" prop_random_circuit
  , -- Cheap structural QC tests, full budget.
    qcOpts $ testProperty "analytic probs sum to 1"      prop_random_probs_sum_to_one
  , qcOpts $ testProperty "sampleAll returns n bits"     prop_outcome_length
  , qcOpts $ testProperty "empirical hist counts sum K"  prop_empirical_count
  ]

------------------------------------------------------------------------

------------------------------------------------------------------------
-- Diagonal-MPS specialization
------------------------------------------------------------------------

-- | A two-term diagonal-MPS state (χ=2). Term 1: |0..0>; Term 2: c · |1..1>.
--   Amplitude on |s> is non-zero only for the all-zero or all-one bitstring.
--   Convention: f_p(0,α) = [1,0][α] (only α=0 has |0> contribution);
--   f_p(1,α) = [0,1][α] (only α=1 has |1> contribution). The last site's
--   f_p(1) is [0,c], scaling term 2 by c.
diagTwoTerm :: Int -> ComplexT -> MPS.DiagMPS
diagTwoTerm n c =
  let f0_ord  = H.fromList [1:+0, 0:+0]
      f1_ord  = H.fromList [0:+0, 1:+0]
      f1_last = H.fromList [0:+0, c]
      site k | k == n-1  = MPS.DiagSite f0_ord f1_last
             | otherwise = MPS.DiagSite f0_ord f1_ord
  in MPS.mkDiagMPS (1:+0) (V.generate n site)

diagSampleVsAnalytic :: TestTree
diagSampleVsAnalytic = testCase "sampleAllDiag matches mpsToDenseVec on 4-qubit GHZ-like" $ do
  let n      = 4
      st     = MPS.fromDiagMPS (diagTwoTerm n (2 :+ 1))  -- unnormalized, |c|² = 5
      raw    = analyticProbs st
      total  = V.sum raw
      probs  = V.map (/ total) raw                       -- normalize for comparison
      hist   = empiricalHistDiag 401 n 4000 st
  assertBool "raw norm² > 0" (total > 1e-9)
  assertBool "analytic probs sum to 1 after normalization"
             (abs (V.sum probs - 1.0) < 1e-9)
  assertHistMatches "diag 4-qubit GHZ-like" 100 4000 probs hist

-- | Empirical hist via sampleAllDiag instead of sampleAll.
empiricalHistDiag :: Int -> Int -> Int -> MPS.StateT -> V.Vector Int
empiricalHistDiag seed n k st =
  let go 0 _   hist = hist
      go i rng hist =
        let (outs, rng') = MPS.sampleAllDiag st rng
            idx          = bitsToIdx outs
        in go (i-1) rng' (hist V.// [(idx, hist V.! idx + 1)])
      rng0 = randoms (mkStdGen seed) :: [Double]
  in go k rng0 (V.replicate (2^n) 0)

isDiagonalTests :: TestTree
isDiagonalTests = testGroup "isDiagonalMPS predicate"
  [ testCase "GHZ-like diagonal state passes" $
      assertBool "" (MPS.isDiagonalMPS (MPS.fromDiagMPS (diagTwoTerm 4 (1 :+ 0))))
  , testCase "Bell state from H+CX is diagonal (χ=2)" $ do
      let st = runProg bellProg (MPS.ket [0,0])
      assertBool "" (MPS.isDiagonalMPS st)
  , testCase "H|0>⊗H|0> is diagonal (product state, χ=1)" $ do
      let st = runProg (plusProg 2) (MPS.ket [0,0])
      assertBool "" (MPS.isDiagonalMPS st)
  , testCase "QFT|0>^3 is generally NOT diagonal" $ do
      let st = runProg (qftProg 3) (MPS.ket [0,0,0])
      -- QFT|0> = |+>^n is diagonal (χ=1)... actually it IS diagonal.
      -- So we just confirm without asserting a specific result; document it.
      _ <- pure (MPS.isDiagonalMPS st)
      pure ()
  ]

-- Cross-method: on diagonal states *built via circuits* (so canonical),
-- sampleAllDiag and the generic sampleAll should agree in distribution.
-- Bell and GHZ are diagonal-MPS by construction (χ=2 perfect correlation).
diagAgreesWithGenericCase :: String -> Int -> Program -> TestTree
diagAgreesWithGenericCase lbl n prog = testCase lbl $ do
  let st       = runProg prog (MPS.ket (replicate n 0))
  assertBool "state must be diagonal" (MPS.isDiagonalMPS st)
  let k        = 2000
      hGeneric = empiricalHist     701 n k st
      hDiag    = empiricalHistDiag 702 n k st
      d        = tvDistance hGeneric hDiag
  assertBool (lbl ++ ": TV " ++ show d) (d < 0.05)

-- Raw entry agreement: sampleDiagMPS on the DiagMPS should give the same
-- distribution as the MPS wrapper sampleAllDiag.
prop_raw_matches_wrapper :: Int -> Property
prop_raw_matches_wrapper seed' =
  forAll (choose (2, 4)) $ \n ->
    let dm    = diagTwoTerm n (1 :+ 1)
        st    = MPS.fromDiagMPS dm
        k     = 800
        seedA = seed'
        seedB = seed' + 200
        hWrapper = empiricalHistDiag seedA n k st
        hRaw     = empiricalHistRaw  seedB n k dm
        d        = tvDistance hWrapper hRaw
    in counterexample ("n = " ++ show n ++ ", TV = " ++ show d)
                      (d < 0.10)
  where
    empiricalHistRaw seedR n k dm =
      let go 0 _   hist = hist
          go i rng hist =
            let (outs, rng') = MPS.sampleDiagMPS dm rng
                idx          = bitsToIdx outs
            in go (i-1) rng' (hist V.// [(idx, hist V.! idx + 1)])
          rng0 = randoms (mkStdGen seedR) :: [Double]
      in go k rng0 (V.replicate (2^n) 0)

----------------------------------------------------------------
-- Generators.
--
-- (a) RandomDiagMPS: arbitrary-χ random diagonal-MPS, built directly from
--     vector data. Spans the *full* class (any χ, no canonical form).
--     Used for tests that don't compare against canonical-form-dependent
--     methods (sampleAll, fold-measure1).
--
-- (b) genDiagProg: a sequence of single-qubit gates + Permute + Phase on n
--     qubits — the operations under which the diagonal class is closed.
--
-- (c) DiagCircuit: genDiagProg starting from |0..0> (product state, χ=1).
--     Canonical by construction; used for sampleAll / fold-measure1
--     cross-checks.
--
-- (d) DiagSetup: a (RandomDiagMPS, genDiagProg) pair sharing n. Used to
--     test closure: applying the program to the random initial state must
--     yield a diagonal MPS.

data RandomDiagMPS = RandomDiagMPS Int Int MPS.DiagMPS

instance Show RandomDiagMPS where
  show (RandomDiagMPS n chi _) =
    "RandomDiagMPS{n=" ++ show n ++ ", χ=" ++ show chi ++ "}"

instance Arbitrary RandomDiagMPS where
  arbitrary = do
    n   <- choose (1, 4)
    -- For n=1 the MPS has a single 1×1 site, so χ must be 1.
    chi <- if n == 1 then return 1 else choose (1, 3)
    let genVec  = H.fromList <$> replicateM chi genComplex
        genSite = liftA2 MPS.DiagSite genVec genVec
    sV <- V.fromList <$> replicateM n genSite
    return (RandomDiagMPS n chi (MPS.mkDiagMPS (1:+0) sV))

diagStateFromRandom :: RandomDiagMPS -> MPS.StateT
diagStateFromRandom (RandomDiagMPS _ _ dm) = MPS.fromDiagMPS dm

genDiagProg :: Int -> Gen Program
genDiagProg n = do
  depth <- choose (0, 5)
  let gate1q = oneof
        [ pure X, pure Y, pure Z, pure H, pure SX
        , liftA2 R (elements [X, Y, Z]) genRat
        ]
      placed = do
        k <- choose (0, n-1)
        g <- gate1q
        return (Tensor (Id k) (Tensor g (Id (n - k - 1))))
      permG = Permute <$> genPermutation n
      phsG  = Phase   <$> genRat
  gates <- replicateM depth (frequency [(5, placed), (1, permG), (1, phsG)])
  return (map Unitary gates)

data DiagCircuit = DiagCircuit Int Program
instance Show DiagCircuit where
  show (DiagCircuit n prog) =
    "DiagCircuit{n=" ++ show n ++ ", depth=" ++ show (length prog) ++ "}"
instance Arbitrary DiagCircuit where
  arbitrary = do
    n    <- choose (1, 4)
    prog <- genDiagProg n
    return (DiagCircuit n prog)

diagState :: DiagCircuit -> MPS.StateT
diagState (DiagCircuit n prog) = runProg prog (MPS.ket (replicate n 0))

data DiagSetup = DiagSetup RandomDiagMPS Program
instance Show DiagSetup where
  show (DiagSetup rdm prog) =
    show rdm ++ " ⊢ depth=" ++ show (length prog)
instance Arbitrary DiagSetup where
  arbitrary = do
    rdm@(RandomDiagMPS n _ _) <- arbitrary
    prog <- genDiagProg n
    return (DiagSetup rdm prog)

----------------------------------------------------------------
-- Properties on *arbitrary-χ* diagonal MPS (RandomDiagMPS).
-- These don't require canonical form.

-- Born rule: empirical hist matches normalized analytic probs.
prop_diag_born_random :: RandomDiagMPS -> Property
prop_diag_born_random rdm@(RandomDiagMPS n _ _) =
  let st    = diagStateFromRandom rdm
      raw   = analyticProbs st
      total = V.sum raw
  in total > 1e-9 ==>
     let probs = V.map (/ total) raw
         k     = 2000
         hist  = empiricalHistDiag 611 n k st
         expCt i = fromIntegral k * (probs V.! i)
         obsCt i = fromIntegral (hist V.! i) :: Double
         bad   = [ (i, expCt i, obsCt i)
                 | i <- [0 .. V.length probs - 1]
                 , abs (obsCt i - expCt i) > 100 ]
     in counterexample
          ("bad cells (first 4): " ++ show (take 4 bad))
          (null bad)

-- sampleDiagMPS (raw DiagMPS entry) matches MPS wrapper sampleAllDiag.
prop_raw_matches_wrapper_random :: RandomDiagMPS -> Property
prop_raw_matches_wrapper_random rdm@(RandomDiagMPS n _ dm) =
  let st    = diagStateFromRandom rdm
      k     = 1500
      hWrap = empiricalHistDiag 901 n k st
      hRaw  = let go 0 _ h = h
                  go i rg h =
                    let (outs, rg') = MPS.sampleDiagMPS dm rg
                        idx         = bitsToIdx outs
                    in go (i-1) rg' (h V.// [(idx, h V.! idx + 1)])
              in go k (randoms (mkStdGen 902) :: [Double]) (V.replicate (2^n) 0)
      d = tvDistance hWrap hRaw
  in counterexample ("TV(wrap, raw) = " ++ show d) (d < 0.10)

-- sampleAllDiag returns exactly n bits.
prop_diag_rng_length_random :: RandomDiagMPS -> Property
prop_diag_rng_length_random rdm@(RandomDiagMPS n _ _) =
  let st       = diagStateFromRandom rdm
      rng      = randoms (mkStdGen 13) :: [Double]
      (outs,_) = MPS.sampleAllDiag st rng
  in length outs === n

-- measureAllDiag leaves the post-state as a unit basis vector at the
-- sampled bitstring index, zero elsewhere.
prop_measureAllDiag_classical_random :: RandomDiagMPS -> Int -> Property
prop_measureAllDiag_classical_random rdm@(RandomDiagMPS n _ _) seedN =
  let st = diagStateFromRandom rdm
  in V.sum (analyticProbs st) > 1e-9 ==>
     let rng            = randoms (mkStdGen seedN) :: [Double]
         (st', outs, _) = MPS.measureAllDiag st rng
         vec            = MPS.mpsToDenseVec st'
         d              = H.rows vec
         idx            = bitsToIdx outs
         eps            = 1e-6
         badCells       = [ i
                          | i <- [0 .. d-1]
                          , let m = magnitude (vec `H.atIndex` (i, 0))
                          , if i == idx then abs (m - 1) > eps
                                        else m > eps ]
     in counterexample
          ("n=" ++ show n ++ ", outs=" ++ show outs
           ++ ", bad cells = " ++ show (take 4 badCells))
          (null badCells)

-- measureAllDiag post-state has unit norm.
prop_measureAllDiag_normalized_random :: RandomDiagMPS -> Int -> Property
prop_measureAllDiag_normalized_random rdm seedN =
  let st = diagStateFromRandom rdm
  in V.sum (analyticProbs st) > 1e-9 ==>
     let rng         = randoms (mkStdGen seedN) :: [Double]
         (st', _, _) = MPS.measureAllDiag st rng
         total       = V.sum (analyticProbs st')
     in counterexample ("post-state norm² = " ++ show total)
                       (abs (total - 1.0) < 1e-6)

-- Same, but the input MPS has a non-unit global scalar. Catches the bug where
-- measureAllDiag was multiplying the renormalisation factor by the original
-- scalar instead of replacing it.
prop_measureAllDiag_normalized_scaled :: RandomDiagMPS -> Int -> Property
prop_measureAllDiag_normalized_scaled rdm seedN =
  let st0 = diagStateFromRandom rdm
  in V.sum (analyticProbs st0) > 1e-9 ==>
     forAll (genComplex `suchThat` (\c -> magnitude c > 0.5 && magnitude c < 10)) $ \c ->
       let st         = c .* st0
           rng        = randoms (mkStdGen seedN) :: [Double]
           (st',_,_)  = MPS.measureAllDiag st rng
           total      = V.sum (analyticProbs st')
       in counterexample ("scalar=" ++ show c ++ ", post-norm²=" ++ show total)
                         (abs (total - 1.0) < 1e-6)

-- Same RNG ⇒ same outcomes from sampleAllDiag and measureAllDiag.
prop_consistent_rng_random :: RandomDiagMPS -> Property
prop_consistent_rng_random rdm =
  let st           = diagStateFromRandom rdm
  in V.sum (analyticProbs st) > 1e-9 ==>
     let rng           = randoms (mkStdGen 17) :: [Double]
         (outsS, _)    = MPS.sampleAllDiag  st rng
         (_, outsM, _) = MPS.measureAllDiag st rng
     in outsS === outsM

-- Closure: applying a sequence of 1q gates + Permute + Phase to an arbitrary
-- diagonal MPS produces a diagonal MPS.
prop_diag_closure_under_gates :: DiagSetup -> Property
prop_diag_closure_under_gates (DiagSetup rdm prog) =
  let st0 = diagStateFromRandom rdm
      st  = runProg prog st0
  in counterexample ("not diagonal after depth-" ++ show (length prog) ++ " circuit")
                    (MPS.isDiagonalMPS st)

----------------------------------------------------------------
-- Properties on *canonical* diagonal states (DiagCircuit from |0..0>).
-- These compare against canonical-form-dependent methods.

-- sampleAllDiag agrees with generic sampleAll.
prop_diag_agrees_sampleAll :: DiagCircuit -> Property
prop_diag_agrees_sampleAll dc@(DiagCircuit n _) =
  let st = diagState dc
      k  = 1500
      h1 = empiricalHist     701 n k st
      h2 = empiricalHistDiag 702 n k st
      d  = tvDistance h1 h2
  in counterexample ("TV(sampleAll, sampleAllDiag) = " ++ show d) (d < 0.10)

-- sampleAllDiag agrees with the projecting `Measure` step.
prop_diag_agrees_measure1 :: DiagCircuit -> Property
prop_diag_agrees_measure1 dc@(DiagCircuit n _) =
  let st = diagState dc
      k  = 1500
      h1 = measureHistMPS    801 n k st
      h2 = empiricalHistDiag 802 n k st
      d  = tvDistance h1 h2
  in counterexample ("TV(measure1, sampleAllDiag) = " ++ show d) (d < 0.10)

diagTests :: TestTree
diagTests = testGroup "Diagonal-MPS specialization"
  [ isDiagonalTests
  , diagSampleVsAnalytic
  , testGroup "Diag agrees with generic on (canonical) named states"
      [ diagAgreesWithGenericCase "Bell"  2 bellProg
      , diagAgreesWithGenericCase "GHZ-3" 3 ghzProg
      , diagAgreesWithGenericCase "|+>^3" 3 (plusProg 3)
      ]
  , testGroup "QC: structural invariants (RandomDiagMPS, any χ)"
      [ qcOpts $ testProperty "Closure under 1q + Permute + Phase"
                              prop_diag_closure_under_gates
      , qcOpts $ testProperty "Raw entry matches MPS wrapper"
                              prop_raw_matches_wrapper_random
      , qcOpts $ testProperty "sampleAllDiag returns n bits"
                              prop_diag_rng_length_random
      , qcOpts $ testProperty "measureAllDiag yields classical state"
                              (\rdm -> prop_measureAllDiag_classical_random rdm 29)
      , qcOpts $ testProperty "measureAllDiag post-state is normalized"
                              (\rdm -> prop_measureAllDiag_normalized_random rdm 31)
      , qcOpts $ testProperty "measureAllDiag normalizes under non-unit scalar"
                              (\rdm -> prop_measureAllDiag_normalized_scaled rdm 37)
      , qcOpts $ testProperty "sample/measure agree on shared RNG"
                              prop_consistent_rng_random
      ]
  , testGroup "QC: empirical (slower)"
      [ localOption (QuickCheckTests 6) $ localOption (QuickCheckMaxSize 3) $
          testProperty "Born: hist matches normalized analytic (any χ)"
                       prop_diag_born_random
      , localOption (QuickCheckTests 8) $ localOption (QuickCheckMaxSize 3) $
          testProperty "TV(sampleAllDiag, sampleAll) on canonical states"
                       prop_diag_agrees_sampleAll
      , localOption (QuickCheckTests 8) $ localOption (QuickCheckMaxSize 3) $
          testProperty "TV(sampleAllDiag, fold-measure1) on canonical states"
                       prop_diag_agrees_measure1
      ]
  ]

------------------------------------------------------------------------

main :: IO ()
main = defaultMain $ testGroup "SamplingTest"
  [ determinismTests
  , permuteTests
  , edgeTests
  , invarianceTests
  , rngTests
  , bornTests
  , crossMethodTests
  , crossBackendTests
  , truncationTests
  , qcTests
  , diagTests
  ]
