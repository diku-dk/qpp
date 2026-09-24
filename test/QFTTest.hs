{-# LANGUAGE ScopedTypeVariables #-}
module Main where

import QPP
import QPP.Semantics (Convertible(..))
import qualified QPP.Semantics.Matrix      as MS
import qualified QPP.Semantics.MPS         as MPS
import qualified QPP.Semantics.Statevector as SV
import QPP.Semantics.Matrix (CMat)
import Programs.QFT (qft)

import Data.Complex (Complex(..), cis)
import qualified Numeric.LinearAlgebra as H

import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck

import Generators

------------------------------------------------------------------------
-- Reference matrices
------------------------------------------------------------------------

-- NB: `qft` uses a qubit-0-as-LSB convention — the integer a basis state
-- represents has qubit 0 contributing the 2^0 bit. The hmatrix tensor layout
-- puts qubit 0 in the 2^(n-1) position, so reference matrices index through
-- `bitRev` to convert layout index ↔ logical integer.

-- | Reverse the n-bit representation of x.
bitRev :: Int -> Int -> Int
bitRev n x = foldl (\acc i -> acc * 2 + (x `div` 2 ^ i) `mod` 2) 0 [0 .. n - 1]

-- | Analytic unitary DFT on N = 2^n points: F[j,k] = e^{2πi jk/N} / √N,
--   with j,k the logical (LSB-convention) integers. This is the *definition*
--   the circuit `qft n` is checked against.
dftMat :: Int -> CMat
dftMat n =
  let bigN = 2 ^ n :: Int
      d    = fromIntegral bigN
      ω a b = let j = bitRev n a; k = bitRev n b
              in cis (2 * pi * fromIntegral (j * k) / d) / (sqrt d :+ 0)
  in (bigN H.>< bigN) [ ω a b | a <- [0 .. bigN - 1], b <- [0 .. bigN - 1] ]

-- | Index-reversal permutation P|j⟩ = |(-j) mod N⟩. The DFT squares to this.
revMat :: Int -> CMat
revMat n =
  let bigN = 2 ^ n :: Int
  in (bigN H.>< bigN)
       [ if bitRev n aOut == ((bigN - bitRev n aIn) `mod` bigN) then 1 else 0
       | aOut <- [0 .. bigN - 1], aIn <- [0 .. bigN - 1] ]

-- | First column of the n-qubit identity-sized uniform vector: 1/√N everywhere.
uniformCol :: Int -> CMat
uniformCol n =
  let bigN = 2 ^ n :: Int
  in H.asColumn (H.fromList (replicate bigN (1 / sqrt (fromIntegral bigN) :+ 0)))

qftMat :: Int -> CMat
qftMat = MS.evalOp . qft

-- | n-fold composition of a QOp with itself.
power :: Int -> QOp -> QOp
power k op = foldr1 (∘) (replicate k op)

------------------------------------------------------------------------
-- Properties
------------------------------------------------------------------------

-- | `qft n` equals the analytic DFT matrix.
prop_qft_is_dft :: QFTSize -> Property
prop_qft_is_dft (QFTSize n) = qftMat n ~~ dftMat n

-- | `qft n` maps |0…0⟩ to the uniform superposition (column 0 is flat).
prop_qft_zero_to_uniform :: QFTSize -> Property
prop_qft_zero_to_uniform (QFTSize n) =
  let col0 = H.asColumn (head (H.toColumns (qftMat n)))
  in col0 ~~ uniformCol n

-- | `qft n` is unitary: M · M† = I.
prop_qft_unitary :: QFTSize -> Property
prop_qft_unitary (QFTSize n) =
  let m = qftMat n
  in (m H.<> H.tr m) ~~ identCMat (2 ^ n)

-- | `qft n ∘ adj (qft n) = I`.
prop_qft_inverse :: QFTSize -> Property
prop_qft_inverse (QFTSize n) =
  MS.evalOp (qft n ∘ Adjoint (qft n)) ~~ identCMat (2 ^ n)

-- | The DFT has order 4: (qft n)⁴ = I.
prop_qft_order4 :: QFTSize -> Property
prop_qft_order4 (QFTSize n) =
  MS.evalOp (power 4 (qft n)) ~~ identCMat (2 ^ n)

-- | (qft n)² is the index-reversal permutation |j⟩ ↦ |(-j) mod N⟩.
prop_qft_squared_reversal :: QFTSize -> Property
prop_qft_squared_reversal (QFTSize n) =
  MS.evalOp (power 2 (qft n)) ~~ revMat n

-- | The Matrix, MPS, and Statevector backends agree on `qft n`.
prop_qft_backends_agree :: QFTSize -> Property
prop_qft_backends_agree (QFTSize n) =
  let mref = qftMat n
      mmps = to (MPS.evalOp (qft n)) :: CMat
      msv  = to (SV.evalOp  (qft n)) :: CMat
  in (mref ~~ mmps) .&&. (mref ~~ msv)

------------------------------------------------------------------------
-- Fixed regression cases
------------------------------------------------------------------------

fixedTests :: TestTree
fixedTests = testGroup "Fixed cases"
  [ testCase "qft 1 = H" $
      assertClose "qft 1" (qftMat 1) (MS.evalOp H)
  , testCase "qft 2 = DFT_4" $
      assertClose "qft 2" (qftMat 2) (dftMat 2)
  , testCase "qft 3 = DFT_8" $
      assertClose "qft 3" (qftMat 3) (dftMat 3)
  ]

------------------------------------------------------------------------
-- Property suite
------------------------------------------------------------------------

qcOpts :: TestTree -> TestTree
qcOpts = localOption (QuickCheckTests 40)

propTests :: TestTree
propTests = qcOpts $ testGroup "Properties"
  [ testProperty "qft n = analytic DFT"        prop_qft_is_dft
  , testProperty "qft n |0⟩ = uniform"         prop_qft_zero_to_uniform
  , testProperty "qft n is unitary"            prop_qft_unitary
  , testProperty "qft n ∘ adj (qft n) = I"     prop_qft_inverse
  , testProperty "(qft n)⁴ = I"                prop_qft_order4
  , testProperty "(qft n)² = index reversal"   prop_qft_squared_reversal
  , testProperty "MS / MPS / SV agree on qft"  prop_qft_backends_agree
  ]

------------------------------------------------------------------------

main :: IO ()
main = defaultMain $ testGroup "QFTTest" [fixedTests, propTests]
