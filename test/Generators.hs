{-# LANGUAGE ScopedTypeVariables #-}

-- | Shared QC generators, matrix-equality helpers, and tolerance for the
--   validation test suites (`BackendValidation`, `AlgebraTest`, etc.).
--   Imported via cabal's `other-modules: Generators` in each test exe.
module Generators
  ( -- * Tolerance & matrix-equality helpers
    tol
  , mDiff
  , identCMat
  , (~~)
  , assertClose
    -- * Random rationals
  , genRat
  , shrinkRat
    -- * Pauli / permutation generators
  , genPauli
  , genPermutation
    -- * QOp generators
  , genQOpAt
  , genQOpBase
    -- * Newtype wrappers with `Arbitrary` instances
  , RandomQOp(..)
  , RotGen(..)
  , Perm(..)
  , QFTSize(..)
    -- * Shrinking
  , shrinkQOp
  ) where

import QPP
import QPP.Semantics.Matrix (CMat)

import Numeric.LinearAlgebra (sumElements, cmap, magnitude, ident, complex)
import qualified Numeric.LinearAlgebra as H

import Data.Ratio ((%), numerator, denominator)

import Test.Tasty.HUnit (Assertion, assertBool)
import Test.Tasty.QuickCheck

------------------------------------------------------------------------
-- Tolerance + matrix equality
------------------------------------------------------------------------

-- | Matrix L1 diff tolerance. Bigger circuits accumulate FP noise; this is
--   the lower bound for distinguishing "same matrix" from a real disagreement.
tol :: Double
tol = 1e-7

-- | Sum of element-wise complex magnitudes of (a − b).
mDiff :: CMat -> CMat -> Double
mDiff a b = sumElements (cmap magnitude (a - b))

identCMat :: Int -> CMat
identCMat n = complex (ident n :: H.Matrix Double)

-- | QuickCheck equality predicate with a counterexample carrying the diff.
(~~) :: CMat -> CMat -> Property
a ~~ b =
  let d = mDiff a b
  in counterexample ("matrix diff = " ++ show d) (d < tol)

infix 4 ~~

assertClose :: String -> CMat -> CMat -> Assertion
assertClose msg a b =
  let d = mDiff a b
  in assertBool (msg ++ " (diff=" ++ show d ++ ")") (d < tol)

------------------------------------------------------------------------
-- Rationals / Pauli / permutations
------------------------------------------------------------------------

-- | Random rational, signed, small numerator and denominator.
genRat :: Gen Rational
genRat = do
  num <- choose ((-7), 7) :: Gen Integer
  den <- choose (1, 8)    :: Gen Integer
  return (num % den)

-- | Shrink the numerator towards 0 and the denominator towards 1. Every
--   candidate has a smaller |numerator| + denominator (or trades a minus sign
--   for a plus), so shrink chains are finite. Halving would never reach 0.
shrinkRat :: Rational -> [Rational]
shrinkRat r =
     [ num' % den | num' <- shrink num ]
  ++ [ num % den' | den' <- shrink den, den' >= 1 ]
  where num = numerator r
        den = denominator r

-- | Pauli string on exactly n qubits: I, X, Y, Z, or tensor of these.
genPauli :: Int -> Gen QOp
genPauli 0 = return (Id 0)
genPauli 1 = elements [Id 1, X, Y, Z]
genPauli n = do
  k <- choose (1, n - 1)
  liftA2 Tensor (genPauli k) (genPauli (n - k))

genPermutation :: Int -> Gen [Int]
genPermutation n = shuffle [0 .. n - 1]

------------------------------------------------------------------------
-- QOp generators
------------------------------------------------------------------------

-- | QOp on exactly n qubits with recursion depth ≤ d.
genQOpAt :: Int -> Int -> Gen QOp
genQOpAt n 0 = genQOpBase n
genQOpAt n d = frequency $
  [ (1, genQOpBase n)
  , (2, do k <- choose (0, n)
           liftA2 Tensor (genQOpAt k (d-1)) (genQOpAt (n-k) (d-1)))
  , (2, liftA2 Compose (genQOpAt n (d-1)) (genQOpAt n (d-1)))
  , (1, fmap Adjoint (genQOpAt n (d-1)))
  ] ++ (if n >= 1 then
          [ (2, fmap C (genQOpAt (n-1) (d-1)))
          , (1, liftA2 DirectSum (genQOpAt (n-1) (d-1)) (genQOpAt (n-1) (d-1)))
          ] else [])

-- | Leaf-level generator (no recursion).
genQOpBase :: Int -> Gen QOp
genQOpBase 0 = oneof [return (Id 0), Phase <$> genRat]
genQOpBase 1 = oneof
  [ return (Id 1)
  , return X, return Y, return Z, return H, return SX
  , liftA2 R (elements [X, Y, Z]) genRat
  ]
genQOpBase n = oneof
  [ return (Id n)
  , Permute <$> genPermutation n
  , liftA2 R (genPauli n) genRat
  ]

------------------------------------------------------------------------
-- Arbitrary wrappers
------------------------------------------------------------------------

-- | Random QOp on a bounded number of qubits and depth (both ≤ 4).
data RandomQOp = RandomQOp Int QOp
instance Show RandomQOp where
  show (RandomQOp n op) = "RandomQOp{n=" ++ show n ++ "} " ++ showOp op

instance Arbitrary RandomQOp where
  arbitrary = do
    n <- choose (1, 4)
    d <- choose (1, 4)
    op <- genQOpAt n d
    return (RandomQOp n op)
  shrink (RandomQOp n op) =
    [ RandomQOp n op' | op' <- shrinkQOp n op ]

-- | Shrink a QOp while keeping its qubit count equal to `n`.
--   A candidate equal to the input (e.g. replacing the `Id` factor of
--   `Id 2 ⊗ b` by `Id 2`) is dropped: QuickCheck would accept it as a
--   successful shrink and loop on it forever.
shrinkQOp :: Int -> QOp -> [QOp]
shrinkQOp n op = filter (/= op) $ case op of
  Id _      -> []
  X         -> [Id n]
  Y         -> [Id n, X]
  Z         -> [Id n, X]
  H         -> [Id n, X]
  SX        -> [Id n, X, H]
  Phase q   -> [Id n] ++ [Phase q' | q' <- shrinkRat q]
  R ax θ    ->
    [Id n] ++ [R ax θ' | θ' <- shrinkRat θ] ++ [R ax' θ | ax' <- shrinkQOp n ax]
  Permute _ -> [Id n]
  Tensor a b ->
    let na = op_qubits a; nb = op_qubits b
    in  [Id n]
     ++ [Tensor (Id na) b, Tensor a (Id nb)]
     ++ [Tensor a' b | a' <- shrinkQOp na a]
     ++ [Tensor a b' | b' <- shrinkQOp nb b]
  Compose a b ->
    [Id n, a, b]
     ++ [Compose a' b | a' <- shrinkQOp n a]
     ++ [Compose a b' | b' <- shrinkQOp n b]
  C a ->
    [Id n]
     ++ [C a' | a' <- shrinkQOp (n-1) a]
  DirectSum a b ->
    let nInner = op_qubits a   -- == op_qubits b
    in  [Id n]
     ++ [DirectSum (Id nInner) b, DirectSum a (Id nInner)]
     ++ [DirectSum a' b | a' <- shrinkQOp nInner a]
     ++ [DirectSum a b' | b' <- shrinkQOp nInner b]
  Adjoint a ->
    [Id n, a]
     ++ [Adjoint a' | a' <- shrinkQOp n a]
  _ -> [Id n]

-- | Random rotation (axis, θ, qubit count) for rotation-specific invariants.
data RotGen = RotGen Int QOp Rational
instance Show RotGen where
  show (RotGen n ax θ) =
    "RotGen{n=" ++ show n ++ ", θ=" ++ show θ ++ "} axis=" ++ showOp ax

instance Arbitrary RotGen where
  arbitrary = do
    n  <- choose (1, 3)
    ax <- genPauli n
    θ  <- genRat
    return (RotGen n ax θ)

-- | Random permutation, for permutation-specific invariants.
newtype Perm = Perm [Int]
instance Show Perm where show (Perm p) = "Perm " ++ show p
instance Arbitrary Perm where
  arbitrary = do
    n <- choose (1, 5)
    Perm <$> genPermutation n

-- | Qubit count for QFT properties. Bounded at 6 so the 2^n × 2^n reference
--   matrices stay cheap under QuickCheck.
newtype QFTSize = QFTSize Int
instance Show QFTSize where show (QFTSize n) = "QFTSize " ++ show n
instance Arbitrary QFTSize where
  arbitrary = QFTSize <$> choose (1, 6)
  shrink (QFTSize n) = [ QFTSize n' | n' <- [1 .. n - 1] ]
