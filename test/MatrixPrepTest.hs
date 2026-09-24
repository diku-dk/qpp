{-# LANGUAGE ScopedTypeVariables #-}
module Main where

import QPP
import QPP.Util (padToPowerOf2)
import QPP.Semantics (Convertible(..))
import Programs.MatrixPreparation

import qualified QPP.Semantics.Matrix      as MS
import qualified QPP.Semantics.MPS         as MPS
import qualified QPP.Semantics.Statevector as SV
import QPP.Semantics.Matrix (CMat)

import Numeric.LinearAlgebra
  ( fromLists, ident, complex, subMatrix, sumElements, cmap, magnitude
  , rows, cols, norm_Frob
  )
import qualified Numeric.LinearAlgebra as H

import qualified Data.Vector as V
import Data.Complex (Complex(..), conjugate, realPart)

import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck

----------------------------------------------------------------------
-- Tolerance + matrix helpers
----------------------------------------------------------------------

tol :: Double
tol = 1e-7

mDiff :: CMat -> CMat -> Double
mDiff a b = sumElements (cmap magnitude (a - b))

identCMat :: Int -> CMat
identCMat n = complex (ident n :: H.Matrix Double)

assertClose :: String -> CMat -> CMat -> Assertion
assertClose msg expected actual =
  let d = mDiff expected actual
  in assertBool
       (msg ++ ": diff = " ++ show d
        ++ "\nexpected:\n" ++ show expected
        ++ "\nactual:\n" ++ show actual)
       (d < tol)

(~~) :: CMat -> CMat -> Property
a ~~ b = counterexample ("matrix diff = " ++ show (mDiff a b)) (mDiff a b < tol)
infix 4 ~~

-- | Upper-left k × k block.
upperLeftBlock :: Int -> CMat -> CMat
upperLeftBlock k m = subMatrix (0, 0) (k, k) m

-- | Next power of 2 ≥ n; nextPow2 1 = 1.
nextPow2 :: Int -> Int
nextPow2 1 = 1
nextPow2 n = 2 ^ (ceiling (logBase 2 (fromIntegral n) :: Double) :: Int)

-- | Pad a complex matrix to k × k with zeros (k ≥ rows, k ≥ cols).
padMatrix :: Int -> CMat -> CMat
padMatrix k m =
  let r = rows m; c = cols m
      extraRows = replicate (k-r) (replicate k (0:+0))
      asLists   = [ row ++ replicate (k-c) (0:+0) | row <- H.toLists m ]
  in fromLists (asLists ++ extraRows)

----------------------------------------------------------------------
-- Pure helpers (splitList, padToPowerOf2)
----------------------------------------------------------------------

cv :: [Int] -> V.Vector ComplexT
cv xs = V.fromList (map (\x -> fromIntegral x :+ 0) xs)

pureHelperTests :: TestTree
pureHelperTests = testGroup "Pure helpers"
  [ testGroup "splitList (lengths only)"
      [ testCase "[a,b]            -> [2]"        $ map V.length (splitList (cv [1,2]))                 @?= [2]
      , testCase "[a,b,c]          -> [2,1]"      $ map V.length (splitList (cv [1,2,3]))               @?= [2,1]
      , testCase "[a,b,c,d]        -> [2,2]"      $ map V.length (splitList (cv [1,2,3,4]))             @?= [2,2]
      , testCase "[a,b,c,d,e]      -> [2,2,1]"    $ map V.length (splitList (cv [1,2,3,4,5]))           @?= [2,2,1]
      , testCase "[a..h] (8 elems) -> [2,2,2,2]"  $ map V.length (splitList (cv [1..8]))                @?= [2,2,2,2]
      , testCase "concat = original (length 5)"   $ V.toList (V.concat (splitList (cv [1..5])))        @?= V.toList (cv [1..5])
      , testCase "concat = original (length 11)"  $ V.toList (V.concat (splitList (cv [1..11])))       @?= V.toList (cv [1..11])
      ]
  , testGroup "padToPowerOf2"
      [ testCase "already 2^n (n=2, len=4)" $
          V.length (padToPowerOf2 2 I (V.replicate 4 I)) @?= 4
      , testCase "pad len 3 to 4 (n=2)" $
          V.length (padToPowerOf2 2 I (V.replicate 3 I)) @?= 4
      , testCase "pad len 5 to 8 (n=3)" $
          V.length (padToPowerOf2 3 I (V.replicate 5 I)) @?= 8
      ]
  ]

----------------------------------------------------------------------
-- buildRowQOp — row preparation
----------------------------------------------------------------------

-- | Apply buildRowQOp to |0..0⟩, return the resulting state-vector
--   (= first column of the unitary).
rowPrepStateMS :: V.Vector ComplexT -> CMat
rowPrepStateMS v =
  let m   = MS.evalOp (buildRowQOp v)
      dim = rows m
  in subMatrix (0, 0) (dim, 1) m

-- | Expected state: v zero-padded to next power of 2, then normalised.
expectedRowPrep :: V.Vector ComplexT -> CMat
expectedRowPrep v =
  let n_dims = V.length v
      n_pad  = nextPow2 (max 2 n_dims)
      v_pad  = V.toList v ++ replicate (n_pad - n_dims) 0
      nrm    = sqrt (realPart (sum [conjugate x * x | x <- v_pad]))
      normC  = nrm :+ 0
  in fromLists (map (\x -> [x / normC]) v_pad)

rowPrepCase :: String -> V.Vector ComplexT -> TestTree
rowPrepCase lbl v = testCase lbl $
  assertClose ("row prep " ++ lbl) (expectedRowPrep v) (rowPrepStateMS v)

rowPrepTests :: TestTree
rowPrepTests = testGroup "buildRowQOp"
  [ rowPrepCase "[3, 4]"               (V.fromList [3 :+ 0, 4 :+ 0])
  , rowPrepCase "[1, 2, 3]"            (V.fromList [1 :+ 0, 2 :+ 0, 3 :+ 0])
  , rowPrepCase "[1, 2, 3, 4]"         (V.fromList [1 :+ 0, 2 :+ 0, 3 :+ 0, 4 :+ 0])
  , rowPrepCase "[1, 2, 3, 4, 5]"      (V.fromList [1 :+ 0, 2 :+ 0, 3 :+ 0, 4 :+ 0, 5 :+ 0])
  , rowPrepCase "[1+i, 2]"             (V.fromList [1 :+ 1, 2 :+ 0])
  , rowPrepCase "[1+2i, 3+4i, 5]"      (V.fromList [1 :+ 2, 3 :+ 4, 5 :+ 0])
  , rowPrepCase "[1+i, 2+2i, 3+3i, 4]" (V.fromList [1:+1, 2:+2, 3:+3, 4:+0])
  ]

----------------------------------------------------------------------
-- matrixPrep — block encoding
----------------------------------------------------------------------

-- | Extract the upper-left block of size matching the original (padded)
--   input from evalOp (matrixPrep M).
blockEncodeMS :: CMat -> CMat
blockEncodeMS m =
  let u     = MS.evalOp (matrixPrep m)
      n_pad = nextPow2 (max (rows m) (cols m))
  in upperLeftBlock n_pad u

-- | Expected block: M (padded to n_pad × n_pad with zeros) / ‖M‖_F.
expectedBlock :: CMat -> CMat
expectedBlock m =
  let n_pad = nextPow2 (max (rows m) (cols m))
      m_pad = padMatrix n_pad m
      nrm   = norm_Frob m :+ 0
  in H.scale (1 / nrm) m_pad

matrixPrepCase :: String -> CMat -> TestTree
matrixPrepCase lbl m = testCase lbl $
  assertClose ("block encoding " ++ lbl) (expectedBlock m) (blockEncodeMS m)

matrixPrepUnitary :: String -> CMat -> TestTree
matrixPrepUnitary lbl m = testCase ("unitary " ++ lbl) $
  let u = MS.evalOp (matrixPrep m)
      n = rows u
  in assertClose ("unitarity " ++ lbl) (identCMat n) (u H.<> H.tr u)

matrixPrepTests :: TestTree
matrixPrepTests = testGroup "matrixPrep (block encoding)"
  [ testGroup "Block extraction"
      [ matrixPrepCase "[[0.6,0],[0,0.3]]"
          (fromLists [[0.6:+0, 0:+0], [0:+0, 0.3:+0]])
      , matrixPrepCase "[[1,2],[3,4]] real"
          (fromLists [[1:+0, 2:+0], [3:+0, 4:+0]])
      , matrixPrepCase "[[1+i,2],[3,4-i]] complex"
          (fromLists [[1:+1, 2:+0], [3:+0, 4:+(-1)]])
      , matrixPrepCase "diag(1,2,3)"
          (fromLists [[1:+0, 0:+0, 0:+0], [0:+0, 2:+0, 0:+0], [0:+0, 0:+0, 3:+0]])
      , matrixPrepCase "4x4 (matrixPrepTester4)"
          (fromLists [[1:+2, 3:+4 , 5:+0 , 6:+0   ]
                     ,[3:+0, 4:+0 , 5:+6 , 7:+0   ]
                     ,[0:+0, 0:+0 , 9:+0 , 0:+0   ]
                     ,[0:+0, 0:+0 , 0:+0 , 10:+0  ]])
      ]
  , testGroup "Unitarity"
      [ matrixPrepUnitary "[[0.6,0],[0,0.3]]"
          (fromLists [[0.6:+0, 0:+0], [0:+0, 0.3:+0]])
      , matrixPrepUnitary "[[1,2],[3,4]]"
          (fromLists [[1:+0, 2:+0], [3:+0, 4:+0]])
      , matrixPrepUnitary "diag(1,2,3)"
          (fromLists [[1:+0, 0:+0, 0:+0], [0:+0, 2:+0, 0:+0], [0:+0, 0:+0, 3:+0]])
      , matrixPrepUnitary "complex 2x2"
          (fromLists [[1:+1, 2:+0], [3:+0, 4:+(-1)]])
      ]
  ]

----------------------------------------------------------------------
-- Backend consistency: MS = MPS = SV for the same QOp
----------------------------------------------------------------------

backendConsistencyTests :: TestTree
backendConsistencyTests = testGroup "Backend consistency (MS = MPS = SV)"
  [ testCase "buildRowQOp [1,2,3,4]" $ checkBackends
      (buildRowQOp (V.fromList [1:+0, 2:+0, 3:+0, 4:+0]))
  , testCase "matrixPrep [[0.6,0],[0,0.3]]" $ checkBackends
      (matrixPrep (fromLists [[0.6:+0, 0:+0], [0:+0, 0.3:+0]]))
  , testCase "matrixPrep [[1,2],[3,4]]" $ checkBackends
      (matrixPrep (fromLists [[1:+0, 2:+0], [3:+0, 4:+0]]))
  , testCase "matrixPrep complex 2x2" $ checkBackends
      (matrixPrep (fromLists [[1:+1, 2:+0], [3:+0, 4:+(-1)]]))
  , testCase "matrixPrep diag(1,2,3)" $ checkBackends
      (matrixPrep (fromLists [[1:+0, 0:+0, 0:+0], [0:+0, 2:+0, 0:+0], [0:+0, 0:+0, 3:+0]]))
  ]
  where
    checkBackends :: QOp -> Assertion
    checkBackends op = do
      let mref = MS.evalOp op
          mmps = to (MPS.evalOp op) :: CMat
          msv  = to (SV.evalOp op)  :: CMat
      assertClose "MPS vs MS" mref mmps
      assertClose "SV  vs MS" mref msv

----------------------------------------------------------------------
-- Property tests (QuickCheck)
----------------------------------------------------------------------

-- | Generator: small integer-valued complex number.
genComplexInt :: Gen ComplexT
genComplexInt = do
  re <- choose ((-3), 3) :: Gen Int
  im <- choose ((-3), 3) :: Gen Int
  return (fromIntegral re :+ fromIntegral im)

-- | Generator: random complex vector with length in 2..6 and small integer entries.
data RandomVec = RandomVec (V.Vector ComplexT)
instance Show RandomVec where show (RandomVec v) = "RandomVec " ++ show (V.toList v)

instance Arbitrary RandomVec where
  arbitrary = do
    n  <- elements [2, 3, 4, 5, 6, 7, 8]
    xs <- vectorOf n genComplexInt
    return (RandomVec (V.fromList xs))
  shrink (RandomVec v)
    | V.length v <= 2 = []
    | otherwise       = [RandomVec (V.take (V.length v - 1) v)]

-- | Generator: random square complex matrix with dim ∈ {2,3,4} and small entries.
data RandomMatrix = RandomMatrix CMat
instance Show RandomMatrix where show (RandomMatrix m) = "RandomMatrix " ++ show m

instance Arbitrary RandomMatrix where
  arbitrary = do
    n  <- elements [2, 3, 4] :: Gen Int
    xs <- vectorOf (n*n) genComplexInt
    return (RandomMatrix (fromLists (chunksOf n xs)))
  shrink (RandomMatrix m) =
    let lists = H.toLists m
        n = length lists
    in if n <= 2 then [] else
       [ RandomMatrix (fromLists (init (map init lists))) ]  -- drop last row + col

chunksOf :: Int -> [a] -> [[a]]
chunksOf _ [] = []
chunksOf n xs = take n xs : chunksOf n (drop n xs)

-- Row-prep correctness for random vectors
prop_rowPrep_correct :: RandomVec -> Property
prop_rowPrep_correct (RandomVec v) =
  let nrm = sqrt (realPart (V.sum (V.map (\x -> conjugate x * x) v)))
  in nrm > tol ==> rowPrepStateMS v ~~ expectedRowPrep v

-- Block-encoding correctness for random matrices
prop_matrixPrep_block :: RandomMatrix -> Property
prop_matrixPrep_block (RandomMatrix m) =
  norm_Frob m > tol ==> blockEncodeMS m ~~ expectedBlock m

-- matrixPrep always produces a unitary
prop_matrixPrep_unitary :: RandomMatrix -> Property
prop_matrixPrep_unitary (RandomMatrix m) =
  norm_Frob m > tol ==>
  let u = MS.evalOp (matrixPrep m)
      n = rows u
  in (u H.<> H.tr u) ~~ identCMat n

-- All three backends agree on matrixPrep
prop_matrixPrep_3backends :: RandomMatrix -> Property
prop_matrixPrep_3backends (RandomMatrix m) =
  norm_Frob m > tol ==>
  let op   = matrixPrep m
      mref = MS.evalOp op
      mmps = to (MPS.evalOp op) :: CMat
      msv  = to (SV.evalOp op)  :: CMat
  in counterexample "MPS != MS" (mref ~~ mmps)
     .&&. counterexample "SV != MS" (mref ~~ msv)

qcOpts :: TestTree -> TestTree
qcOpts = localOption (QuickCheckTests 30) . localOption (QuickCheckMaxSize 4)

propertyTests :: TestTree
propertyTests = qcOpts $ testGroup "Properties"
  [ testProperty "buildRowQOp prepares v/‖v‖"     prop_rowPrep_correct
  , testProperty "matrixPrep block = M/‖M‖_F"       prop_matrixPrep_block
  , testProperty "matrixPrep output is unitary"     prop_matrixPrep_unitary
  , testProperty "matrixPrep matches across backends" prop_matrixPrep_3backends
  ]

----------------------------------------------------------------------

main :: IO ()
main = defaultMain $ testGroup "MatrixPrepTest"
  [ pureHelperTests
  , rowPrepTests
  , matrixPrepTests
  , backendConsistencyTests
  , propertyTests
  ]
