-- | Scalar constants shared by the exact evaluators. Independent expected-value
-- tests check this table; backend agreement alone would not detect a shared sign.
module QPP.Exact.Internal.Gates (primitive, permuteIndex, inversePermutation) where

import Data.Bits (testBit)
import Data.List (sortOn)
import qualified Data.List as List
import Data.Ratio ((%))
import qualified Data.Vector as V
import QPP.Exact.Checked (Gate(..))
import QPP.Exact.Types
import QPP.Exact.Internal.Scalar

-- | Four row-major coefficients of a one-qubit gate.
primitive :: Limits -> Field -> Gate -> Result (V.Vector Coeff)
primitive l f gate = do
  let z = zero f
      o = one f
      no = negC o
  case gate of
    X -> pure (V.fromList [z,o,o,z])
    Y -> do
      i <- rootPower l f 4 1
      pure (V.fromList [z,negC i,i,z])
    Z -> pure (V.fromList [o,z,z,no])
    H -> do
      w <- rootPower l f 8 1
      wi <- rootPower l f 8 (-1)
      h <- addC l w wi >>= scaleC l (1%2)
      pure (V.fromList [h,h,h,negC h])
    SX -> do
      i <- rootPower l f 4 1
      a <- addC l o i >>= scaleC l (1%2)
      b <- subC l o i >>= scaleC l (1%2)
      pure (V.fromList [a,b,b,a])

-- | Scatter convention: output bit j is input bit ks[j].
permuteIndex :: [Int] -> Int -> Int
permuteIndex ks input = List.foldl' step 0 ks
  where step acc k = 2*acc + if testBit input (length ks-1-k) then 1 else 0

inversePermutation :: [Int] -> [Int]
inversePermutation ks = map snd (sortOn fst (zip ks [0..]))
