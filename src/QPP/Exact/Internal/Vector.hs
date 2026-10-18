-- | Statevector operations within one established cyclotomic field.
-- Callers check dimensions and reserve workspace before using these kernels.
-- No operation normalizes its input: intermediate vectors need not have norm 1.
module QPP.Exact.Internal.Vector
  ( Vector, Transform
  , scale, linearCombination, applyGate, directSum, tensor, permute
  ) where

import qualified Data.Vector as V
import QPP.Exact.Types
import QPP.Exact.Internal.Scalar (Field, Coeff)
import qualified QPP.Exact.Internal.Scalar as S
import QPP.Exact.Internal.Gates (permuteIndex)

type Vector = V.Vector Coeff
type Transform = Vector -> Result Vector

scale :: Limits -> Field -> Coeff -> Transform
scale limits field c = V.mapM (S.mulC limits field c)

-- | c*v + d*w, computed entry by entry without two scaled temporary vectors.
linearCombination :: Limits -> Field -> Coeff -> Vector -> Coeff -> Vector -> Result Vector
linearCombination limits field c v d w =
  V.zipWithM (\x y -> weightedSum limits field c x d y) v w

-- | A one-qubit gate stores its four entries in row-major order.
-- The output is (g00*v0 + g01*v1, g10*v0 + g11*v1).
applyGate :: Limits -> Field -> V.Vector Coeff -> Transform
applyGate limits field gate input = V.generateM 2 $ \row ->
  weightedSum limits field
    (gate V.! (2*row))     (input V.! 0)
    (gate V.! (2*row + 1)) (input V.! 1)

weightedSum :: Limits -> Field -> Coeff -> Coeff -> Coeff -> Coeff -> Result Coeff
weightedSum limits field c x d y = do
  cx <- S.mulC limits field c x
  dy <- S.mulC limits field d y
  S.addC limits cx dy

-- | (A direct-sum B)(v0 || v1) = A(v0) || B(v1).
-- The leading qubit selects two equally sized, contiguous blocks.
directSum :: Transform -> Transform -> Transform
directSum applyA applyB input = do
  let half = V.length input `quot` 2
      (v0,v1) = V.splitAt half input
  w0 <- applyA v0
  w1 <- applyB v1
  pure (w0 V.++ w1)

-- | Apply A tensor B to an arbitrary vector, including an entangled one.
-- View v[i*dimB+j] as an array V[i,j]. Applying B to each row, then A to
-- each column, gives A V B^T. The transpose is ordinary, not conjugate:
-- each row is passed to B as a column vector. No operator matrix is built.
tensor :: Int -> Int -> Transform -> Transform -> Transform
tensor dimA dimB applyA applyB input = do
  afterB <- V.generateM dimA $ \i ->
    applyB (V.slice (i*dimB) dimB input)

  afterA <- V.generateM dimB $ \j ->
    let column = V.generate dimA (\i -> (afterB V.! i) V.! j)
    in applyA column

  -- afterA stores columns; flatten back to the original row-major order.
  pure $ V.generate (dimA*dimB) $ \index ->
    let (i,j) = index `quotRem` dimB
    in (afterA V.! j) V.! i

-- | The plan stores a gather map: each output position names its input.
permute :: [Int] -> Vector -> Vector
permute gather input = V.generate (V.length input) $ \output ->
  input V.! permuteIndex gather output
