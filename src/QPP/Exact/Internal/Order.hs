{-# LANGUAGE BangPatterns #-}
-- | Certified ordering in the distinguished complex embedding. This is a
-- sign test on real cyclotomic values, not an approximate scalar conversion.
module QPP.Exact.Internal.Order (signReal) where

import Control.Monad (foldM)
import Data.Bits (shiftL)
import Data.Ratio (numerator, denominator)
import qualified Data.Vector as V
import QPP.Exact.Types
import QPP.Exact.Internal.Scalar (Field, Coeff)
import qualified QPP.Exact.Internal.Scalar as S

-- Bounds lo hi represents [lo/2^precision, hi/2^precision]. Every division
-- rounds outwards. The caller has validated the coefficient and proved it real.
data Bounds = Bounds !Integer !Integer

signReal :: Limits -> Field -> Coeff -> Result Ordering
signReal limits field value
  | Just q <- S.asRationalC value = pure (compare q 0)
  | otherwise = do
      S.reserve limits field S.LinearWork (toInteger (S.degree field) + 64)
      refine 16
  where
    -- Exact zero is caught by the rational case. A nonzero real algebraic
    -- number is eventually separated from zero, unless a resource limit wins.
    refine precision = do
      S.bound "comparison precision bits" (2 * toInteger precision + 16)
        (toInteger (maxIntegerBits limits))
      Bounds lo hi <- realBounds limits precision field value
      if lo > 0 then pure GT
      else if hi < 0 then pure LT
      else refine (2 * precision)

realBounds :: Limits -> Int -> Field -> Coeff -> Result Bounds
realBounds limits precision field value = do
  -- Machin's identity: pi = 16 atan(1/5) - 4 atan(1/239).
  a <- arctanInverse 5
  b <- arctanInverse 239
  a' <- scale 16 1 a
  b' <- scale (-4) 1 b
  piBounds <- add a' b'
  foldM (term piBounds) (Bounds 0 0) (zip [0..] (V.toList (S.coordinates value)))
  where
    unit = 1 `shiftL` precision
    one = Bounds unit unit
    negateBounds (Bounds lo hi) = Bounds (-hi) (-lo)
    ceilDiv a b = negate ((negate a) `div` b)

    add (Bounds a b) (Bounds c d) =
      Bounds <$> S.intAdd limits a c <*> S.intAdd limits b d

    -- Scale by a rational n/d, d > 0, without increasing dyadic precision.
    scale n d (Bounds lo hi) = do
      S.checkInteger limits d
      a <- S.intMul limits n (if n >= 0 then lo else hi)
      b <- S.intMul limits n (if n >= 0 then hi else lo)
      pure (Bounds (a `div` d) (ceilDiv b d))

    multiply (Bounds a b) (Bounds c d) = do
      products <- traverse (uncurry (S.intMul limits)) [(a,c),(a,d),(b,c),(b,d)]
      pure (Bounds (minimum products `div` unit) (ceilDiv (maximum products) unit))

    widen radius (Bounds lo hi) =
      Bounds <$> S.intAdd limits lo (-radius) <*> S.intAdd limits hi radius

    -- Positive arctan terms decrease from the start. Once the next term is
    -- at most one ulp, the alternating remainder fits in [-1,1] ulps.
    arctanInverse q = do
      initial <- scale 1 q one
      let go !k !t@(Bounds _ hi) !acc
            | hi <= 1 = widen hi acc
            | otherwise = do
                acc' <- add acc (if even k then t else negateBounds t)
                denom <- S.intMul limits (q*q) (2*k+3)
                next <- scale (2*k+1) denom t
                go (k+1) next acc'
      go 0 initial (Bounds 0 0)

    -- Arguments are in [0,pi]. Terms x^(2k)/(2k)! decrease for k >= 2;
    -- the first omitted term bounds the alternating cosine remainder.
    cosine angle = do
      square <- multiply angle angle
      let go !k !t@(Bounds _ hi) !acc
            | k >= 2 && hi <= 1 = widen hi acc
            | otherwise = do
                acc' <- add acc (if even k then t else negateBounds t)
                denom <- S.intMul limits (2*k+1) (2*k+2)
                next <- multiply t square >>= scale 1 denom
                go (k+1) next acc'
      go 0 one (Bounds 0 0)

    term _ acc (_,0) = pure acc
    term piBounds acc (j,q) = do
      let n = S.index field
          k = min j (n-j) -- cos(2*pi*j/n) = cos(2*pi*(n-j)/n)
      c <- if j == 0 then pure one else scale (2*k) n piBounds >>= cosine
      weighted <- scale (numerator q) (denominator q) c
      add acc weighted
