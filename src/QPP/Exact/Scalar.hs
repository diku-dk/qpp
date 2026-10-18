-- | Exact complex cyclotomic numbers. There is intentionally no Floating,
-- RealFloat or partial Fractional instance. All arithmetic failures are explicit.
-- Coefficients use the power basis of the distinguished root exp(2*pi*i/N).
module QPP.Exact.Scalar
  ( Cyclo, Limits(..), defaultLimits, ExactError(..), Result
  , rational, rootOfUnity, phase, rotationCoefficients, fromCoordinates
  , fieldIndex, fieldDegree, coefficients, cyclotomicPolynomial, asRational
  , add, sub, mul, divide, negate, conjugate, inverse, pow, embed, equal, isZero
  , compareReal, approximate
  ) where

import Prelude hiding (negate)
import Control.Monad (unless)
import Data.Complex (Complex, cis)
import Data.Ratio (numerator, denominator)
import qualified Data.Vector as V
import QPP.Exact.Types
import QPP.Exact.Internal.Scalar (Cyclo(..), Coeff, Work(..))
import qualified QPP.Exact.Internal.Scalar as S
import qualified QPP.Exact.Internal.Order as Order

fieldIndex :: Cyclo -> Integer
fieldIndex (Cyclo f _) = S.index f

fieldDegree :: Cyclo -> Int
fieldDegree (Cyclo f _) = S.degree f

coefficients :: Cyclo -> [Rational]
coefficients (Cyclo _ c) = V.toList (S.coordinates c)

cyclotomicPolynomial :: Cyclo -> [Integer]
cyclotomicPolynomial (Cyclo f _) = map numerator (V.toList (S.modulus f))

asRational :: Cyclo -> Maybe Rational
asRational (Cyclo _ c) = S.asRationalC c

isZero :: Cyclo -> Bool
isZero (Cyclo _ c) = S.isZeroC c

rational :: Limits -> Rational -> Result Cyclo
rational l q = do
  f <- S.makeField l 1
  r <- S.checkRational l q
  S.storage l f 1
  pure (Cyclo f (S.constant f r))

rootOfUnity :: Limits -> Integer -> Result Cyclo
rootOfUnity l n = do
  f <- S.makeField l n
  S.workspace l f (toInteger (S.degree f))
  Cyclo f <$> S.rootPower l f n 1

-- | exp(i*pi*q), with q a rational angle in QPP's convention.
phase :: Limits -> Rational -> Result Cyclo
phase l angle = do
  q <- S.checkRational l angle
  f <- S.makeField l (2 * denominator q)
  S.workspace l f (toInteger (S.degree f))
  Cyclo f <$> S.phaseC l f q

-- | (cos(pi*q/2), -i*sin(pi*q/2)), computed by roots of unity.
rotationCoefficients :: Limits -> Rational -> Result (Cyclo,Cyclo)
rotationCoefficients l angle = do
  q <- S.checkRational l angle
  f <- S.makeField l (4 * denominator q)
  S.workspace l f (4 * toInteger (S.degree f))
  (c,d) <- S.rotationC l f q
  pure (Cyclo f c, Cyclo f d)

-- | Import exactly phi(N) rational coordinates in the reduced power basis.
fromCoordinates :: Limits -> Integer -> [Rational] -> Result Cyclo
fromCoordinates l n cs = do
  f <- S.makeField l n
  S.storage l f (toInteger (S.degree f))
  Cyclo f <$> S.importCoeff l f cs

validate :: Limits -> Cyclo -> Result ()
validate l (Cyclo f c) = do
  S.checkField l f
  S.storage l f (toInteger (S.degree f))
  S.checkCoeff l f c

unary :: Limits -> Work -> (S.Field -> Coeff -> Result Coeff) -> Cyclo -> Result Cyclo
unary l work action x@(Cyclo f c) = do
  validate l x
  S.reserve l f work (2 * toInteger (S.degree f))
  Cyclo f <$> action f c

-- Reuse a context when possible; equality always compares in a common field.
binary :: Limits -> Work -> (S.Field -> Coeff -> Coeff -> Result a) -> Cyclo -> Cyclo -> Result a
binary l work action x@(Cyclo f a) y@(Cyclo g b) = do
  validate l x
  validate l y
  k <- S.commonField l f g
  let same = S.index f == S.index g
  S.reserve l k (if same then work else PolynomialWork)
    ((if same then 3 else 5) * toInteger (S.degree k))
  a' <- S.embedC l k f a
  b' <- S.embedC l k g b
  action k a' b'

add, sub, mul, divide :: Limits -> Cyclo -> Cyclo -> Result Cyclo
add l = binary l LinearWork (\f a b -> Cyclo f <$> S.addC l a b)
sub l = binary l LinearWork (\f a b -> Cyclo f <$> S.subC l a b)
mul l = binary l PolynomialWork (\f a b -> Cyclo f <$> S.mulC l f a b)
divide l = binary l PolynomialWork (\f a b -> S.inverseC l f b >>= \bi -> Cyclo f <$> S.mulC l f a bi)

negate, conjugate, inverse :: Limits -> Cyclo -> Result Cyclo
negate l = unary l LinearWork (\_ -> pure . S.negC)
conjugate l = unary l PolynomialWork (S.conjugateC l)
inverse l = unary l PolynomialWork (S.inverseC l)

pow :: Limits -> Cyclo -> Integer -> Result Cyclo
pow l a n = S.checkInteger l n >> unary l PolynomialWork (\f c -> S.powC l f c n) a

-- | Promote to a multiple of the stored index. No subfield descent is implied.
embed :: Limits -> Integer -> Cyclo -> Result Cyclo
embed l n x@(Cyclo f c) = do
  validate l x
  if n == S.index f then pure x else do
    g <- S.fieldFor l n [f]
    S.workspace l g (2 * toInteger (S.degree g))
    Cyclo g <$> S.embedC l g f c

equal :: Limits -> Cyclo -> Cyclo -> Result Bool
equal l x@(Cyclo f a) y@(Cyclo g b)
  | S.index f == S.index g = validate l x >> validate l y >> pure (a == b)
  | otherwise = binary l PolynomialWork (\_ a' b' -> pure (a' == b')) x y

-- | Exact ordering of real values in the distinguished embedding. Nonreal
-- inputs are rejected. Refinement may fail with ResourceLimit, never a guess.
compareReal :: Limits -> Cyclo -> Cyclo -> Result Ordering
compareReal l = binary l PolynomialWork $ \f a b -> do
  ac <- S.conjugateC l f a
  bc <- S.conjugateC l f b
  unless (a == ac && b == bc) (Left (InvalidInput "compareReal expects real scalars"))
  difference <- S.subC l a b
  Order.signReal l f difference

-- | Lossy observation only; never used by equality, evaluation or zero tests.
-- Large rational coordinates can overflow or suffer catastrophic cancellation.
approximate :: Cyclo -> Complex Double
approximate (Cyclo f c) = V.ifoldl' step 0 (S.coordinates c)
  where
    step acc j q = acc + fromRational q * cis (2 * pi * fromIntegral j / fromInteger (S.index f))
