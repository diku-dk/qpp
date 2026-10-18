{-# LANGUAGE BangPatterns #-}
-- | Private power-basis representation and checked polynomial kernels.
-- Constructors never escape through an exposed module.
module QPP.Exact.Internal.Scalar
  ( Field, Coeff, Cyclo(..), Result, Work(..)
  , index, degree, modulus, coordinates, importCoeff
  , validateLimits, bound, slots, bits, checkInteger, checkRational
  , intAdd, intMul
  , storage, reserve, workspace, makeField, fieldDegree, fieldFor, commonField
  , checkField, checkCoeff, jointIndex, dimension, widthOfLength
  , constant, zero, one, isZeroC, asRationalC
  , addC, subC, negC, scaleC, mulC, powC, inverseC
  , rootPower, embedC, conjugateC, phaseC, rotationC
  ) where

import Control.Monad (foldM, unless, when)
import Data.Bits (shiftL)
import Data.List (group)
import qualified Data.List as List
import Data.Ratio ((%), numerator, denominator)
import qualified Data.Vector as V
import Math.NumberTheory.Logarithms (integerLog2)
import QPP.Exact.Types

type Poly = V.Vector Rational -- ascending powers; no trailing zeros

-- The modulus is monic with integer coefficients. The embedding is always
-- zeta_N = exp(2*pi*i/N), including when two indices describe the same field.
data Field = Field !Integer !Poly deriving (Eq, Show)
newtype Coeff = Coeff (V.Vector Rational) deriving (Eq, Show)
data Cyclo = Cyclo !Field !Coeff

instance Show Cyclo where
  showsPrec _ (Cyclo f (Coeff xs)) =
    showString "Cyclo " . shows (index f) . showChar ' ' . shows (V.toList xs)

index :: Field -> Integer
index (Field n _) = n

degree :: Field -> Int
degree (Field _ p) = V.length p - 1

modulus :: Field -> Poly
modulus (Field _ p) = p

coordinates :: Coeff -> V.Vector Rational
coordinates (Coeff p) = p

importCoeff :: Limits -> Field -> [Rational] -> Result Coeff
importCoeff l f xs = do
  unless (length xs == degree f)
    (Left (InvalidInput "coordinate count must equal the field degree"))
  Coeff . V.fromList <$> traverse (checkRational l) xs

validateLimits :: Limits -> Result ()
validateLimits l = unless (maxRootIndex l >= 1 && maxFieldDegree l >= 1
                       && maxCoefficientSlots l >= 1 && maxIntegerBits l >= 1) $
  Left (InvalidInput "exact limits must all be positive")

bound :: String -> Integer -> Integer -> Result ()
bound label requested allowed = when (requested > allowed) $
  Left (ResourceLimit label requested allowed)

slots :: Limits -> Integer -> Result ()
slots l n = do
  bound "coefficient slots" n (maxCoefficientSlots l)
  bound "array length" n (toInteger (maxBound :: Int))

bits :: Integer -> Int
bits 0 = 0
bits n = 1 + integerLog2 (abs n)

checkInteger :: Limits -> Integer -> Result ()
checkInteger l n = bound "integer bits" (toInteger (bits n)) (toInteger (maxIntegerBits l))

-- Conservative prechecks also bound unreduced intermediate integers.
intMul :: Limits -> Integer -> Integer -> Result Integer
intMul l a b
  | a == 0 || b == 0 = Right 0
  | abs a == 1 = checked (a*b)
  | abs b == 1 = checked (a*b)
  | otherwise = do
      bound "integer bits (product)" (toInteger (bits a) + toInteger (bits b))
        (toInteger (maxIntegerBits l))
      checked (a*b)
  where checked x = checkInteger l x >> pure x

intAdd :: Limits -> Integer -> Integer -> Result Integer
intAdd l a b = do
  let extra = if a /= 0 && b /= 0 && signum a == signum b then 1 else 0
  bound "integer bits (sum)" (toInteger (max (bits a) (bits b)) + extra)
    (toInteger (maxIntegerBits l))
  let !c = a+b
  checkInteger l c
  pure c

checkRational :: Limits -> Rational -> Result Rational
checkRational l q = do
  unless (denominator q > 0) (Left (InvalidInput "rational denominator must be positive"))
  checkInteger l (numerator q)
  checkInteger l (denominator q)
  pure (numerator q % denominator q)

ratAdd :: Limits -> Rational -> Rational -> Result Rational
ratAdd l a b
  | a == 0 = pure b
  | b == 0 = pure a
  | otherwise = do
      let g = gcd (denominator a) (denominator b)
          u = denominator a `quot` g
          v = denominator b `quot` g
      x <- intMul l (numerator a) v
      y <- intMul l (numerator b) u
      z <- intAdd l x y
      d <- intMul l (denominator a) v
      checkRational l (z % d)

ratMul :: Limits -> Rational -> Rational -> Result Rational
ratMul l a b
  | a == 0 || b == 0 = pure 0
  | otherwise = do
      let g = gcd (numerator a) (denominator b)
          h = gcd (numerator b) (denominator a)
      n <- intMul l (numerator a `quot` g) (numerator b `quot` h)
      d <- intMul l (denominator a `quot` h) (denominator b `quot` g)
      checkRational l (n % d)

ratDiv :: Limits -> Rational -> Rational -> Result Rational
ratDiv _ _ 0 = Left DivisionByZero
ratDiv l a b = ratMul l a (denominator b % numerator b)

trim :: Poly -> Poly
trim = V.reverse . V.dropWhile (==0) . V.reverse

at :: Poly -> Int -> Rational
at p i = if i < V.length p then p V.! i else 0

polyAdd :: Limits -> Poly -> Poly -> Result Poly
polyAdd l a b = do
  let n = max (V.length a) (V.length b)
  slots l (3 * toInteger n)
  trim <$> V.generateM n (\i -> ratAdd l (at a i) (at b i))

polyMul :: Limits -> Poly -> Poly -> Result Poly
polyMul l a b
  | V.null a || V.null b = pure V.empty
  | otherwise = do
      let n = V.length a + V.length b - 1
      slots l (4 * toInteger n)
      trim <$> V.generateM n (\k -> foldM (term k) 0
        [max 0 (k - V.length b + 1) .. min k (V.length a - 1)])
  where
    term k acc i = ratMul l (a V.! i) (b V.! (k-i)) >>= ratAdd l acc

polyDivMod :: Limits -> Poly -> Poly -> Result (Poly, Poly)
polyDivMod _ _ b | V.null b = Left DivisionByZero
polyDivMod l a b = do
  slots l (8 * toInteger (max (V.length a) (V.length b)))
  go (V.replicate (max 0 (V.length a - V.length b + 1)) 0) (trim a)
  where
    go !q !r
      | V.length r < V.length b = pure (trim q, r)
      | otherwise = do
          let k = V.length r - V.length b
          c <- ratDiv l (V.last r) (V.last b)
          r' <- V.generateM (V.length r) $ \i ->
            if i < k then pure (r V.! i) else do
              x <- ratMul l c (b V.! (i-k))
              ratAdd l (r V.! i) (-x)
          go (q V.// [(k,c)]) (trim r')

primeFactors :: Integer -> [Integer]
primeFactors = go 2
  where
    go p n
      | n == 1 = []
      | p*p > n = [n]
      | n `mod` p == 0 = p : go p (n `quot` p)
      | otherwise = go (if p == 2 then 3 else p+2) n

-- Validation charges only the stored coordinates and the field modulus.
-- Arithmetic reserves logical coordinate slots, not boxed Rational heap bytes.
storage :: Limits -> Field -> Integer -> Result ()
storage l f resident = slots l (resident + toInteger (degree f) + 1)

data Work = LinearWork | PolynomialWork deriving (Eq)

reserve :: Limits -> Field -> Work -> Integer -> Result ()
reserve l f work resident = storage l f (resident + scratch)
  where
    -- Linear kernels need only a constant number of Rational accumulators.
    -- Polynomial kernels reserve sixteen vectors of length at most 2*d+1:
    -- the four EEA iterates, quotient/remainder and Bezout updates, normalized
    -- copies, and division/reduction temporaries. Multiplication, conjugation
    -- and embedding use a subset of this allowance. No syntax-depth factor is
    -- needed here; evaluator-owned live buffers are charged by their callers.
    scratch = case work of
      LinearWork -> 4
      PolynomialWork -> 16 * (2 * toInteger (degree f) + 1)

workspace :: Limits -> Field -> Integer -> Result ()
workspace l f = reserve l f PolynomialWork

fieldDegree :: Limits -> Integer -> Result Int
fieldDegree l n = fst <$> describeField l n

describeField :: Limits -> Integer -> Result (Int, [Integer])
describeField l n = do
  validateLimits l
  unless (n > 0) (Left (InvalidInput "root index must be positive"))
  bound "root index" n (maxRootIndex l)
  let fs = primeFactors n
      d = List.foldl' (\acc p -> acc `quot` p * (p-1)) n [p | p:_ <- group fs]
  bound "field degree" d (toInteger (maxFieldDegree l))
  pure (fromInteger d,fs)

makeField :: Limits -> Integer -> Result Field
makeField l n = do
  (d,fs) <- describeField l n
  -- Each recursive level retains a modulus and a lifted polynomial of at most
  -- 2*d+1 coordinates. Eight such buffers cover division and its copies;
  -- the factor count bounds the recursion depth. Q itself needs just x-1.
  slots l (if n == 1 then 2 else 8 * (2 * toInteger d + 1) * toInteger (length fs + 1))
  p <- cyclotomic n fs
  pure (Field n p)
  where
    cyclotomic 1 _ = pure (V.fromList [-1,1])
    cyclotomic k (p:ps) = do
      let m = k `quot` p
      a <- cyclotomic m ps
      let len = (toInteger (V.length a)-1)*p+1
      slots l (8*len)
      let lifted = V.generate (fromInteger len) $ \j ->
            if toInteger j `mod` p == 0 then a V.! fromInteger (toInteger j `quot` p) else 0
      if m `mod` p == 0 then pure lifted else do
        (q,r) <- polyDivMod l lifted a
        unless (V.null r) (Left (InvalidInput "internal cyclotomic polynomial division"))
        pure q
    cyclotomic _ [] = Left (InvalidInput "internal field factorization")

checkField :: Limits -> Field -> Result ()
checkField l f = do
  validateLimits l
  bound "root index" (index f) (maxRootIndex l)
  bound "field degree" (toInteger (degree f)) (toInteger (maxFieldDegree l))
  mapM_ (checkRational l) (modulus f)

checkCoeff :: Limits -> Field -> Coeff -> Result ()
checkCoeff l f (Coeff a) = do
  unless (V.length a == degree f) (Left (InvalidInput "internal coefficient degree mismatch"))
  mapM_ (checkRational l) a

-- Context selection belongs in one place. Candidates have already passed
-- checkField at the public boundary; immutable matching contexts can be reused.
fieldFor :: Limits -> Integer -> [Field] -> Result Field
fieldFor l n candidates = case filter ((== n) . index) candidates of
  f:_ -> pure f
  [] -> makeField l n

commonField :: Limits -> Field -> Field -> Result Field
commonField l a b = do
  n <- jointIndex l (index a) (index b)
  fieldFor l n [a,b]

jointIndex :: Limits -> Integer -> Integer -> Result Integer
jointIndex l n m = do
  let k = (n `quot` gcd n m)*m
  bound "root index" k (maxRootIndex l)
  pure k

constant :: Field -> Rational -> Coeff
constant f x = Coeff (V.generate (degree f) (\i -> if i == 0 then x else 0))

zero :: Field -> Coeff
zero f = constant f 0

one :: Field -> Coeff
one f = constant f 1

isZeroC :: Coeff -> Bool
isZeroC (Coeff p) = V.all (==0) p

asRationalC :: Coeff -> Maybe Rational
asRationalC (Coeff p)
  | V.all (==0) (V.tail p) = Just (V.head p)
  | otherwise = Nothing

reduce :: Limits -> Field -> Poly -> Result Coeff
reduce l f p = do
  (_,r) <- polyDivMod l (trim p) (modulus f)
  pure (Coeff (V.generate (degree f) (at r)))

addC :: Limits -> Coeff -> Coeff -> Result Coeff
addC l = zipCoefficients (ratAdd l)

-- Fuse subtraction into the coordinate loop: a separately negated vector
-- would need degree(f) temporary slots, exceeding the linear-work allowance.
subC :: Limits -> Coeff -> Coeff -> Result Coeff
subC l = zipCoefficients (\a b -> ratAdd l a (-b))

zipCoefficients :: (Rational -> Rational -> Result Rational) -> Coeff -> Coeff -> Result Coeff
zipCoefficients action (Coeff a) (Coeff b) = do
  unless (V.length a == V.length b)
    (Left (InvalidInput "internal coefficient degree mismatch"))
  Coeff <$> V.zipWithM action a b

negC :: Coeff -> Coeff
negC (Coeff a) = Coeff (V.map negate a)

scaleC :: Limits -> Rational -> Coeff -> Result Coeff
scaleC l c (Coeff a) = Coeff <$> V.mapM (ratMul l c) a

mulC :: Limits -> Field -> Coeff -> Coeff -> Result Coeff
mulC l f a@(Coeff p) b@(Coeff q)
  | Just c <- asRationalC a = scaleC l c b
  | Just c <- asRationalC b = scaleC l c a
  | otherwise = polyMul l (trim p) (trim q) >>= reduce l f

powC :: Limits -> Field -> Coeff -> Integer -> Result Coeff
powC l f a n
  | n < 0 = inverseC l f a >>= \b -> powC l f b (-n)
  | otherwise = go (one f) a n
  where
    go !acc _ 0 = pure acc
    go !acc !b k = do
      acc' <- if odd k then mulC l f acc b else pure acc
      if k == 1 then pure acc' else do
        b' <- mulC l f b b
        go acc' b' (k `quot` 2)

rootPower :: Limits -> Field -> Integer -> Integer -> Result Coeff
rootPower l f order powerIndex = do
  unless (order > 0 && index f `mod` order == 0)
    (Left (InvalidInput "root order must divide the field index"))
  let k = ((powerIndex `mod` order) * (index f `quot` order)) `mod` index f
  if k < toInteger (degree f)
    then pure (Coeff (V.generate (degree f) (\j -> if toInteger j == k then 1 else 0)))
    else do
      x <- reduce l f (V.fromList [0,1])
      powC l f x k

-- Horner substitution avoids constructing degree N*d intermediate polynomials.
substitute :: Limits -> Field -> Coeff -> Coeff -> Result Coeff
substitute l f x (Coeff p) = V.foldM' step (zero f) (V.reverse p)
  where step acc c = mulC l f acc x >>= \a -> addC l a (constant f c)

embedC :: Limits -> Field -> Field -> Coeff -> Result Coeff
embedC l target source c
  | index target == index source = pure c
  | otherwise = rootPower l target (index source) 1 >>= \x -> substitute l target x c

conjugateC :: Limits -> Field -> Coeff -> Result Coeff
conjugateC l f c = rootPower l f (index f) (-1) >>= \x -> substitute l f x c

inverseC :: Limits -> Field -> Coeff -> Result Coeff
inverseC _ _ a | isZeroC a = Left DivisionByZero
inverseC l f (Coeff a) = do
  (r,t) <- monic (trim a) (V.singleton 1)
  go (modulus f) r V.empty t
  where
    -- Maintain r_j = s_j*Phi + t_j*a. Dividing BOTH r_j and t_j by
    -- the leading coefficient keeps this identity while removing the scalar
    -- content that otherwise swells through successive Euclidean remainders.
    monic r t = do
      let leading = V.last r -- only called for a nonzero remainder
      r' <- V.mapM (\x -> ratDiv l x leading) r
      t' <- V.mapM (\x -> ratDiv l x leading) t
      pure (r',t')
    go r0 r1 t0 t1
      | V.length r1 == 1 = reduce l f t1 -- the monic constant gcd is 1
      | otherwise = do
          (q,r2) <- polyDivMod l r0 r1
          when (V.null r2) (Left (InvalidInput "internal non-field modulus"))
          qt <- polyMul l q t1
          t2 <- polyAdd l t0 (V.map negate qt)
          (r,t) <- monic r2 t2
          go r1 r t1 t

phaseC :: Limits -> Field -> Rational -> Result Coeff
phaseC l f angle = do
  q <- checkRational l angle
  rootPower l f (2*denominator q) (numerator q)

rotationC :: Limits -> Field -> Rational -> Result (Coeff,Coeff)
rotationC l f angle = do
  q <- checkRational l angle
  w <- rootPower l f (4*denominator q) (numerator q)
  wi <- rootPower l f (4*denominator q) (-numerator q)
  c <- addC l w wi >>= scaleC l (1%2)
  d <- subC l wi w >>= scaleC l (1%2)
  pure (c,d)

-- Check before shifting, including Id maxBound: never construct 2^n merely
-- to discover that its index or allocation cannot fit.
dimension :: Limits -> Int -> Result Int
dimension l n = do
  validateLimits l
  when (n < 0) (Left (InvalidInput "negative qubit count"))
  let cap = min (maxCoefficientSlots l) (toInteger (maxBound :: Int))
      maxN = integerLog2 cap
  bound "qubits for vector dimension" (toInteger n) (toInteger maxN)
  pure (1 `shiftL` n)

widthOfLength :: Int -> Result Int
widthOfLength n
  | n <= 0 = Left (InvalidInput "amplitude count must be a nonzero power of two")
  | otherwise = go n 0
  where
    go 1 k = pure k
    go m k | even m = go (m `quot` 2) (k+1)
           | otherwise = Left (InvalidInput "amplitude count must be a nonzero power of two")
