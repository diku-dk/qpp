-- | Field-owned coefficient buffers, independent of quantum shapes. Public
-- boundaries validate once; promotion and kernels then share a context.
module QPP.Exact.Internal.Buffer
  ( Buffer(..), field, values, count, validate, fromScalars, scalars
  , promote, preparePair, scale, add, sub, inner, equal, tensor
  ) where

import Control.Monad (foldM, unless)
import qualified Data.Vector as V
import QPP.Exact.Types
import QPP.Exact.Internal.Scalar (Field, Coeff, Cyclo(..), Work(..))
import qualified QPP.Exact.Internal.Scalar as S

-- Invariant: every coefficient belongs to this field and has degree(f) slots.
-- The constructor is shared only by the library's storage/evaluator modules.
-- This module knows array lengths, never qubit widths or matrix dimensions.
data Buffer = Buffer !Field !(V.Vector Coeff) deriving (Show)

field :: Buffer -> Field
field (Buffer f _) = f

values :: Buffer -> V.Vector Coeff
values (Buffer _ xs) = xs

count :: Buffer -> Integer
count = toInteger . V.length . values

scalars :: Buffer -> [Cyclo]
scalars (Buffer f xs) = map (Cyclo f) (V.toList xs)

validate :: Limits -> Buffer -> Result ()
validate l b@(Buffer f xs) = do
  S.checkField l f
  S.storage l f (count b * toInteger (S.degree f))
  V.mapM_ (S.checkCoeff l f) xs

fromScalars :: Limits -> [Cyclo] -> Result Buffer
fromScalars l xs = do
  let fields = [f | Cyclo f _ <- xs]
  mapM_ (\(Cyclo f c) -> S.checkField l f >> S.checkCoeff l f c) xs
  idx <- foldM (S.jointIndex l) 1 (map S.index fields)
  f <- S.fieldFor l idx fields
  let work = if all ((== idx) . S.index) fields then LinearWork else PolynomialWork
  S.reserve l f work (2 * toInteger (length xs) * toInteger (S.degree f))
  Buffer f . V.fromList <$> traverse (\(Cyclo g c) -> S.embedC l f g c) xs

-- Inputs have been validated and the caller has reserved source/target storage.
promote :: Limits -> Field -> Buffer -> Result Buffer
promote l target b@(Buffer source xs)
  | S.index source == S.index target = pure b
  | otherwise = Buffer target <$> V.mapM (S.embedC l target source) xs

-- Charge inputs, possible promoted copies and the caller's output before any
-- promotion. Shape checks belong to State/Matrix, before entering this kernel.
preparePair :: Limits -> Work -> Integer -> Buffer -> Buffer -> Result (Buffer, Buffer)
preparePair l work output a b = do
  validate l a
  validate l b
  prepareValidated l work output a b

prepareValidated :: Limits -> Work -> Integer -> Buffer -> Buffer -> Result (Buffer, Buffer)
prepareValidated l work output a b = do
  f <- S.commonField l (field a) (field b)
  let changes = S.index (field a) /= S.index f || S.index (field b) /= S.index f
      copies = if changes then 2 else 1
      actualWork = if changes then PolynomialWork else work
  S.reserve l f actualWork ((copies * (count a + count b) + output) * toInteger (S.degree f))
  (,) <$> promote l f a <*> promote l f b

sameLength :: Buffer -> Buffer -> Result ()
sameLength a b = unless (count a == count b)
  (Left (InvalidInput "internal buffer length mismatch"))

scale :: Limits -> Cyclo -> Buffer -> Result Buffer
scale l (Cyclo scalarField c) b = do
  validate l b
  S.checkField l scalarField
  S.checkCoeff l scalarField c
  f <- S.commonField l scalarField (field b)
  -- Original/promoted/output arrays plus the original/promoted scalar.
  S.workspace l f ((3 * count b + 2) * toInteger (S.degree f))
  c' <- S.embedC l f scalarField c
  Buffer _ xs <- promote l f b
  Buffer f <$> V.mapM (S.mulC l f c') xs

add, sub :: Limits -> Buffer -> Buffer -> Result Buffer
add l = combine l (S.addC l)
sub l = combine l (S.subC l)

combine :: Limits -> (Coeff -> Coeff -> Result Coeff) -> Buffer -> Buffer -> Result Buffer
combine l action a b = do
  sameLength a b
  (Buffer f xs, Buffer _ ys) <- preparePair l LinearWork (count a) a b
  Buffer f <$> V.zipWithM action xs ys

inner :: Limits -> Buffer -> Buffer -> Result Cyclo
inner l a b = do
  sameLength a b
  (Buffer f xs, Buffer _ ys) <- preparePair l PolynomialWork 3 a b
  Cyclo f <$> V.ifoldM' (\acc j x -> do
    xc <- S.conjugateC l f x
    p <- S.mulC l f xc (ys V.! j)
    S.addC l acc p) (S.zero f) xs

equal :: Limits -> Buffer -> Buffer -> Result Bool
equal l a b = do
  validate l a
  validate l b
  if count a /= count b then pure False
  else if S.index (field a) == S.index (field b) then pure (values a == values b)
  else do
    (Buffer _ xs, Buffer _ ys) <- prepareValidated l PolynomialWork 0 a b
    pure (xs == ys)

tensor :: Limits -> Buffer -> Buffer -> Result Buffer
tensor l a b = do
  let size = count a * count b
  (Buffer f xs, Buffer _ ys) <- preparePair l PolynomialWork size a b
  let r = V.length ys
  Buffer f <$> V.generateM (fromInteger size)
    (\j -> S.mulC l f (xs V.! (j `quot` r)) (ys V.! (j `rem` r)))
