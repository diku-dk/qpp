-- | Dense square-matrix algebra in an established cyclotomic field.
-- The evaluator owns validation, field promotion and workspace reservations;
-- these kernels own the row-major representation and matrix algorithms.
-- Dimensions are positive and compatible, and every entry belongs to the
-- caller's field. The constructor is available only inside the library.
module QPP.Exact.Internal.Matrix
  ( Matrix(..), at, toRows
  , identity, diagonal, scalar, permutation
  , add, scale, multiply, tensor, directSum, adjoint, apply
  ) where

import Control.Monad (foldM)
import qualified Data.Vector as V
import QPP.Exact.Types
import QPP.Exact.Internal.Scalar (Field, Coeff)
import qualified QPP.Exact.Internal.Scalar as S
import QPP.Exact.Internal.Gates (permuteIndex)

data Matrix = Matrix
  { dimension :: !Int
  , entries :: !(V.Vector Coeff)
  }

at :: Matrix -> Int -> Int -> Coeff
at (Matrix dim xs) row column = xs V.! (row * dim + column)

toRows :: Matrix -> [[Coeff]]
toRows mat = [[at mat row column | column <- indices] | row <- indices]
  where indices = [0 .. dimension mat - 1]

generate :: Int -> (Int -> Int -> Coeff) -> Matrix
generate dim entry = Matrix dim $ V.generate (dim * dim) $ \index ->
  let (row,column) = index `quotRem` dim
  in entry row column

generateM :: Int -> (Int -> Int -> Result Coeff) -> Result Matrix
generateM dim entry = Matrix dim <$> V.generateM (dim * dim) (\index ->
  let (row,column) = index `quotRem` dim
  in entry row column)

identity :: Field -> Int -> Matrix
identity field dim = diagonal field dim (S.one field)

diagonal :: Field -> Int -> Coeff -> Matrix
diagonal field dim value = generate dim $ \row column ->
  if row == column then value else S.zero field

scalar :: Coeff -> Matrix
scalar value = Matrix 1 (V.singleton value)

permutation :: Field -> [Int] -> Matrix
permutation field wires = generate (2 ^ length wires) $ \row column ->
  if row == permuteIndex wires column then S.one field else S.zero field

add :: Limits -> Matrix -> Matrix -> Result Matrix
add limits a b = Matrix (dimension a) <$>
  V.zipWithM (S.addC limits) (entries a) (entries b)

scale :: Limits -> Field -> Coeff -> Matrix -> Result Matrix
scale limits field value mat = Matrix (dimension mat) <$>
  V.mapM (S.mulC limits field value) (entries mat)

-- (AB)[i,j] = sum_k A[i,k] B[k,j].
multiply :: Limits -> Field -> Matrix -> Matrix -> Result Matrix
multiply limits field a b = generateM (dimension a) $ \row column ->
  dot limits field (dimension a) (at a row) (\k -> at b k column)

-- (A tensor B)[i*r+k,j*r+l] = A[i,j] B[k,l], where r = dimension B.
tensor :: Limits -> Field -> Matrix -> Matrix -> Result Matrix
tensor limits field a b = generateM (dimension a * dimension b) $ \row column ->
  let (rowA,rowB) = row `quotRem` dimension b
      (columnA,columnB) = column `quotRem` dimension b
  in S.mulC limits field (at a rowA columnA) (at b rowB columnB)

-- The off-diagonal blocks are zero.
directSum :: Field -> Matrix -> Matrix -> Matrix
directSum field a b = generate (dimension a + dimension b) $ \row column ->
  let split = dimension a
  in if row < split && column < split
       then at a row column
     else if row >= split && column >= split
       then at b (row-split) (column-split)
     else S.zero field

adjoint :: Limits -> Field -> Matrix -> Result Matrix
adjoint limits field mat = generateM (dimension mat) $ \row column ->
  S.conjugateC limits field (at mat column row)

apply :: Limits -> Field -> Matrix -> V.Vector Coeff -> Result (V.Vector Coeff)
apply limits field mat vector = V.generateM (dimension mat) $ \row ->
  dot limits field (dimension mat) (at mat row) (vector V.!)

dot :: Limits -> Field -> Int -> (Int -> Coeff) -> (Int -> Coeff) -> Result Coeff
dot limits field dim a b = foldM accumulate (S.zero field) [0 .. dim-1]
  where
    accumulate total k = do
      product' <- S.mulC limits field (a k) (b k)
      S.addC limits total product'
