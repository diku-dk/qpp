-- | Read-only views of validated QOp trees. The checked constructor is private;
-- clients cannot forge widths or bypass the existing symbolic type checker.
module QPP.Exact.Checked
  ( CheckedOp, View(..), Gate(..), checkOp, qubits, dimension, fieldIndex
  , size, depth, view, checkProgram, programFieldIndex
  ) where

import Control.Monad (foldM, unless)
import Data.Bifunctor (first)
import Data.Ratio (denominator)
import QPP.Syntax (QOp(..), Step(..), Program)
import qualified QPP.TypeCheck as T
import QPP.Exact.Types
import qualified QPP.Exact.Internal.Scalar as S

data Gate = GX | GY | GZ | GH | GSX deriving (Eq, Show)
data View
  = Identity | ScalarPhase Rational | Primitive Gate
  | Rotation CheckedOp Rational | Control CheckedOp | Permutation [Int]
  | TensorProduct CheckedOp CheckedOp | BlockSum CheckedOp CheckedOp
  | Composition CheckedOp CheckedOp | Dagger CheckedOp
  deriving (Show)

-- Keep record selectors private: exporting them would allow record updates
-- that forge metadata even though the constructor itself is hidden.
data CheckedOp = CheckedOp
  { checkedQubits :: !Int
  , checkedDimension :: !Int
  , checkedFieldIndex :: !Integer
  , checkedSize :: !Integer
  , checkedDepth :: !Int
  , checkedView :: View
  } deriving (Show)

qubits, dimension, depth :: CheckedOp -> Int
qubits = checkedQubits
dimension = checkedDimension
depth = checkedDepth

fieldIndex, size :: CheckedOp -> Integer
fieldIndex = checkedFieldIndex
size = checkedSize

view :: CheckedOp -> View
view = checkedView

checkOp :: Limits -> QOp -> Either ExactError CheckedOp
checkOp l original = do
  S.validateLimits l
  -- Validate raw ratios before the type checker's rational phase calculations.
  op <- normalizeAngles original
  _ <- first InvalidOperator (T.typeCheck op)
  annotate op
  where
    normalizeAngles = \case
      Phase q -> Phase <$> S.checkRational l q
      R u q -> R <$> normalizeAngles u <*> S.checkRational l q
      C u -> C <$> normalizeAngles u
      Tensor u t -> Tensor <$> normalizeAngles u <*> normalizeAngles t
      DirectSum u t -> DirectSum <$> normalizeAngles u <*> normalizeAngles t
      Compose u t -> Compose <$> normalizeAngles u <*> normalizeAngles t
      Adjoint u -> Adjoint <$> normalizeAngles u
      op -> pure op
    node n own children v = do
      d <- S.dimension l n
      idx <- foldM (S.jointIndex l) own (map fieldIndex children)
      S.bound "root index" idx (maxRootIndex l)
      pure CheckedOp
        { checkedQubits = n, checkedDimension = d, checkedFieldIndex = idx
        , checkedSize = 1 + sum (map size children)
        , checkedDepth = 1 + maximum (0 : map depth children)
        , checkedView = v
        }
    leaf n = node n 8 []
    annotate = \case
      Id n -> leaf n Identity
      Phase q -> node 0 (lcm 8 (2 * denominator q)) [] (ScalarPhase q)
      X -> leaf 1 (Primitive GX)
      Y -> leaf 1 (Primitive GY)
      Z -> leaf 1 (Primitive GZ)
      H -> leaf 1 (Primitive GH)
      SX -> leaf 1 (Primitive GSX)
      R u q -> do
        a <- annotate u
        node (qubits a) (lcm 8 (4 * denominator q)) [a] (Rotation a q)
      C u -> do
        a <- annotate u
        node (qubits a + 1) 8 [a] (Control a)
      Permute ks -> leaf (length ks) (Permutation ks)
      Tensor u t -> pair (\a b -> qubits a + qubits b) TensorProduct u t
      DirectSum u t -> pair (\a _ -> qubits a + 1) BlockSum u t
      Compose u t -> pair (\a _ -> qubits a) Composition u t
      Adjoint u -> do
        a <- annotate u
        node (qubits a) 8 [a] (Dagger a)
    pair width constructor u t = do
      a <- annotate u
      b <- annotate t
      node (width a b) 8 [a,b] (constructor a b)

-- | Preflight the entire program before executing any unitary step. Step indices
-- are zero-based. Even unsupported operations on an empty wire list are rejected.
checkProgram :: Limits -> Int -> Program -> Either ExactError [CheckedOp]
checkProgram l width program = do
  _ <- S.dimension l width
  traverse check (zip [0..] program)
  where
    check (j, Unitary op) = first (AtStep j) $ do
      c <- checkOp l op
      unless (qubits c == width) (Left (WidthMismatch (qubits c) width))
      pure c
    check (j, Measure _) = Left (UnsupportedStep j "Measure")
    check (j, Initialize _ _) = Left (UnsupportedStep j "Initialize")

programFieldIndex :: Limits -> Integer -> [CheckedOp] -> Either ExactError Integer
programFieldIndex l start ops = do
  S.validateLimits l
  unless (start > 0) (Left (InvalidInput "root index must be positive"))
  S.bound "root index" start (maxRootIndex l)
  foldM (\n (j,c) -> first (AtStep j) $ do
    idx <- S.jointIndex l n (fieldIndex c)
    _ <- S.fieldDegree l idx
    pure idx)
    start (zip [0..] ops)
