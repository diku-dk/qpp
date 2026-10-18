-- | Read-only views of validated QOp trees, with zero rotations eliminated.
-- The checked constructor is private; clients cannot forge widths or bypass
-- the existing symbolic type checker.
module QPP.Exact.Checked
  ( CheckedOp, View(..), Gate(..), checkOp, qubits, dimension, fieldIndex
  , size, depth, view, toQOp, CheckedStep(..), checkStep, checkProgram, programFieldIndex
  ) where

import Control.Monad (foldM, unless)
import Data.Bifunctor (first)
import Data.Ratio (denominator)
import qualified Data.Set as Set
import QPP.Syntax (QOp, Step, Program, op_qubits)
import qualified QPP.Syntax as Q
import qualified QPP.TypeCheck as T
import QPP.Exact.Types
import qualified QPP.Exact.Internal.Scalar as S

data Gate = X | Y | Z | H | SX deriving (Eq, Show)
-- Constructor names follow QOp. The checked node carries the width of Id;
-- the five one-qubit gates share a Primitive case.
data View
  = Id | Phase Rational | Primitive Gate
  | R CheckedOp Rational | C CheckedOp | Permute [Int]
  | Tensor CheckedOp CheckedOp | DirectSum CheckedOp CheckedOp
  | Compose CheckedOp CheckedOp | Adjoint CheckedOp
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

-- | Recover the validated syntax. Angles are normalized and zero rotations
-- have already become identities; the discarded axes cannot be reintroduced.
toQOp :: CheckedOp -> QOp
toQOp op = case view op of
  Id -> Q.Id (qubits op)
  Phase q -> Q.Phase q
  Primitive X -> Q.X
  Primitive Y -> Q.Y
  Primitive Z -> Q.Z
  Primitive H -> Q.H
  Primitive SX -> Q.SX
  R a q -> Q.R (toQOp a) q
  C a -> Q.C (toQOp a)
  Permute ks -> Q.Permute ks
  Tensor a b -> Q.Tensor (toQOp a) (toQOp b)
  DirectSum a b -> Q.DirectSum (toQOp a) (toQOp b)
  Compose a b -> Q.Compose (toQOp a) (toQOp b)
  Adjoint a -> Q.Adjoint (toQOp a)

checkOp :: Limits -> QOp -> Either ExactError CheckedOp
checkOp l original = do
  S.validateLimits l
  -- Validate raw ratios before the type checker's rational phase calculations.
  op <- normalizeAngles original
  _ <- first InvalidOperator (T.typeCheck op)
  annotate op
  where
    normalizeAngles = \case
      Q.Phase q -> Q.Phase <$> S.checkRational l q
      Q.R u q -> Q.R <$> normalizeAngles u <*> S.checkRational l q
      Q.C u -> Q.C <$> normalizeAngles u
      Q.Tensor u t -> Q.Tensor <$> normalizeAngles u <*> normalizeAngles t
      Q.DirectSum u t -> Q.DirectSum <$> normalizeAngles u <*> normalizeAngles t
      Q.Compose u t -> Q.Compose <$> normalizeAngles u <*> normalizeAngles t
      Q.Adjoint u -> Q.Adjoint <$> normalizeAngles u
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
      Q.Id n -> leaf n Id
      Q.Phase q -> node 0 (lcm 8 (2 * denominator q)) [] (Phase q)
      Q.X -> leaf 1 (Primitive X)
      Q.Y -> leaf 1 (Primitive Y)
      Q.Z -> leaf 1 (Primitive Z)
      Q.H -> leaf 1 (Primitive H)
      Q.SX -> leaf 1 (Primitive SX)
      -- The complete term has passed typeCheck, including this axis. Its
      -- width is valid, but its angles and fields play no part in evaluation.
      Q.R u 0 -> leaf (op_qubits u) Id
      Q.R u q -> do
        a <- annotate u
        node (qubits a) (lcm 8 (4 * denominator q)) [a] (R a q)
      Q.C u -> do
        a <- annotate u
        node (qubits a + 1) 8 [a] (C a)
      Q.Permute ks -> leaf (length ks) (Permute ks)
      Q.Tensor u t -> pair (\a b -> qubits a + qubits b) Tensor u t
      Q.DirectSum u t -> pair (\a _ -> qubits a + 1) DirectSum u t
      Q.Compose u t -> pair (\a _ -> qubits a) Compose u t
      Q.Adjoint u -> do
        a <- annotate u
        node (qubits a) 8 [a] (Adjoint a)
    pair width constructor u t = do
      a <- annotate u
      b <- annotate t
      node (width a b) 8 [a,b] (constructor a b)

data CheckedStep
  = Unitary CheckedOp
  | Measure [Int]
  | Initialize [Int] [Bool]
  deriving (Show)

-- | Validate quantum widths and wire lists without consuming randomness.
-- Repeated measurements are allowed; initialization needs distinct targets.
checkStep :: Limits -> Int -> Step -> Result CheckedStep
checkStep l width step = do
  _ <- S.dimension l width
  case step of
    Q.Unitary op -> do
      c <- checkOp l op
      unless (qubits c == width) (Left (WidthMismatch (qubits c) width))
      pure (Unitary c)
    Q.Measure ks -> wires ks >> pure (Measure ks)
    Q.Initialize ks vs -> do
      wires ks
      unless (length ks == length vs)
        (Left (InvalidInput "Initialize needs one value per wire"))
      unless (Set.size (Set.fromList ks) == length ks)
        (Left (InvalidInput "Initialize wires must be distinct"))
      pure (Initialize ks vs)
  where
    wires ks = unless (all (\k -> k >= 0 && k < width) ks)
      (Left (InvalidInput "program wire out of bounds"))

-- | Preflight the entire program. Errors carry zero-based step indices.
checkProgram :: Limits -> Int -> Program -> Result [CheckedStep]
checkProgram l width program = do
  _ <- S.dimension l width
  traverse (\(j,step) -> first (AtStep j) (checkStep l width step)) (zip [0..] program)

programFieldIndex :: Limits -> Integer -> [CheckedStep] -> Result Integer
programFieldIndex l start ops = do
  S.validateLimits l
  unless (start > 0) (Left (InvalidInput "root index must be positive"))
  S.bound "root index" start (maxRootIndex l)
  foldM (\n (j,step) -> first (AtStep j) $ case step of
    Unitary c -> do
      idx <- S.jointIndex l n (fieldIndex c)
      _ <- S.fieldDegree l idx
      pure idx
    _ -> pure n)
    start (zip [0..] ops)
