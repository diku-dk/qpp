-- | Exact dense QOp semantics. Entries include global phase. OpT can also hold
-- general square matrices via fromRows/linear combinations; only evalOp's
-- checked QOp construction guarantees unitarity.
module QPP.Semantics.Exact.Matrix
  ( OpT, evalOp, apply, evalStep, evalProg
  , evalOpWithLimits, applyWithLimits, evalStepWithLimits, evalProgWithLimits
  , opQubits, opFieldIndex, entry, entries, fromRows
  , addOp, scaleOp, innerOp, normSquaredOp, equalOp, adjointOp, composeOp
  , module QPP.Exact.State
  ) where

import qualified Control.Applicative as A
import Control.Monad (unless)
import GHC.Stack (HasCallStack, withFrozenCallStack)
import QPP.Exact.Types (orFail)
import QPP.Syntax (QOp(..), Program, Step, RNG, Outcomes, op_qubits, HasQubits(..))
import QPP.Exact.State
import qualified QPP.Exact.Checked as C
import QPP.Exact.Internal.Scalar (Field, Cyclo(..))
import qualified QPP.Exact.Internal.Scalar as S
import qualified QPP.Exact.Internal.State as St
import QPP.Exact.Internal.Buffer (Buffer(..))
import qualified QPP.Exact.Internal.Buffer as B
import qualified QPP.Exact.Internal.Evaluation as E
import QPP.Exact.Internal.Gates (primitive)
import qualified QPP.Exact.Internal.Matrix as M

data OpT = OpT !Int !Buffer deriving (Show)

instance HasQubits OpT where n_qubits = opQubits

-- | Check the complete term and choose one field before interpreting it.
evalOp :: HasCallStack => QOp -> OpT
evalOp op = withFrozenCallStack $ orFail "evalOp" (evalOpWithLimits defaultLimits op)

evalOpWithLimits :: Limits -> QOp -> Result OpT
evalOpWithLimits limits op = do
  checked <- C.checkOp limits op
  field <- S.makeField limits (C.fieldIndex checked)
  matrixBudget limits field checked
  evaluate limits field checked

-- Interpret the validated QOp directly, using the numerical backend's notation.
-- The local matrix operations propagate Result failures in one common field;
-- validation, promotion and storage layout stay outside the semantic clauses.
evaluate :: Limits -> Field -> C.CheckedOp -> Result OpT
evaluate limits field checked = pack (C.qubits checked) field <$> eval (C.toQOp checked)
  where
    eval op = case op of
      Id n    -> ident (2^n)
      Phase q -> phase q

      X  -> gate C.X
      Y  -> gate C.Y
      Z  -> gate C.Z
      H  -> gate C.H
      SX -> gate C.SX

      -- R_A(q) = c I + d A, where c = cos(pi*q/2), d = -i*sin(pi*q/2).
      R axis q -> do
        (c,d) <- S.rotationC limits field q
        if S.isZeroC d
          then c .* ident (2 ^ op_qubits axis) -- no axis evaluation for scalar rotations
          else do
            a <- eval axis
            (c .* ident (M.dimension a)) .+ (d .* pure a)

      Permute ks -> pure (M.permutation field ks)

      C a -> do
        mat <- eval a
        ident (M.dimension mat) <+> pure mat

      Tensor a b    -> eval a ⊗ eval b
      Compose a b   -> eval a ∘ eval b
      DirectSum a b -> eval a <+> eval b
      Adjoint a     -> adj (eval a)

    ident dim = pure (M.identity field dim)
    phase q = M.scalar <$> S.phaseC limits field q
    gate g = M.Matrix 2 <$> primitive limits field g

    -- The same matrix notation as Semantics.Matrix, lifted through explicit
    -- arithmetic failures. These are local operations, not partial instances.
    (⊗) = binary (M.tensor limits field)
    (∘) = binary (M.multiply limits field)
    (<+>) = A.liftA2 (M.directSum field)
    (.+) = binary (M.add limits)
    value .* matrix = matrix >>= M.scale limits field value
    adj matrix = matrix >>= M.adjoint limits field

    infixr 8 ⊗, .*
    infixr 7 <+>, .+
    infixr 6 ∘

    binary operation left right = do
      a <- left
      b <- right
      operation a b

apply :: HasCallStack => OpT -> StateT -> StateT
apply op state = withFrozenCallStack $ orFail "apply" (applyWithLimits defaultLimits op state)

applyWithLimits :: Limits -> OpT -> StateT -> Result StateT
applyWithLimits limits (OpT n matrix) state@(St.StateT _ input) = do
  unless (n == qubits state) (Left (WidthMismatch n (qubits state)))
  (mat,vector) <- B.preparePair limits S.PolynomialWork (B.count input) matrix input
  applyIn limits (OpT n mat) (St.StateT n vector)

-- Inputs already share a field; preparation and checks belong to the caller.
applyIn :: Limits -> OpT -> StateT -> Result StateT
applyIn limits op@(OpT n (Buffer field _)) (St.StateT _ (Buffer _ vector)) = do
  result <- M.apply limits field (dense op) vector
  pure (St.StateT n (Buffer field result))

evalProg :: HasCallStack => Program -> StateT -> RNG -> (StateT, Outcomes, RNG)
evalProg program state rng = withFrozenCallStack $
  orFail "evalProg" (evalProgWithLimits defaultLimits program state rng)

evalProgWithLimits :: Limits -> Program -> StateT -> RNG -> Result (StateT, Outcomes, RNG)
evalProgWithLimits limits =
  E.evalProgWith limits (matrixBudget limits) (applyChecked limits)

evalStep :: HasCallStack => (StateT, Outcomes, RNG) -> Step -> (StateT, Outcomes, RNG)
evalStep input step = withFrozenCallStack $
  orFail "evalStep" (evalStepWithLimits defaultLimits input step)

evalStepWithLimits :: Limits -> (StateT, Outcomes, RNG) -> Step -> Result (StateT, Outcomes, RNG)
evalStepWithLimits limits =
  E.evalStepWith limits (matrixBudget limits) (applyChecked limits)

applyChecked :: Limits -> Field -> C.CheckedOp -> StateT -> Result StateT
applyChecked limits field op state = do
  matrix <- evaluate limits field op
  applyIn limits matrix state

-- Checked linear algebra also accepts nonunitary matrices imported by fromRows.
-- Buffer operations handle common-field promotion for addition and scaling.
addOp :: Limits -> OpT -> OpT -> Result OpT
addOp limits a@(OpT n xs) b@(OpT _ ys) = do
  sameWidth a b
  OpT n <$> B.add limits xs ys

scaleOp :: Limits -> Cyclo -> OpT -> Result OpT
scaleOp limits value (OpT n buffer) = OpT n <$> B.scale limits value buffer

composeOp :: Limits -> OpT -> OpT -> Result OpT
composeOp limits a@(OpT n xs) b@(OpT _ ys) = do
  sameWidth a b
  (left,right) <- B.preparePair limits S.PolynomialWork (B.count xs) xs ys
  let field = B.field left
  pack n field <$> M.multiply limits field (dense (OpT n left)) (dense (OpT n right))

adjointOp :: Limits -> OpT -> Result OpT
adjointOp limits op@(OpT n buffer) = do
  B.validate limits buffer
  let field = B.field buffer
  S.workspace limits field (2 * B.count buffer * toInteger (S.degree field))
  pack n field <$> M.adjoint limits field (dense op)

innerOp :: Limits -> OpT -> OpT -> Result Cyclo
innerOp limits a@(OpT _ xs) b@(OpT _ ys) = do
  sameWidth a b
  B.inner limits xs ys

normSquaredOp :: Limits -> OpT -> Result Cyclo
normSquaredOp limits a = innerOp limits a a

equalOp :: Limits -> OpT -> OpT -> Result Bool
equalOp limits (OpT _ a) (OpT _ b) = B.equal limits a b

-- Import and observation retain the common field and the quantum shape.
fromRows :: Limits -> [[Cyclo]] -> Result OpT
fromRows limits rows = do
  n <- S.widthOfLength (length rows)
  dim <- S.dimension limits n
  S.slots limits (toInteger dim * toInteger dim)
  unless (all ((== dim) . length) rows)
    (Left (InvalidInput "matrix must be square with power-of-two dimensions"))
  OpT n <$> B.fromScalars limits (concat rows)

opQubits :: OpT -> Int
opQubits (OpT n _) = n

opFieldIndex :: OpT -> Integer
opFieldIndex (OpT _ buffer) = S.index (B.field buffer)

-- | Complete exact rows, with field context attached to each observed scalar.
entries :: OpT -> [[Cyclo]]
entries op@(OpT _ buffer) = map (map (Cyclo (B.field buffer))) (M.toRows (dense op))

entry :: OpT -> Int -> Int -> Result Cyclo
entry op@(OpT _ buffer) row column
  | row >= 0 && row < dim && column >= 0 && column < dim =
      pure (Cyclo (B.field buffer) (M.at matrix row column))
  | otherwise = Left (InvalidInput "matrix index out of bounds")
  where
    matrix = dense op
    dim = M.dimension matrix

-- Views of the same coefficient array; neither conversion copies its entries.
dense :: OpT -> M.Matrix
dense (OpT n buffer) = M.Matrix (2^n) (B.values buffer)

pack :: Int -> Field -> M.Matrix -> OpT
pack n field matrix = OpT n (Buffer field (M.entries matrix))

sameWidth :: OpT -> OpT -> Result ()
sameWidth a b = unless (opQubits a == opQubits b)
  (Left (WidthMismatch (opQubits a) (opQubits b)))

-- At each recursion level retain up to two child matrices and one result;
-- eight further buffers cover input promotion, application and block assembly.
matrixBudget :: Limits -> Field -> C.CheckedOp -> Result ()
matrixBudget limits field op = S.workspace limits field $
  (3 * toInteger (C.depth op) + 8) * dim * dim * toInteger (S.degree field)
  where dim = toInteger (C.dimension op)
