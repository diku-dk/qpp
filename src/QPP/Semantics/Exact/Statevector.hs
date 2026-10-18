-- | Direct exact QOp semantics: apply operators without building their matrices.
-- evalOp prepares scalar constants once; apply interprets the resulting plan
-- on a vector. The recursive statevector equations are in runPlan below.
-- Plan preparation and array algorithms live in Internal.Plan and Internal.Vector.
module QPP.Semantics.Exact.Statevector
  ( OpT, evalOp, apply, evalStep, evalProg, opQubits, opFieldIndex
  , evalOpWithLimits, applyWithLimits, evalStepWithLimits, evalProgWithLimits
  , module QPP.Exact.State
  ) where

import Control.Monad (unless)
import GHC.Stack (HasCallStack, withFrozenCallStack)
import QPP.Exact.Types (orFail)
import QPP.Syntax (QOp, Program, Step, RNG, Outcomes, HasQubits(..))
import QPP.Exact.State
import qualified QPP.Exact.Checked as C
import QPP.Exact.Internal.Scalar (Field, Coeff)
import qualified QPP.Exact.Internal.Scalar as S
import qualified QPP.Exact.Internal.State as St
import QPP.Exact.Internal.Buffer (Buffer(..))
import qualified QPP.Exact.Internal.Buffer as B
import qualified QPP.Exact.Internal.Evaluation as E
import qualified QPP.Exact.Internal.Plan as P
import qualified QPP.Exact.Internal.Vector as Vector

-- Keep only the compiled plan and its execution requirements, not a second
-- copy of the source tree. Selectors stay private to prevent record updates.
data OpT = OpT
  { operatorQubits :: !Int
  , operatorRequirements :: !P.Requirements
  , operatorField :: !Field
  , operatorPlan :: !(P.Plan Coeff)
  } deriving (Show)

instance HasQubits OpT where n_qubits = opQubits

opQubits :: OpT -> Int
opQubits = operatorQubits

opFieldIndex :: OpT -> Integer
opFieldIndex = S.index . operatorField

-- | Validate the complete term, select a field and prepare a reusable operator.
evalOp :: HasCallStack => QOp -> OpT
evalOp op = withFrozenCallStack $ orFail "evalOp" (evalOpWithLimits defaultLimits op)

evalOpWithLimits :: Limits -> QOp -> Result OpT
evalOpWithLimits limits op = do
  checked <- C.checkOp limits op
  field <- S.makeField limits (C.fieldIndex checked)
  let required = P.requirements checked
  P.budget limits field (C.dimension checked) required
  plan <- P.compile limits field checked
  pure OpT
    { operatorQubits = C.qubits checked
    , operatorRequirements = required
    , operatorField = field
    , operatorPlan = plan
    }

-- Interpret each prepared operator as a function on vectors. Plan actions
-- mirror QOp: rational angles have become constants, and adjoints have been
-- pushed to the leaves. All coefficients already belong to the same field.
runPlan :: Limits -> Field -> P.Plan Coeff -> Vector.Transform
runPlan limits field = eval
  where
    eval :: P.Plan Coeff -> Vector.Transform
    eval (P.Plan _ action) = case action of
      P.Id -> pure
      P.Scalar c -> Vector.scale limits field c
      P.Gate g   -> Vector.applyGate limits field g

      -- (cI + dA)v = cv + d(Av). Keep v while recursively computing Av.
      P.R axis c d -> \v -> do
        av <- eval axis v
        Vector.linearCombination limits field c v d av

      -- Control is the direct sum I ⊕ A, selected by the leading qubit.
      P.C a           -> pure <+> eval a
      P.DirectSum a b -> eval a <+> eval b
      P.Compose a b   -> eval a ∘ eval b

      -- Apply B to rows, then A to columns: the array rule is V -> A V B^T.
      P.Tensor a b -> Vector.tensor (P.dimension a) (P.dimension b) (eval a) (eval b)
      P.Permute gather -> pure . Vector.permute gather

    -- The matrix interpreter's notation, now acting on vector transformations.
    (<+>) = Vector.directSum
    (∘) a b v = b v >>= a

    infixr 7 <+>
    infixr 6 ∘

-- | Apply a prepared operator, embedding its constants and the input into a
-- common field if needed. Input vectors may be unnormalized or entangled.
apply :: HasCallStack => OpT -> StateT -> StateT
apply op state = withFrozenCallStack $ orFail "apply" (applyWithLimits defaultLimits op state)

applyWithLimits :: Limits -> OpT -> StateT -> Result StateT
applyWithLimits limits op state@(St.StateT n input) = do
  unless (opQubits op == n) (Left (WidthMismatch (opQubits op) n))
  let sourceField = operatorField op
      plan = operatorPlan op
  P.validate limits sourceField plan
  St.validate limits state
  field <- S.commonField limits sourceField (B.field input)
  P.budget limits field (P.dimension plan) (operatorRequirements op)
  promotedPlan <- P.promote limits field sourceField plan
  promotedState <- St.StateT n <$> B.promote limits field input
  applyIn limits field promotedPlan promotedState

-- The plan and input already share a field; preparation belongs to the caller.
applyIn :: Limits -> Field -> P.Plan Coeff -> StateT -> Result StateT
applyIn limits field plan (St.StateT n (Buffer _ input)) = do
  output <- runPlan limits field plan input
  pure (St.StateT n (Buffer field output))

-- | Check every step before execution and use one field for the whole program.
evalProg :: HasCallStack => Program -> StateT -> RNG -> (StateT, Outcomes, RNG)
evalProg program state rng = withFrozenCallStack $
  orFail "evalProg" (evalProgWithLimits defaultLimits program state rng)

evalProgWithLimits :: Limits -> Program -> StateT -> RNG -> Result (StateT, Outcomes, RNG)
evalProgWithLimits limits =
  E.evalProgWith limits (planBudget limits) (applyChecked limits)

evalStep :: HasCallStack => (StateT, Outcomes, RNG) -> Step -> (StateT, Outcomes, RNG)
evalStep input step = withFrozenCallStack $
  orFail "evalStep" (evalStepWithLimits defaultLimits input step)

evalStepWithLimits :: Limits -> (StateT, Outcomes, RNG) -> Step -> Result (StateT, Outcomes, RNG)
evalStepWithLimits limits =
  E.evalStepWith limits (planBudget limits) (applyChecked limits)

applyChecked :: Limits -> Field -> C.CheckedOp -> StateT -> Result StateT
applyChecked limits field op state = do
  plan <- P.compile limits field op
  applyIn limits field plan state

planBudget :: Limits -> Field -> C.CheckedOp -> Result ()
planBudget limits field op =
  P.budget limits field (C.dimension op) (P.requirements op)
