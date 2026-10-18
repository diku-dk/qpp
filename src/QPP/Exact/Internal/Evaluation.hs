-- | Shared checked program interpreter. Backends supply only the unitary
-- workspace bound and action; measurement and reset use the common state type.
module QPP.Exact.Internal.Evaluation (evalProgWith, evalStepWith) where

import Control.Monad (foldM)
import Data.Bifunctor (first)
import QPP.Exact.Types
import QPP.Syntax (Program, Step, Outcomes, RNG)
import qualified QPP.Exact.Checked as C
import QPP.Exact.Internal.Scalar (Field)
import qualified QPP.Exact.Internal.Scalar as S
import QPP.Exact.Internal.State (StateT(..))
import qualified QPP.Exact.Internal.State as St
import qualified QPP.Exact.Internal.Buffer as B
import qualified QPP.Exact.Internal.Measurement as Measurement

type Run = (StateT, Outcomes, RNG)
type Reserve = Field -> C.CheckedOp -> Result ()
type ApplyUnitary = Field -> C.CheckedOp -> StateT -> Result StateT

evalProgWith :: Limits -> Reserve -> ApplyUnitary -> Program -> StateT -> RNG -> Result Run
evalProgWith limits reserve applyUnitary program state@(StateT n input) rng = do
  St.validate limits state
  (field,steps) <- prepareProgram limits n (B.field input) program
  mapM_ (\(j,step) -> first (AtStep j) (reserveStep field step)) (zip [0..] steps)
  start <- StateT n <$> B.promote limits field input
  foldM (\acc (j,step) -> first (AtStep j) (runStep field acc step))
    (start, [], rng) (zip [0..] steps)
  where
    reserveStep field = \case
      C.Unitary op -> reserve field op
      C.Measure [] -> pure ()
      C.Initialize [] [] -> pure ()
      _ -> Measurement.budget limits field n

    runStep field (vector, outcomes, randoms) = \case
      C.Unitary op -> do
        output <- applyUnitary field op vector
        pure (output, outcomes, randoms)
      C.Measure ks -> measure (vector, outcomes, randoms) ks
      C.Initialize ks targets -> do
        (projected, observed, remaining) <- measure (vector, [], randoms) ks
        corrected <- Measurement.correctBits limits (zip ks (zipWith (/=) observed targets)) projected
        pure (corrected, outcomes, remaining)

    -- Match QPP.Semantics: reverse the wires, prepend each result. This leaves
    -- the new outcomes in the supplied wire order, ahead of earlier outcomes.
    measure acc ks = foldM (Measurement.measure1WithLimits limits) acc (reverse ks)

-- | A single step uses the same validation, field selection and error context
-- as a one-step program, while preserving previously recorded outcomes.
evalStepWith :: Limits -> Reserve -> ApplyUnitary -> Run -> Step -> Result Run
evalStepWith limits reserve applyUnitary (state, outcomes, rng) step = do
  (output, observed, remaining) <- evalProgWith limits reserve applyUnitary [step] state rng
  pure (output, observed ++ outcomes, remaining)

-- The input has already been validated at the public boundary. Keep its field
-- for an empty program; otherwise construct the common context only once.
prepareProgram :: Limits -> Int -> Field -> Program -> Result (Field, [C.CheckedStep])
prepareProgram l width inputField program = do
  ops <- C.checkProgram l width program
  idx <- C.programFieldIndex l (S.index inputField) ops
  f <- S.fieldFor l idx [inputField]
  pure (f,ops)
