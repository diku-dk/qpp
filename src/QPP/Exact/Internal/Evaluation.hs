-- | Shared whole-program preflight. Execution stays backend-specific.
module QPP.Exact.Internal.Evaluation (prepareProgram) where

import QPP.Exact.Types
import QPP.Syntax (Program)
import qualified QPP.Exact.Checked as C
import QPP.Exact.Internal.Scalar (Field)
import qualified QPP.Exact.Internal.Scalar as S

-- The input has already been validated at the public boundary. Keep its field
-- for an empty program; otherwise construct the common context only once.
prepareProgram :: Limits -> Int -> Field -> Program -> Result (Field, [C.CheckedOp])
prepareProgram l width inputField program = do
  ops <- C.checkProgram l width program
  idx <- C.programFieldIndex l (S.index inputField) ops
  f <- S.fieldFor l idx [inputField]
  pure (f,ops)
