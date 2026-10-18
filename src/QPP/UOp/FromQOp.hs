-- | Translate a checked qubit operator into general unitary syntax.
module QPP.UOp.FromQOp (fromQOp) where

import qualified QPP.Syntax as Q
import QPP.TypeCheck (TypeError, typeCheck)
import qualified QPP.UOp.Syntax as U

-- | Validate every QOp subterm, then lower it without constructing matrices.
-- The result has dimension 2^n, where n is the input's qubit count.
-- Global phases are preserved, including inside controlled operators.
fromQOp :: Q.QOp -> Either TypeError U.UOp
fromQOp op = do
  _ <- typeCheck op
  pure (lower op)
  where
    lower q = case q of
      Q.Id n       -> U.identity (2 ^ n)
      Q.Phase t    -> U.phase t
      Q.X          -> x
      Q.Y          -> U.tensor [U.phase (1/2), U.compose 2 [x,z]]  -- Y = i XZ.
      Q.Z          -> z
      Q.H          -> U.tensor
        [ U.phase (1/2)
        , U.compose 2 [U.R z (1/2), U.R x (1/2), U.R z (1/2)]
        ]  -- The three rotations give -i H; retain the correcting phase.
      Q.SX         -> U.tensor [U.phase (1/4), U.R x (1/2)]
      Q.R a t      -> U.R (lower a) t
      Q.C a        -> let b = lower a in U.directSum [U.identity (U.dim b), b]
      Q.Permute p  -> U.PermuteTensor (replicate (length p) 2) p
      Q.Tensor a b -> U.tensor [lower a, lower b]
      Q.DirectSum a b -> U.directSum [lower a, lower b]
      Q.Compose a b -> let a' = lower a in U.compose (U.dim a') [a', lower b]
      Q.Adjoint a  -> U.dagger (lower a)

    x = U.PermuteSum [1,1] [1,0]
    z = U.directSum [U.identity 1, U.phase 1]
