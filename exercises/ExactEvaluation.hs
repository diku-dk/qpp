-- Run with: cabal run ExactEvaluation
module Main where

import Data.Ratio ((%))
import QPP.Syntax (QOp(..))
import QPP.Exact.Types
import qualified QPP.Exact.Scalar as S
import qualified QPP.Exact.State as St
import qualified QPP.Semantics.Exact.Matrix as M
import qualified QPP.Semantics.Exact.Statevector as V

demo :: Either ExactError (Bool, S.Cyclo, [S.Cyclo])
demo = do
  let limits = defaultLimits
      bell = Compose (C X) (Tensor H (Id 1))
      op = Compose (Tensor (Id 1) (C (Phase (1%3)))) bell
  input <- St.ket limits [False,False]
  direct <- V.evalOp limits op
  output <- V.apply limits direct input
  matrix <- M.evalOp limits op
  reference <- M.apply limits matrix input
  agree <- St.equal limits output reference
  norm2 <- St.normSquared limits output
  unbell <- V.evalOp limits (Adjoint bell)
  interfered <- V.apply limits unbell output
  probabilities <- traverse (\a -> S.conjugate limits a >>= S.mul limits a)
    (St.amplitudes interfered)
  pure (agree,norm2,probabilities)

main :: IO ()
main = case demo of
  Left err -> ioError (userError (show err))
  Right (agree,norm2,probabilities) -> do
    putStrLn ("Matrix/direct agreement: " ++ show agree)
    putStrLn ("Squared norm: " ++ show (S.asRational norm2))
    putStrLn ("Probabilities after inverse Bell circuit: " ++ show (map S.asRational probabilities))
