-- Run with: cabal run ExactEvaluation
module Main where

import Data.Ratio ((%))
import QPP.Syntax (QOp(..), Step(..), Outcomes)
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
      input = St.ket [0,0]
      output = V.apply (V.evalOp op) input
      reference = M.apply (M.evalOp op) input
  agree <- St.equal limits output reference
  norm2 <- St.normSquared limits output
  let interfered = V.apply (V.evalOp (Adjoint bell)) output
  probabilities <- traverse (\a -> S.conjugate limits a >>= S.mul limits a)
    (St.amplitudes interfered)
  pure (agree,norm2,probabilities)

-- The stored branch has squared norm 1/2; conditional probability is still 1.
measurementDemo :: Result (Outcomes, S.Cyclo, S.Cyclo)
measurementDemo = do
  let bell = Compose (C X) (Tensor H (Id 1))
      program = [Unitary bell, Measure [0], Initialize [1] [False]]
      (state,outcomes,_) = V.evalProg program (St.ket [0,0]) [0.8,0.2]
  weight <- St.normSquared defaultLimits state
  probability <- St.probability defaultLimits state 0 True
  pure (outcomes,weight,probability)

main :: IO ()
main = case (,) <$> demo <*> measurementDemo of
  Left err -> ioError (userError (show err))
  Right ((agree,norm2,probabilities),(outcomes,weight,probability)) -> do
    putStrLn ("Matrix/direct agreement: " ++ show agree)
    putStrLn ("Squared norm: " ++ show (S.asRational norm2))
    putStrLn ("Probabilities after inverse Bell circuit: " ++ show (map S.asRational probabilities))
    putStrLn ("Recorded outcomes (reset samples hidden): " ++ show outcomes)
    putStrLn ("Stored branch squared norm: " ++ show (S.asRational weight))
    putStrLn ("Probability of measuring wire 0 as one again: " ++ show (S.asRational probability))
