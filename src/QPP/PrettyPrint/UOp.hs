-- | Human-readable UOp expressions, with the same infix notation as QOp.
-- Identities are printed as I[d], where d is a dimension, not a qubit count.
module QPP.PrettyPrint.UOp (showUOp, printUOp) where

import Data.List (intercalate)
import Data.Ratio (numerator, denominator)
import QPP.UOp.Syntax

showUOp :: UOp -> String
showUOp op = case op of
  R a t             -> "R(" ++ showUOp a ++ ", " ++ showAngle t ++ ")"
  Tensor xs         -> showProduct (showIdentity 1) " ⊗ " xs
  DirectSum xs      -> showProduct (showIdentity 0) " ⊕ " xs
  Compose d xs      -> showProduct (showIdentity d) " ∘ " xs
  PermuteTensor ds p -> "PermuteTensor " ++ show ds ++ " " ++ show p
  PermuteSum ds p    -> "PermuteSum " ++ show ds ++ " " ++ show p

showProduct :: String -> String -> [UOp] -> String
showProduct unit _ [] = unit
showProduct _ _ [a] = showUOp a
showProduct _ separator xs =
  "(" ++ intercalate separator (map showUOp xs) ++ ")"

showIdentity :: Dim -> String
showIdentity d = "I[" ++ show d ++ "]"

showAngle :: Rational -> String
showAngle t
  | denominator t == 1 = show (numerator t)
  | otherwise = show (numerator t) ++ "/" ++ show (denominator t)

printUOp :: UOp -> IO ()
printUOp = putStrLn . showUOp
