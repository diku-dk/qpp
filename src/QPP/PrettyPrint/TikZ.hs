-- | Operator syntax trees, laid out by Forest (a TikZ package).
--
-- Load @latex/qpp-tikz.sty@ in the surrounding LaTeX document. The output is
-- a fragment suitable for @\input@. Children retain their syntactic order;
-- the rightmost operand of a composition acts first.
-- Printing neither normalizes nor evaluates the operator.
module QPP.PrettyPrint.TikZ
  ( showTikZQOp, showTikZUOp
  , showTikZQOpWithoutQubits, showTikZUOpWithoutQubits
  ) where

import Data.List (intercalate)
import Data.Ratio (numerator, denominator)
import Data.Tree (Tree (Node))
import QPP.PrettyPrint.SyntaxTree
import qualified QPP.Syntax as Q
import qualified QPP.UOp.Syntax as U

-- | QOp as an ordered syntax tree with @: n@ outside every node for its
-- qubit count. Identity subscripts remain matrix dimensions. As with
-- 'Q.n_qubits', widths assume well-formed input; this is not a type-check.
showTikZQOp :: Q.QOp -> String
showTikZQOp = render qopNode (annotation . show) . qopTree

-- | UOp as an ordered syntax tree with exact qubit counts when the dimension
-- is a power of two. Other dimensions, including zero, use @: d = ...@.
showTikZUOp :: U.UOp -> String
showTikZUOp = render uopNode (annotation . qubitLabel) . uopTree

-- | Omit the outside qubit-count annotations.
showTikZQOpWithoutQubits :: Q.QOp -> String
showTikZQOpWithoutQubits = render qopNode (const "") . qopTree

-- | Omit the outside qubit-count and dimension annotations.
showTikZUOpWithoutQubits :: U.UOp -> String
showTikZUOpWithoutQubits = render uopNode (const "") . uopTree

-- The trees come from "QPP.PrettyPrint.SyntaxTree", which describes structure
-- only. These two functions give a node's style and math label from the
-- outermost constructor of its term, so adding a constructor does not require
-- writing layout code.
qopNode :: Q.QOp -> (String, String)
qopNode op = case op of
  Q.Id n          -> ("identity", identityLabel (qubitDimension n))
  Q.Phase t       -> ("phase", application "\\mathrm{Phase}" (rational t))
  Q.X             -> gate "X"
  Q.Y             -> gate "Y"
  Q.Z             -> gate "Z"
  Q.H             -> gate "H"
  Q.SX            -> gate "\\sqrt{X}"
  Q.R _ t         -> ("rotation", application "R" (rational t))
  Q.C _           -> ("control", "\\mathrm{C}")
  Q.Adjoint _     -> ("adjoint", "\\dagger")
  Q.Permute p     -> ("permutation", application "P" (list p))
  Q.Tensor _ _    -> ("tensor", "\\otimes")
  Q.DirectSum _ _ -> ("sum", "\\oplus")
  Q.Compose _ _   -> ("compose", "\\circ")
  where
    gate label = ("gate", label)
    qubitDimension :: Q.Nat -> String
    qubitDimension 0 = "1"
    qubitDimension 1 = "2"
    qubitDimension n = "2^{" ++ show n ++ "}"

uopNode :: U.UOp -> (String, String)
uopNode op = case op of
  U.R _ t          -> ("rotation", application "R" (rational t))
  U.Tensor []      -> ("tensor", "\\otimes[\\,]=" ++ identityLabel "1")
  U.DirectSum []   -> ("sum", "\\oplus[\\,]=" ++ identityLabel "0")
  U.Compose d []   -> ("identity", identityLabel (show d))
  U.Tensor _       -> ("tensor", "\\otimes")
  U.DirectSum _    -> ("sum", "\\oplus")
  U.Compose d _    -> ("compose", "\\circ_{" ++ show d ++ "}")
  U.PermuteTensor ds p -> permutation "\\otimes" ds p
  U.PermuteSum ds p    -> permutation "\\oplus" ds p
  where
    permutation kind ds p =
      ("permutation", application ("P_{" ++ kind ++ "}") (list ds ++ ";\\," ++ list p))

annotation :: String -> String
annotation n = ", qpp annotation={" ++ n ++ "}"

-- Forest's bracket notation mirrors Data.Tree directly. Braces around each
-- label protect permutation brackets and commas from Forest's parser.
render :: (term -> (String, String)) -> (size -> String) -> SyntaxTree term size -> String
render describe annotate tree = "\\begin{forest}\nqpp tree\n" ++ draw "" tree ++ "\\end{forest}%\n"
  where
    draw indent (Node (term, n) children) =
      let (style, label) = describe term
      in indent ++ "[{$" ++ label ++ "$}, node options+={qpp/" ++ style ++ "}" ++ annotate n ++ "\n" ++
         concatMap (draw (indent ++ "  ")) children ++ indent ++ "]\n"

command :: String -> [String] -> String
command name args = "\\" ++ name ++ concatMap (\a -> "{" ++ a ++ "}") args

identityLabel :: String -> String
identityLabel d = "I_{" ++ d ++ "}"

application :: String -> String -> String
application f a = f ++ "\\left(" ++ a ++ "\\right)"

list :: Show a => [a] -> String
list xs = "[" ++ intercalate "," (map show xs) ++ "]"

rational :: Rational -> String
rational t
  | denominator t == 1 = show (numerator t)
  | t < 0 = "-" ++ rational (-t)
  | otherwise = command "frac" [show (numerator t), show (denominator t)]
