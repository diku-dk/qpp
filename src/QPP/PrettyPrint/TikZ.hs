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

import Data.Bits ((.&.))
import Data.List (intercalate)
import Data.Ratio (numerator, denominator)
import Data.Tree (Tree (Node))
import Math.NumberTheory.Logarithms (integerLog2)
import qualified QPP.Syntax as Q
import qualified QPP.UOp.Syntax as U

-- Sizes are computed from the children once: qubit counts for QOp, dimensions
-- for UOp. The unannotated renderer leaves these lazy fields unevaluated.
type OperatorTree size = Tree (String, String, size)

node :: String -> String -> size -> [OperatorTree size] -> OperatorTree size
node style label n = Node (style, label, n)

leaf :: String -> String -> size -> OperatorTree size
leaf style label n = node style label n []

size :: OperatorTree size -> size
size (Node (_, _, n) _) = n

unary :: String -> String -> (size -> size) -> OperatorTree size -> OperatorTree size
unary style label measure a = node style label (measure (size a)) [a]

-- | QOp as an ordered syntax tree with @: n@ outside every node for its
-- qubit count. Identity subscripts remain matrix dimensions. As with
-- 'Q.n_qubits', widths assume well-formed input; this is not a type-check.
showTikZQOp :: Q.QOp -> String
showTikZQOp = render (annotation . show) . qopTree

-- | UOp as an ordered syntax tree with exact qubit counts when the dimension
-- is a power of two. Other dimensions, including zero, use @: d = ...@.
showTikZUOp :: U.UOp -> String
showTikZUOp = render (annotation . qubitLabel) . uopTree

-- | Omit the outside qubit-count annotations.
showTikZQOpWithoutQubits :: Q.QOp -> String
showTikZQOpWithoutQubits = render (const "") . qopTree

-- | Omit the outside qubit-count and dimension annotations.
showTikZUOpWithoutQubits :: U.UOp -> String
showTikZUOpWithoutQubits = render (const "") . uopTree

-- The two translations describe structure only. Rendering and styling are
-- shared, so adding a constructor does not require writing layout code.
qopTree :: Q.QOp -> OperatorTree Integer
qopTree op = case op of
  Q.Id n          -> leaf "identity" (identityLabel (qubitDimension n)) (toInteger n)
  Q.Phase t       -> leaf "phase" (application "\\mathrm{Phase}" (rational t)) 0
  Q.X             -> gate "X"
  Q.Y             -> gate "Y"
  Q.Z             -> gate "Z"
  Q.H             -> gate "H"
  Q.SX            -> gate "\\sqrt{X}"
  Q.R a t         -> unary "rotation" (application "R" (rational t)) id (qopTree a)
  Q.C a           -> unary "control" "\\mathrm{C}" (+1) (qopTree a)
  Q.Adjoint a     -> unary "adjoint" "\\dagger" id (qopTree a)
  Q.Permute p     -> leaf "permutation" (application "P" (list p)) (toInteger (length p))
  Q.Tensor a b    -> binary "tensor" "\\otimes" (+) a b
  Q.DirectSum a b -> binary "sum" "\\oplus" (\n _ -> n + 1) a b
  Q.Compose a b   -> binary "compose" "\\circ" const a b
  where
    gate label = leaf "gate" label 1
    binary style label measure a b =
      let left = qopTree a
          right = qopTree b
      in node style label (measure (size left) (size right)) [left, right]
    qubitDimension 0 = "1"
    qubitDimension 1 = "2"
    qubitDimension n = "2^{" ++ show n ++ "}"

uopTree :: U.UOp -> OperatorTree U.Dim
uopTree op = case op of
  U.R a t          -> unary "rotation" (application "R" (rational t)) id (uopTree a)
  U.Tensor []      -> leaf "tensor" ("\\otimes[\\,]=" ++ identityLabel "1") 1
  U.DirectSum []   -> leaf "sum" ("\\oplus[\\,]=" ++ identityLabel "0") 0
  U.Compose d []   -> leaf "identity" (identityLabel (show d)) d
  U.Tensor xs      -> branch "tensor" "\\otimes" product xs
  U.DirectSum xs   -> branch "sum" "\\oplus" sum xs
  U.Compose d xs   -> branch "compose" ("\\circ_{" ++ show d ++ "}") (const d) xs
  U.PermuteTensor ds p -> permutation "\\otimes" (product ds) ds p
  U.PermuteSum ds p    -> permutation "\\oplus" (sum ds) ds p
  where
    branch style label measure xs =
      let children = map uopTree xs
      in node style label (measure (map size children)) children
    permutation kind d ds p = leaf "permutation"
      (application ("P_{" ++ kind ++ "}") (list ds ++ ";\\," ++ list p)) d

-- A positive integer is a power of two iff exactly one bit is set.
-- Use integer arithmetic: floating logarithms can round a nearby dimension.
qubitLabel :: U.Dim -> String
qubitLabel d
  | d > 0 && d .&. (d - 1) == 0 = show (integerLog2 (toInteger d))
  | otherwise = "d=" ++ show d

annotation :: String -> String
annotation n = ", qpp annotation={" ++ n ++ "}"

-- Forest's bracket notation mirrors Data.Tree directly. Braces around each
-- label protect permutation brackets and commas from Forest's parser.
render :: (size -> String) -> OperatorTree size -> String
render annotate tree = "\\begin{forest}\nqpp tree\n" ++ draw "" tree ++ "\\end{forest}%\n"
  where
    draw indent (Node (style, label, n) children) =
      indent ++ "[{$" ++ label ++ "$}, node options+={qpp/" ++ style ++ "}" ++ annotate n ++ "\n" ++
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
