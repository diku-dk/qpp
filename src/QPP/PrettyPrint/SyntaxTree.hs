-- | The syntax trees that the tree printers draw, as text in
--   "QPP.PrettyPrint.Tree" and through Forest in "QPP.PrettyPrint.TikZ".
--
--   A tree has one node per constructor application. The node holds the term
--   rooted there, and a printer labels the node by that term's outermost
--   constructor. The children are the operands, in the order they are written.
--   The node also holds the term's size, a qubit count for QOp and a dimension
--   for UOp. Building the tree neither normalizes nor evaluates the operator.
module QPP.PrettyPrint.SyntaxTree
  ( SyntaxTree, qopTree, uopTree, qubitLabel
  ) where

import Data.Bits ((.&.))
import Data.Tree (Tree (Node))
import Math.NumberTheory.Logarithms (integerLog2)
import qualified QPP.Syntax as Q
import qualified QPP.UOp.Syntax as U

-- The sizes are computed from the children once. A printer that shows no
-- sizes leaves these lazy fields unevaluated.
type SyntaxTree term size = Tree (term, size)

size :: SyntaxTree term size -> size
size (Node (_, n) _) = n

-- | A QOp with the qubit count of every subterm. As with 'Q.n_qubits', the
--   counts assume well-formed input; this is not a type-check.
qopTree :: Q.QOp -> SyntaxTree Q.QOp Integer
qopTree op = case op of
  Q.Id n          -> leaf (toInteger n)
  Q.Phase _       -> leaf 0
  Q.X             -> leaf 1
  Q.Y             -> leaf 1
  Q.Z             -> leaf 1
  Q.H             -> leaf 1
  Q.SX            -> leaf 1
  Q.R a _         -> unary id a
  Q.C a           -> unary (+ 1) a
  Q.Adjoint a     -> unary id a
  Q.Permute p     -> leaf (toInteger (length p))
  Q.Tensor a b    -> binary (+) a b
  Q.DirectSum a b -> binary (\n _ -> n + 1) a b
  Q.Compose a b   -> binary const a b
  where
    leaf n = Node (op, n) []
    unary measure a =
      let child = qopTree a
      in Node (op, measure (size child)) [child]
    binary measure a b =
      let left = qopTree a
          right = qopTree b
      in Node (op, measure (size left) (size right)) [left, right]

-- | A UOp with the dimension of every subterm.
uopTree :: U.UOp -> SyntaxTree U.UOp U.Dim
uopTree op = case op of
  U.R a _            -> let child = uopTree a in Node (op, size child) [child]
  U.Tensor xs        -> branch product xs
  U.DirectSum xs     -> branch sum xs
  U.Compose d xs     -> branch (const d) xs
  U.PermuteTensor ds _ -> Node (op, product ds) []
  U.PermuteSum ds _    -> Node (op, sum ds) []
  where
    branch measure xs =
      let children = map uopTree xs
      in Node (op, measure (map size children)) children

-- | A dimension as a size annotation: the qubit count @n@ of dimension
--   @2^n@, and @d=...@ for any other dimension, zero included.
--   A positive integer is a power of two iff exactly one bit is set.
--   Use integer arithmetic: floating logarithms can round a nearby dimension.
qubitLabel :: U.Dim -> String
qubitLabel d
  | d > 0 && d .&. (d - 1) == 0 = show (integerLog2 (toInteger d))
  | otherwise = "d=" ++ show d
