-- | Operator syntax trees drawn as text, for the terminal.
--
--   These are the trees that "QPP.PrettyPrint.TikZ" typesets: one node per
--   constructor application, with its operands below it from left to right,
--   so the rightmost operand of a composition acts first. A rotation node
--   shows its angle and has its axis as its child.
--
--   > ghci> printOpTree (C X ∘ (H ⊗ I))
--   >    ∘
--   > ┌──┴──┐
--   > C     ⊗
--   > │   ┌─┴─┐
--   > X   H   I
--
--   Nodes are spelled as 'showOp' and 'showUOp' spell them, with angles as
--   they are written in the syntax, in units of π. Printing neither normalizes
--   nor evaluates the operator. A drawing is as wide as its leaves put side by
--   side, so a large term needs a wide terminal.
module QPP.PrettyPrint.Tree
  ( showOpTree, printOpTree
  , showOpTreeWithQubits, printOpTreeWithQubits
  , showUOpTree, printUOpTree
  , showUOpTreeWithQubits, printUOpTreeWithQubits
  ) where

import Data.Char (isSpace)
import Data.List (dropWhileEnd, intercalate, transpose)
import Data.Ratio (denominator, numerator)
import Data.Tree (Tree (Node))
import QPP.PrettyPrint.Op (showOp)
import QPP.PrettyPrint.SyntaxTree
import QPP.PrettyPrint.UOp (showUOp)
import qualified QPP.Syntax as Q
import qualified QPP.UOp.Syntax as U

-- | The syntax tree of a QOp, root at the top.
showOpTree :: Q.QOp -> String
showOpTree = render qopLabel (const "") . qopTree

-- | The same tree with @: n@ after every node for the qubit count of the term
--   rooted there. As with 'Q.n_qubits', the counts assume well-formed input;
--   this is not a type-check.
showOpTreeWithQubits :: Q.QOp -> String
showOpTreeWithQubits = render qopLabel (annotation . show) . qopTree

-- | The syntax tree of a UOp, root at the top. Identities and compositions
--   show their dimension in brackets, as in @I[3]@ and @∘[8]@.
showUOpTree :: U.UOp -> String
showUOpTree = render uopLabel (const "") . uopTree

-- | The same tree with @: n@ after every node of dimension @2^n@, and
--   @: d=...@ after a node of any other dimension, zero included.
showUOpTreeWithQubits :: U.UOp -> String
showUOpTreeWithQubits = render uopLabel (annotation . qubitLabel) . uopTree

printOpTree, printOpTreeWithQubits :: Q.QOp -> IO ()
printOpTree = putStrLn . showOpTree
printOpTreeWithQubits = putStrLn . showOpTreeWithQubits

printUOpTree, printUOpTreeWithQubits :: U.UOp -> IO ()
printUOpTree = putStrLn . showUOpTree
printUOpTreeWithQubits = putStrLn . showUOpTreeWithQubits

--------------------------------------------------------------------------------
-- Labels
--------------------------------------------------------------------------------

-- The trees come from "QPP.PrettyPrint.SyntaxTree". These two functions name a
-- node by the outermost constructor of its term.
qopLabel :: Q.QOp -> String
qopLabel op = case op of
  Q.Id _          -> showOp op
  Q.Phase t       -> "Phase " ++ angle t
  Q.X             -> showOp op
  Q.Y             -> showOp op
  Q.Z             -> showOp op
  Q.H             -> showOp op
  Q.SX            -> showOp op
  Q.R _ t         -> "R " ++ angle t
  Q.C _           -> "C"
  Q.Adjoint _     -> "adj"
  Q.Permute _     -> showOp op
  Q.Tensor _ _    -> "⊗"
  Q.DirectSum _ _ -> "⊕"
  Q.Compose _ _   -> "∘"

-- The two empty lists are different identities, and the drawing keeps them
-- apart from the empty composition @I[d]@.
uopLabel :: U.UOp -> String
uopLabel op = case op of
  U.R _ t             -> "R " ++ angle t
  U.Tensor []         -> "⊗[] = " ++ showUOp op
  U.DirectSum []      -> "⊕[] = " ++ showUOp op
  U.Compose _ []      -> showUOp op
  U.Tensor _          -> "⊗"
  U.DirectSum _       -> "⊕"
  U.Compose d _       -> "∘[" ++ show d ++ "]"
  U.PermuteTensor _ _ -> showUOp op
  U.PermuteSum _ _    -> showUOp op

-- | An angle as an argument in the syntax: @2@, @(1/4)@, @(-1/3)@.
angle :: Rational -> String
angle t
  | d == 1 && n >= 0 = show n
  | d == 1           = "(" ++ show n ++ ")"
  | otherwise        = "(" ++ show n ++ "/" ++ show d ++ ")"
  where
    n = numerator t
    d = denominator t

annotation :: String -> String
annotation n = " : " ++ n

--------------------------------------------------------------------------------
-- Layout
--------------------------------------------------------------------------------

-- | A drawn tree: its rows, root first, the width of its widest row, and the
--   column in which a branch from above meets the root.
data Drawing = Drawing { stem :: Int, width :: Int, rows :: [String] }

render :: (term -> String) -> (size -> String) -> SyntaxTree term size -> String
render label annotate =
  intercalate "\n" . map (dropWhileEnd isSpace) . rows . draw . fmap (\(term, n) -> (label term, annotate n))

-- | Siblings, and the subtrees below them, are this many columns apart.
gap :: Int
gap = 3

-- | A node is its label followed by its annotation. Its children are drawn
--   side by side below it, each subtree in a rectangle of its own. A row of
--   branches joins the middle of each child's label to the middle of the
--   node's label, which sits halfway between its first and its last child.
draw :: Tree (String, String) -> Drawing
draw (Node (label, note) children) = case map draw children of
  []    -> Drawing middle (length text) [text]
  drawn -> Drawing (shift + fork) (shift + max total (fork - middle + length text))
             (indent (shift + fork - middle) text : map (indent shift) (branches : below))
    where
      -- The children's left edges, and one more column past the last gap.
      offsets = scanl (\column d -> column + width d + gap) 0 drawn
      total   = last offsets - gap
      -- The children's stems, in increasing order.
      stems   = zipWith (\column d -> column + stem d) offsets drawn
      first   = minimum stems
      final   = maximum stems
      fork    = (first + final) `div` 2
      -- A label wider than its children moves them to the right.
      shift   = max 0 (middle - fork)
      branches = map branch [0 .. final]
      branch column
        | first == final        = if column == fork then '│' else ' '
        | column == first       = '┌'
        | column == final       = '┐'
        | column `elem` stems   = if column == fork then '┼' else '┬'
        | column == fork        = '┴'
        | column > first        = '─'
        | otherwise             = ' '
      height = maximum (map (length . rows) drawn)
      below  = map concat (transpose
        [ map (pad (width d + gap)) (take height (rows d ++ repeat "")) | d <- drawn ])
  where
    text   = label ++ note
    middle = (length label - 1) `div` 2

indent :: Int -> String -> String
indent n s = replicate n ' ' ++ s

pad :: Int -> String -> String
pad n s = s ++ replicate (n - length s) ' '
