-- | Text drawings of operator syntax trees (QPP.PrettyPrint.Tree).
--
--   Small drawings are pinned. Over random terms, each drawing is read back
--   into a tree of labels and compared with the constructor tree, which is
--   written out here independently of the printer, with 'n_qubits' and 'U.dim'
--   for the sizes. The Forest printer must give the same nodes at the same
--   depths.
module Main where

import Data.List (intercalate, isInfixOf)
import Data.Ratio (denominator, numerator)
import Data.Tree (Tree (Node), flatten, levels)
import Generators (RandomQOp (..))
import QPP
import qualified QPP.UOp as U
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck

main :: IO ()
main = defaultMain $ localOption (QuickCheckTests 200) $ testGroup "TreeTest"
  [ qopDrawings, qopLabels, uopDrawings, readBack ]

block :: [String] -> String
block = intercalate "\n"

------------------------------------------------------------------------
-- Pinned drawings
------------------------------------------------------------------------

qopDrawings :: TestTree
qopDrawings = testGroup "QOp drawings"
  [ testCase "operands hang below their constructor, left to right" $
      showOpTree (C X ∘ (H ⊗ I)) @?= block
        [ "   ∘"
        , "┌──┴──┐"
        , "C     ⊗"
        , "│   ┌─┴─┐"
        , "X   H   I" ]
  , testCase "qubit counts follow the labels, and branches still meet the labels" $
      showOpTreeWithQubits (C X ∘ (H ⊗ I)) @?= block
        [ "      ∘ : 2"
        , "┌─────┴─────┐"
        , "C : 2       ⊗ : 2"
        , "│       ┌───┴───┐"
        , "X : 1   H : 1   I : 1" ]
  , testCase "a rotation shows its angle and has its axis as its child" $
      showOpTree (R Z (-1/3)) @?= block
        [ "R (-1/3)"
        , "   │"
        , "   Z" ]
  , testCase "the two-qubit Fourier transform of the lecture notes" $
      showOpTree (Permute [1,0] ∘ (I ⊗ H) ∘ C (Phase (1/4) ⊗ R Z (1/2)) ∘ (H ⊗ I)) @?= block
        [ "                  ∘"
        , "      ┌───────────┴───────────┐"
        , "Permute [1,0]                 ∘"
        , "                  ┌───────────┴───────────┐"
        , "                  ⊗                       ∘"
        , "                ┌─┴─┐              ┌──────┴───────┐"
        , "                I   H              C              ⊗"
        , "                                   │            ┌─┴─┐"
        , "                                   ⊗            H   I"
        , "                             ┌─────┴─────┐"
        , "                        Phase (1/4)   R (1/2)"
        , "                                         │"
        , "                                         Z" ]
  , testCase "printing does not normalize" $
      showOpTree (Adjoint (Adjoint H) ∘ Id 1) @?= block
        [ "   ∘"
        , " ┌─┴──┐"
        , "adj   I"
        , " │"
        , "adj"
        , " │"
        , " H" ]
  ]

qopLabels :: TestTree
qopLabels = testGroup "QOp labels"
  [ testCase "leaves are spelled as showOp spells them" $
      mapM_ (\op -> showOpTree op @?= showOp op) [Id 0, Id 1, Id 3, X, Y, Z, H, SX, Permute [2,0,1]]
  , testCase "angles are exact arguments of the syntax" $
      map showOpTree [Phase 0, Phase 2, Phase (-3), Phase (1/4), Phase (-1/3)]
        @?= ["Phase 0", "Phase 2", "Phase (-3)", "Phase (1/4)", "Phase (-1/3)"]
  , testCase "qubit counts follow the width rules" $
      mapM_ (\(op, n) -> head (lines (showOpTreeWithQubits op)) `endsWith` (" : " ++ n))
        [ (Phase (1/3), "0"), (X, "1"), (SX, "1"), (Id 1000000, "1000000")
        , (Id 3 ⊗ X, "4"), (C (Id 3), "4"), (Id 3 ⊕ Id 3, "4")
        , (Id 3 ∘ Id 3, "3"), (Adjoint (Id 3), "3"), (R (Id 3) 0, "3"), (Permute [2,0,1], "3") ]
  , testCase "a wide identity is not enumerated" $
      showOpTreeWithQubits (Id 1000000) @?= "Id 1000000 : 1000000"
  ]
  where
    endsWith s suffix = assertBool (s ++ " should end with" ++ suffix) (suffix == drop (length s - length suffix) s)

uopDrawings :: TestTree
uopDrawings = testGroup "UOp drawings"
  [ testCase "a list constructor has all its operands, with a branch for each" $
      U.showUOpTree (U.Tensor [U.identity 2, U.identity 3, U.identity 5]) @?= block
        [ "        ⊗"
        , " ┌──────┼──────┐"
        , "I[2]   I[3]   I[5]" ]
  , testCase "empty lists keep their different dimensions" $
      U.showUOpTree (U.DirectSum [U.Tensor [], U.DirectSum [], U.identity 5]) @?= block
        [ "               ⊕"
        , "    ┌──────────┴─┬─────────┐"
        , "⊗[] = I[1]   ⊕[] = I[0]   I[5]" ]
  , testCase "a singleton list keeps its constructor" $
      U.showUOpTree (U.Tensor [U.identity 2]) @?= block [ " ⊗", " │", "I[2]" ]
  , testCase "sizes are qubit counts for powers of two and dimensions otherwise" $
      U.showUOpTreeWithQubits
        (U.Compose 3 [ U.PermuteSum [1,2] [1,0]
                     , U.DirectSum [U.phase (1/3), U.R (U.PermuteSum [1,1] [1,0]) (1/2)] ]) @?= block
        [ "                          ∘[3] : d=3"
        , "          ┌────────────────┴─────────────────┐"
        , "PermuteSum [1,2] [1,0] : d=3                 ⊕ : d=3"
        , "                                  ┌──────────┴──────────┐"
        , "                               R (-2/3) : 0          R (1/2) : 1"
        , "                                  │                     │"
        , "                                 I[1] : 0     PermuteSum [1,1] [1,0] : 1" ]
  , testCase "dimensions do not overflow machine integers" $ do
      let d = 2 ^ (100 :: Int)
      U.showUOpTreeWithQubits (U.identity d) @?= "I[" ++ show d ++ "] : 100"
      U.showUOpTreeWithQubits (U.identity (d + 1)) @?= "I[" ++ show (d + 1) ++ "] : d=" ++ show (d + 1)
      U.showUOpTreeWithQubits (U.identity 0) @?= "I[0] : d=0"
  ]

------------------------------------------------------------------------
-- Reading drawings back
------------------------------------------------------------------------

readBack :: TestTree
readBack = testGroup "A drawing read back is the constructor tree"
  [ testProperty "QOp, with and without qubit counts" prop_qop
  , testProperty "UOp as lowered and as normalized, with and without sizes" prop_uop
  , testProperty "the Forest printer gives the same nodes at the same depths" prop_forest
  ]

prop_qop :: RandomQOp -> Property
prop_qop (RandomQOp _ op) = conjoin
  [ drawn (showOpTree op) (qopSpec (const "") op)
  , drawn (showOpTreeWithQubits op) (qopSpec (\a -> " : " ++ show (n_qubits a)) op) ]

prop_uop :: RandomQOp -> Property
prop_uop (RandomQOp _ op) = case U.fromQOp op of
  Left _  -> property Discard
  Right u -> conjoin (map check [u, U.normalize u])
  where
    check u = conjoin
      [ drawn (U.showUOpTree u) (uopSpec (const "") u)
      , drawn (U.showUOpTreeWithQubits u) (uopSpec (\a -> " : " ++ size (U.dim a)) u) ]
    size d = case [ n | n <- [0 .. 16 :: Int], 2 ^ n == d ] of
      [n] -> show n
      _   -> "d=" ++ show d

prop_forest :: RandomQOp -> Property
prop_forest (RandomQOp _ op) =
  depths 0 (readDrawing (showOpTree op))
    === [ length (takeWhile (== ' ') line) `div` 2 | line <- lines (showTikZQOp op), "[{$" `isInfixOf` line ]
  where
    depths :: Int -> Tree String -> [Int]
    depths d (Node _ children) = d : concatMap (depths (d + 1)) children

-- | The drawing reads back as the expected tree and holds nothing else: a row
--   of labels per level with a row of branches between levels, every label
--   and every run of branches in the tree, and no space at the end of a row.
drawn :: String -> Tree String -> Property
drawn drawing expected = counterexample drawing $ conjoin
  [ readDrawing drawing === expected
  , length rows === 2 * length (levels expected) - 1
  , sum (map (length . fields) labelRows) === length (flatten expected)
  , sum (map (length . fields) branchRows) === length (filter (not . null) (flatten (parents expected)))
  , counterexample "a row ends in a space" (all (\row -> not (null row) && last row /= ' ') rows)
  ]
  where
    rows = lines drawing
    labelRows = [ row | (i, row) <- zip [0 :: Int ..] rows, even i ]
    branchRows = [ row | (i, row) <- zip [0 :: Int ..] rows, odd i ]
    parents (Node _ children) = Node children (map parents children)

-- | Read a drawing back into its tree of labels. Rows alternate between
--   labels and branches. Under a label, the run of branch characters that
--   meets it (at ┴, at ┼ or at a lone │) leads to its children, one below each
--   ┌, ┬, ┼, ┐ or │.
readDrawing :: String -> Tree String
readDrawing drawing = case rows of
  top : _ | [root] <- fields top -> grow rows root
  _ -> Node "no single root" []
  where
    rows = lines drawing

grow :: [String] -> (Int, String) -> Tree String
grow (_ : branches : below@(next : _)) parent =
  Node (snd parent)
    [ grow below child
    | (start, run) <- fields branches
    , any (covers parent . (start +)) (marks "┴┼│" run)
    , stem <- marks "┌┬┼┐│" run
    , child <- fields next
    , covers child (start + stem) ]
grow _ leaf = Node (snd leaf) []

covers :: (Int, String) -> Int -> Bool
covers (start, text) column = start <= column && column < start + length text

marks :: String -> String -> [Int]
marks set run = [ i | (i, c) <- zip [0 ..] run, c `elem` set ]

-- | The texts in a row with their first columns. A text has no two spaces in
--   a row, and two texts are at least two spaces apart.
fields :: String -> [(Int, String)]
fields = go 0
  where
    go _ "" = []
    go i (' ' : row) = go (i + 1) row
    go i row = (i, text) : go (i + length text) rest
      where (text, rest) = field row
    field (c : ' ' : d : row) | c /= ' ', d /= ' ' = (c : ' ' : text, rest)
      where (text, rest) = field (d : row)
    field (c : row) | c /= ' ' = (c : text, rest)
      where (text, rest) = field row
    field row = ("", row)

------------------------------------------------------------------------
-- The constructor trees, written out
------------------------------------------------------------------------

qopSpec :: (QOp -> String) -> QOp -> Tree String
qopSpec note op = Node (name ++ note op) (map (qopSpec note) operands)
  where
    (name, operands) = case op of
      Tensor a b    -> ("⊗", [a, b])
      DirectSum a b -> ("⊕", [a, b])
      Compose a b   -> ("∘", [a, b])
      C a           -> ("C", [a])
      Adjoint a     -> ("adj", [a])
      R a t         -> ("R " ++ argument t, [a])
      Phase t       -> ("Phase " ++ argument t, [])
      Id 0          -> ("One", [])
      Id 1          -> ("I", [])
      _             -> (show op, [])

uopSpec :: (U.UOp -> String) -> U.UOp -> Tree String
uopSpec note op = Node (name ++ note op) (map (uopSpec note) operands)
  where
    (name, operands) = case op of
      U.R a t           -> ("R " ++ argument t, [a])
      U.Tensor []       -> ("⊗[] = I[1]", [])
      U.DirectSum []    -> ("⊕[] = I[0]", [])
      U.Compose d []    -> ("I[" ++ show d ++ "]", [])
      U.Tensor xs       -> ("⊗", xs)
      U.DirectSum xs    -> ("⊕", xs)
      U.Compose d xs    -> ("∘[" ++ show d ++ "]", xs)
      U.PermuteTensor{} -> (show op, [])
      U.PermuteSum{}    -> (show op, [])

-- | A rational as an argument in Haskell source.
argument :: Rational -> String
argument t
  | t >= 0, denominator t == 1 = show (numerator t)
  | otherwise = "(" ++ show (numerator t) ++ concat [ "/" ++ show (denominator t) | denominator t /= 1 ] ++ ")"
