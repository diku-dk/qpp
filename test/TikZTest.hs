module Main where

import Data.List (isInfixOf)
import qualified QPP as Q
import qualified QPP.UOp as U
import Test.Tasty
import Test.Tasty.HUnit

main :: IO ()
main = defaultMain $ testGroup "TikZ printer contracts"
  [ testCase "QOp composition keeps mathematical operand order" $
      nodes (Q.showTikZQOpWithoutQubits (Q.Compose Q.H Q.X)) @?=
        [ "[{$\\circ$}, node options+={qpp/compose}"
        , "  [{$H$}, node options+={qpp/gate}"
        , "  [{$X$}, node options+={qpp/gate}"
        ]
  , testCase "Nested operators keep their parent-child structure" $
      nodes (Q.showTikZQOpWithoutQubits (Q.Tensor Q.H (Q.DirectSum Q.X Q.Y))) @?=
        [ "[{$\\otimes$}, node options+={qpp/tensor}"
        , "  [{$H$}, node options+={qpp/gate}"
        , "  [{$\\oplus$}, node options+={qpp/sum}"
        , "    [{$X$}, node options+={qpp/gate}"
        , "    [{$Y$}, node options+={qpp/gate}"
        ]
  , testCase "UOp composition retains its dimension and all operands" $
      nodes (U.showTikZUOpWithoutQubits (U.Compose 3 (replicate 3 (U.identity 3)))) @?=
        "[{$\\circ_{3}$}, node options+={qpp/compose}" : replicate 3 "  [{$I_{3}$}, node options+={qpp/identity}"
  , testCase "Printing retains singleton constructor structure" $
      nodes (U.showTikZUOpWithoutQubits (U.Tensor [U.identity 2])) @?=
        ["[{$\\otimes$}, node options+={qpp/tensor}", "  [{$I_{2}$}, node options+={qpp/identity}"]
  , testCase "Empty lists retain their different dimensions" $ do
      contains (U.showTikZUOp (U.Tensor [])) "\\otimes[\\,]=I_{1}"
      contains (U.showTikZUOp (U.DirectSum [])) "\\oplus[\\,]=I_{0}"
      contains (U.showTikZUOp (U.Compose 7 [])) "I_{7}"
  , testCase "Rotation angles are exact and signed" $
      contains (Q.showTikZQOp (Q.R Q.Z (-1/3))) "R\\left(-\\frac{1}{3}\\right)"
  , testCase "A UOp rotation prints its axis even at angle zero" $ do
      let tex = U.showTikZUOp (U.R (U.PermuteSum [1,2] [1,0]) 0)
      contains tex "R\\left(0\\right)"
      contains tex "qpp/rotation"
      contains tex "P_{\\oplus}"
  , testCase "Braces protect permutation brackets and commas from Forest" $
      contains (U.showTikZUOp (U.PermuteTensor [2,3,4] [2,0,1]))
        "[{$P_{\\otimes}\\left([2,3,4];\\,[2,0,1]\\right)$}, node options+={qpp/permutation}"
  , testCase "QOp identity width is shown symbolically" $ do
      let tex = Q.showTikZQOp (Q.Id 1000000)
      contains tex "I_{2^{1000000}}"
      assertBool "printing must not enumerate dimensions" (length tex < 200)
  , testCase "UOp dimensions do not overflow machine integers" $
      contains (U.showTikZUOp (U.identity (2 ^ (100 :: Int)))) "1267650600228229401496703205376"
  , testCase "Styles stay in LaTeX, not hard-coded on individual nodes" $ do
      let tex = Q.showTikZQOp (Q.C (Q.Adjoint (Q.DirectSum Q.X Q.Y)))
      mapM_ (contains tex) ["qpp/control", "qpp/adjoint", "qpp/sum"]
      assertBool "fragment should not override the document theme" (not ("\\tikzset" `isInfixOf` tex))
  , testCase "Fragments use Forest without paragraph breaks" $ do
      let tex = Q.showTikZQOp (Q.Tensor Q.H Q.X)
      assertBool "an empty line would break the surrounding paragraph" (not ("\n\n" `isInfixOf` tex))
      contains tex "\\begin{forest}\nqpp tree\n"
      contains tex "\\end{forest}%\n"
  , testCase "Qubit counts annotate every node outside its math label" $
      nodes (Q.showTikZQOp (Q.DirectSum (Q.Id 3) (Q.Id 3))) @?=
        [ "[{$\\oplus$}, node options+={qpp/sum}, qpp annotation={4}"
        , "  [{$I_{2^{3}}$}, node options+={qpp/identity}, qpp annotation={3}"
        , "  [{$I_{2^{3}}$}, node options+={qpp/identity}, qpp annotation={3}"
        ]
  , testCase "QOp counts follow the width rules" $
      mapM_ (\(op, n) -> rootCount (Q.showTikZQOp op) n)
        [ (Q.Phase (1/3), "0"), (Q.X, "1"), (Q.Y, "1"), (Q.Z, "1")
        , (Q.H, "1"), (Q.SX, "1"), (Q.Id 1000000, "1000000")
        , (Q.Tensor (Q.Id 3) Q.X, "4"), (Q.C (Q.Id 3), "4")
        , (Q.Compose (Q.Id 3) (Q.Id 3), "3"), (Q.Adjoint (Q.Id 3), "3")
        , (Q.R (Q.Id 3) 0, "3"), (Q.Permute [2,0,1], "3")
        ]
  , testCase "UOp counts use total dimensions, including non-qubit operands" $
      mapM_ (\(op, n) -> rootCount (U.showTikZUOp op) n)
        [ (U.Tensor [], "0"), (U.DirectSum [], "d=0"), (U.identity 0, "d=0")
        , (U.phase (1/3), "0"), (U.identity 3, "d=3")
        , (U.Tensor [U.identity 2, U.identity 8], "4")
        , (U.DirectSum [U.identity 3, U.identity 5], "3")
        , (U.Compose 8 [U.identity 8, U.identity 8], "3")
        , (U.PermuteTensor [2,4] [1,0], "3")
        , (U.PermuteSum [3,5] [1,0], "3"), (U.R (U.identity 16) 0, "4")
        ]
  , testCase "Power-of-two detection is exact for large dimensions" $ do
      let d = 2 ^ (100 :: Int)
      rootCount (U.showTikZUOp (U.identity d)) "100"
      rootCount (U.showTikZUOp (U.identity (d + 1))) ("d=" ++ show (d + 1))
  , testCase "Annotations can be explicitly omitted" $
      mapM_ (\tex -> assertBool "explicit opt-out should omit counts"
        (not ("qpp annotation" `isInfixOf` tex)))
        [Q.showTikZQOpWithoutQubits (Q.C Q.X), U.showTikZUOpWithoutQubits (U.identity 2)]
  ]

rootCount :: String -> String -> Assertion
rootCount tex n = contains (head (nodes tex)) ("qpp annotation={" ++ n ++ "}")

-- Indentation records depth, so these lines expose both child order and scope.
nodes :: String -> [String]
nodes = filter (isInfixOf "[{$") . lines

contains :: String -> String -> Assertion
contains tex fragment = assertBool ("Missing LaTeX fragment: " ++ fragment) (fragment `isInfixOf` tex)
