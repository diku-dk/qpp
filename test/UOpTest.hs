module Main where

import Data.Ratio ((%))
import Programs.QFT (qft)
import QPP (adj, (<.>), (<+>), (∘), (>:), showUOp)
import qualified QPP as Q
import qualified QPP.UOp as U
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck
import UOpNormalizationTests (normalizationTests)

main :: IO ()
main = defaultMain $ testGroup "UOp"
  [dimensionTests, constructionTests, adjointTests, prettyTests, loweringTests, normalizationTests]

dimensionTests :: TestTree
dimensionTests = testGroup "Dimensions count basis vectors"
  [ testCase (show op) $ U.dim op @?= d
  | (op, d) <-
      [ (U.identity 0, 0), (U.identity 1, 1), (U.identity 3, 3)
      , (U.phase (1 % 3), 1)
      , (U.R (U.identity 7) (2 % 3), 7)
      , (U.Tensor [], 1), (U.DirectSum [], 0)
      , (U.Tensor [U.identity 2, U.identity 3], 6)
      , (U.DirectSum [U.identity 2, U.identity 0, U.identity 3], 5)
      , (U.Tensor [U.identity 2, U.identity 0, U.identity 3], 0)
      , (U.Compose 3 [U.identity 3, U.identity 3], 3)
      , (U.PermuteTensor [2,3,4] [2,0,1], 24)
      , (U.PermuteSum [2,3,4] [2,0,1], 9)
      , (U.PermuteTensor [] [], 1), (U.PermuteSum [] [], 0)
      , (U.PermuteTensor [0,3] [1,0], 0)
      , (U.PermuteSum [0,3] [1,0], 3)
      , (U.Tensor (replicate 100 (U.identity 2)), 2 ^ (100 :: Int))
      ]
  ]

constructionTests :: TestTree
constructionTests = testGroup "Construction and shared notation"
  [ testCase "Tensor associations build the same list" $ do
      (a <.> b) <.> c @?= a <.> (b <.> c)
      a <.> b <.> c @?= U.Tensor [a,b,c]
  , testCase "Direct-sum associations build the same list" $ do
      (a <+> b) <+> c @?= a <+> (b <+> c)
      a <+> b <+> c @?= U.DirectSum [a,b,c]
  , testCase "Composition preserves mathematical order" $ do
      (a ∘ b) ∘ c @?= a ∘ (b ∘ c)
      a ∘ b ∘ c @?= U.Compose 1 [a,b,c]
      a >: b @?= U.Compose 1 [b,a]
  , testCase "Nested tensors flatten recursively" $
      U.tensor [a, U.Tensor [U.Tensor [], U.Tensor [b,c]]]
        @?= U.Tensor [a,b,c]
  , testCase "Nested sums flatten recursively" $
      U.directSum [U.DirectSum [], U.DirectSum [a, U.DirectSum [b,c]]]
        @?= U.DirectSum [a,b,c]
  , testCase "Composition identities retain the outer dimension" $
      U.compose 3 [U.identity 3, U.Compose 3 [U.identity 3]] @?= U.identity 3
  , testCase "A mismatched empty composition is not erased" $
      U.compose 2 [U.identity 3] @?= U.Compose 2 [U.identity 3]
  , testCase "Scalar phases keep the exact angle" $
      U.phase (-7 % 3) @?= U.R (U.identity 1) (14 % 3)
  ]
  where
    a = U.phase (1 % 3)
    b = U.phase (1 % 4)
    c = U.phase (1 % 5)

adjointTests :: TestTree
adjointTests = testGroup "Structural adjoints"
  [ testCase "Negates the angle, not the Hermitian axis" $
      U.dagger (U.R x (1 % 3)) @?= U.R x (-1 % 3)
  , testCase "Only composition reverses operand order" $ do
      U.dagger (U.Tensor [x,z]) @?= U.Tensor [U.dagger x, U.dagger z]
      U.dagger (U.DirectSum [x,z]) @?= U.DirectSum [U.dagger x, U.dagger z]
      U.dagger (U.Compose 2 [x,z]) @?= U.Compose 2 [U.dagger z, U.dagger x]
  , testCase "Unequal block swap has a different adjoint" $
      U.dagger (U.PermuteSum [1,2] [1,0]) @?= U.PermuteSum [2,1] [1,0]
  , testCase "Tensor inverse starts from permuted dimensions" $
      U.dagger (U.PermuteTensor [2,3,4] [1,2,0])
        @?= U.PermuteTensor [3,4,2] [2,0,1]
  , testCase "Adjoint supports zero-sized blocks" $
      U.dagger (U.PermuteSum [0,2,1] [2,0,1])
        @?= U.PermuteSum [1,0,2] [1,2,0]
  , testCase "Shared adj notation" $ adj x @?= U.dagger x
  , testProperty "Permutation adjoint is involutive and preserves dimension" $
      forAll permutationCase $ \(ds,p) -> conjoin
        [ let op = con ds p
          in conjoin [U.dagger (U.dagger op) === op, U.dim (U.dagger op) === U.dim op]
        | con <- [U.PermuteTensor, U.PermuteSum]
        ]
  , testCase "Adjoint is involutive for a composite term" $
      let op = U.Tensor [U.Compose 2 [U.R x (1 % 3),z], U.DirectSum [U.identity 0,x]]
      in U.dagger (U.dagger op) @?= op
  ]
  where
    x = U.PermuteSum [1,1] [1,0]
    z = U.DirectSum [U.identity 1, U.phase 1]
    permutationCase = do
      k <- chooseInt (0,6)
      ds <- vectorOf k (fromIntegral <$> chooseInt (0,5))
      p <- shuffle [0..k-1]
      pure (ds,p)

prettyTests :: TestTree
prettyTests = testGroup "Pretty printing"
  [ testCase (show op) $ showUOp op @?= text
  | (op, text) <-
      [ (U.identity 0, "I[0]")
      , (U.identity 1, "I[1]")
      , (U.identity 3, "I[3]")
      , (U.Tensor [], "I[1]")
      , (U.DirectSum [], "I[0]")
      , (U.Tensor [U.identity 3], "I[3]")
      , (U.DirectSum [U.identity 3], "I[3]")
      , (U.Compose 3 [U.identity 3], "I[3]")
      , (U.Tensor [U.identity 2,U.identity 3,U.identity 4], "(I[2] ⊗ I[3] ⊗ I[4])")
      , (U.DirectSum [U.identity 2,U.identity 3], "(I[2] ⊕ I[3])")
      , (U.Compose 2 [U.identity 2,U.identity 2], "(I[2] ∘ I[2])")
      , (U.R (U.identity 2) (1 % 3), "R(I[2], 1/3)")
      , (U.R (U.identity 2) (-2 % 3), "R(I[2], -2/3)")
      , (U.phase (-2), "R(I[1], 4)")
      , (U.phase 0, "R(I[1], 0)")
      , (U.PermuteTensor [2,3] [1,0], "PermuteTensor [2,3] [1,0]")
      , (U.PermuteSum [1,2] [1,0], "PermuteSum [1,2] [1,0]")
      , (U.PermuteTensor [] [], "PermuteTensor [] []")
      , (U.PermuteSum [] [], "PermuteSum [] []")
      , (U.R (U.DirectSum [U.identity 1, U.phase 1]) (1 % 2),
         "R((I[1] ⊕ R(I[1], -2)), 1/2)")
      , (U.Tensor [U.DirectSum [U.identity 1,U.identity 2],U.identity 3],
         "((I[1] ⊕ I[2]) ⊗ I[3])")
      ]
  ]

loweringTests :: TestTree
loweringTests = testGroup "Symbolic QOp lowering"
  [ testCase "Rotation lowering preserves QOp's angle" $
      U.fromQOp (Q.R Q.X (1 % 3)) @?= Right (U.R (U.PermuteSum [1,1] [1,0]) (1 % 3))
  , testGroup "Preserves dimension and produces consistent composition shapes"
      [ testCase (show q) $ case (Q.typeCheck q, U.fromQOp q) of
          (Right n, Right u) -> do
            U.dim u @?= 2 ^ n
            assertBool "well-formed composition dimensions" (consistent u)
          result -> assertFailure (show result)
      | q <- [Q.Id 0, Q.Id 3, Q.Phase (-3 % 7), Q.X, Q.Y, Q.Z, Q.H, Q.SX,
              Q.R Q.H (2 % 5), Q.R Q.SX 0, Q.C Q.H, Q.C Q.Y, Q.C Q.SX,
              Q.C (Q.Phase (1 % 3)), Q.Permute [], Q.Permute [1,2,0],
              Q.Tensor Q.H Q.SX, Q.DirectSum Q.H Q.SX, Q.Compose Q.X Q.Y,
              Q.Adjoint (Q.Permute [1,2,0])] ++ map qft [0..16]
      ]
  , testGroup "Rejects malformed source terms before lowering"
      [ testCase (show q) $ case Q.typeCheck q of
          Left err -> U.fromQOp q @?= Left err
          Right _ -> assertFailure "test input should be rejected"
      | q <- [Q.Id (-1), Q.Permute [0,0], Q.DirectSum Q.X (Q.Id 2),
              Q.Compose Q.X (Q.Id 2), Q.R Q.SX 1, Q.R (Q.Id (-1)) 0]
      ]
  ]
  where
    consistent u = case u of
      U.R a _ -> consistent a
      U.Tensor xs -> all consistent xs
      U.DirectSum xs -> all consistent xs
      U.Compose d xs -> all (\a -> U.dim a == d && consistent a) xs
      _ -> True
