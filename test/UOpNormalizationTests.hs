module UOpNormalizationTests (normalizationTests) where

import Data.Complex (Complex((:+)), magnitude)
import Data.List (sortOn, transpose)
import Data.Ratio ((%))
import Programs.QFT (qft)
import qualified QPP as Q
import qualified QPP.UOp as U
import qualified QPP.Semantics.Matrix as M
import qualified Numeric.LinearAlgebra as LA
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck

normalizationTests :: TestTree
normalizationTests = testGroup "Symbolic normalization"
  [ testGroup "Algebraic identities"
      [ testCase name $ U.normalize input @?= expected
      | (name, input, expected) <- examples
      ]
  , testGroup "Small examples preserve vector action, including phase"
      [ testCase name $ assertBool "vector action changed" (sameAction input)
      | (name, input, _) <- examples
      ]
  , testGroup "Lowering matches QOp matrix semantics, including phase"
      [ testCase (show q) $ case U.fromQOp q of
          Left err -> assertFailure (show err)
          Right u -> do
            let matrix = M.evalOp q
                basis j = [if k == j then 1 else 0 | k <- [0 .. LA.cols matrix - 1]]
            mapM_ (\j -> do
              let v = basis j
                  expected = LA.toList (matrix LA.#> LA.fromList v)
              assertBool "raw lowering changed an amplitude" (closeVectors expected (apply u v))
              assertBool "normalized lowering changed an amplitude"
                (closeVectors expected (apply (U.normalize u) v))) [0 .. LA.cols matrix - 1]
      | q <- [Q.Id 0, Q.Phase (1/3), Q.X, Q.Y, Q.Z, Q.H, Q.SX,
              Q.C Q.H, Q.C Q.SX, Q.C (Q.Phase (1/3)),
              Q.Adjoint (Q.R Q.Y (1/3)), Q.Permute [1,2,0]]
          ++ [Q.R a t | a <- [Q.Id 0, Q.X, Q.Y, Q.Z, Q.H, Q.C Q.X],
                        t <- [-4,-2,-1,-1/3,0,1/3,1,2,4]]
          ++ map qft [0..3]
      ]
  , localOption (QuickCheckTests 500) $
      testProperty "Preserves dimension and vector action; is idempotent" $
        forAll (chooseInt (0,8) >>= genOp 3) $ \a ->
          let b = U.normalize a
          in conjoin [U.dim b === U.dim a, U.normalize b === b,
                      property (compositionsOutside b), property (sameAction a)]
  , testGroup "QFT normalization stays symbolic"
      [ testCase ("qft " ++ show n) $ case U.fromQOp (qft n) of
          Left err -> assertFailure (show err)
          Right a -> do
            let b = U.normalize a
            U.dim b @?= U.dim a
            U.normalize b @?= b
            assertBool "composition remains inside a parallel layer" (compositionsOutside b)
            if n >= 2 then case b of
              U.Compose _ _ -> pure ()
              _ -> assertFailure "QFT should be an outer composition"
            else pure ()
            if n <= 4 then assertBool "QFT vector action changed" (sameAction a) else pure ()
      | n <- [0..16]
      ]
  , testCase "Controlled relative phase becomes a diagonal block" $
      fmap U.normalize (U.fromQOp (Q.C (Q.Tensor (Q.Phase (1/4)) (Q.R Q.Z (1/2)))))
        @?= Right (U.DirectSum [U.identity 3, U.phase (1/2)])
  ]

examples :: [(String, U.UOp, U.UOp)]
examples =
  [ ("Empty tensor", U.Tensor [], i 1)
  , ("Empty sum", U.DirectSum [], i 0)
  , ("Nested units and singletons", U.Tensor [i 1, U.Tensor [U.DirectSum [i 0,x]]], x)
  , ("Zero-dimensional tensor", U.Tensor [x,i 0], i 0)
  , ("Merge adjacent tensor identities", U.Tensor [i 2,i 3,x,i 4], U.Tensor [i 6,x,i 4])
  , ("Merge adjacent sum identities", U.DirectSum [i 2,i 3,x,i 4], U.DirectSum [i 5,x,i 4])
  , ("Composition units", U.Compose 2 [i 2,U.Compose 2 [x,i 2]], x)
  , ("Tensor lifts composition with identity padding", U.Tensor [U.Compose 2 [x,z],s],
        U.Compose 6 [U.Tensor [x,s],U.Tensor [z,i 3]])
  , ("Tensor zips unequal composition lengths", U.Tensor [U.Compose 2 [x,z],U.Compose 3 [s,s,s]],
        U.Compose 6 [U.Tensor [x,s],U.Tensor [z,s],U.Tensor [i 2,s]])
  , ("An empty composition pads with its own dimension", U.Tensor [i 3,U.Compose 2 [x,z]],
        U.Compose 6 [U.Tensor [i 3,x],U.Tensor [i 3,z]])
  , ("Direct sum lifts composition with identity padding", U.DirectSum [U.Compose 2 [x,z],s],
        U.Compose 5 [U.DirectSum [x,s],U.DirectSum [i 1,U.phase 1,i 3]])
  , ("Direct sum zips unequal composition lengths", U.DirectSum [U.Compose 2 [x,z],U.Compose 3 [s,s,s]],
        U.Compose 5 [U.DirectSum [x,s],U.DirectSum [i 1,U.phase 1,s],U.DirectSum [i 2,s]])
  , ("A global phase multiplies one composition factor", U.Tensor [U.phase (1/3),U.Compose 2 [x,z]],
        U.Compose 2 [U.Tensor [U.phase (1/3),x],z])
  , ("A relative phase occurs in only one layer", U.DirectSum [U.Compose 2 [x,z],U.phase (1/3)],
        U.Compose 3 [U.DirectSum [x,U.phase (1/3)],U.DirectSum [i 1,U.phase 1,i 1]])
  , ("Signed phase representative", U.phase (7/4), U.phase (-1/4))
  , ("Tensor scalar phases", U.Tensor [U.phase (1/3),x,U.phase (5/3)], x)
  , ("Composition scalar phases", U.Compose 2 [U.Tensor [U.phase (1/3),x],
                                               U.Tensor [U.phase (5/3),x]], i 2)
  , ("Relative phases are retained", U.DirectSum [U.phase (1/3),U.phase (2/3)],
                                     U.DirectSum [U.phase (1/3),U.phase (2/3)])
  , ("Absorb a scalar through tensor padding", U.Tensor [U.phase (1/3),i 5,
        U.DirectSum [U.phase (-1/3),U.phase (1/3)]],
        U.Tensor [i 5,U.DirectSum [i 1,U.phase (2/3)]])
  , ("Absorb into only one direct-sum factor", U.Tensor [U.phase 1,z,z],
        U.Tensor [U.DirectSum [U.phase 1,i 1],z])
  , ("Zero angle allows a non-Hermitian axis", U.R (U.phase (1/3)) 0, i 1)
  , ("Rotation period is four", U.R x (13/3), U.R x (1/3))
  , ("Two full turns give identity", U.R x 4, i 2)
  , ("One full turn retains minus identity", U.R x 2, U.Tensor [U.phase 1,i 2])
  , ("Positive half turn gives minus i times the axis", U.R x 1, U.Tensor [U.phase (-1/2),x])
  , ("Negative half turn gives i times the axis", U.R x (-1), U.Tensor [U.phase (1/2),x])
  , ("Three half turns retain their phase", U.R x 3, U.Tensor [U.phase (1/2),x])
  , ("Negative axis", U.R (U.Tensor [U.phase 1,x]) (1/3), U.R x (-1/3))
  , ("Direct-sum rotation", U.R z (1/3), U.DirectSum [U.phase (-1/6),U.phase (1/6)])
  , ("Identity-axis rotation", U.R (i 5) (1/3), U.Tensor [U.phase (-1/6),i 5])
  , ("Same-axis fusion", U.Compose 2 [U.R x (1/3),U.R x (1/3)], U.R x (2/3))
  , ("Cancellation exposes another pair", U.Compose 2 [U.R x (1/3),x,x,U.R x (-1/3)], i 2)
  , ("Identity tensor permutation", U.PermuteTensor [2,3] [0,1], i 6)
  , ("Zero-dimensional permutation", U.PermuteTensor [0,3] [1,0], i 0)
  , ("Unequal tensor permutation inverse", U.Compose 6
        [U.PermuteTensor [3,2] [1,0],U.PermuteTensor [2,3] [1,0]], i 6)
  , ("Unequal block permutation inverse", U.Compose 3
        [U.PermuteSum [2,1] [1,0],U.PermuteSum [1,2] [1,0]], i 3)
  , ("Unequal block swap is not self-inverse", U.Compose 3 [s,s], U.Compose 3 [s,s])
  , ("Sum permutation with a zero block", U.PermuteSum [0,2,1] [2,0,1],
                                         U.PermuteSum [0,2,1] [2,0,1])
  ]
  where
    i = U.identity
    x = U.PermuteSum [1,1] [1,0]
    z = U.DirectSum [i 1,U.phase 1]
    s = U.PermuteSum [1,2] [1,0]

-- Identity is an empty Compose, so only nonempty compositions are forbidden
-- immediately inside a tensor or sum. Rotation axes may still compose.
compositionsOutside :: U.UOp -> Bool
compositionsOutside op = case op of
  U.Tensor xs -> all parallelChild xs
  U.DirectSum xs -> all parallelChild xs
  U.Compose _ xs -> all compositionsOutside xs
  U.R a _ -> compositionsOutside a
  _ -> True
  where
    parallelChild (U.Compose _ (_:_)) = False
    parallelChild a = compositionsOutside a

-- Generate well-formed terms of small arbitrary dimension. Conjugated sign
-- diagonals give Hermitian axes that need not themselves be diagonal.
genOp :: Int -> Int -> Gen U.UOp
genOp depth d
  | depth == 0 = leaf
  | otherwise = oneof
      [ leaf
      , U.Compose (fromIntegral d) <$> vectorOf 3 child
      , do k <- chooseInt (0,d)
           a <- genOp (depth-1) k
           b <- genOp (depth-1) (d-k)
           pure (U.DirectSum [a,b])
      , do (a,b) <- elements [(a,b) | a <- [0..d+1], b <- [0..d+1], a*b == d]
           U.Tensor <$> sequence [genOp (depth-1) a, genOp (depth-1) b]
      , do u <- child
           signs <- vectorOf d (U.phase <$> elements [0,1])
           t <- rational
           pure (U.R (U.Compose (fromIntegral d) [u,U.DirectSum signs,U.dagger u]) t)
      ]
  where
    child = genOp (depth-1) d
    leaf = oneof $
      [ pure (U.identity (fromIntegral d))
      , U.DirectSum <$> vectorOf d (U.phase <$> rational)
      , U.PermuteSum (replicate d 1) <$> shuffle [0..d-1]
      ] ++ [pure (U.PermuteTensor [2,3] [1,0]) | d == 6]
    rational = (%) <$> chooseInteger (-8,8) <*> chooseInteger (1,8)

sameAction :: U.UOp -> Bool
sameAction a = closeVectors before result
  where
    -- Complex, nonuniform amplitudes distinguish reordered basis positions.
    v = [fromIntegral (k+1) :+ fromIntegral (k*k-3) | k <- [0 .. fromIntegral (U.dim a)-1 :: Int]]
    before = apply a v
    result = apply (U.normalize a) v

closeVectors :: [Complex Double] -> [Complex Double] -> Bool
closeVectors xs ys = length xs == length ys && and (zipWith close xs ys)
  where close x y = magnitude (x-y) < 1e-8 * (1 + magnitude x)

-- Independent test-only vector interpreter. Rotations act via
-- cos(pi*t/2) v - i sin(pi*t/2) A v; no operator matrix is constructed.
apply :: U.UOp -> [Complex Double] -> [Complex Double]
apply op _
  | U.dim op == 0 = []
apply (U.R _ 0) v = v
apply (U.R a t) v = zipWith (+) (map (c *) v) (map (s *) (apply a v))
  where
    theta = pi * fromRational t / 2
    c = cos theta :+ 0
    s = 0 :+ (-sin theta)
apply (U.Compose _ xs) v = foldr apply v xs
apply (U.DirectSum xs) v = concat (zipWith apply xs (blocks (map (fromIntegral . U.dim) xs) v))
apply (U.Tensor []) v = v
apply (U.Tensor (a:as)) v =
  concat (transpose (map (apply a) (transpose (map (apply rest) rows))))
  where
    rest = U.Tensor as
    rows = blocks (replicate (fromIntegral (U.dim a)) (fromIntegral (U.dim rest))) v
apply (U.PermuteSum ds p) v = concat [chunks !! k | k <- p]
  where chunks = blocks (map fromIntegral ds) v
apply (U.PermuteTensor ds p) v = map snd (sortOn fst indexed)
  where
    sizes = map fromIntegral ds
    reorder xs = [xs !! k | k <- p]
    indexed = [(encode (reorder sizes) (reorder (decode sizes k)), a) | (k,a) <- zip [0..] v]
    decode [] _ = []
    decode (_:rest) k = let (q,r) = k `divMod` product rest in q : decode rest r
    encode sizes' digits = foldl (\k (d,j) -> k*d+j) (0 :: Int) (zip sizes' digits)

blocks :: [Int] -> [a] -> [[a]]
blocks [] _ = []
blocks (d:ds) xs = let (a,rest) = splitAt d xs in a : blocks ds rest
