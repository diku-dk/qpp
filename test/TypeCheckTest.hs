module Main where

import Data.Ratio ((%))
import QPP
import qualified QPP.Semantics.Matrix as MS
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck
import Generators

arityTests :: TestTree
arityTests = testGroup "Valid operators and inferred widths"
  [ testCase (show op) $ typeCheck op @?= Right n
  | (op, n) <-
      [ (Id 0, 0), (Id 5, 5), (Phase (1 % 3), 0)
      , (X, 1), (Y, 1), (Z, 1), (H, 1), (SX, 1)
      , (R (Tensor X Z) (1 % 3), 2)
      , (C (C X), 3), (C (Phase (1 % 2)), 1)
      , (Tensor (Id 2) H, 3), (Tensor (Phase (1 % 2)) X, 1)
      , (DirectSum X H, 2), (DirectSum (Id 0) (Phase 1), 1)
      , (Compose H SX, 1), (Adjoint (C X), 2)
      , (Permute [], 0), (Permute [0], 1)
      , (Permute [1,0], 2), (Permute [1,2,0], 3)
      , (2 <@ X @> 3, 6)
      , (Id maxBound, maxBound)
      ]
  ]

errorTests :: TestTree
errorTests = testGroup "Invalid operators"
  [ testGroup "Negative widths, including padding"
      [ testCase (show op) $ typeCheck op @?= Left (NegativeQubitCount (-1))
      | op <- [Id (-1), C (Id (-1)), Adjoint (Id (-1))
              , Tensor X (Id (-1)), Compose X (Id (-1))
              , DirectSum X (Id (-1)), R (Id (-1)) 0
              , (-1) <@ X, X @> (-1)]
      ]
  , testGroup "Invalid permutations"
      [ testCase (show ks) $ typeCheck (Permute ks) @?= Left (InvalidPermutation ks)
      | ks <- [[-1], [1], [0,0], [0,2], [1,1,0], [2,0], [-1,0], [maxBound]]
      ]
  , testGroup "Arity mismatches"
      [ testCase (show op) $ typeCheck op @?= Left (ArityMismatch op m n)
      | (op, m, n) <-
          [ (Compose X (Id 2), 1, 2), (Compose (Id 2) X, 2, 1)
          , (DirectSum X (Id 2), 1, 2), (DirectSum (Id 2) X, 2, 1)
          , (Compose (Phase 1) X, 0, 1)
          ]
      ]
  , testCase "Reports a nested mismatch at the offending subterm" $
      let bad = DirectSum X (Id 2)
      in typeCheck (Adjoint (Tensor H (C bad))) @?= Left (ArityMismatch bad 1 2)
  , testCase "Checks both children even if the first determines the width" $
      let bad = Compose X (Id 2)
      in typeCheck (DirectSum X bad) @?= Left (ArityMismatch bad 1 2)
  , testCase "Zero angle still checks the axis permutation" $
      typeCheck (R (Permute [0,0]) 0) @?= Left (InvalidPermutation [0,0])
  , testCase "Zero angle still checks the axis arity" $
      let bad = Compose X (Id 2)
      in typeCheck (R bad 0) @?= Left (ArityMismatch bad 1 2)
  , testCase "Zero angle still checks nested nonzero-angle rotations" $
      typeCheck (R (R SX (1 % 3)) 0) @?= Left (UnprovenHermitian SX)
  , testCase "Zero angle still checks for width overflow" $
      let bad = C (Id maxBound)
      in typeCheck (R bad 0) @?= Left (QubitCountOverflow bad)
  , testGroup "Widths cannot overflow Int"
      [ testCase (show op) $ typeCheck op @?= Left (QubitCountOverflow op)
      | op <- [C (Id maxBound), Tensor (Id maxBound) X
              , DirectSum (Id maxBound) (Id maxBound)]
      ]
  ]

hermitianAxes :: [QOp]
hermitianAxes =
  [ Id 0, Id 2, Phase (-3), X, Y, Z, H, Compose SX SX
  , Tensor (Phase (-1)) (Tensor X Y)
  , Adjoint (Tensor H Z), C H, DirectSum X Z
  , Permute [], Permute [1,0], Permute [2,1,0]
  , R X 2, R Y (-4), R Z 0
  , Compose X X, Compose (Id 1) H, Compose H (Id 1)
  , Compose (R X (1 % 3)) (Adjoint (R X (1 % 3)))
  , Compose (Adjoint SX) SX
  , Compose H (Compose Z H), Compose (Compose H Z) H
  , Compose SX (Compose Z (Adjoint SX))
  , Compose (Compose SX Z) (Adjoint SX)
  , Compose (Tensor X (Id 1)) (Tensor (Id 1) Z)
  , Compose (Tensor X X) (Tensor Z Z)
  , Tensor (Tensor (Phase (1 % 3)) H) (Tensor (Phase (2 % 3)) Z)
  , Tensor (DirectSum (Tensor (Phase (1 % 2)) X) (Tensor (Phase (-1 % 2)) Z))
           (Tensor (Phase (1 % 2)) H)
  , Compose (Tensor (Phase (1 % 2)) H) (Tensor (Phase (1 % 2)) H)
  , Compose (R H (1 % 3)) (R H (5 % 3))
  , Compose X (R Z (1 % 3)), Compose (R Z (1 % 3)) X
  , Compose (Tensor X X) (R (Tensor Z (Id 1)) (2 % 7))
  , DirectSum (Compose (Tensor X X) (Tensor Z Z)) (Tensor H Y)
  ]

nonHermitianAxes :: [QOp]
nonHermitianAxes =
  [ SX, Phase (1 % 2), R Z 1, R Z (1 % 2)
  , Permute [1,2,0], Compose X Z
  , Tensor (Phase (1 % 2)) X, C SX, DirectSum X SX, Adjoint SX
  ]

rotationTests :: TestTree
rotationTests = testGroup "Rotation axes"
  [ testGroup "Proved Hermitian"
      [ testCase (show axis) $ do
          typeCheck (R axis (1 % 3)) @?= Right (op_qubits axis)
          let mat = MS.evalOp axis
          assertClose "axis equals its adjoint" mat (adj mat)
      | axis <- hermitianAxes
      ]
  , testGroup "Non-Hermitian axes are rejected at nonzero angles"
      [ testCase (show (R axis q)) $
          typeCheck (R axis q) @?= Left (UnprovenHermitian axis)
      | axis <- nonHermitianAxes
      , q <- [1 % 3, 2, -4] -- Only literal zero is exempt, not a nonzero multiple of a period.
      ]
  , testGroup "Zero angle accepts axes without a Hermiticity proof"
      [ testCase (show axis) $ do
          let op = R axis 0
              n = op_qubits axis
          typeCheck op @?= Right n
          adjointPhase op @?= Right (Just 0)
          assertClose "zero-angle rotation is identity" (MS.evalOp op) (identCMat (2^n))
      | axis <- nonHermitianAxes ++ [Compose (Tensor H H) (Tensor Y Y)]
      ]
  , testCase "The proof remains conservative outside the Pauli fragment" $ do
      let axis = Compose (Tensor H H) (Tensor Y Y)
          mat = MS.evalOp axis
      assertClose "this axis is Hermitian" mat (adj mat)
      typeCheck (R axis (1 % 3)) @?= Left (UnprovenHermitian axis)
  ]

phaseTests :: TestTree
phaseTests = testGroup "Adjoint phases"
  [ testGroup "Exact phases, checked against matrix semantics"
      [ testCase (show op) $ do
          adjointPhase op @?= Right (Just p)
          assertClose "A dagger = exp(i*pi*p) A" (adj (MS.evalOp op))
            (MS.evalOp (Tensor (Phase p) op))
      | (op, p) <-
          [ (Phase (-7 % 3), 2 % 3)
          , (Tensor (Phase (1 % 3)) X, 4 % 3)
          , (Adjoint (Tensor (Phase (1 % 3)) H), 2 % 3)
          , (Compose X Y, 1), (Compose Y X, 1), (Compose X Z, 1)
          , (Compose (Tensor X X) (Tensor Z Z), 0)
          , (DirectSum (Tensor (Phase (1 % 3)) H)
                       (Tensor (Phase (4 % 3)) Z), 4 % 3)
          , (DirectSum (Tensor (Phase (1 % 2)) X)
                       (Tensor (Phase (-1 % 2)) Z), 1)
          , (R H 1, 1), (R H (-3), 1), (R H 2, 0)
          , (R SX 0, 0), (R (R SX 0) (1 % 3), 1 % 3)
          , (Compose (R SX 0) (R SX 0), 0)
          , (R (Id 2) (1 % 3), 1 % 3)
          , (R (Tensor (Phase 1) (Id 2)) (1 % 3), 5 % 3)
          , (R (R H 2) (1 % 3), 5 % 3)
          , (R (R H 4) (1 % 3), 1 % 3)
          , (R (Compose X (Compose Z (Compose X Z))) (1 % 3), 5 % 3)
          , (Compose (R H (1 % 3)) (R H (2 % 3)), 1)
          , (Compose SX (Compose (Tensor (Phase (1 % 3)) Z) (Adjoint SX)), 4 % 3)
          , (Compose (Compose SX (Tensor (Phase (1 % 3)) Z)) (Adjoint SX), 4 % 3)
          , (C (Phase 1), 0), (C (Phase 2), 0)
          , (DirectSum X (Tensor (Phase 1) X), 0)
          , (Compose (DirectSum X (Tensor (Phase 1) X)) (Tensor X X), 1)
          ]
      ]
  , testGroup "Anticommuting Pauli operators reverse rotations"
      [ testCase (show op) $ do
          adjointPhase op @?= Right (Just (2 - 2 * s))
          assertClose "reflection phase" (adj (MS.evalOp op))
            (MS.evalOp (Tensor (Phase (-2 * s)) op))
      | a <- [X,Y,Z], k <- [X,Y,Z], a /= k
      , t <- [1 % 3, -3 % 7]
      , s <- [1 % 5, 1 % 2]
      , let labelled = Tensor (Phase s) a
      , op <- [Compose labelled (R k t), Compose (R k t) labelled]
      ]
  , testGroup "Unknown phases are not proofs"
      [ testCase (show op) $ adjointPhase op @?= Right Nothing
      | op <- [ SX, R H (1 % 3), C (Tensor (Phase (1 % 2)) X)
              , DirectSum (Tensor (Phase (1 % 3)) H) (Tensor (Phase (2 % 3)) Z)
              , Compose X (R X (1 % 3))
              ]
      ]
  , testCase "Malformed subterms cannot be hidden by cancellation" $
      let bad = R SX (1 % 3)
      in adjointPhase (Compose bad (Adjoint bad)) @?= Left (UnprovenHermitian SX)
  , testCase "Phase queries validate widths" $
      adjointPhase (Tensor (Id (-1)) X) @?= Left (NegativeQubitCount (-1))
  , testCase "Sparse forms do not allocate one entry per identity qubit" $ do
      let x = Tensor (Id (maxBound - 1)) X
          z = Tensor (Id (maxBound - 1)) Z
      adjointPhase (Compose x z) @?= Right (Just 1)
      typeCheck (R (Compose x x) (1 % 3)) @?= Right maxBound
  ]

-- The complete fragment, with arbitrary association and scalar factors.
-- Eighth-integer phases keep the matrix oracle well away from roundoff
-- boundaries: a non-Hermitian value cannot be arbitrarily close to Hermitian.
genPauliExpr :: Int -> Int -> Gen QOp
genPauliExpr n depth
  | depth <= 0 = oneof [genPauli n, Tensor . Phase <$> angle <*> genPauli n]
  | otherwise = oneof
      [ genPauliExpr n 0
      , Compose <$> sub n <*> sub n
      , Adjoint <$> sub n
      , Tensor . Phase <$> angle <*> sub n
      , do k <- chooseInt (0, n)
           Tensor <$> sub k <*> sub (n - k)
      ]
  where
    sub k = genPauliExpr k (depth - 1)
    angle = (% 8) <$> chooseInteger (-16, 16)

prop_pauliComplete :: Property
prop_pauliComplete = forAll (chooseInt (0, 4) >>= \n -> genPauliExpr n 4) $ \op ->
  let mat = MS.evalOp op
      isHermitian = mDiff mat (adj mat) < tol
      result = typeCheck (R op (1 % 3))
  in case adjointPhase op of
    Right (Just p) -> conjoin
      [ counterexample "Hermiticity decision is complete" $
          (result == Right (op_qubits op)) === isHermitian
      , adj mat ~~ MS.evalOp (Tensor (Phase p) op)
      ]
    other -> counterexample ("missing exact Pauli phase: " ++ show other) False

prop_phaseSound :: RandomQOp -> Property
prop_phaseSound (RandomQOp _ op) = case adjointPhase op of
  Right (Just p) -> classify True "proved phase" $
    conjoin [property (0 <= p && p < 2), adj (MS.evalOp op) ~~ MS.evalOp (Tensor (Phase p) op)]
  Right Nothing -> classify False "proved phase" $ property True
  Left err -> counterexample (show err) False

prop_conjugationPhase :: RandomQOp -> Property
prop_conjugationPhase (RandomQOp n u) = forAll ((,) <$> genRat <*> genPauli n) $ \(q, p) ->
  let a = Tensor (Phase q) p
      expected = adjointPhase a
  in conjoin
      [ adjointPhase (Compose u (Compose a (Adjoint u))) === expected
      , adjointPhase (Compose (Compose u a) (Adjoint u)) === expected
      ]

prop_generatedOperators :: RandomQOp -> Property
prop_generatedOperators (RandomQOp n op) = typeCheck op === Right n

prop_zeroRotation :: RandomQOp -> Property
prop_zeroRotation (RandomQOp n axis) =
  let op = R axis 0
  in conjoin
      [ typeCheck op === Right n
      , adjointPhase op === Right (Just 0)
      , MS.evalOp op ~~ identCMat (2^n)
      ]

-- Use an independent matrix oracle to check the soundness of accepted axes.
prop_rotationSound :: RandomQOp -> Property
prop_rotationSound (RandomQOp n axis) = case typeCheck (R axis (1 % 3)) of
  Left _ -> classify False "accepted axis" $ property True
  Right m -> classify True "accepted axis" $
    let a = MS.evalOp axis
        r = MS.evalOp (R axis (1 % 3))
    in conjoin [m === n, a ~~ adj a, (r <> adj r) ~~ identCMat (2^n)]

main :: IO ()
main = defaultMain $ localOption (QuickCheckTests 100) $
  testGroup "TypeCheckTest"
    [ arityTests, errorTests, rotationTests, phaseTests
    , testProperty "Pauli fragment is complete, with exact adjoint phases" prop_pauliComplete
    , testProperty "Every inferred adjoint phase agrees with matrix semantics" prop_phaseSound
    , testProperty "Conjugation preserves arbitrary adjoint phases" prop_conjugationPhase
    , testProperty "Existing well-formed generator is accepted" prop_generatedOperators
    , testProperty "Zero-angle rotations of well-formed axes are identity" prop_zeroRotation
    , testProperty "Accepted rotation axes give Hermitian matrices and unitary rotations"
        prop_rotationSound
    ]
