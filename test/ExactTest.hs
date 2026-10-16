module Main where

import Control.Monad (forM_, replicateM)
import Data.Bits (testBit)
import Data.Either (isLeft)
import Data.Ratio ((%), numerator, denominator)
import qualified QPP.Exact.Scalar as S
import qualified QPP.Exact.Checked as Checked
import qualified QPP.Exact.State as St
import QPP.Exact.Types
import QPP.Syntax (QOp(..), Step(..))
import QPP.Rewrite (cleanOnes)
import qualified QPP.Semantics.Exact.Matrix as M
import qualified QPP.Semantics.Exact.Statevector as V
import Programs.QFT (qft)
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck

l :: Limits
l = defaultLimits

-- Fixture failures are test failures, never an arithmetic fallback.
fixture :: Either ExactError a -> a
fixture = either (error . show) id

right :: Either ExactError a -> IO a
right = either (\e -> assertFailure (show e) >> fail "unreachable") pure

rat :: Rational -> S.Cyclo
rat = fixture . S.rational l

root :: Integer -> S.Cyclo
root = fixture . S.rootOfUnity l

plus, times :: S.Cyclo -> S.Cyclo -> S.Cyclo
plus a b = fixture (S.add l a b)
times a b = fixture (S.mul l a b)

negative :: S.Cyclo -> S.Cyclo
negative = fixture . S.negate l

power :: S.Cyclo -> Integer -> S.Cyclo
power a = fixture . S.pow l a

h, imaginary :: S.Cyclo
h = fixture (S.fromCoordinates l 8 [0,1%2,0,(-1)%2])
imaginary = root 4

assertScalar :: S.Cyclo -> S.Cyclo -> Assertion
assertScalar expected actual = do
  same <- right (S.equal l expected actual)
  assertBool ("expected " ++ show expected ++ ", got " ++ show actual) same

assertValues :: [S.Cyclo] -> St.StateT -> Assertion
assertValues expected actual = do
  length (St.amplitudes actual) @?= length expected
  sequence_ (zipWith assertScalar expected (St.amplitudes actual))

assertMatrix :: [[S.Cyclo]] -> M.OpT -> Assertion
assertMatrix expected actual = do
  map length (M.entries actual) @?= map length expected
  sequence_ (zipWith assertScalar (concat expected) (concat (M.entries actual)))

assertError :: (ExactError -> Bool) -> Either ExactError a -> Assertion
assertError predicate result = case result of
  Left e -> assertBool ("unexpected error: " ++ show e) (predicate e)
  Right _ -> assertFailure "expected an explicit error"

isResource :: ExactError -> Bool
isResource (ResourceLimit _ _ _) = True
isResource _ = False

runDirect :: QOp -> [S.Cyclo] -> Either ExactError St.StateT
runDirect op xs = do
  p <- V.evalOp l op
  v <- St.fromAmplitudes l xs
  V.apply l p v

scalarTests :: TestTree
scalarTests = testGroup "Cyclotomic arithmetic"
  [ testGroup "Known cyclotomic polynomials"
      [ testCase (show n) $ do
          x <- right (S.rootOfUnity l n)
          S.cyclotomicPolynomial x @?= cs
          assertScalar (rat 1) (power x n)
          forM_ [d | d <- [1..n-1], n `mod` d == 0] $ \d ->
            S.equal l (power x d) (rat 1) @?= Right False
      | (n,cs) <-
          [ (1,[-1,1]), (2,[1,1]), (3,[1,1,1]), (4,[1,0,1])
          , (5,[1,1,1,1,1]), (6,[1,-1,1]), (8,[1,0,0,0,1])
          , (9,[1,0,0,1,0,0,1]), (12,[1,0,-1,0,1])
          , (15,[1,-1,0,1,-1,1,0,-1,1])
          , (24,[1,0,0,0,-1,0,0,0,1])
          , (30,[1,1,0,-1,-1,-1,0,1,1]) ]
      ]
  , testCase "Hadamard coordinates and field embedding" $ do
      assertScalar (rat (1%2)) (times h h)
      promoted <- right (S.embed l 24 h)
      S.coefficients promoted @?= [0,1%2,0,1%2,0,(-1)%2,0,0]
      assertScalar h promoted
      assertScalar (root 6) (plus (rat 1) (root 3))
  , testCase "Polynomial inversion and negative powers" $ do
      actual <- right (S.inverse l (plus (rat 1) imaginary))
      assertScalar (times (rat (1%2)) (plus (rat 1) (negative imaginary))) actual
      assertScalar (times (rat 2) h) (power h (-1))
      assertScalar (rat 1) (power (rat 0) 0)
      assertError (== DivisionByZero) (S.inverse l (rat 0))
      assertError (== DivisionByZero) (S.divide l h (rat 0))
  , testCase "Phase and rotation sign conventions" $ do
      right (S.phase l (1%2)) >>= assertScalar imaginary
      right (S.phase l ((-1)%2)) >>= assertScalar (negative imaginary)
      (c,d) <- right (S.rotationCoefficients l 1)
      assertScalar (rat 0) c
      assertScalar (negative imaginary) d
      (c2,d2) <- right (S.rotationCoefficients l 2)
      assertScalar (rat (-1)) c2
      assertScalar (rat 0) d2
  , testCase "Equal values in different contexts" $ do
      zero24 <- right (S.fromCoordinates l 24 (replicate 8 0))
      assertBool "exact zero" (S.isZero zero24)
      assertScalar (rat 0) zero24
      S.asRational zero24 @?= Just 0
      assertScalar (rat (-1)) (root 2)
      S.asRational imaginary @?= Nothing
  ]

gateTests :: TestTree
gateTests = testGroup "Independent gate values"
  [ testGroup "Primitive matrices"
      [ testCase (show op) $ right (M.evalOp l op) >>= assertMatrix expected
      | (op,expected) <-
          [ (X,[[rat 0,rat 1],[rat 1,rat 0]])
          , (Y,[[rat 0,negative imaginary],[imaginary,rat 0]])
          , (Z,[[rat 1,rat 0],[rat 0,rat (-1)]])
          , (H,[[h,h],[h,negative h]])
          , (SX,let a = times (rat (1%2)) (plus (rat 1) imaginary)
                    b = times (rat (1%2)) (plus (rat 1) (negative imaginary))
                in [[a,b],[b,a]]) ]
      ]
  , testCase "R Z 1 = -i Z and controlled identity rotation = Z" $ do
      right (M.evalOp l (R Z 1)) >>= assertMatrix
        [[negative imaginary,rat 0],[rat 0,imaginary]]
      right (M.evalOp l (C (R (Id 0) 2))) >>= assertMatrix
        [[rat 1,rat 0],[rat 0,rat (-1)]]
      right (runDirect (C (R (Id 0) 2)) [rat 2,rat 3]) >>= assertValues [rat 2,rat (-3)]
  , testCase "S and T preserve their specified global phases" $ do
      right (M.evalOp l (C (Phase (1%2)))) >>= assertMatrix
        [[rat 1,rat 0],[rat 0,imaginary]]
      right (M.evalOp l (Tensor (Phase (1%8)) (R Z (1%4)))) >>= assertMatrix
        [[rat 1,rat 0],[rat 0,root 8]]
  , testCase "Zero-qubit scalar evolution on an unnormalized value" $ do
      let input = plus (rat 2) imaginary
          expected = plus (rat (-1)) (times (rat 2) imaginary)
      v <- right (runDirect (Phase (1%2)) [input])
      assertValues [expected] v
      right (St.normSquared l v) >>= assertScalar (rat 5)
      p <- right (M.evalOp l (Phase (1%2)))
      s <- right (St.fromAmplitudes l [input])
      right (M.apply l p s) >>= assertValues [expected]
      right (M.evalOp l (Permute [])) >>= assertMatrix [[rat 1]]
  , testCase "Zero-angle non-Hermitian axis follows the committed checker" $ do
      right (M.evalOp l (R SX 0)) >>= assertMatrix [[rat 1,rat 0],[rat 0,rat 1]]
      right (runDirect (R SX 0) [rat 2,imaginary]) >>= assertValues [rat 2,imaginary]
  , testCase "Hadamard-axis rotation has the paper's exact coordinates" $ do
      let a = fixture (S.fromCoordinates l 24 [0,1%4,1,(-1)%4,0,(-1)%4,(-1)%2,0])
          b = fixture (S.fromCoordinates l 24 [0,1%4,0,(-1)%4,0,(-1)%4,0,0])
      v <- right (runDirect (R H (1%3)) [rat 1,rat 0])
      assertValues [a,b] v
      right (St.normSquared l v) >>= assertScalar (rat 1)
      forM_ (zip (St.amplitudes v) [7%8,1%8]) $ \(x,p) -> do
        xc <- right (S.conjugate l x)
        right (S.mul l xc x) >>= assertScalar (rat p)
      m <- right (M.evalOp l (R H (1%3)))
      right (M.entry m 0 0) >>= assertScalar a
      right (M.entry m 1 0) >>= assertScalar b
  , testCase "Phased Bell state and interference probabilities" $ do
      let bell = Compose (C X) (Tensor H (Id 1))
          op = Compose (Tensor (Id 1) (C (Phase (1%3)))) bell
          z = root 6
      v <- right (runDirect op [rat 1,rat 0,rat 0,rat 0])
      assertValues [h,rat 0,rat 0,times h z] v
      p <- right (V.evalOp l (Adjoint bell))
      w <- right (V.apply l p v)
      assertValues [times (rat (1%2)) (plus (rat 1) z),rat 0,
                    times (rat (1%2)) (plus (rat 1) (negative z)),rat 0] w
      forM_ (zip (St.amplitudes w) [3%4,0,1%4,0]) $ \(x,prob) -> do
        xc <- right (S.conjugate l x)
        right (S.mul l xc x) >>= assertScalar (rat prob)
  , testCase "Complex tensor row action uses transpose, not adjoint" $ do
      let z = root 8
      right (runDirect (Tensor H (C (Phase (1%4)))) (map rat [1,2,3,5]))
        >>= assertValues [times h (rat 4),times (times h z) (rat 7),
                          times h (rat (-2)),times (times h z) (rat (-3))]
  , testCase "Non-self-inverse permutation scatter and dagger" $ do
      let ks = [1,2,0]
      p <- right (M.evalOp l (Permute ks))
      right (M.entry p 1 4) >>= assertScalar (rat 1)
      v <- right (St.ket l [True,False,False])
      expected <- right (St.ket l [False,False,True])
      op <- right (V.evalOp l (Permute ks))
      actual <- right (V.apply l op v)
      St.equal l expected actual @?= Right True
      inv <- right (V.evalOp l (Adjoint (Permute ks)))
      restored <- right (V.apply l inv actual)
      St.equal l v restored @?= Right True
  ]

linearTests :: TestTree
linearTests = testGroup "Exact linear algebra and programs"
  [ testCase "Linear combinations and conjugate-linear first argument" $ do
      v <- right (St.fromAmplitudes l [rat 1,imaginary])
      w <- right (St.fromAmplitudes l [imaginary,rat 2])
      a <- right (St.scale l (root 6) v)
      b <- right (St.scale l h w)
      c <- right (St.add l a b)
      assertValues [plus (root 6) (times h imaginary),
                    plus (times (root 6) imaginary) (times (rat 2) h)] c
      x <- right (St.inner l a w)
      y <- right (St.inner l v w)
      zc <- right (S.conjugate l (root 6))
      assertScalar (times zc y) x
      right (St.sub l c b) >>= \a' -> St.equal l a a' @?= Right True
  , testCase "Squared norm needs no noncyclotomic square root" $ do
      let sqrtTwo = times (rat 2) h
      v <- right (St.fromAmplitudes l [rat 1,plus (rat 2) sqrtTwo])
      right (St.normSquared l v) >>= assertScalar (plus (rat 7) (times (rat 4) sqrtTwo))
  , testCase "Matrix linear algebra is Frobenius, without a unitarity assertion" $ do
      a <- right (M.evalOp l H)
      right (M.normSquaredOp l a) >>= assertScalar (rat 2)
      b <- right (M.scaleOp l (rat 2) a)
      c <- right (M.addOp l a a)
      M.equalOp l b c @?= Right True
      d <- right (M.composeOp l a a)
      right (M.evalOp l (Id 1)) >>= \ident -> M.equalOp l d ident @?= Right True
      right (M.fromRows l (M.entries b)) >>= \copy -> M.equalOp l b copy @?= Right True
  , testCase "Zero-qubit and ordinary state tensor products" $ do
      s <- right (St.fromAmplitudes l [imaginary])
      v <- right (St.fromAmplitudes l [rat 2,rat 3])
      right (St.tensor l s v) >>= assertValues [times imaginary (rat 2),times imaginary (rat 3)]
      right (St.tensor l v s) >>= assertValues [times imaginary (rat 2),times imaginary (rat 3)]
      right (St.tensor l s s) >>= assertValues [rat (-1)]
  , testCase "Program order and common-field preflight" $ do
      initial <- right (St.ket l [False])
      let steps = [Unitary H, Unitary Z, Unitary (C (Phase (1%3)))]
          expected = [h,negative (times h (root 6))]
      forM_ [M.evalProg,V.evalProg] $ \run -> right (run l steps initial) >>= assertValues expected
      v <- right (V.evalProg l steps initial)
      St.fieldIndex v @?= 24
  , testCase "An empty program preserves the input context" $ do
      v <- right (St.ket l [])
      let tight = l {maxRootIndex=1}
      forM_ [M.evalProg,V.evalProg] $ \run -> do
        x <- right (run tight [] v)
        St.fieldIndex x @?= 1
        assertValues [rat 1] x
  , testCase "External amplitude field promotes the operator once" $ do
      v <- right (St.fromAmplitudes l [root 5,rat 0])
      p <- right (V.evalOp l H)
      m <- right (M.evalOp l H)
      x <- right (V.apply l p v)
      y <- right (M.apply l m v)
      St.fieldIndex x @?= 40
      St.equal l x y @?= Right True
      assertValues [times h (root 5),times h (root 5)] x
  , testCase "Direct evaluation fits when a dense matrix cannot" $ do
      let limited = l {maxCoefficientSlots=256000}
          op = Tensor H (Id 9)
      p <- right (V.evalOp limited op)
      s <- right (St.ket limited (replicate 10 False))
      v <- right (V.apply limited p s)
      right (St.amplitude v 0) >>= assertScalar h
      right (St.amplitude v 512) >>= assertScalar h
      right (St.normSquared limited v) >>= assertScalar (rat 1)
      assertError isResource (M.evalOp limited op)
  ]

errorTests :: TestTree
errorTests = testGroup "Explicit failures"
  [ testCase "Invalid syntax and non-Hermitian axes" $ forM_
      [Id (-1),Permute [0,0],Compose X (Id 2),R SX (1%3),R (Id (-1)) 0,
       R (Compose X (Id 2)) 0,Tensor (Id maxBound) X] $ \op -> do
      assertBool "dense rejected" (isLeft (M.evalOp l op))
      assertBool "direct rejected" (isLeft (V.evalOp l op))
  , testCase "Huge widths are rejected before shifting or allocating" $ do
      assertError isResource (M.evalOp l (Id maxBound))
      assertError isResource (V.evalOp l (Id maxBound))
  , testCase "Field-degree explosion is rejected before coordinate allocation" $ do
      let limits = l {maxRootIndex=100000}
          op = Compose (R Z (1%101)) (R Z (1%103))
          degreeError (ResourceLimit "field degree" 40800 256) = True
          degreeError _ = False
      assertError degreeError (M.evalOp limits op)
      assertError degreeError (V.evalOp limits op)
  , testCase "Import validation, field limits, and division by zero" $ do
      assertBool "empty state" (isLeft (St.fromAmplitudes l []))
      assertBool "non-power-of-two state" (isLeft (St.fromAmplitudes l (replicate 3 (rat 1))))
      assertBool "bad coordinates" (isLeft (S.fromCoordinates l 8 [1]))
      assertBool "bad root" (isLeft (S.rootOfUnity l 0))
      assertBool "invalid embedding" (isLeft (S.embed l 8 (root 3)))
      assertBool "bad matrix" (isLeft (M.fromRows l [[rat 1,rat 0]]))
      assertBool "bad limits" (isLeft (S.rootOfUnity (l {maxFieldDegree=0}) 8))
      assertBool "empty preflight still checks width" (isLeft (Checked.checkProgram l (-1) []))
      assertBool "empty preflight still checks limits"
        (isLeft (Checked.checkProgram (l {maxFieldDegree=0}) 0 []))
      assertBool "invalid starting field" (isLeft (Checked.programFieldIndex l 0 []))
      assertError isResource (Checked.programFieldIndex l 8192 [])
  , testCase "Tighter limits recheck previously constructed inputs" $ do
      let tight = l {maxIntegerBits=8}
      assertError isResource (S.mul tight (rat 128) (rat 128))
      assertError isResource (S.add tight (rat 256) (rat 0))
      v <- right (St.fromAmplitudes l [rat 256])
      p <- right (V.evalOp l (Id 0))
      assertError isResource (V.apply tight p v)
      assertError isResource (St.scale tight (rat 1) v)
      assertError isResource (S.equal (l {maxRootIndex=8}) (root 3) (root 8))
  , testCase "Program preflight reports step and rejects measurement/reset" $ do
      v <- right (St.ket l [False])
      forM_ [M.evalProg,V.evalProg] $ \run -> do
        assertError (== UnsupportedStep 1 "Measure") (run l [Unitary H,Measure []] v)
        assertError (== UnsupportedStep 0 "Initialize") (run l [Initialize [] []] v)
        assertError (== AtStep 1 (WidthMismatch 2 1)) (run l [Unitary H,Unitary (C X)] v)
  , testCase "Shape checks and observation indices" $ do
      v <- right (St.ket l [False])
      w <- right (St.ket l [])
      assertError (== WidthMismatch 1 0) (St.add l v w)
      p <- right (V.evalOp l (Id 0))
      assertError (== WidthMismatch 0 1) (V.apply l p v)
      assertBool "amplitude index" (isLeft (St.amplitude v 2))
      m <- right (M.evalOp l H)
      assertBool "matrix index" (isLeft (M.entry m (-1) 0))
  ]

newtype SmallCyclo = SmallCyclo S.Cyclo deriving (Show)
instance Arbitrary SmallCyclo where
  arbitrary = do
    -- Every index divides 120, keeping triple common fields bounded while
    -- exercising both odd primes 3 and 5 and nontrivial embeddings.
    n <- elements [1,2,3,4,5,6,8,10,12,15,20,24,30,40,60,120]
    SmallCyclo <$> cycloIn n
  shrink (SmallCyclo a) = map SmallCyclo (shrinkCyclo a)

cycloIn :: Integer -> Gen S.Cyclo
cycloIn n = do
  let d = S.fieldDegree (root n)
  cs <- replicateM d ((%) <$> chooseInteger (-3,3) <*> chooseInteger (1,4))
  pure (fixture (S.fromCoordinates l n cs))

shrinkAngle :: Rational -> [Rational]
shrinkAngle q = [n % denominator q | n <- shrink (numerator q)] ++
  [numerator q % d | d <- shrink (denominator q), d > 0]

-- Dropping coefficients pads with zeros to preserve the field's basis length.
-- Exclude unchanged padded values so repeated shrinking always makes progress.
shrinkCyclo :: S.Cyclo -> [S.Cyclo]
shrinkCyclo a =
  [fixture (S.fromCoordinates l (S.fieldIndex a) ys)
  | zs <- shrinkList shrinkAngle xs
  , let ys = take (length xs) (zs ++ repeat 0)
  , ys /= xs]
  where xs = S.coefficients a

newtype LargerCyclo = LargerCyclo S.Cyclo deriving (Show)
instance Arbitrary LargerCyclo where
  arbitrary = LargerCyclo <$> (elements [35,45,60,84,105,120,168,210] >>= cycloIn)
  shrink (LargerCyclo a) = map LargerCyclo (shrinkCyclo a)

scalarProperties :: TestTree
scalarProperties = localOption (QuickCheckTests 60) $ testGroup "Field laws across contexts"
  [ testProperty "distributivity" $ \(SmallCyclo a) (SmallCyclo b) (SmallCyclo c) ->
      S.equal l (times a (plus b c)) (plus (times a b) (times a c)) === Right True
  , testProperty "inverse and conjugation" $ \(SmallCyclo a) ->
      let ac = fixture (S.conjugate l a)
          acc = fixture (S.conjugate l ac)
      in S.equal l a acc === Right True .&&.
         (if S.isZero a then property True else S.equal l (times a (power a (-1))) (rat 1) === Right True)
  , testProperty "normalized inversion in larger fields" $ \(LargerCyclo a) ->
      if S.isZero a then property True else
        let result = S.inverse l a >>= S.mul l a >>= S.equal l (rat 1)
        in counterexample (show result) (result === Right True)
  , testProperty "embedding composition and homomorphism" $
      forAll (elements [(3,15,120),(5,10,120),(8,24,120),(7,35,210),(9,45,90)]) $ \(n,m,k) ->
      forAllShrink (cycloIn n) shrinkCyclo $ \a ->
      forAllShrink (cycloIn n) shrinkCyclo $ \b ->
        let result = do
              direct <- S.embed l k a
              via <- S.embed l m a >>= S.embed l k
              eb <- S.embed l k b
              esum <- S.add l a b >>= S.embed l k
              eproduct <- S.mul l a b >>= S.embed l k
              sum' <- S.add l direct eb
              product' <- S.mul l direct eb
              and <$> sequence [S.equal l direct via, S.equal l esum sum', S.equal l eproduct product']
        in counterexample (show result) (result === Right True)
  ]

angle :: Gen Rational
angle = (%) <$> chooseInteger (-4,4) <*> elements [1,2,3,4]

axis :: Int -> Gen QOp
axis 0 = elements [Id 0,Phase 1]
axis 1 = elements [X,Y,Z,H,Compose H (Compose Z H),Compose X (R Z (1%3))]
axis n = Tensor <$> axis 1 <*> axis (n-1)

term :: Int -> Int -> Gen QOp
term n depth
  | depth <= 0 = leaf
  | otherwise = frequency $
      [ (3,leaf)
      , (2,Compose <$> sub n <*> sub n)
      , (1,Adjoint <$> sub n)
      , (2,do k <- chooseInt (0,n); Tensor <$> sub k <*> sub (n-k))
      ] ++ if n == 0 then [] else
      [(1,C <$> sub (n-1)),(1,DirectSum <$> sub (n-1) <*> sub (n-1))]
  where
    sub k = term k (depth-1)
    leaf = oneof $
      [pure (Id n), Permute <$> shuffle [0..n-1], R <$> axis n <*> angle] ++
      if n == 0 then [Phase <$> angle] else
      if n == 1 then map pure [X,Y,Z,H,SX] else []

data CircuitCase = CircuitCase Int QOp [S.Cyclo] deriving (Show)
instance Arbitrary CircuitCase where
  arbitrary = do
    n <- chooseInt (0,3)
    d <- chooseInt (0,3)
    op <- term n d
    cs <- replicateM (2^n) $ do
      a <- chooseInteger (-2,2)
      b <- chooseInteger (-2,2)
      pure (plus (rat (fromInteger a)) (times imaginary (rat (fromInteger b))))
    pure (CircuitCase n op cs)
  shrink (CircuitCase n op cs) =
    [CircuitCase n smaller cs | smaller <- shrinkTerm n op] ++
    [CircuitCase n op ys
    | zs <- shrinkList shrinkCyclo cs
    , let ys = take (length cs) (zs ++ repeat (rat 0))
    , map S.coefficients ys /= map S.coefficients cs]

-- Preserve width and well-formedness, including the Hermitian-axis obligation.
shrinkTerm :: Int -> QOp -> [QOp]
shrinkTerm n op = [Id n | op /= Id n] ++ case op of
  Phase q -> map Phase (shrinkAngle q)
  R a q -> [R a q' | q' <- shrinkAngle q]
  C a -> map C (shrinkTerm (n-1) a)
  Compose a b -> [a,b] ++ [Compose a' b | a' <- shrinkTerm n a] ++ [Compose a b' | b' <- shrinkTerm n b]
  Adjoint a -> a : map Adjoint (shrinkTerm n a)
  Tensor a b ->
    let na = Checked.qubits (fixture (Checked.checkOp l a))
    in [Tensor a' b | a' <- shrinkTerm na a] ++ [Tensor a b' | b' <- shrinkTerm (n-na) b]
  DirectSum a b -> [DirectSum a' b | a' <- shrinkTerm (n-1) a] ++ [DirectSum a b' | b' <- shrinkTerm (n-1) b]
  _ -> []

reviewRegressions :: TestTree
reviewRegressions = testGroup "Review regressions"
  [ testCase "Cyclotomic product identity over integer polynomials, n=1..120" $
      forM_ [1..120] $ \n -> do
        let factors = [S.cyclotomicPolynomial (root d) | d <- [1..n], n `mod` d == 0]
            expected = [-1] ++ replicate (fromInteger n-1) 0 ++ [1]
        foldl multiplyPolynomials [1] factors @?= expected
  , testCase "Degree-48 inverses fit the default bit limit" $
      forM_ [105,168,210] $ \n -> do
        let d = S.fieldDegree (root n)
            a = fixture (S.fromCoordinates l n [((toInteger j `mod` 7)-3) % 3 | j <- [0..d-1]])
        ai <- right (S.inverse l a)
        right (S.mul l a ai) >>= assertScalar (rat 1)
        assertBool "result coordinates remain small"
          (all (\q -> abs (numerator q) < 2^(140::Int) && denominator q < 2^(140::Int)) (S.coefficients ai))
  , testCase "Matrix operations report actual qubit widths" $ do
      a <- right (M.evalOp l (Id 2))
      b <- right (M.evalOp l (Id 3))
      assertError (== WidthMismatch 2 3) (M.innerOp l a b)
      assertError (== WidthMismatch 2 3) (M.addOp l a b)
      assertError (== WidthMismatch 2 3) (M.composeOp l a b)
  , testCase "An empty program needs only the input's stored coordinates" $
      forM_ [3,257] $ \slots -> do
        let tiny = l {maxCoefficientSlots=slots}
        v <- right (St.ket tiny [])
        St.equal tiny v v @?= Right True
        a <- right (S.rational tiny 1)
        S.equal tiny a a @?= Right True
        forM_ [M.evalProg,V.evalProg] $ \run -> do
          w <- right (run tiny [] v)
          assertValues [rat 1] w
        assertError isResource (V.evalProg (tiny {maxCoefficientSlots=2}) [] v)
  , testCase "Linear arithmetic does not reserve polynomial workspace" $ do
      let tiny = l {maxCoefficientSlots=12}
      a <- right (S.rational tiny 2)
      b <- right (S.rational tiny 3)
      right (S.add tiny a b) >>= assertScalar (rat 5)
      right (S.negate tiny a) >>= assertScalar (rat (-2))
      S.equal tiny a a @?= Right True
  , testCase "Existing field contexts do not require construction workspace" $ do
      let tight = l {maxCoefficientSlots=2000}
      a <- right (S.rootOfUnity l 120)
      v <- right (St.fromAmplitudes l [a])
      assertError isResource (S.rootOfUnity tight 120)
      right (St.scale tight (rat 1) v) >>= assertValues [a]
      right (St.fromAmplitudes tight [a]) >>= assertValues [a]
      right (S.embed tight 120 a) >>= assertScalar a
  , testCase "Prepared direct plans recheck constants under tighter limits" $ do
      p <- right (V.evalOp l H)
      v <- right (St.ket l [False])
      assertError isResource (V.apply (l {maxIntegerBits=1}) p v)
  , testCase "Program field-degree errors identify the introducing step" $ do
      let limits = l {maxRootIndex=100000,maxFieldDegree=512}
          program = [Unitary (R Z (1%101)),Unitary (R Z (1%103))]
      v <- right (St.ket l [False])
      forM_ [M.evalProg,V.evalProg] $ \run ->
        assertError (== AtStep 1 (ResourceLimit "field degree" 40800 512)) (run limits program v)
  ]
  where
    -- Independent integer convolution: no cyclotomic scalar operations.
    multiplyPolynomials :: [Integer] -> [Integer] -> [Integer]
    multiplyPolynomials as bs =
      [sum [a*b | (i,a) <- zip [0..] as, (j,b) <- zip [0..] bs, i+j == k]
      | k <- [0..length as+length bs-2]]

evaluatorProperties :: TestTree
evaluatorProperties = localOption (QuickCheckTests 60) $ testGroup "Exact evaluator properties"
  [ testProperty "shrinking preserves well-formed operators and input widths" $ \example@(CircuitCase n _ _) ->
      conjoin
        [ counterexample (show candidate) $
            (Checked.qubits <$> Checked.checkOp l op) === Right n .&&. length cs === 2^n
        | candidate@(CircuitCase _ op cs) <- take 50 (shrink example)]
  , testProperty "dense/direct agreement on arbitrary unnormalized vectors" $ \(CircuitCase _ op cs) ->
      let result = do
            v <- St.fromAmplitudes l cs
            m <- M.evalOp l op
            p <- V.evalOp l op
            a <- M.apply l m v
            b <- V.apply l p v
            St.equal l a b
      in counterexample (show result) (result === Right True)
  , testProperty "direct inverse and squared norm" $ \(CircuitCase _ op cs) ->
      let result = do
            v <- St.fromAmplitudes l cs
            p <- V.evalOp l op
            pi' <- V.evalOp l (Adjoint op)
            w <- V.apply l p v
            restored <- V.apply l pi' w
            a <- St.equal l v restored
            nv <- St.normSquared l v
            nw <- St.normSquared l w
            b <- S.equal l nv nw
            pure (a && b)
      in counterexample (show result) (result === Right True)
  ]

qftTests :: TestTree
qftTests = testGroup "QFT equals the analytic DFT, exactly"
  [ testCase (show n ++ " qubits") $ do
      let d = 2^n
          bitReverse j = sum [if testBit j k then 2^(n-1-k) else 0 | k <- [0..n-1]] :: Int
          omega = root (toInteger d)
          normalization = power h (toInteger n)
          expected = [[times normalization (power omega (toInteger (bitReverse i * bitReverse j)))
                      | j <- [0..d-1]] | i <- [0..d-1]]
      m <- right (M.evalOp l (qft n))
      assertMatrix expected m
      p <- right (V.evalOp l (qft n))
      forM_ [0..d-1] $ \j -> do
        v <- right (St.ket l [testBit j (n-1-k) | k <- [0..n-1]])
        actual <- right (V.apply l p v)
        assertValues [row !! j | row <- expected] actual
  | n <- [0..3]
  ]

constructorTests :: TestTree
constructorTests = testGroup "Every column of composite operators"
  [ testCase (show op) $ do
      m <- right (M.evalOp l op)
      p <- right (V.evalOp l op)
      let n = M.opQubits m
          d = 2^n
      forM_ [0..d-1] $ \j -> do
        input <- right (St.ket l [testBit j (n-1-k) | k <- [0..n-1]])
        expected <- traverse (\i -> right (M.entry m i j)) [0..d-1]
        right (V.apply l p input) >>= assertValues expected
      md <- right (M.adjointOp l m)
      fromSyntax <- right (M.evalOp l (Adjoint op))
      M.equalOp l md fromSyntax @?= Right True
  | op <-
      [ Tensor (Phase (1%3)) (Phase ((-1)%4))
      , Tensor (Phase (1%3)) (Tensor H SX)
      , Tensor (Tensor H SX) (Phase (1%3))
      , Tensor H (C X), Tensor (C X) H
      , DirectSum (Phase (1%3)) (Phase ((-1)%4))
      , DirectSum H (R Z (1%3)), C (C (Phase (1%3)))
      , R (C X) (1%3), R (DirectSum H Z) ((-1)%4)
      , R (Compose X (R Z (1%3))) (1%4)
      , Adjoint (Compose (R H (1%3)) SX)
      , Adjoint (Tensor (Phase (1%3)) SX)
      , Adjoint (DirectSum SX (R H (1%3)))
      , C (Permute [1,2,0]), Adjoint (C (Permute [1,2,0]))
      ]
  ]

rewriteTests :: TestTree
rewriteTests = testGroup "Identity-axis rewrite preserves exact phase"
  [ testCase (show op) $ do
      before <- right (M.evalOp l op)
      rewritten <- right (M.evalOp l (cleanOnes op))
      M.equalOp l before rewritten @?= Right True
  | op <- [R (Id 0) (1%3), C (R (Id 0) 2), Tensor H (R (Id 0) 2),
           DirectSum (R (Id 1) (1%3)) H, C (R (Tensor (Id 0) (Id 1)) 2)]
  ]

main :: IO ()
main = defaultMain (testGroup "ExactTest"
  [scalarTests, gateTests, linearTests, errorTests, scalarProperties, evaluatorProperties,
   constructorTests, rewriteTests, qftTests, reviewRegressions])
