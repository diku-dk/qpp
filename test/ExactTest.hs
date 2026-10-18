module Main where

import Control.Monad (forM_, replicateM)
import Control.Exception (ErrorCall, evaluate, try)
import Data.Bits (testBit)
import Data.Either (isLeft)
import Data.List (isInfixOf)
import Data.Ratio ((%), numerator, denominator)
import qualified QPP.Exact.Scalar as S
import qualified QPP.Exact.Checked as Checked
import qualified QPP.Exact.State as St
import QPP.Exact.Types
import QPP.Syntax (QOp(..), Step(..), Program, RNG, Outcomes, n_qubits)
import QPP.Rewrite (cleanOnes)
import qualified QPP.Semantics.Exact.Matrix as M
import qualified QPP.Semantics.Exact.Statevector as V
import qualified QPP.Semantics.Matrix as NumericalMatrix
import qualified QPP.Semantics.Statevector as NumericalVector
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
  p <- V.evalOpWithLimits l op
  v <- St.fromAmplitudes l xs
  V.applyWithLimits l p v

interfaceTests :: TestTree
interfaceTests = testGroup "Conventional backend interface"
  [ testCase "Default entry points give the exact Bell state" $ do
      let bell = Compose (C X) (Tensor H (Id 1))
          matrix = M.evalOp bell
          direct = V.evalOp bell
      n_qubits matrix @?= 2
      n_qubits direct @?= 2
      forM_ [M.apply matrix (M.ket [0,0]), V.apply direct (V.ket [0,0])] $ \state -> do
        n_qubits state @?= 2
        assertValues [h,rat 0,rat 0,h] state
  , testCase "Default program interface preserves the RNG for unitaries" $
      forM_ [M.evalProg,V.evalProg] $ \run -> do
        let rng = [0.2,0.8]
            (state, outcomes, remaining) = run [Unitary H,Unitary Z] (St.ket [0]) rng
        assertValues [h,negative h] state
        outcomes @?= []
        remaining @?= rng
  , testCase "Checked kets reject non-binary integers" $
      forM_ [[2],[-1],[0,1,2]] $ \bits ->
        assertError (== InvalidInput "ket expects bits equal to 0 or 1")
          (St.ketWithLimits l bits)
  , testCase "Default errors identify the operation" $ do
      result <- try (evaluate (St.ket [2])) :: IO (Either ErrorCall St.StateT)
      case result of
        Left err -> assertBool "ket diagnostic" ("QPP.Exact.ket" `isInfixOf` show err)
        Right _ -> assertFailure "expected an exception from the conventional interface"
  ]

-- Retain concise assertions for the unitary-only regression cases.
unitaryProgram :: (Limits -> Program -> St.StateT -> RNG -> Result (St.StateT, Outcomes, RNG))
               -> Limits -> Program -> St.StateT -> Result St.StateT
unitaryProgram run limits program state = do
  (output, _, _) <- run limits program state []
  pure output

orderingTests :: TestTree
orderingTests = testGroup "Certified real ordering"
  [ testCase "Real order, equality, and nonreal rejection" $ do
      let sqrtTwo = times (rat 2) h
      S.compareReal l sqrtTwo (rat (7%5)) @?= Right GT
      S.compareReal l sqrtTwo (rat (3%2)) @?= Right LT
      promoted <- right (S.embed l 24 sqrtTwo)
      S.compareReal l sqrtTwo promoted @?= Right EQ
      assertError (== InvalidInput "compareReal expects real scalars")
        (S.compareReal l imaginary imaginary)
  , testCase "Conjugate root pairs have the signs of their quadrants" $
      forM_ [3,5,7,8,9,12,15,16,24,30] $ \n ->
        forM_ [0..n-1] $ \j -> do
          let z = power (root n) j
              zc = fixture (S.conjugate l z)
              expected | 4*j == n || 4*j == 3*n = EQ
                       | 4*j < n || 4*j > 3*n = GT
                       | otherwise = LT
          S.compareReal l (plus z zc) (rat 0) @?= Right expected
  , testCase "Pell approximants require refinement far beyond Double precision" $ do
      let sqrtTwo = times (rat 2) h
          pell = iterate (\(p,q) -> (p+2*q,p+q)) (1,0)
      forM_ [30,100] $ \j -> do
        let (p,q) = pell !! j
            expected = if p*p > 2*q*q then LT else GT
        S.compareReal l sqrtTwo (rat (p%q)) @?= Right expected
      let (p,q) = pell !! 30
      assertError isResource (S.compareReal (l {maxIntegerBits=96}) sqrtTwo (rat (p%q)))
  , testProperty "Nonzero squared magnitudes are positive" $ \(SmallCyclo a) ->
      let result = do
            ac <- S.conjugate l a
            square <- S.mul l ac a
            S.compareReal l square (rat 0)
      in result === Right (if S.isZero a then EQ else GT)
  ]

measurementTests :: TestTree
measurementTests = testGroup "Exact measurement and initialization"
  [ testCase "Bell collapse retains its amplitude and repeated measurements consume RNG" $
      forM_ [M.evalProgWithLimits,V.evalProgWithLimits] $ \run -> do
        let bell = Compose (C X) (Tensor H (Id 1))
        (state,outs,remaining) <- right $
          run l [Unitary bell,Measure [0],Measure [1]] (St.ket [0,0]) [0.8,0,0.3]
        assertValues [rat 0,rat 0,rat 0,h] state
        outs @?= [True,True]
        remaining @?= [0.3]
        right (St.normSquared l state) >>= assertScalar (rat (1%2))
        right (St.probability l state 0 True) >>= assertScalar (rat 1)
  , testCase "Exact rational threshold does not round to the RNG Double" $ do
      let sqrtTwo = times (rat 2) h
      state <- right (St.fromAmplitudes l [rat 1,sqrtTwo,rat 1,rat 1])
      right (St.probability l state 0 False) >>= assertScalar (rat (3%5))
      -- The Double literal 0.6 is strictly below the exact rational 3/5.
      assertBool "dyadic threshold is below 3/5" (toRational (0.6 :: Double) < 3%5)
      (projected,outs,_) <- right (St.measure1WithLimits l (state,[],[0.6]) 0)
      outs @?= [False]
      assertValues [rat 1,sqrtTwo,rat 0,rat 0] projected
  , testCase "Threshold equality selects one; zero-probability branches are never chosen" $ do
      let plusState = V.apply (V.evalOp H) (St.ket [0])
          (_,outs,_) = St.measure1 (plusState,[],[0.5]) 0
      outs @?= [True]
      forM_ [(0,0),(1,0),(0,0.999),(1,0.999)] $ \(bit,r) -> do
        (state,result,remaining) <- right (St.measure1WithLimits l (St.ket [bit],[],[r]) 0)
        assertValues (map rat (if bit == 0 then [1,0] else [0,1])) state
        result @?= [bit == 1]
        remaining @?= []
  , testCase "Irrational probabilities and rescaling preserve outcomes" $ do
      let circuit = Compose H (Compose (C (Phase (1%4))) H)
          state = V.apply (V.evalOp circuit) (St.ket [0])
          p0 = times (rat (1%4)) (plus (rat 2) (times (rat 2) h))
      right (St.probability l state 0 False) >>= assertScalar p0
      scaled <- right (St.scale l (plus (rat 2) imaginary) state)
      forM_ [(0.85,False),(0.86,True)] $ \(r,expected) ->
        forM_ [state,scaled] $ \input -> do
          (_,outs,_) <- right (St.measure1WithLimits l (input,[],[r]) 0)
          outs @?= [expected]
          right (St.probability l input 0 False) >>= assertScalar p0
  , testCase "A branch whose normalization leaves cyclotomics remains usable" $ do
      let sqrtTwo = times (rat 2) h
          b = plus (rat 2) sqrtTwo
          weight = plus (rat 7) (times (rat 4) sqrtTwo)
      state <- right (St.fromAmplitudes l [rat 1,b,rat 1,rat 0])
      (projected,_,_) <- right (St.measure1WithLimits l (state,[],[0]) 0)
      assertValues [rat 1,b,rat 0,rat 0] projected
      right (St.normSquared l projected) >>= assertScalar weight
      expected <- right (S.inverse l weight)
      right (St.probability l projected 1 False) >>= assertScalar expected
      (_,outs,_) <- right (St.measure1WithLimits l (projected,[],[0]) 1)
      outs @?= [False]
  , testCase "Programs match numerical outcome order, RNG consumption and subset reset" $ do
      let program = [Unitary (Tensor H (Tensor H H)),Measure [1],
                     Initialize [2,0] [True,False],Measure [0,1,2]]
          rng = [0.2,0.8,0.2,0.8,0.2,0.8,0.3]
          (_,matrixOuts,matrixRng) = NumericalMatrix.evalProg program (NumericalMatrix.ket [0,0,0]) rng
          (_,vectorOuts,vectorRng) = NumericalVector.evalProg program (NumericalVector.ket [0,0,0]) rng
      matrixOuts @?= [False,False,True,False]
      vectorOuts @?= matrixOuts
      vectorRng @?= matrixRng
      forM_ [M.evalProg,V.evalProg] $ \run -> do
        let (state,outs,remaining) = run program (St.ket [0,0,0]) rng
        outs @?= matrixOuts
        remaining @?= matrixRng
        assertValues [rat 0,times (rat (1%2)) h,rat 0,rat 0,rat 0,rat 0,rat 0,rat 0] state
  , testCase "Standalone step preserves previous outcomes and reset hides its samples" $
      forM_ [M.evalStepWithLimits,V.evalStepWithLimits] $ \run -> do
        let state = V.apply (V.evalOp H) (St.ket [0])
        (reset,outs,rng) <- right (run l (state,[True],[0.8,0.2]) (Initialize [0] [False]))
        outs @?= [True]
        rng @?= [0.2]
        assertValues [h,rat 0] reset
        (_,outs',rest) <- right (run l (reset,outs,[0.2,0.9]) (Measure [0,0]))
        outs' @?= [False,False,True]
        rest @?= []
  , testCase "Empty steps preserve state, outcomes, field and RNG" $
      forM_ [M.evalProgWithLimits,V.evalProgWithLimits] $ \run -> do
        (state,outs,rng) <- right (run l [Measure [],Initialize [] []] (St.ket []) [0.4])
        assertValues [rat 1] state
        St.fieldIndex state @?= 1
        outs @?= []
        rng @?= [0.4]
  , testCase "Malformed measurement inputs report explicit failures" $ do
      zero <- right (St.fromAmplitudes l [rat 0,rat 0])
      assertError (== InvalidInput "cannot measure the zero vector")
        (St.measure1WithLimits l (zero,[],[0]) 0)
      assertError (== InvalidInput "cannot measure the zero vector") (St.probability l zero 0 False)
      assertError (== InvalidInput "measurement RNG is empty")
        (St.measure1WithLimits l (St.ket [0],[],[]) 0)
      forM_ [0/0,1/0,-1/0,-0.1,1] $ \r ->
        assertError (== InvalidInput "measurement RNG value must be finite and in [0,1)")
          (St.measure1WithLimits l (St.ket [0],[],[r]) 0)
      forM_ [-1,1] $ \k ->
        assertError (== InvalidInput "measurement wire out of bounds")
          (St.measure1WithLimits l (St.ket [0],[],[0]) k)
      forM_ [M.evalProgWithLimits,V.evalProgWithLimits] $ \run -> do
        assertError (== AtStep 0 (InvalidInput "Initialize wires must be distinct"))
          (run l [Initialize [0,0] [False,True]] (St.ket [0]) [])
        assertError (== AtStep 1 (InvalidInput "measurement RNG is empty"))
          (run l [Measure [0],Measure [0]] (St.ket [0]) [0.2])
        assertError (\err -> case err of AtStep 0 e -> isResource e; _ -> False)
          (run (l {maxCoefficientSlots=20}) [Measure [0]] (St.ket [0]) [0])
  ]

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
      [ testCase (show op) $ right (M.evalOpWithLimits l op) >>= assertMatrix expected
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
      right (M.evalOpWithLimits l (R Z 1)) >>= assertMatrix
        [[negative imaginary,rat 0],[rat 0,imaginary]]
      right (M.evalOpWithLimits l (C (R (Id 0) 2))) >>= assertMatrix
        [[rat 1,rat 0],[rat 0,rat (-1)]]
      right (runDirect (C (R (Id 0) 2)) [rat 2,rat 3]) >>= assertValues [rat 2,rat (-3)]
  , testCase "S and T preserve their specified global phases" $ do
      right (M.evalOpWithLimits l (C (Phase (1%2)))) >>= assertMatrix
        [[rat 1,rat 0],[rat 0,imaginary]]
      right (M.evalOpWithLimits l (Tensor (Phase (1%8)) (R Z (1%4)))) >>= assertMatrix
        [[rat 1,rat 0],[rat 0,root 8]]
  , testCase "Zero-qubit scalar evolution on an unnormalized value" $ do
      let input = plus (rat 2) imaginary
          expected = plus (rat (-1)) (times (rat 2) imaginary)
      v <- right (runDirect (Phase (1%2)) [input])
      assertValues [expected] v
      right (St.normSquared l v) >>= assertScalar (rat 5)
      p <- right (M.evalOpWithLimits l (Phase (1%2)))
      s <- right (St.fromAmplitudes l [input])
      right (M.applyWithLimits l p s) >>= assertValues [expected]
      right (M.evalOpWithLimits l (Permute [])) >>= assertMatrix [[rat 1]]
  , testCase "Zero-angle non-Hermitian axis follows the committed checker" $ do
      right (M.evalOpWithLimits l (R SX 0)) >>= assertMatrix [[rat 1,rat 0],[rat 0,rat 1]]
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
      m <- right (M.evalOpWithLimits l (R H (1%3)))
      right (M.entry m 0 0) >>= assertScalar a
      right (M.entry m 1 0) >>= assertScalar b
  , testCase "Phased Bell state and interference probabilities" $ do
      let bell = Compose (C X) (Tensor H (Id 1))
          op = Compose (Tensor (Id 1) (C (Phase (1%3)))) bell
          z = root 6
      v <- right (runDirect op [rat 1,rat 0,rat 0,rat 0])
      assertValues [h,rat 0,rat 0,times h z] v
      p <- right (V.evalOpWithLimits l (Adjoint bell))
      w <- right (V.applyWithLimits l p v)
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
      p <- right (M.evalOpWithLimits l (Permute ks))
      right (M.entry p 1 4) >>= assertScalar (rat 1)
      v <- right (St.ketWithLimits l [1,0,0])
      expected <- right (St.ketWithLimits l [0,0,1])
      op <- right (V.evalOpWithLimits l (Permute ks))
      actual <- right (V.applyWithLimits l op v)
      St.equal l expected actual @?= Right True
      inv <- right (V.evalOpWithLimits l (Adjoint (Permute ks)))
      restored <- right (V.applyWithLimits l inv actual)
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
      a <- right (M.evalOpWithLimits l H)
      right (M.normSquaredOp l a) >>= assertScalar (rat 2)
      b <- right (M.scaleOp l (rat 2) a)
      c <- right (M.addOp l a a)
      M.equalOp l b c @?= Right True
      d <- right (M.composeOp l a a)
      right (M.evalOpWithLimits l (Id 1)) >>= \ident -> M.equalOp l d ident @?= Right True
      right (M.fromRows l (M.entries b)) >>= \copy -> M.equalOp l b copy @?= Right True
  , testCase "Zero-qubit and ordinary state tensor products" $ do
      s <- right (St.fromAmplitudes l [imaginary])
      v <- right (St.fromAmplitudes l [rat 2,rat 3])
      right (St.tensor l s v) >>= assertValues [times imaginary (rat 2),times imaginary (rat 3)]
      right (St.tensor l v s) >>= assertValues [times imaginary (rat 2),times imaginary (rat 3)]
      right (St.tensor l s s) >>= assertValues [rat (-1)]
  , testCase "Program order and common-field preflight" $ do
      initial <- right (St.ketWithLimits l [0])
      let steps = [Unitary H, Unitary Z, Unitary (C (Phase (1%3)))]
          expected = [h,negative (times h (root 6))]
      forM_ [(unitaryProgram M.evalProgWithLimits),(unitaryProgram V.evalProgWithLimits)] $ \run -> right (run l steps initial) >>= assertValues expected
      v <- right ((unitaryProgram V.evalProgWithLimits) l steps initial)
      St.fieldIndex v @?= 24
  , testCase "An empty program preserves the input context" $ do
      v <- right (St.ketWithLimits l [])
      let tight = l {maxRootIndex=1}
      forM_ [(unitaryProgram M.evalProgWithLimits),(unitaryProgram V.evalProgWithLimits)] $ \run -> do
        x <- right (run tight [] v)
        St.fieldIndex x @?= 1
        assertValues [rat 1] x
  , testCase "External amplitude field promotes the operator once" $ do
      v <- right (St.fromAmplitudes l [root 5,rat 0])
      p <- right (V.evalOpWithLimits l H)
      m <- right (M.evalOpWithLimits l H)
      x <- right (V.applyWithLimits l p v)
      y <- right (M.applyWithLimits l m v)
      St.fieldIndex x @?= 40
      St.equal l x y @?= Right True
      assertValues [times h (root 5),times h (root 5)] x
  , testCase "Direct evaluation fits when a dense matrix cannot" $ do
      let limited = l {maxCoefficientSlots=256000}
          op = Tensor H (Id 9)
      p <- right (V.evalOpWithLimits limited op)
      s <- right (St.ketWithLimits limited (replicate 10 0))
      v <- right (V.applyWithLimits limited p s)
      right (St.amplitude v 0) >>= assertScalar h
      right (St.amplitude v 512) >>= assertScalar h
      right (St.normSquared limited v) >>= assertScalar (rat 1)
      assertError isResource (M.evalOpWithLimits limited op)
  ]

errorTests :: TestTree
errorTests = testGroup "Explicit failures"
  [ testCase "Invalid syntax and non-Hermitian axes" $ forM_
      [Id (-1),Permute [0,0],Compose X (Id 2),R SX (1%3),R (Id (-1)) 0,
       R (Compose X (Id 2)) 0,Tensor (Id maxBound) X] $ \op -> do
      assertBool "dense rejected" (isLeft (M.evalOpWithLimits l op))
      assertBool "direct rejected" (isLeft (V.evalOpWithLimits l op))
  , testCase "Huge widths are rejected before shifting or allocating" $ do
      assertError isResource (M.evalOpWithLimits l (Id maxBound))
      assertError isResource (V.evalOpWithLimits l (Id maxBound))
  , testCase "Field-degree explosion is rejected before coordinate allocation" $ do
      let limits = l {maxRootIndex=100000}
          op = Compose (R Z (1%101)) (R Z (1%103))
          degreeError (ResourceLimit "field degree" 40800 256) = True
          degreeError _ = False
      assertError degreeError (M.evalOpWithLimits limits op)
      assertError degreeError (V.evalOpWithLimits limits op)
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
      p <- right (V.evalOpWithLimits l (Id 0))
      assertError isResource (V.applyWithLimits tight p v)
      assertError isResource (St.scale tight (rat 1) v)
      assertError isResource (S.equal (l {maxRootIndex=8}) (root 3) (root 8))
  , testCase "Program preflight reports malformed steps before execution" $ do
      v <- right (St.ketWithLimits l [0])
      forM_ [(unitaryProgram M.evalProgWithLimits),(unitaryProgram V.evalProgWithLimits)] $ \run -> do
        assertError (== AtStep 1 (InvalidInput "program wire out of bounds"))
          (run l [Measure [0],Measure [1]] v)
        assertError (== AtStep 0 (InvalidInput "Initialize needs one value per wire"))
          (run l [Initialize [0] []] v)
        assertError (== AtStep 1 (WidthMismatch 2 1)) (run l [Unitary H,Unitary (C X)] v)
  , testCase "Shape checks and observation indices" $ do
      v <- right (St.ketWithLimits l [0])
      w <- right (St.ketWithLimits l [])
      assertError (== WidthMismatch 1 0) (St.add l v w)
      p <- right (V.evalOpWithLimits l (Id 0))
      assertError (== WidthMismatch 0 1) (V.applyWithLimits l p v)
      assertBool "amplitude index" (isLeft (St.amplitude v 2))
      m <- right (M.evalOpWithLimits l H)
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
  , testProperty "subtraction undoes addition across fields" $ \(SmallCyclo a) (SmallCyclo b) ->
      let result = S.add l a b >>= \sum' -> S.sub l sum' b >>= S.equal l a
      in counterexample (show result) (result === Right True)
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
      a <- right (M.evalOpWithLimits l (Id 2))
      b <- right (M.evalOpWithLimits l (Id 3))
      assertError (== WidthMismatch 2 3) (M.innerOp l a b)
      assertError (== WidthMismatch 2 3) (M.addOp l a b)
      assertError (== WidthMismatch 2 3) (M.composeOp l a b)
  , testCase "An empty program needs only the input's stored coordinates" $
      forM_ [3,257] $ \slots -> do
        let tiny = l {maxCoefficientSlots=slots}
        v <- right (St.ketWithLimits tiny [])
        St.equal tiny v v @?= Right True
        a <- right (S.rational tiny 1)
        S.equal tiny a a @?= Right True
        forM_ [(unitaryProgram M.evalProgWithLimits),(unitaryProgram V.evalProgWithLimits)] $ \run -> do
          w <- right (run tiny [] v)
          assertValues [rat 1] w
        assertError isResource ((unitaryProgram V.evalProgWithLimits) (tiny {maxCoefficientSlots=2}) [] v)
  , testCase "Linear arithmetic does not reserve polynomial workspace" $ do
      let tiny = l {maxCoefficientSlots=12}
      a <- right (S.rational tiny 2)
      b <- right (S.rational tiny 3)
      right (S.add tiny a b) >>= assertScalar (rat 5)
      right (S.negate tiny a) >>= assertScalar (rat (-2))
      S.equal tiny a a @?= Right True
  , testCase "Subtraction fits linear workspace even in a degree-32 field" $ do
      let xs = [toInteger j % 3 | j <- [0..31 :: Int]]
          ys = [toInteger (31-j) % 5 | j <- [0..31 :: Int]]
          expected = zipWith (-) xs ys
      a <- right (S.fromCoordinates l 120 xs)
      b <- right (S.fromCoordinates l 120 ys)
      -- Two inputs, one output, the modulus, and four scalar temporaries.
      c <- right (S.sub (l {maxCoefficientSlots=133}) a b)
      S.coefficients c @?= expected
      v <- right (St.fromAmplitudes l [a,b])
      w <- right (St.fromAmplitudes l [b,a])
      result <- right (St.sub (l {maxCoefficientSlots=229}) v w)
      map S.coefficients (St.amplitudes result) @?= [expected,map negate expected]
  , testCase "Long compositions reuse vector workspace in either association" $ do
      let gates = map (\g -> Tensor g (Id 9)) (replicate 127 X ++ [Y])
          rightAssociated = foldr Compose (Id 10) gates
          leftAssociated = foldl Compose (Id 10) gates
      input <- right (St.ketWithLimits l (replicate 10 0))
      -- X^127 Y = i Z on the first wire, so this is observably nonidentity.
      expected <- right (St.scale l imaginary input)
      sequential <- right ((unitaryProgram V.evalProgWithLimits) l (map Unitary (reverse gates)) input)
      St.equal l sequential expected @?= Right True
      forM_ [rightAssociated,leftAssociated] $ \op -> do
        p <- right (V.evalOpWithLimits l op)
        output <- right (V.applyWithLimits l p input)
        St.equal l output expected @?= Right True
        p' <- right (V.evalOpWithLimits l (Adjoint op))
        restored <- right (V.applyWithLimits l p' output)
        St.equal l restored input @?= Right True
      -- Whole-program preflight must use the same structural budget.
      right ((unitaryProgram V.evalProgWithLimits) l [Unitary rightAssociated] input) >>= \output ->
        St.equal l output expected @?= Right True
  , testCase "Long compositions still charge all retained plan constants" $ do
      let op = foldr Compose (Id 0) (replicate 500 (Phase (1%2)))
      assertError isResource (V.evalOpWithLimits (l {maxCoefficientSlots=2000}) op)
  , testCase "Zero rotations discard unused fields after validating their axes" $ do
      forM_ [103,1009] $ \den -> do
        let op = R (Phase (1%den)) 0
        checked <- right (Checked.checkOp l op)
        Checked.fieldIndex checked @?= 8
        Checked.qubits checked @?= 0
        m <- right (M.evalOpWithLimits l op)
        assertMatrix [[rat 1]] m
        p <- right (V.evalOpWithLimits l op)
        v <- right (St.fromAmplitudes l [rat 2])
        right (V.applyWithLimits l p v) >>= assertValues [rat 2]
        forM_ [(unitaryProgram M.evalProgWithLimits),(unitaryProgram V.evalProgWithLimits)] $ \run ->
          right (run l [Unitary op] v) >>= assertValues [rat 2]
      forM_ [R (Compose (Phase (1%1009)) X) 0, R (R SX (1%3)) 0] $ \op -> do
        assertError isInvalidOperator (M.evalOpWithLimits l op)
        assertError isInvalidOperator (V.evalOpWithLimits l op)
      assertError isResource (V.evalOpWithLimits l (R (Id maxBound) 0))
  , testCase "Existing field contexts do not require construction workspace" $ do
      let tight = l {maxCoefficientSlots=2000}
      a <- right (S.rootOfUnity l 120)
      v <- right (St.fromAmplitudes l [a])
      assertError isResource (S.rootOfUnity tight 120)
      right (St.scale tight (rat 1) v) >>= assertValues [a]
      right (St.fromAmplitudes tight [a]) >>= assertValues [a]
      right (S.embed tight 120 a) >>= assertScalar a
  , testCase "Prepared direct plans recheck constants under tighter limits" $ do
      p <- right (V.evalOpWithLimits l H)
      v <- right (St.ketWithLimits l [0])
      assertError isResource (V.applyWithLimits (l {maxIntegerBits=1}) p v)
  , testCase "Program field-degree errors identify the introducing step" $ do
      let limits = l {maxRootIndex=100000,maxFieldDegree=512}
          program = [Unitary (R Z (1%101)),Unitary (R Z (1%103))]
      v <- right (St.ketWithLimits l [0])
      forM_ [(unitaryProgram M.evalProgWithLimits),(unitaryProgram V.evalProgWithLimits)] $ \run ->
        assertError (== AtStep 1 (ResourceLimit "field degree" 40800 512)) (run limits program v)
  ]
  where
    isInvalidOperator (InvalidOperator _) = True
    isInvalidOperator _ = False
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
            m <- M.evalOpWithLimits l op
            p <- V.evalOpWithLimits l op
            a <- M.applyWithLimits l m v
            b <- V.applyWithLimits l p v
            St.equal l a b
      in counterexample (show result) (result === Right True)
  , testProperty "direct inverse and squared norm" $ \(CircuitCase _ op cs) ->
      let result = do
            v <- St.fromAmplitudes l cs
            p <- V.evalOpWithLimits l op
            pi' <- V.evalOpWithLimits l (Adjoint op)
            w <- V.applyWithLimits l p v
            restored <- V.applyWithLimits l pi' w
            a <- St.equal l v restored
            nv <- St.normSquared l v
            nw <- St.normSquared l w
            b <- S.equal l nv nw
            pure (a && b)
      in counterexample (show result) (result === Right True)
  , testProperty "compiled plan promotion preserves scalar linearity" $ \(CircuitCase _ op cs) ->
      let result = do
            v <- St.fromAmplitudes l cs
            p <- V.evalOpWithLimits l op
            w <- V.applyWithLimits l p v
            -- The generated operators use indices dividing 48; prime 5
            -- forces every compiled plan into a larger common field.
            promoted <- St.scale l (root 5) v
            actual <- V.applyWithLimits l p promoted
            expected <- St.scale l (root 5) w
            St.equal l actual expected
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
      m <- right (M.evalOpWithLimits l (qft n))
      assertMatrix expected m
      p <- right (V.evalOpWithLimits l (qft n))
      forM_ [0..d-1] $ \j -> do
        v <- right (St.ketWithLimits l [fromEnum (testBit j (n-1-k)) | k <- [0..n-1]])
        actual <- right (V.applyWithLimits l p v)
        assertValues [row !! j | row <- expected] actual
  | n <- [0..3]
  ]

constructorTests :: TestTree
constructorTests = testGroup "Every column of composite operators"
  [ testCase (show op) $ do
      m <- right (M.evalOpWithLimits l op)
      p <- right (V.evalOpWithLimits l op)
      let n = M.opQubits m
          d = 2^n
      forM_ [0..d-1] $ \j -> do
        input <- right (St.ketWithLimits l [fromEnum (testBit j (n-1-k)) | k <- [0..n-1]])
        expected <- traverse (\i -> right (M.entry m i j)) [0..d-1]
        right (V.applyWithLimits l p input) >>= assertValues expected
      md <- right (M.adjointOp l m)
      fromSyntax <- right (M.evalOpWithLimits l (Adjoint op))
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
      before <- right (M.evalOpWithLimits l op)
      rewritten <- right (M.evalOpWithLimits l (cleanOnes op))
      M.equalOp l before rewritten @?= Right True
  | op <- [R (Id 0) (1%3), C (R (Id 0) 2), Tensor H (R (Id 0) 2),
           DirectSum (R (Id 1) (1%3)) H, C (R (Tensor (Id 0) (Id 1)) 2)]
  ]

main :: IO ()
main = defaultMain (testGroup "ExactTest"
  [interfaceTests, orderingTests, measurementTests, scalarTests, gateTests, linearTests, errorTests, scalarProperties, evaluatorProperties,
   constructorTests, rewriteTests, qftTests, reviewRegressions])
