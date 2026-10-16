-- | Exact dense QOp semantics. Entries include global phase. OpT can also hold
-- general square matrices via fromRows/linear combinations; only evalOp's
-- checked QOp construction guarantees unitarity.
module QPP.Semantics.Exact.Matrix
  ( OpT, evalOp, apply, evalProg
  , opQubits, opFieldIndex, entry, entries, fromRows
  , addOp, scaleOp, innerOp, normSquaredOp, equalOp, adjointOp, composeOp
  , module QPP.Exact.State
  ) where

import Control.Monad (foldM, unless)
import Data.Bifunctor (first)
import qualified Data.Vector as V
import QPP.Syntax (QOp, Program)
import QPP.Exact.State
import qualified QPP.Exact.Checked as C
import QPP.Exact.Internal.Scalar (Field, Coeff, Cyclo(..))
import qualified QPP.Exact.Internal.Scalar as S
import qualified QPP.Exact.Internal.State as St
import QPP.Exact.Internal.Buffer (Buffer(..))
import qualified QPP.Exact.Internal.Buffer as B
import QPP.Exact.Internal.Evaluation (prepareProgram)
import QPP.Exact.Internal.Gates (primitive, permuteIndex)

data OpT = OpT !Int !Buffer deriving (Show)

opQubits :: OpT -> Int
opQubits (OpT n _) = n

opFieldIndex :: OpT -> Integer
opFieldIndex (OpT _ b) = S.index (B.field b)

-- | Complete exact rows, with field context attached to each observed scalar.
entries :: OpT -> [[Cyclo]]
entries (OpT n (Buffer f xs)) =
  [[Cyclo f (xs V.! (i*d+j)) | j <- [0..d-1]] | i <- [0..d-1]]
  where d = 2^n

entry :: OpT -> Int -> Int -> Result Cyclo
entry (OpT n (Buffer f xs)) i j
  | i >= 0 && i < d && j >= 0 && j < d = pure (Cyclo f (xs V.! (i*d+j)))
  | otherwise = Left (InvalidInput "matrix index out of bounds")
  where d = 2^n

fromRows :: Limits -> [[Cyclo]] -> Result OpT
fromRows l rows = do
  n <- S.widthOfLength (length rows)
  d <- S.dimension l n
  S.slots l (toInteger d * toInteger d)
  unless (all ((== length rows) . length) rows)
    (Left (InvalidInput "matrix must be square with power-of-two dimensions"))
  OpT n <$> B.fromScalars l (concat rows)

matrixBudget :: Limits -> Field -> C.CheckedOp -> Result ()
-- At each recursion level retain up to two child matrices and one result;
-- eight further buffers cover input promotion, application and block assembly.
matrixBudget l f op = S.workspace l f $
  (3 * toInteger (C.depth op) + 8) * d * d * toInteger (S.degree f)
  where d = toInteger (C.dimension op)

evalOp :: Limits -> QOp -> Result OpT
evalOp l op = do
  c <- C.checkOp l op
  f <- S.makeField l (C.fieldIndex c)
  matrixBudget l f c
  evaluate l f c

-- Inputs and workspace have been preflighted. All recursive matrices share f.
evaluate :: Limits -> Field -> C.CheckedOp -> Result OpT
evaluate l f c = OpT (C.qubits c) . Buffer f <$> go c
  where
    ident d = V.generate (d*d) $ \k -> if k `quot` d == k `rem` d then S.one f else S.zero f
    go op = let d = C.dimension op in case C.view op of
      C.Identity -> pure (ident d)
      C.ScalarPhase q -> V.singleton <$> S.phaseC l f q
      C.Primitive g -> primitive l f g
      C.Rotation axis q -> do
        (a,b) <- S.rotationC l f q
        if S.isZeroC b
          then pure (V.generate (d*d) (\k -> if k `quot` d == k `rem` d then a else S.zero f))
          else do
            xs <- go axis
            V.imapM (\k x -> do
              y <- S.mulC l f b x
              if k `quot` d == k `rem` d then S.addC l a y else pure y) xs
      C.Control child -> do
        b <- go child
        pure (block (d `quot` 2) (ident (d `quot` 2)) b)
      C.BlockSum a b -> block (C.dimension a) <$> go a <*> go b
      C.Composition a b -> do
        xs <- go a
        ys <- go b
        multiply l f d xs ys
      C.TensorProduct a b -> do
        xs <- go a
        ys <- go b
        let da = C.dimension a
            db = C.dimension b
        V.generateM (d*d) $ \k -> do
          let (i,j) = k `quotRem` d
          S.mulC l f (xs V.! ((i `quot` db)*da+j `quot` db))
                     (ys V.! ((i `rem` db)*db+j `rem` db))
      C.Dagger a -> go a >>= transposeConjugate l f d
      C.Permutation ks -> pure $ V.generate (d*d) $ \k ->
        let (i,j) = k `quotRem` d
        in if i == permuteIndex ks j then S.one f else S.zero f
    block r a b = V.generate (4*r*r) $ \k ->
      let (i,j) = k `quotRem` (2*r)
      in if i < r && j < r then a V.! (i*r+j)
         else if i >= r && j >= r then b V.! ((i-r)*r+j-r)
         else S.zero f

multiply :: Limits -> Field -> Int -> V.Vector Coeff -> V.Vector Coeff -> Result (V.Vector Coeff)
multiply l f d a b = V.generateM (d*d) $ \k ->
  let (i,j) = k `quotRem` d
  in foldM (\acc t -> S.mulC l f (a V.! (i*d+t)) (b V.! (t*d+j)) >>= S.addC l acc)
           (S.zero f) [0..d-1]

transposeConjugate :: Limits -> Field -> Int -> V.Vector Coeff -> Result (V.Vector Coeff)
transposeConjugate l f d xs = V.generateM (d*d) $ \k ->
  let (i,j) = k `quotRem` d in S.conjugateC l f (xs V.! (j*d+i))

apply :: Limits -> OpT -> StateT -> Result StateT
apply l (OpT n matrix) state@(St.StateT _ input) = do
  unless (n == qubits state) (Left (WidthMismatch n (qubits state)))
  (mat,v) <- B.preparePair l S.PolynomialWork (B.count input) matrix input
  applyIn l (OpT n mat) (St.StateT n v)

-- Established common context; public apply and program preflight own checks.
applyIn :: Limits -> OpT -> StateT -> Result StateT
applyIn l (OpT n (Buffer f mat)) (St.StateT _ (Buffer _ v)) = do
  let d = 2^n
  out <- V.generateM d $ \i -> foldM
    (\acc j -> S.mulC l f (mat V.! (i*d+j)) (v V.! j) >>= S.addC l acc)
    (S.zero f) [0..d-1]
  pure (St.StateT n (Buffer f out))

evalProg :: Limits -> Program -> StateT -> Result StateT
evalProg l program state@(St.StateT n input) = do
  St.validate l state
  (f,ops) <- prepareProgram l n (B.field input) program
  mapM_ (\(j,op) -> first (AtStep j) (matrixBudget l f op)) (zip [0..] ops)
  start <- St.StateT n <$> B.promote l f input
  foldM (\v (j,op) -> first (AtStep j) (evaluate l f op >>= \m -> applyIn l m v)) start (zip [0..] ops)

sameWidth :: OpT -> OpT -> Result ()
sameWidth a b = unless (opQubits a == opQubits b)
  (Left (WidthMismatch (opQubits a) (opQubits b)))

addOp :: Limits -> OpT -> OpT -> Result OpT
addOp l a@(OpT n xs) b@(OpT _ ys) = sameWidth a b >> (OpT n <$> B.add l xs ys)

scaleOp :: Limits -> Cyclo -> OpT -> Result OpT
scaleOp l c (OpT n b) = OpT n <$> B.scale l c b

innerOp :: Limits -> OpT -> OpT -> Result Cyclo
innerOp l a@(OpT _ xs) b@(OpT _ ys) = sameWidth a b >> B.inner l xs ys

normSquaredOp :: Limits -> OpT -> Result Cyclo
normSquaredOp l a = innerOp l a a

equalOp :: Limits -> OpT -> OpT -> Result Bool
equalOp l (OpT _ a) (OpT _ b) = B.equal l a b

adjointOp :: Limits -> OpT -> Result OpT
adjointOp l (OpT n b@(Buffer f xs)) = do
  B.validate l b
  S.workspace l f (2 * B.count b * toInteger (S.degree f))
  OpT n . Buffer f <$> transposeConjugate l f (2^n) xs

composeOp :: Limits -> OpT -> OpT -> Result OpT
composeOp l a@(OpT n av) b@(OpT _ bv) = do
  sameWidth a b
  (Buffer f xs, Buffer _ ys) <- B.preparePair l S.PolynomialWork (B.count av) av bv
  OpT n . Buffer f <$> multiply l f (2^n) xs ys
