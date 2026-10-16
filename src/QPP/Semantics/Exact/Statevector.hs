-- | Direct exact application of QOp: no dense operator matrices, no numerical
-- trig/exponential, no branch normalization. The checked plan is compiled to
-- scalar constants once by evalOp, then reused or embedded when applied.
module QPP.Semantics.Exact.Statevector
  ( OpT, evalOp, apply, evalProg, opQubits, opFieldIndex
  , module QPP.Exact.State
  ) where

import Control.Monad (foldM, unless)
import Data.Bifunctor (first)
import qualified Data.Vector as V
import QPP.Syntax (QOp, Program)
import QPP.Exact.State
import qualified QPP.Exact.Checked as C
import QPP.Exact.Internal.Scalar (Field, Coeff)
import qualified QPP.Exact.Internal.Scalar as S
import qualified QPP.Exact.Internal.State as St
import QPP.Exact.Internal.Buffer (Buffer(..))
import qualified QPP.Exact.Internal.Buffer as B
import QPP.Exact.Internal.Evaluation (prepareProgram)
import QPP.Exact.Internal.Gates (primitive, permuteIndex, inversePermutation)

-- Keep only the compiled plan and its execution requirements, not a second
-- copy of the source tree. Selectors stay private to prevent record updates.
data OpT = OpT
  { operatorQubits :: !Int
  , operatorSlots :: !Integer
  , operatorField :: !Field
  , operatorPlan :: !Plan
  } deriving (Show)

opQubits :: OpT -> Int
opQubits = operatorQubits

opFieldIndex :: OpT -> Integer
opFieldIndex = S.index . operatorField

-- Four vectors per active level cover saved rotation inputs and tensor
-- rows/columns; eight cover input promotion and final assembly. Each plan
-- node has at most four constants, with old/new copies during embedding.
planSlots :: C.CheckedOp -> Integer
planSlots op = (4 * toInteger (C.depth op) + 8) * toInteger (C.dimension op) + 8*C.size op

planBudget :: Limits -> Field -> Integer -> Result ()
planBudget l f count = S.workspace l f (count * toInteger (S.degree f))

evalOp :: Limits -> QOp -> Result OpT
evalOp l op = do
  c <- C.checkOp l op
  f <- S.makeField l (C.fieldIndex c)
  let slots = planSlots c
  planBudget l f slots
  plan <- compile l f False c
  pure OpT {operatorQubits = C.qubits c, operatorSlots = slots,
            operatorField = f, operatorPlan = plan}

data Plan = Plan !Int Action deriving (Show)
data Action
  = Identity | Scalar Coeff | Gate (V.Vector Coeff)
  | Rotation Coeff Coeff Plan | Control Plan | BlockSum Plan Plan
  | Tensor Plan Plan | Compose Plan Plan | Permutation [Int]
  deriving (Show)

dimension :: Plan -> Int
dimension (Plan d _) = d

compile :: Limits -> Field -> Bool -> C.CheckedOp -> Result Plan
compile l f dagger op = Plan (C.dimension op) <$> case C.view op of
  C.Identity -> pure Identity
  C.ScalarPhase q -> Scalar <$> S.phaseC l f (if dagger then -q else q)
  C.Primitive g -> do
    xs <- primitive l f g
    ys <- if dagger then V.generateM 4 (\j -> S.conjugateC l f (xs V.! ((j `rem` 2)*2+j `quot` 2))) else pure xs
    pure (Gate ys)
  C.Rotation a q -> do
    (c,d) <- S.rotationC l f (if dagger then -q else q)
    if S.isZeroC d then pure (Scalar c)
      else Rotation c d <$> compile l f False a -- valid nonzero axis is Hermitian
  C.Control a -> Control <$> recur a
  C.BlockSum a b -> BlockSum <$> recur a <*> recur b
  C.TensorProduct a b -> Tensor <$> recur a <*> recur b
  C.Composition a b -> if dagger then Compose <$> recur b <*> recur a else Compose <$> recur a <*> recur b
  C.Dagger a -> do
    Plan _ action <- compile l f (not dagger) a
    pure action
  C.Permutation ks ->
    -- Store the gather map, the inverse of the denotation's scatter map.
    pure (Permutation (if dagger then ks else inversePermutation ks))
  where recur = compile l f dagger

-- Visit constants without rebuilding syntax or recomputing roots of unity.
planCoefficients :: Plan -> [Coeff]
planCoefficients p = collect p []
  where
    collect (Plan _ action) rest = case action of
      Identity -> rest
      Scalar c -> c:rest
      Gate xs -> V.foldr (:) rest xs
      Rotation c d a -> c:d:collect a rest
      Control a -> collect a rest
      BlockSum a b -> collect a (collect b rest)
      Tensor a b -> collect a (collect b rest)
      Compose a b -> collect a (collect b rest)
      Permutation _ -> rest

promotePlan :: Limits -> Field -> Field -> Plan -> Result Plan
promotePlan l target source p
  | S.index target == S.index source = pure p
  | otherwise = walk p
  where
    scalar = S.embedC l target source
    walk (Plan d action) = Plan d <$> case action of
      Identity -> pure Identity
      Scalar c -> Scalar <$> scalar c
      Gate xs -> Gate <$> V.mapM scalar xs
      Rotation c e a -> Rotation <$> scalar c <*> scalar e <*> walk a
      Control a -> Control <$> walk a
      BlockSum a b -> BlockSum <$> walk a <*> walk b
      Tensor a b -> Tensor <$> walk a <*> walk b
      Compose a b -> Compose <$> walk a <*> walk b
      Permutation ks -> pure (Permutation ks)

runPlan :: Limits -> Field -> Plan -> V.Vector Coeff -> Result (V.Vector Coeff)
runPlan l f (Plan dim action) v = case action of
  Identity -> pure v
  Scalar c -> V.mapM (S.mulC l f c) v
  Gate g -> V.generateM 2 $ \i -> do
    a <- S.mulC l f (g V.! (2*i)) (v V.! 0)
    b <- S.mulC l f (g V.! (2*i+1)) (v V.! 1)
    S.addC l a b
  Rotation c d axis -> do
    u <- recur axis v -- v must remain intact until both contributions exist
    V.zipWithM (\x y -> do
      a <- S.mulC l f c x
      b <- S.mulC l f d y
      S.addC l a b) v u
  Control p -> do
    let half = dim `quot` 2
    b <- recur p (V.slice half half v)
    pure (V.take half v V.++ b)
  BlockSum a b -> do
    let half = dim `quot` 2
    x <- recur a (V.take half v)
    y <- recur b (V.slice half half v)
    pure (x V.++ y)
  Compose a b -> recur b v >>= recur a
  Tensor a b -> do
    let da = dimension a
        db = dimension b
    -- Rows receive B as column vectors; columns then receive A. This is
    -- A V B^T, not A V B^dagger, including for entangled input vectors.
    rows <- V.generateM da (\i -> recur b (V.slice (i*db) db v))
    columns <- V.generateM db $ \j ->
      recur a (V.generate da (\i -> (rows V.! i) V.! j))
    pure (V.generate dim (\k -> (columns V.! (k `rem` db)) V.! (k `quot` db)))
  Permutation gather -> pure (V.generate dim (\i -> v V.! permuteIndex gather i))
  where recur = runPlan l f

-- The plan and input already share f; preparation and checks are external.
applyIn :: Limits -> Field -> Plan -> StateT -> Result StateT
applyIn l f p (St.StateT n (Buffer _ xs)) =
  St.StateT n . Buffer f <$> runPlan l f p xs

apply :: Limits -> OpT -> StateT -> Result StateT
apply l op state@(St.StateT n input) = do
  unless (opQubits op == n) (Left (WidthMismatch (opQubits op) n))
  let sourceField = operatorField op
      plan = operatorPlan op
  S.checkField l sourceField
  mapM_ (S.checkCoeff l sourceField) (planCoefficients plan)
  St.validate l state
  f <- S.commonField l sourceField (B.field input)
  planBudget l f (operatorSlots op)
  p <- promotePlan l f sourceField plan
  v <- St.StateT n <$> B.promote l f input
  applyIn l f p v

evalProg :: Limits -> Program -> StateT -> Result StateT
evalProg l program state@(St.StateT n input) = do
  St.validate l state
  (f,ops) <- prepareProgram l n (B.field input) program
  mapM_ (\(j,c) -> first (AtStep j) (planBudget l f (planSlots c))) (zip [0..] ops)
  start <- St.StateT n <$> B.promote l f input
  foldM (\v (j,c) -> first (AtStep j) (compile l f False c >>= \p -> applyIn l f p v))
    start (zip [0..] ops)
