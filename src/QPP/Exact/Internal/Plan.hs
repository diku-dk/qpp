{-# LANGUAGE DeriveTraversable #-}
-- | Prepared operators for direct statevector evaluation.
-- Preparation replaces angles by exact constants and pushes adjoints to the
-- leaves. Execution is defined in QPP.Semantics.Exact.Statevector.
module QPP.Exact.Internal.Plan
  ( Plan(..), Action(..), dimension, compile, promote, validate
  , Requirements, requirements, budget
  ) where

import qualified Data.Vector as V
import QPP.Exact.Types
import qualified QPP.Exact.Checked as C
import QPP.Exact.Internal.Scalar (Field, Coeff)
import qualified QPP.Exact.Internal.Scalar as S
import QPP.Exact.Internal.Gates (primitive, inversePermutation)

-- Each node records its vector dimension. Structural constructors retain their
-- QOp names. R axis c d stores the coefficients of cI + dA in place of an angle;
-- adjoints have already been pushed to the leaves.
-- Scalar stores a multiplier for a phase or a scalar rotation on any width.
-- Gate stores a prepared 2-by-2 table, including adjointed primitive gates.
-- Derived traversals visit exactly the constants, preserving the structure.
data Plan a = Plan !Int (Action a)
  deriving (Show, Functor, Foldable, Traversable)

data Action a
  = Id
  | Scalar a
  | Gate (V.Vector a)
  | R (Plan a) a a
  | C (Plan a)
  | DirectSum (Plan a) (Plan a)
  | Tensor (Plan a) (Plan a)
  | Compose (Plan a) (Plan a)
  | Permute [Int]
  deriving (Show, Functor, Foldable, Traversable)

dimension :: Plan a -> Int
dimension (Plan dim _) = dim

-- Backward prepares the adjoint. Carrying this direction through the tree
-- avoids constructing and then conjugating an intermediate forward plan.
data Direction = Forward | Backward

opposite :: Direction -> Direction
opposite Forward = Backward
opposite Backward = Forward

directedAngle :: Direction -> Rational -> Rational
directedAngle Forward q = q
directedAngle Backward q = -q

compile :: Limits -> Field -> C.CheckedOp -> Result (Plan Coeff)
compile limits field = eval Forward
  where
    eval direction op = Plan (C.dimension op) <$> case C.view op of
      C.Id -> pure Id
      C.Phase q -> Scalar <$> S.phaseC limits field (directedAngle direction q)
      C.Primitive g -> Gate <$> gate direction g

      C.R axis q -> do
        (c,d) <- S.rotationC limits field (directedAngle direction q)
        if S.isZeroC d
          then pure (Scalar c)
          -- The checked axis is Hermitian: (cI + dA)^dagger = cI - dA.
          else do
            a <- eval Forward axis
            pure (R a c d)

      C.C a -> C <$> eval direction a
      C.DirectSum a b -> DirectSum <$> eval direction a <*> eval direction b
      C.Tensor a b -> Tensor <$> eval direction a <*> eval direction b

      -- (AB)^dagger = B^dagger A^dagger.
      C.Compose a b -> case direction of
        Forward -> Compose <$> eval Forward a <*> eval Forward b
        Backward -> Compose <$> eval Backward b <*> eval Backward a

      C.Adjoint a -> do
        Plan _ action <- eval (opposite direction) a
        pure action

      -- QOp specifies a scatter map; execution needs its inverse gather map.
      -- Taking the adjoint inverts the scatter map, so its gather map is ks.
      C.Permute ks -> pure $ Permute $ case direction of
        Forward -> inversePermutation ks
        Backward -> ks

    gate Forward g = primitive limits field g
    gate Backward g = do
      entries <- primitive limits field g
      V.generateM 4 $ \index ->
        let (row,column) = index `quotRem` 2
        in S.conjugateC limits field (entries V.! (2*column + row))

-- A prepared plan can be reused in a larger field without recompiling gates.
promote :: Limits -> Field -> Field -> Plan Coeff -> Result (Plan Coeff)
promote limits target source plan
  | S.index target == S.index source = pure plan
  | otherwise = traverse (S.embedC limits target source) plan

-- Recheck stored constants because apply may receive tighter limits than evalOp.
validate :: Limits -> Field -> Plan Coeff -> Result ()
validate limits field plan = do
  S.checkField limits field
  mapM_ (S.checkCoeff limits field) plan

-- Counts are field elements, converted to rational-coordinate slots once the
-- common field is known. Live values include active input and output, but not
-- constants or caller-retained input. These are conservative structural bounds:
-- slices, sharing and eliminated scalar-rotation axes can use less storage.
data Requirements = Requirements
  { constantCount :: !Integer
  , liveValues :: !Integer
  } deriving (Show)

requirements :: C.CheckedOp -> Requirements
requirements op = case C.view op of
  C.Id -> Requirements 0 dim
  C.Phase _ -> Requirements 1 (2*dim)
  C.Primitive _ -> Requirements 4 (2*dim)
  C.Permute _ -> Requirements 0 (2*dim)

  C.R axis _ ->
    let Requirements constants live = requirements axis
    -- Retain v while evaluating A*v, then combine v, A*v and the output.
    in Requirements (2+constants) (max (dim+live) (3*dim))

  C.C a ->
    let Requirements constants live = requirements a
    in Requirements constants (max (dim+live) (3*dim))

  C.DirectSum a b ->
    -- Retain the input and one completed half while evaluating the other.
    binary (\x y -> max (3*dim) (dim + dim `quot` 2 + max x y)) a b

  C.Tensor a b ->
    -- Input, rows, columns and final assembly, plus one active child.
    binary (\x y -> 4*dim + max x y) a b

  C.Compose a b -> binary max a b -- b's output becomes a's input
  C.Adjoint a -> requirements a
  where
    dim = toInteger (C.dimension op)
    binary peak a b =
      let Requirements ca wa = requirements a
          Requirements cb wb = requirements b
      in Requirements (ca+cb) (peak wa wb)

-- Two further vectors cover caller-retained input and promotion. Two copies
-- of constants cover old/new plans during embedding. Polynomial scratch also
-- covers the scalar temporaries used while compiling gates.
budget :: Limits -> Field -> Int -> Requirements -> Result ()
budget limits field dim required = S.workspace limits field $
  (2 * toInteger dim + 2 * constantCount required + liveValues required)
    * toInteger (S.degree field)
