-- | Shared infrastructure for the semantics backends.
--
--   == The backend module convention
--
--   A semantics backend is a plain module — not a class instance. Every
--   complete backend module (@QPP.Semantics.Matrix@, @.Statevector@, @.MPS@;
--   @.Stabilizer@ is an exercise hole still to be filled) defines the same set
--   of names as free functions and types:
--
--   [@StateT@]   the backend's state representation.
--   [@OpT@]      the backend's operator representation (which may be a
--                function/closure rather than an evaluated object, as long as
--                @apply@ can consume it).
--   [@ket@]      @:: [Int] -> StateT@ — build a computational-basis state.
--   [@evalOp@]   @:: QOp -> OpT@ — evaluate a symbolic operator.
--   [@apply@]    @:: OpT -> StateT -> StateT@.
--   [@measure1@] @:: (StateT, Outcomes, RNG) -> Nat -> (StateT, Outcomes, RNG)@.
--   [@evalStep@] @:: (StateT, Outcomes, RNG) -> Step -> (StateT, Outcomes, RNG)@.
--   [@evalProg@] @:: Program -> StateT -> RNG -> (StateT, Outcomes, RNG)@.
--
--   A user picks a backend by importing exactly one @QPP.Semantics.*@ module;
--   swapping that single import line switches the whole program to another
--   backend. This is why the "QPP" umbrella re-exports syntax, rewriting and
--   pretty-printing but never a backend — the choice has to be explicit.
--
--   `evalStep` and `evalProg` are the same fold over `Step`s for every
--   backend, so they are defined once here as `evalStepWith`/`evalProgWith`
--   and instantiated by each backend with its own unitary application and
--   single-qubit measurement.
--
--   Backends also give their @OpT@ a @Convertible OpT CMat@ instance (see
--   `Convertible` below), which is how the test suites compare them against
--   the dense reference "QPP.Semantics.Matrix". This module additionally holds
--   `HasWork` (for backends whose working representation differs from the
--   stored one) and `SparseMat`.
module QPP.Semantics where

import QPP.Syntax
import Data.Bits(xor)
import Data.Array(accumArray,elems)
import Data.List (foldl')


--------------------------------------------------------------------------------
-- Numerical tolerance
--------------------------------------------------------------------------------

-- | The tolerance every backend uses for its floating-point sanity checks:
--   measurement probabilities must sum to 1 within `tol`, a state whose norm is
--   below `tol` is not renormalised, and sparse conversions drop entries whose
--   magnitude is at most `tol`. Defined once here so that swapping backends
--   never changes a threshold. (The MPS backend's per-state `EvalCfg` carries
--   its own `tol` field, defaulting to the same value.)
tol :: RealT
tol = 1e-12


--------------------------------------------------------------------------------
-- Conversions between representations
--------------------------------------------------------------------------------

-- | Conversion between two representations of the same object, e.g. a
--   backend's `StateT` and a dense matrix or a sparse amplitude list. Used to
--   compare backends against each other and to print states uniformly.
class Convertible a b where
  to   :: a -> b
  from :: b -> a

-- | TODO: 1) SparseMat t, 2) Useful functions for SparseMat
data SparseMat = SparseMat ((Integer,Integer), [((Integer,Integer), ComplexT)])
   deriving (Show,Eq)

-- | Backends that compute in a different representation from the one they
--   store states in (e.g. a delayed array while working, a manifest one at
--   rest) provide the conversion pair here. `Work t` is the backend's working
--   representation for the stored type `t`.
class HasWork t where
  type Work t
  toWork   :: t -> Work t
  fromWork :: Work t -> t


--------------------------------------------------------------------------------
-- The shared Step/Program interpreter
--------------------------------------------------------------------------------

-- | Shared `Step` interpreter, generalized over a monoidal annotation `w`
--   (e.g. the MPS backend's truncation `Profile`; `w ~ ()` for backends with
--   nothing to report). Backends supply unitary application (including any
--   arity checking they want) and single-qubit measurement; `Measure` and
--   `Initialize` handling is identical across backends and lives here.
--   `Initialize` measures the qubits and applies X-corrections where the
--   outcome differs from the requested classical value.
evalStepWithW :: (Monoid w, HasQubits st)
              => (QOp -> st -> (st, w))                                    -- ^ apply a unitary
              -> ((st, Outcomes, RNG) -> Nat -> ((st, Outcomes, RNG), w))  -- ^ measure one qubit
              -> (st, Outcomes, RNG) -> Step -> ((st, Outcomes, RNG), w)
evalStepWithW applyU meas (st, outs, rng) step = case step of
  Unitary op -> let (st', w) = applyU op st in ((st', outs, rng), w)

  -- outcomes are latest-first, so ks is reversed on input
  Measure ks ->
    foldl' (\(acc, w) k -> let (acc', w') = meas acc k in (acc', w <> w'))
           ((st, outs, rng), mempty) (reverse ks)

  Initialize ks vs ->
    let n = n_qubits st
        ((st', os, rng'), w1) = evalStepWithW applyU meas (st, [], rng) (Measure ks)
        -- List of outcomes xor values for each initialized qubit
        corrections     = zipWith xor os vs
        -- Now we build the full list, including unaffected qubits
        corrFull        = accumArray xor False (0,n-1) (zip ks corrections)
        corrOp          = foldl (⊗) One [ if c then X else Id 1 | c <- elems corrFull ]
        (res, w2)       = evalStepWithW applyU meas (st', outs, rng') (Unitary corrOp)
    in (res, w1 <> w2)

-- | Shared program evaluator: fold `evalStepWithW` over the steps.
evalProgWithW :: (Monoid w, HasQubits st)
              => (QOp -> st -> (st, w))
              -> ((st, Outcomes, RNG) -> Nat -> ((st, Outcomes, RNG), w))
              -> Program -> st -> RNG -> ((st, Outcomes, RNG), w)
evalProgWithW applyU meas prog st rng =
  foldl' (\(acc, w) s -> let (acc', w') = evalStepWithW applyU meas acc s in (acc', w <> w'))
         ((st, [], rng), mempty) prog

-- | Annotation-free instantiations (w ~ ()).
evalStepWith :: HasQubits st
             => (QOp -> st -> st)                                    -- ^ apply a unitary
             -> ((st, Outcomes, RNG) -> Nat -> (st, Outcomes, RNG))  -- ^ measure one qubit
             -> (st, Outcomes, RNG) -> Step -> (st, Outcomes, RNG)
evalStepWith applyU meas acc step =
  fst $ evalStepWithW (\op st -> (applyU op st, ()))
                      (\a k -> (meas a k, ()))
                      acc step

evalProgWith :: HasQubits st
             => (QOp -> st -> st)
             -> ((st, Outcomes, RNG) -> Nat -> (st, Outcomes, RNG))
             -> Program -> st -> RNG -> (st, Outcomes, RNG)
evalProgWith applyU meas prog st rng =
  fst $ evalProgWithW (\op s -> (applyU op s, ()))
                      (\a k -> (meas a k, ()))
                      prog st rng
