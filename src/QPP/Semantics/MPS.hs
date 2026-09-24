-- | The MPS semantics backend: evaluates `QOp` syntax trees and `Program`
--   steps on the matrix-product-state representation.
--
--   The MPS data type and its canonical-form machinery live in
--   "QPP.MPS"; the diagonal-core specialization in "QPP.MPS.Diagonal".
--   Both are re-exported here so importing this module is enough to use the
--   backend.
--
--   Truncation statistics: every entry point exists in a plain and a
--   `Profiled` variant. The plain ones discard the run's `Profile`; the
--   profiled ones return it. Profiling is always on internally (the full
--   singular spectrum is computed by every SVD anyway), so the variants cost
--   the same.
module QPP.Semantics.MPS
  ( Site(..), MPS(..), StateT, OpT(..)
  , Interval(..)
  , Trunc(..), EvalCfg(..), defaultCfg
  , Profile(..)
  , ket, toSparseMat, mpsToDenseVec, opToDenseMat
  , measureProjection, measure1, measure1Profiled, sampleAll
    -- Diagonal-core MPS
  , DiagSite(..), DiagMPS(..)
  , mkDiagMPS, dNSites, siteToDiag
  , isDiagonalSite, isDiagonalMPS
  , toDiagMPS, fromDiagMPS
  , sampleDiagMPS, measureDiagMPS
  , sampleAllDiag, measureAllDiag
  , apply, applyProfiled, evalStep, evalProg, evalProgProfiled
  , evalOp, evalOpAt
  , dagger
  , approxCfg, compressRange, moveCenterToPhys, bondDimensions, maxBondDimension
  ) where
-- | TODO: 1. Implement "physical" permutation of MPS sites
---        2. Optimized multi-qubit rotation and multi-controlled gates
---        3. Visualization tools for MPS

import Data.Complex (Complex(..), conjugate)
import Data.List (foldl')
import Data.Vector ((!))
import qualified Data.Set as S
import Control.Monad.Writer.CPS (Writer, runWriter, writer)
import QPP.Syntax
import QPP.Util
import QPP.Semantics
import QPP.MPS
import QPP.MPS.Diagonal
import QPP.Semantics.Matrix (CMat)
import qualified Numeric.LinearAlgebra as H
import GHC.Stack (HasCallStack)

type StateT = MPS
data OpT    = OpT { opQubits :: !Int, runOp :: MPS -> (MPS, Profile) }

apply :: OpT -> StateT -> StateT
apply (OpT _ f) x = fst (f x)

-- | `apply` with the run's truncation statistics.
applyProfiled :: OpT -> StateT -> (StateT, Profile)
applyProfiled (OpT _ f) x = f x

instance HasQubits OpT where n_qubits = opQubits

pauliGate :: QOp -> Gate1
pauliGate = \case
  X -> (0,1,1,0)
  Y -> (0, 0:+(-1), 0:+1, 0)
  Z -> (1,0,0,-1)
  _ -> error "pauliGate"

-- | exp(−iπθ/2 · P) as a 2×2 gate, for a bare 1-qubit Pauli axis P.
--   Nothing for any other axis (handled by the general LCU path).
pauliRotGate :: QOp -> Rational -> Maybe Gate1
pauliRotGate axis θ =
  let t  = pi * fromRational θ / 2
      c  = cos t :+ 0
      s  = sin t
      is = 0 :+ s
  in case axis of
       X -> Just (c, -is, -is, c)
       Y -> Just (c, (-s):+0, s:+0, c)
       Z -> Just (c-is, 0, 0, c+is)
       _ -> Nothing

-- | Decompose a rotation axis (Phase · Pauli-string in Tensor form,
--   type-checked upstream) into (global phase φ, per-qubit operator list).
axisPaulis :: Int -> QOp -> (ComplexT, [QOp])
axisPaulis n op = case op of
  Id m      -> (1, replicate m (Id 1))
  Phase q   -> (cisPi q, replicate n (Id 1))
  X         -> (1, [X])
  Y         -> (1, [Y])
  Z         -> (1, [Z])

  Tensor a b ->
    let (ϕ1,p1) = axisPaulis (op_qubits a) a
        (ϕ2,p2) = axisPaulis (op_qubits b) b
    in (ϕ1*ϕ2, p1++p2)

  -- TODO: Allow Compose

  -- (ϕ1·P)† = ϕ1*·P for a Pauli string P (Hermitian).
  Adjoint a -> let (ϕ1,ps) = axisPaulis n a in (conjugate ϕ1, ps)

  _         -> error "axisPaulis: axis not Pauli string"

applyPauliString :: Int -> QOp -> MPS -> (ComplexT, MPS, Interval)
applyPauliString base axis m =
  let n = op_qubits axis
      (ϕ, ops) = axisPaulis n axis
      actIdx = op_support axis
      phys   = [ log2phys m ! (base+i) | i <- S.toList actIdx ]
      iSupp = case phys of
                [] -> singletonIval (log2phys m ! base)
                _  -> Ival (minimum phys) (maximum phys)
      m' = foldl' (\acc (i,g) -> case g of
                                   Id _ -> acc
                                   _    -> apply1Logical base i (pauliGate g) acc
                   ) m (zip [0..n-1] ops)
  in (ϕ, m', iSupp)

-- support interval (physical hull) from op_support
supportInterval :: MPS -> Int -> QOp -> Interval
supportInterval st base op =
  let supL = S.toList (op_support op)
      ps   = [ log2phys st ! (base+q) | q <- supL ]
  in case ps of
       [] -> singletonIval (log2phys st ! base)
       _  -> Ival (minimum ps) (maximum ps)

-- | Projection onto qubit k = out, as an operator. Unlike `measure1` the
--   result is *not* renormalized.
measureProjection :: HasCallStack => Int -> Int -> Int -> OpT
measureProjection arity k out = OpT arity $ \st ->
  projectCtrl (log2phys st ! k) (out==1) st

-- | Materialize an `OpT` as a dense matrix by applying it to each canonical basis ket
--   prepared with the given EvalCfg, then assembling the resulting column vectors. The cfg
--   propagates into the lambda's internal MPS operations, so the matrix reflects the
--   operator's actual truncation behavior under that cfg. Intended for small n.
opToDenseMat :: EvalCfg -> OpT -> CMat
opToDenseMat cfg' (OpT n f) =
  let mkKet bits = (ket bits) { cfg = cfg' }
      columns    = [ mpsToDenseVec (fst (f (mkKet (toBits' n j)))) | j <- [0 .. pow2 n - 1] ]
  in H.fromBlocks [columns]

instance Convertible OpT CMat where
  to     = opToDenseMat defaultCfg
  from _ = error "Convertible OpT CMat: from is not implemented; build OpT via evalOp"

-- evaluator
evalOp :: QOp -> OpT
evalOp op = OpT (op_qubits op) (runWriter . evalOpAt 0 op)

-- | Evaluate `op` with its own qubit window starting at logical qubit `base`.
--   Runs in @Writer Profile@: every SVD performed anywhere below is emitted
--   exactly once, so the accumulated profile covers *all* branches of
--   controlled/direct-sum evaluation.
evalOpAt :: HasCallStack => Int -> QOp -> MPS -> Writer Profile MPS
evalOpAt base op st = case op of
  Id _      -> pure st
  Phase q   -> pure ((cisPi q) .* st)
  Permute π -> pure (applyPermute base π st)
  X  -> pure (apply1Logical base 0 (pauliGate X) st)
  Y  -> pure (apply1Logical base 0 (pauliGate Y) st)
  Z  -> pure (apply1Logical base 0 (pauliGate Z) st)
  H  -> let s = (1/sqrt 2):+0 in pure (apply1Logical base 0 (s,s,s,-s) st)
  SX ->
    let p = 0.5:+0.5
        m = 0.5:+(-0.5)
    in pure (apply1Logical base 0 (p, m, m, p) st)

  Tensor a b -> evalOpAt base a st >>= evalOpAt (base + op_qubits a) b

  Compose a b -> evalOpAt base b st >>= evalOpAt base a

  Adjoint a -> evalOpAt base (dagger a) st

  R axis θ
    | θ == 0 -> pure st
    -- Bare 1-qubit Pauli axis: exp(−iπθ/2·P) is a single-site gate. Applying
    -- it directly (no branch-add, no SVD) keeps the state's canonical and
    -- diagonal-core structure intact.
    | Just g <- pauliRotGate axis θ -> pure (apply1Logical base 0 g st)
    | otherwise ->
        let t = pi * fromRational θ / 2
            (phi, pst, iSupp) = applyPauliString base axis st
            -- The axis is φ·P with P a bare Pauli string, so (axis)² = φ²·I
            -- and exp(−it·φP) = cos(tφ)·I − i·sin(tφ)·P — exact for any
            -- complex φ (complex-argument cos/sin), matching QPP.Semantics.Matrix's
            -- matrix exponential. For the meaningful axes φ = ±1 this reduces
            -- to the familiar c·I − i·s·φ·P.
            tphi = (t :+ 0) * phi
            iC = 0 :+ 1
            ψ1 = cos tphi .* st
            ψ2 = ((-iC) * sin tphi) .* pst
        in writer (addLocal iSupp ψ1 ψ2)

  C a -> do
    let ctrlP  = log2phys st ! base
        ctrlIv = singletonIval ctrlP
        iA     = if S.null (op_support a)
                   then ctrlIv                          -- empty support: stay at control bit
                   else supportInterval st (base+1) a
        iHull  = hull ctrlIv iA
    psi0 <- writer (projectCtrl ctrlP False st)
    psi1 <- writer (projectCtrl ctrlP True st) >>= evalOpAt (base+1) a
    writer (addLocal iHull psi0 psi1)

  DirectSum a b -> do
    let ctrlP  = log2phys st ! base
        ctrlIv = singletonIval ctrlP
        iA     = if S.null (op_support a) then ctrlIv else supportInterval st (base+1) a
        iB     = if S.null (op_support b) then ctrlIv else supportInterval st (base+1) b
        iHull  = hull ctrlIv (hull iA iB)
    psi0 <- writer (projectCtrl ctrlP False st) >>= evalOpAt (base+1) a
    psi1 <- writer (projectCtrl ctrlP True  st) >>= evalOpAt (base+1) b
    writer (addLocal iHull psi0 psi1)

-- Steps / programs: shared interpreter from QPP.Semantics.
evalStep :: (StateT, Outcomes, RNG) -> Step -> (StateT, Outcomes, RNG)
evalStep = evalStepWith (apply . evalOp) measure1

evalProg :: Program -> StateT -> RNG -> (StateT, Outcomes, RNG)
evalProg = evalProgWith (apply . evalOp) measure1

-- | `evalProg` with the run's truncation statistics (unitaries *and* the SVD
--   work inside measurements).
evalProgProfiled :: Program -> StateT -> RNG -> ((StateT, Outcomes, RNG), Profile)
evalProgProfiled = evalProgWithW (applyProfiled . evalOp) measure1Profiled
