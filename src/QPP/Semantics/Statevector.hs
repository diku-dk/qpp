-- | The state-vector semantics: what an operator /does to a state/.
--
--   The matrix semantics ("QPP.Semantics.Matrix") evaluates an n-qubit operator
--   to its 2^n × 2^n matrix. This backend applies the operator to a state
--   instead and never builds that matrix: the largest matrix of an operator
--   that it forms is 2×2.
--
--   `applyOp` first regroups the amplitudes into one row for each value of the
--   operator's own qubits. Each rule is then a formula on the rows:
--
--   >  u ∘ v     apply v, then u
--   >  u ⊕ v     u on the upper half of the rows, v on the lower half
--   >  C u       u on the lower half of the rows
--   >  R u θ     cos(πθ/2) ψ − i sin(πθ/2) (u applied to ψ)
--   >  u ⊗ v     with one row for each value of u's qubits: v on every row, then u
--
--   A one-qubit gate multiplies the two rows by its 2×2 matrix, and `Permute`
--   permutes the rows. The gate rule is written out entry by entry, not handed
--   to a matrix product: each product is then rounded on its own, so that
--   amplitudes which cancel do so exactly. A measurement with the sample 0.0
--   relies on that (see the Deutsch case in @StatevectorTest@).
--
--   A state may have more qubits than the operator: `applyOp` acts on the
--   first ones and carries the rest along, so that for an operator @e@ with
--   matrix M(e) it computes (M(e) ⊗ I) ψ. The recursive calls rely on this.
--   The public `apply` asks for a state of exactly the operator's qubits.
--
--   Evaluated operators (`OpT`) have the operator interface of "QPP.Syntax":
--   `<>` and `∘`, `⊗`, `⊕`, and the vector-space operations. These combine
--   the operators' functions by the rules above. Only `adj`, `inner` and
--   `norm` of an evaluated operator form its matrix. To take an adjoint
--   without one, take it on the expression: @evalOp (adj e)@.
module QPP.Semantics.Statevector where

import QPP.Syntax
import QPP.Semantics
import QPP.Semantics.Matrix (CMat)
import qualified QPP.Semantics.Matrix as M
import Numeric.LinearAlgebra hiding (normalize, step, (<>))
import Data.Bits (bit, (.&.))

-- | States are shared with the matrix semantics: a column of 2^n amplitudes.
type StateT = CMat

-- | An operator is the function that applies it, together with its arity. The
--   function accepts a state of at least that many qubits, in any shape: it
--   acts on the first qubits and returns the shape it was given.
data OpT = OpT { opQubits :: Int, runOp :: CMat -> CMat }

instance HasQubits OpT where n_qubits = opQubits

ket :: [Bit] -> StateT
ket = M.ket

evalOp :: QOp -> OpT
evalOp e = OpT (op_qubits e) (applyOp e)

-- | Apply an operator to a state, or to several states side by side.
apply :: OpT -> StateT -> StateT
apply (OpT n f) ψ
  | rows ψ == 2^n = f ψ
  | otherwise     = error $ "A " ++ show n ++ "-qubit operator cannot act on "
                         ++ show (rows ψ) ++ " amplitudes"

-- | Applied to the identity, an operator returns its matrix; a matrix is
--   applied by multiplication.
instance Convertible OpT CMat where
  to (OpT n f) = f (ident (2^n))
  from m       = onRows (n_qubits m) (m <>)

-- | @applyOp e φ@ applies @e@ to the amplitudes in φ, and returns them in the
--   shape of φ. The rules work on ψ, which holds the same amplitudes with one
--   row for each value of the qubits of @e@.
applyOp :: QOp -> CMat -> CMat
applyOp e φ = rowsAs (rows φ) $ case e of
    Id _          -> ψ
    Phase θ       -> scale (cis (pi * fromRational θ)) ψ
    Permute ks    -> ψ ? map (sourceRow ks) [0 .. rows ψ - 1]
    R u θ         -> scale (cos a :+ 0) ψ - scale (0 :+ sin a) (applyOp u ψ)
      where a = pi * fromRational θ / 2
    C u           -> ψ0 `above` applyOp u ψ1
    DirectSum u v -> applyOp u ψ0 `above` applyOp v ψ1
    Compose u v   -> applyOp u (applyOp v ψ)
    Adjoint u     -> applyOp (dagger u) ψ
    Tensor u v    -> applyOp u (tr' (applyOp v (tr' (rowsAs d ψ))))
      where d = 2 ^ op_qubits u
    _             -> (g 0 0 ψ0 + g 0 1 ψ1) `above` (g 1 0 ψ0 + g 1 1 ψ1)
      where g i j = scale (M.evalOp e `atIndex` (i, j))     -- X, Y, Z, H, SX
  where
    ψ  = rowsAs (2 ^ op_qubits e) φ
    h  = rows ψ `div` 2
    ψ0 = takeRows h ψ
    ψ1 = dropRows h ψ

-- | The same entries, in the same row-by-row order, in d rows.
rowsAs :: Int -> CMat -> CMat
rowsAs d m = reshape (rows m * cols m `div` d) (flatten m)

-- | One matrix above another.
above :: CMat -> CMat -> CMat
above top bottom = reshape (cols top) (vjoin [flatten top, flatten bottom])

-- | @Permute ks@ sends |b₀…bₙ₋₁⟩ to |b_{k₀}…b_{kₙ₋₁}⟩. Row x of the result is
--   therefore row @sourceRow ks x@ of ψ: the one that has qubit kⱼ set exactly
--   when x has qubit j set.
sourceRow :: [Int] -> Int -> Int
sourceRow ks x = sum [ qubit k | (j, k) <- zip [0..] ks, x .&. qubit j /= 0 ]
  where qubit i = bit (length ks - 1 - i)     -- qubit i of a row index

-- | The operator on n qubits that maps the 2^n rows of a state by f.
onRows :: Int -> (CMat -> CMat) -> OpT
onRows n f = OpT n (\φ -> rowsAs (rows φ) (f (rowsAs (2^n) φ)))

-- | The number of qubits of two operators that are composed, summed or added.
sameQubits :: Int -> Int -> Int
sameQubits n m
  | n == m    = n
  | otherwise = error $ "Operators on " ++ show n ++ " and " ++ show m
                     ++ " qubits cannot be composed, summed or added"

-- | @u <> v@ is u ∘ v: v, then u.
instance Semigroup OpT where
  OpT n f <> OpT m g = OpT (sameQubits n m) (f . g)

-- | u ⊕ v: u on the upper half of the rows, v on the lower half, as in `applyOp`.
instance HasDirectSum OpT where
  OpT n f ⊕ OpT m g = onRows (k + 1) (\ψ -> f (takeRows (2^k) ψ) `above` g (dropRows (2^k) ψ))
    where k = sameQubits n m

-- | u ⊗ v: with one row for each value of u's qubits, v on every row, then
--   u, as in `applyOp`.
instance HasTensorProduct OpT where
  OpT n f ⊗ OpT m g = onRows (n + m) (\ψ -> f (tr' (g (tr' (rowsAs (2^n) ψ)))))

-- | The adjoint of an evaluated operator is taken through its matrix. For an
--   expression, evaluate its adjoint instead: @evalOp (adj e)@ forms no matrix.
instance HasAdjoint OpT where
  adj u = from (adj (to u :: CMat))

instance Operator OpT

-- | Operators also form a vector space. The inner product and the norm are
--   those of the operators' matrices, which are formed for the purpose.
instance HilbertSpace OpT where
  type Scalar  OpT = ComplexT
  type Realnum OpT = RealT
  c .* OpT n f       = OpT n (scale c . f)
  OpT n f .+ OpT m g = OpT (sameQubits n m) (\ψ -> f ψ + g ψ)
  OpT n f .- OpT m g = OpT (sameQubits n m) (\ψ -> f ψ - g ψ)
  inner u v          = inner (to u :: CMat) (to v)
  normalize u        = (1 / (norm u :+ 0)) .* u

-- | The projection onto "qubit k reads c" is (I + (−1)^c Z_k) / 2.
measureProjection :: Int -> Int -> Bit -> OpT
measureProjection n k c = OpT n $ \ψ -> scale 0.5 (ψ + scale sign (zk ψ))
  where sign = (-1) ^ fromEnum c
        zk   = applyOp (k <@ Z @> (n-k-1))     -- Z on qubit k of n

measure1 :: (StateT, Outcomes, RNG) -> Nat -> (StateT, Outcomes, RNG)
measure1 (_, _, []) _ = error "Random numbers exhausted."
measure1 (state, outcomes, r:rng) k =
  let n       = n_qubits state
      s0      = apply (measureProjection n k 0) state
      s1      = apply (measureProjection n k 1) state
      prob0   = realPart (inner state s0)
      prob1   = realPart (inner state s1)
      outcome = Bit (r >= prob0)
  in if abs (prob0 + prob1 - 1) > tol
       then error $ "Probabilities don't sum to 1: " ++ show (prob0, prob1)
       else (normalize (if outcome == 1 then s1 else s0), outcome : outcomes, rng)

evalStep :: (StateT, Outcomes, RNG) -> Step -> (StateT, Outcomes, RNG)
evalStep = evalStepWith (apply . evalOp) measure1

evalProg :: Program -> StateT -> RNG -> (StateT, Outcomes, RNG)
evalProg = evalProgWith (apply . evalOp) measure1
