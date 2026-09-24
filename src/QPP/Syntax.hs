-- | The symbolic layer: the operator interface classes, the `QOp` AST of
--   unitary operators, the qubit-placement sugar, quantum programs with
--   measurement (`Step`/`Program`), and the purely syntactic analysis
--   functions over them.
--
--   Nothing here evaluates anything — semantics live in "QPP.Semantics" and
--   the backend modules below it, rewriting in "QPP.Rewrite".
module QPP.Syntax where
import Data.Complex
import Data.Bits(shiftL)
import qualified Data.Set as S
import QPP.Util(permInvert, permSupport)

type Nat = Int
type RealT = Double  -- Can be replaced by e.g. exact fractions or constructive reals
type ComplexT = Complex RealT

--------------------------------------------------------------------------------
-- Interface classes
--------------------------------------------------------------------------------
-- Shared by the symbolic QOp and by every concrete semantics (matrices,
-- state vectors, tensor networks, stabilizers), so the same program text
-- works against any of them.

class HasTensorProduct o where
  (⊗) :: o -> o -> o

  (<.>) :: o -> o -> o
  (<.>) = (⊗)

class HasDirectSum o where
  (⊕) :: o -> o -> o

  (<+>) :: o -> o -> o
  (<+>) = (⊕)

class HasAdjoint o where adj :: o -> o
class HasQubits o where n_qubits :: o -> Nat

{-| The Operator type class allows us to work with both the Op symbolic operators, and concrete semantics (e.g. matrices, tensor networks, stabilizers) using the same syntax. I.e., no matter which representation we're working with, we can use the same code to compose, tensor, take adjoints, etc.

The operators form a Semigroup under composition, so we inherit Haskell's standard composition operator <>. Note that a<>b means "first apply b, then a".
-}
class (Semigroup o, HasTensorProduct o, HasDirectSum o, HasAdjoint o) => Operator o where
  -- Compose is semigroup operator <> with synonym ∘ for math order and >: for left-to-right.
  -- Direct sum: UTF ⊕, ASCII <+>
  -- Tensor product: UTF ⊗, ASCII <.>

  -- Syntactic sugar in Unicode and ASCII
  (∘),(>:) :: o -> o -> o
  (∘)   = (<>)         -- right-to-left composition (math operator order)
  (>:) a b = b ∘ a   -- left-to-right composition


{-| We define a HilbertSpace typeclass, which we will use for states.
    Tensor product, direct sum, adjoint, and composition are inherited from Operator.
 -}
class (Scalar v ~ Complex (Realnum v), Floating (Realnum v), HasTensorProduct v)
    => HilbertSpace v where
  type Realnum v
  type Scalar  v

  (.*)  :: Scalar v -> v -> v -- Scalar-vector multiplication
  (.+)  :: v -> v -> v        -- Vector-vector addition
  (.-)  :: v -> v -> v        -- Vector-vector subtraction

  inner     :: v -> v -> Scalar v -- Inner product
  normalize :: v -> v

  norm  :: v -> Realnum v     -- Vector 2-norm
  norm x = sqrt(realPart $ inner x x)

infixr 8 ⊗, <.>, .*
infixr 7 ⊕, <+>, .+, .-
infixr 6 ∘, >:


--------------------------------------------------------------------------------
-- QOp algebra
--------------------------------------------------------------------------------

{-|
    The Op type is a symbolic unitary operator, which just builds an abstract syntax tree (AST).
    It provides building blocks for building any n-qubit unitary operator. Explanation of constructors is given below.
 -}
data QOp
  = Id Nat -- Identity n: C^{2^n} -> C^{2^n} is the n-qubit identity operator.
                 -- Identity 0: C^1 -> C^1 scalar multiplication by 1, unit for ⊗.
                 -- Identity 1 = I: C^2 -> C^2
                 -- Identity n is the family of units for ∘.
  | Phase Rational -- Global phase e^{i π θ} (scalar multiplication)
  | X | Y | Z | H | SX
  | R QOp Rational  -- Rotation around (possibly multi-qubit) axis defined by QOp by angle (in units of π)
  | C QOp           -- Controlled (possibly multi-qubit) operator
  | Permute [Int]
  | Tensor QOp QOp
  | DirectSum QOp QOp -- Direct sum of operators with same arity.
  | Compose QOp QOp
  | Adjoint QOp
  deriving (Show,Eq)

-- | The 0-qubit and 1-qubit identities, spelled out.
pattern One, I :: QOp
pattern One <- Id 0
  where One  = Id 0

pattern I <- Id 1
  where I  = Id 1

instance Semigroup QOp where
  (<>) = Compose

instance HasTensorProduct QOp where (⊗) = Tensor
instance HasDirectSum QOp     where (⊕) = DirectSum
instance HasAdjoint QOp       where adj = Adjoint
instance HasQubits QOp where  n_qubits op = op_qubits op
instance Operator QOp


--------------------------------------------------------------------------------
-- Qubit-placement sugar
--------------------------------------------------------------------------------
-- `k <@ op` puts op at qubit k; `op @> l` pads l idle qubits below it, so
-- `k <@ op @> l` places op at qubit k of a (k + op_qubits op + l)-qubit register.

-- | Syntactic sugar patterns
pattern AtQubit :: QOp -> Nat -> QOp
pattern AtQubit op n <- Tensor (Id n) op
  where AtQubit op n  = Tensor (Id n) op

(@>) :: QOp -> Nat -> QOp
(@>) op n = Tensor op (Id n)
(<@) n op = AtQubit op n
(<@) :: Nat -> QOp -> QOp

infixr 5 @>, <@


--------------------------------------------------------------------------------
-- Programs: unitaries, initialization and measurement
--------------------------------------------------------------------------------

{- Quantum programs including measurement. -}
data Step
  = Unitary QOp             -- A unitary quantum program
  | Initialize [Nat] [Bool] -- Initialize qubits qs to classical values vs.
  | Measure    [Nat] -- Measurement of qubits ks (stochastic non-reversible process)
  deriving (Show, Eq)

type Program = [Step]

type Outcomes = [Bool]     -- head = most recent
type RNG      = [Double]   -- infinite steam in [0,1)


--------------------------------------------------------------------------------
-- Analysis
--------------------------------------------------------------------------------
-- Purely syntactic queries and rewrites: how wide an operator or program is,
-- which qubits it touches, and its structural adjoint.

op_qubits :: QOp -> Nat
op_qubits op = case op of
    Id n          -> n
    Phase _       -> 0
    R a _         -> op_qubits a
    C a           -> 1 + op_qubits a
    Tensor    a b -> op_qubits a + op_qubits b
    DirectSum a _ -> 1 + op_qubits a -- Assume op_qubits a == op_qubits b is type checked
    Compose   a _ -> op_qubits a     -- Assume op_qubits a == op_qubits b is type checked
    Adjoint   a   -> op_qubits a
    Permute   ks  -> length ks
    _             -> 1 -- 1-qubit gates

-- | Signature of an operator a: C^{2^m} -> C^{2^n} is (m,n) = (op_domain a, op_range a)
op_dimension :: QOp -> Nat
op_dimension op = 1 `shiftL` (op_qubits op)

-- | Support of an operator: the list of qubits it acts non-trivially on.
op_support :: QOp -> S.Set Nat
op_support op = let
    shift ns k   = S.map (+k) ns
    union xs ys  = S.union xs ys
  in case op of
  Id _          -> S.empty
  Phase _       -> S.empty
  R _ 0         -> S.empty
  R a _         -> op_support a
  C a           -> S.insert 0 ((op_support a) `shift` 1)
  Tensor a b    -> (op_support a) `union` ((op_support b) `shift` (op_qubits a))
  DirectSum a b -> S.insert 0 ((op_support a `union` op_support b) `shift` 1)
  -- Compose's support is the union of its operands'. (Earlier special cases for
  -- Permute were unsound — e.g. `Compose (Permute ks) (Id n)` returned ∅ instead of
  -- permSupport ks, because the bits Permute moves count toward support even when the
  -- other side has empty support.)
  Compose a b   -> union (op_support a) (op_support b)
  Adjoint a     -> op_support a
  Permute ks    -> S.fromList $ permSupport ks
  _             -> S.singleton 0 -- 1-qubit gates

step_qubits :: Step -> Nat
step_qubits step = case step of
  Unitary op -> op_qubits op
  Measure ks      -> 1 + foldr max 0 ks
  Initialize ks _ -> 1 + foldr max 0 ks

prog_qubits :: Program -> Nat
prog_qubits program = maximum $ map step_qubits program

-- | Structural adjoint: rewrite a QOp to an equivalent one with no `Adjoint`
--   at the head. Shared by the backends that evaluate `Adjoint a` by
--   recursing on `dagger a`.
dagger :: QOp -> QOp
dagger = \case
  Id n          -> Id n
  Phase q       -> Phase (-q)
  X             -> X
  Y             -> Y
  Z             -> Z
  H             -> H
  SX            -> Compose SX X       -- SX² = X, so SX⁻¹ = SX³ = SX·X

  R a t         -> R a (-t)
  C a           -> C (dagger a)
  Permute ks    -> Permute (permInvert ks)
  Tensor a b    -> Tensor (dagger a) (dagger b)
  DirectSum a b -> DirectSum (dagger a) (dagger b)
  Compose a b   -> Compose (dagger b) (dagger a)
  Adjoint a     -> a
