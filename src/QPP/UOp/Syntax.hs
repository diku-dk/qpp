-- | Symbolic unitary operators on arbitrary finite-dimensional spaces.
--
-- Dimensions count basis vectors, not qubits; dimension zero is allowed.
-- Like "QPP.Syntax", this module builds syntax, without evaluating or
-- type-checking it. Composition dimensions, permutations, and rotation
-- axes must satisfy the conditions documented on the constructors.
module QPP.UOp.Syntax
  ( UOp(..), Dim, Permutation
  , HasTensorProduct(..), HasDirectSum(..), HasAdjoint(..), Operator(..)
  , identity, phase, tensor, directSum, compose
  , dim, dagger
  ) where

import Numeric.Natural (Natural)
import QPP.Syntax
  ( HasTensorProduct(..), HasDirectSum(..), HasAdjoint(..), Operator(..) )
import QPP.Util (permApply, permInvert)

--------------------------------------------------------------------------------
-- Core syntax
--------------------------------------------------------------------------------

-- | Matrix dimension, with no machine-integer overflow.
type Dim = Natural

-- | Output position j selects input position p[j]. The list must be a
-- permutation of [0 .. k-1], where k is the number of factors or blocks.
type Permutation = [Int]

data UOp
  = R UOp Rational                   -- ^ exp(-i*pi*t*A/2), as in QOp; A must be Hermitian unless t = 0.
  | Tensor [UOp]                     -- ^ Kronecker product; the empty list denotes I_1.
  | DirectSum [UOp]                  -- ^ Block diagonal; the empty list denotes I_0.
  | Compose Dim [UOp]                -- ^ Matrix product, rightmost first; all operands have this dimension.
  | PermuteTensor [Dim] Permutation  -- ^ Reorder tensor factors with the given input dimensions.
  | PermuteSum [Dim] Permutation     -- ^ Reorder direct-sum blocks with the given input dimensions.
  deriving (Show, Eq)

--------------------------------------------------------------------------------
-- Shared operator notation
--------------------------------------------------------------------------------

instance HasTensorProduct UOp where a ⊗ b = tensor [a, b]
instance HasDirectSum UOp where a ⊕ b = directSum [a, b]
instance HasAdjoint UOp where adj = dagger
instance Semigroup UOp where a <> b = compose (dim a) [a, b]
instance Operator UOp

--------------------------------------------------------------------------------
-- Convenient construction
--------------------------------------------------------------------------------

-- | Identity on C^d. In particular, I_0 is the direct-sum unit and I_1
-- is the tensor unit. The empty composition needs an explicit dimension.
identity :: Dim -> UOp
identity d = Compose d []

-- | Scalar phase exp(i*pi*t), represented using the one-dimensional identity.
-- The phase angle is kept exactly; construction does not normalize rotations.
phase :: Rational -> UOp
phase t = R (identity 1) (-2*t)

-- | Flatten nested tensors while preserving operand order.
tensor :: [UOp] -> UOp
tensor = Tensor . concatMap factors
  where
    factors (Tensor xs) = concatMap factors xs
    factors x           = [x]

-- | Flatten nested direct sums, including unequal and zero-sized blocks.
directSum :: [UOp] -> UOp
directSum = DirectSum . concatMap blocks
  where
    blocks (DirectSum xs) = concatMap blocks xs
    blocks x              = [x]

-- | Flatten compositions of the same dimension. All operands must have
-- dimension d; a mismatched annotation is retained for later validation.
compose :: Dim -> [UOp] -> UOp
compose d = Compose d . concatMap factors
  where
    factors (Compose e xs)
      | d == e = concatMap factors xs
    factors x = [x]

--------------------------------------------------------------------------------
-- Syntactic queries
--------------------------------------------------------------------------------

-- | Matrix dimension, assuming the term is well formed. As with QOp's
-- width query, this does not check composition operands or permutations.
dim :: UOp -> Dim
dim op = case op of
  R a _              -> dim a
  Tensor xs          -> product (map dim xs)
  DirectSum xs       -> sum (map dim xs)
  Compose d _        -> d
  PermuteTensor ds _ -> product ds
  PermuteSum ds _    -> sum ds

--------------------------------------------------------------------------------
-- Structural adjoint
--------------------------------------------------------------------------------

-- | Adjoint of a well-formed operator. Rotation axes satisfy the
-- Hermiticity condition above, and permutation indices must be valid.
dagger :: UOp -> UOp
dagger op = case op of
  R a t       -> R a (-t)  -- The Hermitian axis stays unchanged.
  Tensor xs     -> Tensor (map dagger xs)
  DirectSum xs  -> DirectSum (map dagger xs)
  Compose d xs  -> Compose d (map dagger (reverse xs))
  PermuteTensor ds p ->
    PermuteTensor (permApply p ds) (permInvert p)
  PermuteSum ds p ->
    PermuteSum (permApply p ds) (permInvert p)
    -- The inverse starts from the reordered dimensions, not the original list.
