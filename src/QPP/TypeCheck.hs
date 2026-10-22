-- | Symbolic well-formedness and Hermiticity proofs for QOp.
--   Complete for phase-labelled Pauli expressions; conservative elsewhere.
--   No matrices or exponentially large Pauli expansions are constructed.
module QPP.TypeCheck (TypeError(..), typeCheck, adjointPhase) where

import Control.Applicative ((<|>))
import Control.Monad (guard, when)
import Data.List (sort)
import Data.Maybe (isJust)
import Data.Ratio (denominator)
import qualified Data.Set as S
import QPP.Rewrite (fuseOps)
import QPP.Syntax
import qualified QPP.TypeCheck.Pauli as P
import QPP.Util (permInvert)

data TypeError
  = NegativeQubitCount Nat
  | InvalidPermutation [Nat]
  | ArityMismatch QOp Nat Nat
  | UnprovenHermitian QOp
  | QubitCountOverflow QOp
  deriving (Show, Eq)

-- | Check every subterm, returning its number of qubits or the first error.
--   A rotation needs a valid axis, and a Hermiticity proof when its angle is nonzero.
--   This checks QOp terms, not Step/Program or backend capabilities.
--
--   >>> typeCheck (C X)
--   Right 2
--   >>> typeCheck (DirectSum X (Id 2))
--   Left (ArityMismatch (DirectSum X (Id 2)) 1 2)
typeCheck :: QOp -> Either TypeError Nat
typeCheck = fmap qubits . check -- Extract the width from Right; propagate any Left error.

-- | Check a term and try to prove A† = exp(i*pi*p) A, with 0 <= p < 2.
--   Just 0 proves Hermiticity; any other Just p proves non-Hermiticity.
--   Nothing means unknown, not non-Hermitian. Rational arithmetic is exact.
--
--   >>> adjointPhase (Compose X Y)
--   Right (Just (1 % 1))
adjointPhase :: QOp -> Either TypeError (Maybe Rational)
adjointPhase = fmap (phaseProof . evidence) . check

-- Check children, enforce the typing premises, then construct the result.
-- The checked... helpers infer properties from the already checked children.
check :: QOp -> Either TypeError Checked
check op = case op of
  Id n
    | n >= 0    -> leaf n
    | otherwise -> Left (NegativeQubitCount n)
  Phase _ -> leaf 0 -- A global phase is a scalar, so it has no qubit wires.
  X       -> leaf 1
  Y       -> leaf 1
  Z       -> leaf 1
  H       -> leaf 1
  SX      -> leaf 1

  R axis angle -> do
    a <- check axis -- Validate every subterm before using its algebraic properties.
    when (angle /= 0) (requireHermitian axis a) -- At exactly zero, exp(0 A) = I for any valid axis.
    pure (checkedRotation a angle)

  C a -> do
    a' <- check a
    n <- addQubits op 1 (qubits a') -- Add the control wire.
    pure (checkedControl n a')

  Permute ks
    | sort ks == [0 .. length ks - 1] -> leaf (length ks) -- Every wire index occurs exactly once.
    | otherwise -> Left (InvalidPermutation ks)

  Tensor a b -> do
    a' <- check a
    b' <- check b
    n <- addQubits op (qubits a') (qubits b')
    pure (checkedTensor n a' b')

  DirectSum a b -> do
    a' <- check a
    b' <- check b
    requireSameArity op a' b'
    n <- addQubits op 1 (qubits a') -- Selecting between two equal blocks adds one qubit.
    pure (checkedDirectSum n a' b')

  Compose a b -> do
    a' <- check a
    b' <- check b
    requireSameArity op a' b'
    pure (checkedCompose a' b')

  Adjoint a -> do
    a' <- check a
    pure (checkedAdjoint a')
  where
    leaf n = pure (checkedLeaf n op)

requireSameArity :: QOp -> Checked -> Checked -> Either TypeError ()
requireSameArity op a b
  | qubits a == qubits b = Right () -- Success; this constraint produces no new value.
  | otherwise = Left (ArityMismatch op (qubits a) (qubits b))

-- Used for nonzero-angle rotations, after checking the axis's well-formedness.
requireHermitian :: QOp -> Checked -> Either TypeError ()
requireHermitian axis a
  | phaseProof (evidence a) == Just 0 = Right ()
  | otherwise = Left (UnprovenHermitian axis) -- Includes an unknown phase, not just a nonzero one.

-- The child widths are nonnegative, but their sum must also fit in Int.
addQubits :: QOp -> Nat -> Nat -> Either TypeError Nat
addQubits op m n
  | m > maxBound - n = Left (QubitCountOverflow op) -- Compare before adding, to avoid overflow.
  | otherwise        = Right (m + n)

-- Internal analysis: evidence and construction of checked results.

-- ExactPauli knows the whole operator exp(i*pi*q) P. AdjointPhase only
-- knows A† = exp(i*pi*p) A; H, for example, has AdjointPhase 0.
-- Unknown means no proof was found, not that the operator is non-Hermitian.
data Evidence = Unknown | AdjointPhase Rational | ExactPauli P.Form

data Checked = Checked
  { qubits :: Nat
  , evidence :: Evidence
  , view :: View
  }

-- Product rules inspect nested compositions and rotations. Keep their
-- checked children so proofs can be reused; other nodes need only their
-- original syntax. Each view has exactly the children its constructor needs.
data View = Other QOp | Rotation Checked Rational | Product Checked Checked

-- These constructors use facts from already checked operands. Validation,
-- including width overflow and Hermiticity of nonzero-angle axes, belongs to check.
checkedLeaf :: Nat -> QOp -> Checked
checkedLeaf n op = Checked n proof (Other op)
  where
    proof = case op of
      Id _    -> ExactPauli (P.scalar 0) -- Phase zero represents +I; width is stored separately.
      Phase q -> ExactPauli (P.scalar q)
      X       -> ExactPauli (P.letter P.PX)
      Y       -> ExactPauli (P.letter P.PY)
      Z       -> ExactPauli (P.letter P.PZ)
      H       -> AdjointPhase 0 -- Hermitian, although it is not a phase times a Pauli string.
      Permute ks
        | ks == [0 .. n - 1] -> ExactPauli (P.scalar 0) -- No wires move: the operator is I.
        | ks == permInvert ks -> AdjointPhase 0 -- For a permutation matrix, P† = P⁻¹.
      _ -> Unknown -- SX and non-involutive permutations have no phase proof.

checkedRotation :: Checked -> Rational -> Checked
checkedRotation a angle =
  Checked (qubits a) (rotationEvidence (evidence a) angle) (Rotation a angle)

-- C A = I ⊕ A, so control uses the direct-sum evidence rule.
checkedControl :: Nat -> Checked -> Checked
checkedControl n a = Checked n proof (Other (C (term a)))
  where
    proof = sumEvidence (ExactPauli (P.scalar 0)) (evidence a)

checkedTensor :: Nat -> Checked -> Checked -> Checked
checkedTensor n a b = Checked n (preferPauli form phase) (Other (Tensor (term a) (term b)))
  where
    ea = evidence a
    eb = evidence b
    form = P.tensor (qubits a) <$> pauliForm ea <*> pauliForm eb -- <*> requires both forms to be known.
    phase = (+) <$> phaseProof ea <*> phaseProof eb -- Adjoint phases add under tensor product.

checkedDirectSum :: Nat -> Checked -> Checked -> Checked
checkedDirectSum n a b =
  Checked n (sumEvidence (evidence a) (evidence b)) (Other (DirectSum (term a) (term b)))

checkedCompose :: Checked -> Checked -> Checked
checkedCompose a b = Checked (qubits a) (composeEvidence a b) (Product a b)

checkedAdjoint :: Checked -> Checked
checkedAdjoint a = Checked (qubits a) proof (Other (Adjoint (term a)))
  where
    proof = case evidence a of
      Unknown        -> Unknown
      AdjointPhase p -> AdjointPhase (P.normalPhase (-p)) -- Taking the adjoint inverts its phase.
      ExactPauli p   -> ExactPauli (P.adjoint p) -- Conjugate the scalar; the Pauli string is Hermitian.

-- Prefer an exact Pauli form whenever one is found. Otherwise keep the
-- weaker phase proof, if any. A nonzero phase is a proof, not a failure.
-- This is the common boundary between exact and conservative analysis.
preferPauli :: Maybe P.Form -> Maybe Rational -> Evidence
preferPauli (Just p) _       = ExactPauli p -- Exact information already determines the adjoint phase.
preferPauli Nothing (Just p) = AdjointPhase (P.normalPhase p) -- Keep only the phase, reduced modulo 2.
preferPauli Nothing Nothing  = Unknown -- Neither route supplied a proof.

pauliForm :: Evidence -> Maybe P.Form
pauliForm (ExactPauli p) = Just p
pauliForm _             = Nothing

-- For A = exp(i*pi*q) P, A† = exp(-2*i*pi*q) A. Derive this phase from
-- the Pauli form rather than storing a second, redundant piece of evidence.
phaseProof :: Evidence -> Maybe Rational
phaseProof Unknown          = Nothing
phaseProof (AdjointPhase p) = Just p
phaseProof (ExactPauli p)   = Just (P.adjointPhase p) -- Scalar phase q gives adjoint phase -2q.

sumEvidence :: Evidence -> Evidence -> Evidence
sumEvidence a b = preferPauli form phase
  where
    -- Equal or opposite Pauli blocks give I tensor P or Z tensor P.
    form = do
      p <- pauliForm a
      q <- pauliForm b
      P.directSum p q -- Succeeds only when the two Pauli blocks are equal or opposite.
    -- Blocks need the same adjoint phase; their phases cannot cancel.
    -- For control, the identity block therefore forces phase 0.
    phase = do
      p <- phaseProof a
      q <- phaseProof b
      guard (p == q) -- In Maybe, a failed guard returns Nothing and ends this proof attempt.
      pure p

-- Zero angle gives I for any valid axis. For nonzero q, A is checked Hermitian,
-- so R_A(q) = cos(pi*q/2) I - i sin(pi*q/2) A.
-- A general angle gives a sum of I and A; we do not expand such sums.
rotationEvidence :: Evidence -> Rational -> Evidence
rotationEvidence a q
  | q == 0 = ExactPauli (P.scalar 0) -- This case needs no Hermiticity proof for A.
  -- Even q gives (-1)^(q/2) I, regardless of the Hermitian axis.
  | integral (q / 2) = ExactPauli (P.scalar (-q / 2))
  -- A scalar Hermitian axis is +/-I, so s is 0 or 1. Any angle is then scalar.
  | Just s <- pauliForm a >>= P.scalarPhase = -- First find a Pauli form, then test whether it is scalar.
      ExactPauli (P.scalar (if s == 0 then -q / 2 else q / 2))
  -- The remaining integral angles are odd: R_A(q) = exp(-i*pi*q/2) A.
  | integral q = case a of
      ExactPauli p -> ExactPauli (P.scale (-q / 2) p)
      _            -> AdjointPhase (P.normalPhase q) -- Odd q still proves R† = -R.
  | otherwise = Unknown -- No exact form or adjoint-phase proof was found.

composeEvidence :: Checked -> Checked -> Evidence
composeEvidence a b
  -- Multiplying known Pauli forms includes phases such as XY = iZ.
  | Just p <- pauliForm ea, Just q <- pauliForm eb = ExactPauli (P.compose p q)
  -- Even when the operands are not Pauli, their product may be.
  | otherwise = preferPauli (pauliForm collapsed) phase -- Try the collapsed node's exact form, then a phase-only proof.
  where
    ea = evidence a
    eb = evidence b
    -- <|> tries another proof only on Nothing. Some rules do not need a
    -- phase proof for each operand: for example, X R_Z(t) is Hermitian.
    phase = conjugationPhase a b <|> phaseProof collapsed <|> commuting
        <|> reflection a b <|> reflection b a -- Cover both A R_K(t) and R_K(t) A.
    -- A product that a gate identity collapses to one node (A A† = I,
    -- SX SX = X, same-axis rotations, permutations) takes that node's
    -- evidence. A fused rotation reuses the checked axis: check builds every
    -- R with checkedRotation, so a's view is Rotation whenever its term is an
    -- R. The other results are leaves, with the ordinary leaf proofs.
    collapsed = case fuseOps (term a) (term b) of
      Just (R _ t) | Rotation axis _ <- view a -> rotationEvidence (evidence axis) t
      Just c  -> evidence (checkedLeaf (qubits a) c)
      Nothing -> Unknown
    -- In the commuting case, both operand phases must be known:
    -- (AB)† = exp(i*pi*(p_a+p_b+w)) AB.
    commuting = do
      p <- phaseProof ea
      q <- phaseProof eb
      w <- commutation a b -- 0 for commuting, 1 for anticommuting.
      pure (p + q + w) -- Reversing the factors in (AB)† contributes w.

-- U A U† has A's adjoint phase, even when U has no phase proof of its own.
-- These patterns cover U (A U†) and (U A) U† without reassociating the tree.
conjugationPhase :: Checked -> Checked -> Maybe Rational
conjugationPhase a b = case (view a, view b) of
  (_, Product middle u) | adjointPair (term a) (term u) -> phaseProof (evidence middle) -- U (A U†).
  (Product u middle, _) | adjointPair (term u) (term b) -> phaseProof (evidence middle) -- (U A) U†.
  _ -> Nothing

-- If A anticommutes with K, then A R_K(t) A† = R_K(t)†.
-- Thus A R_K(t), and R_K(t) A, inherit A's adjoint phase, for any t.
reflection :: Checked -> Checked -> Maybe Rational
reflection a b = case view b of
  Rotation axis _ -> do
    guard (commutation a axis == Just 1) -- Phase 1 means AK = -KA.
    phaseProof (evidence a)
  _ -> Nothing

-- BA = exp(i*pi*w) AB. The rules here only return w = 0 or 1:
-- exact Pauli commutation, or sufficient structural proofs of commutation.
commutation :: Checked -> Checked -> Maybe Rational
commutation a b
  | Just p <- pauliForm (evidence a), Just q <- pauliForm (evidence b) =
      Just (P.commutationPhase p q)
  | term a == term b = Just 0
  | isScalar a || isScalar b = Just 0 -- Scalars commute with every operator.
  | S.disjoint (op_support (term a)) (op_support (term b)) = Just 0 -- Disjoint wires commute.
  | otherwise = Nothing
  where
    isScalar x = isJust (pauliForm (evidence x) >>= P.scalarPhase)

-- Recover the source term for equality, adjoint and support checks. Keeping
-- syntax in the view avoids a separate term/children pair that must agree.
term :: Checked -> QOp
term a = case view a of
  Other op       -> op
  Rotation b q   -> R (term b) q
  Product b c    -> Compose (term b) (term c)

integral :: Rational -> Bool
integral q = denominator q == 1 -- Rational values are stored as reduced fractions.
