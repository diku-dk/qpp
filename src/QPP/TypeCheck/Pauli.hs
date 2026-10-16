-- | Exact arithmetic for phase-labelled Pauli strings, used by the checker.
--   A form denotes exp(i*pi*q) times a tensor product of I, X, Y and Z.
--   Identity factors are implicit: storage depends on the expression, not
--   on the number of qubits in an Id node. Widths are checked by TypeCheck.
module QPP.TypeCheck.Pauli
  ( Form(..), Letter(..), normalPhase, scalar, letter, scalarPhase
  , scale, adjoint, tensor, compose, directSum, adjointPhase, commutationPhase
  ) where

import qualified Data.IntMap.Strict as IM

data Letter = PX | PY | PZ deriving (Eq, Show)

data Form = Form
  { phase :: Rational
  , letters :: IM.IntMap Letter
  } deriving (Eq, Show)

-- | Angles are in units of pi, so phases have period 2.
normalPhase :: Rational -> Rational
normalPhase q = q - 2 * fromInteger (floor (q / 2))

scalar :: Rational -> Form
scalar q = Form (normalPhase q) IM.empty

letter :: Letter -> Form
letter p = Form 0 (IM.singleton 0 p)

scalarPhase :: Form -> Maybe Rational
scalarPhase (Form q ps)
  | IM.null ps = Just q
  | otherwise  = Nothing

scale :: Rational -> Form -> Form
scale q (Form r ps) = Form (normalPhase (q + r)) ps

adjoint :: Form -> Form
adjoint (Form q ps) = Form (normalPhase (-q)) ps

-- | The left width gives the starting index of the right operand.
tensor :: Int -> Form -> Form -> Form
tensor n (Form q ps) (Form r qs) =
  Form (normalPhase (q + r)) (IM.union ps (IM.mapKeysMonotonic (+n) qs))

compose :: Form -> Form -> Form
compose (Form q ps) (Form r qs) = IM.foldlWithKey' multiplyAt (Form (normalPhase (q + r)) ps) qs
  where
    multiplyAt (Form s rs) k b = case IM.lookup k rs of
      Nothing -> Form s (IM.insert k b rs)
      Just a  -> let (t, c) = multiply a b
                 in Form (normalPhase (s + t))
                      (maybe (IM.delete k rs) (\p -> IM.insert k p rs) c)

-- XY = iZ, YZ = iX, ZX = iY; reversing the order negates the phase.
multiply :: Letter -> Letter -> (Rational, Maybe Letter)
multiply PX PX = (0, Nothing)
multiply PY PY = (0, Nothing)
multiply PZ PZ = (0, Nothing)
multiply PX PY = ( 1/2, Just PZ)
multiply PY PZ = ( 1/2, Just PX)
multiply PZ PX = ( 1/2, Just PY)
multiply PY PX = (-1/2, Just PZ)
multiply PZ PY = (-1/2, Just PX)
multiply PX PZ = (-1/2, Just PY)

-- | A direct sum is a Pauli string when the blocks are equal or opposite.
--   diag(P,P) = I tensor P; diag(P,-P) = Z tensor P.
directSum :: Form -> Form -> Maybe Form
directSum (Form q ps) (Form r qs)
  | ps /= qs = Nothing
  | d == 0   = Just (Form q shifted)
  | d == 1   = Just (Form q (IM.insert 0 PZ shifted))
  | otherwise = Nothing
  where
    d = normalPhase (r - q)
    shifted = IM.mapKeysMonotonic (+1) ps

-- | A† = exp(i*pi*p) A, with p = -2q because every bare Pauli string is Hermitian.
adjointPhase :: Form -> Rational
adjointPhase = normalPhase . (* (-2)) . phase

-- | AB = exp(i*pi*w) BA: w is 0 for commuting and 1 for anticommuting strings.
--   Each shared position with different nonidentity letters contributes -1.
commutationPhase :: Form -> Form -> Rational
commutationPhase (Form _ ps) (Form _ qs) =
  if odd (IM.size (IM.filter id (IM.intersectionWith (/=) ps qs))) then 1 else 0
