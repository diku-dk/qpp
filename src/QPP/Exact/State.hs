-- | Checked linear algebra on exact, possibly unnormalized statevectors.
-- These are the field-compatible operations of HilbertSpace. Norm and
-- normalization are deliberately absent: they need not remain cyclotomic.
module QPP.Exact.State
  ( StateT, Limits(..), defaultLimits, ExactError(..), Result
  , ket, fromAmplitudes, qubits, fieldIndex, amplitudes, amplitude
  , scale, add, sub, tensor, inner, normSquared, equal
  ) where

import QPP.Exact.Types
import QPP.Exact.Internal.State
