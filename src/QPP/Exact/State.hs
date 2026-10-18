-- | Checked linear algebra on exact, possibly unnormalized statevectors.
-- These are the field-compatible operations of HilbertSpace. Norm and
-- normalization are deliberately absent: they need not remain cyclotomic.
-- Measurement retains an unnormalized projection; probability divides by the
-- stored squared norm. The linear operations still describe the stored vector.
module QPP.Exact.State
  ( StateT, Limits(..), defaultLimits, ExactError(..), Result
  , ket, ketWithLimits, fromAmplitudes, qubits, fieldIndex, amplitudes, amplitude
  , scale, add, sub, tensor, inner, normSquared, equal
  , probability, measure1, measure1WithLimits
  ) where

import QPP.Exact.Types
import GHC.Stack (HasCallStack, withFrozenCallStack)
import QPP.Exact.Internal.State hiding (ket)
import qualified QPP.Exact.Internal.State as Internal
import QPP.Exact.Internal.Measurement (probability, measure1, measure1WithLimits)

-- | Computational-basis state, using the same bit convention as every backend.
ket :: HasCallStack => [Int] -> StateT
ket bits = withFrozenCallStack $ orFail "ket" (ketWithLimits defaultLimits bits)

-- | Checked construction with explicit representation limits.
ketWithLimits :: Limits -> [Int] -> Result StateT
ketWithLimits = Internal.ket
