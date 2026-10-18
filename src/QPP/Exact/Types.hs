-- | Configuration and explicit failures for the exact prototype.
module QPP.Exact.Types (Limits(..), defaultLimits, ExactError(..), Result, orFail) where

import GHC.Stack (HasCallStack, withFrozenCallStack)
import QPP.TypeCheck (TypeError)

type Result a = Either ExactError a

-- | The conventional backend API reports failures at the caller. Call the
-- corresponding checked function when the program needs to recover instead.
orFail :: HasCallStack => String -> Result a -> a
orFail operation = either
  (\err -> withFrozenCallStack (error ("QPP.Exact." ++ operation ++ ": " ++ show err)))
  id

-- | Representation limits, not wall-time or GHC heap guarantees. Workspace
-- estimates include temporary scalar and evaluator buffers and are deliberately
-- conservative. Validation checks stored coordinates only; operations reserve
-- their own workspace. Raise the limits explicitly for larger fields or circuits.
data Limits = Limits
  { maxRootIndex :: Integer
  , maxFieldDegree :: Int
  , maxCoefficientSlots :: Integer
  , maxIntegerBits :: Int
  } deriving (Eq, Show)

defaultLimits :: Limits
defaultLimits = Limits
  { maxRootIndex = 4096
  , maxFieldDegree = 256
  , maxCoefficientSlots = 2000000
  , maxIntegerBits = 4096
  }

data ExactError
  = InvalidInput String
  | InvalidOperator TypeError
  | WidthMismatch Int Int
  | DivisionByZero
  | ResourceLimit String Integer Integer -- ^ resource, requested, allowed
  | AtStep Int ExactError
  deriving (Eq, Show)
