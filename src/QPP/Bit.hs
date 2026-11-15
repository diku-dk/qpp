-- | Classical bits.
--
--   A bit is written 0 or 1, as in the mathematics: @ket [0,0,1]@,
--   @Initialize [0,2] [0,1]@, and measurement outcomes that print as @[0,1]@.
--   This module depends on nothing else in QPP.
module QPP.Bit (Bit(..)) where

import Control.DeepSeq (NFData)

-- | A classical bit. The literals 0 and 1 are bits, and arithmetic is modulo 2:
--   `+` is exclusive or and `*` is and. Any other literal is an error.
--
--   The constructor converts from `Bool` and the pattern @Bit b@ converts back;
--   `fromEnum` gives the `Int` 0 or 1, and `fromIntegral` turns an integer that
--   is 0 or 1 into a bit.
--
--   The sum of a list of bits is its parity, not the number of ones. Count the
--   ones with @length (filter (== 1) bits)@.
--
--   `NFData` is there so that `force` and `deepseq` work on lists of outcomes,
--   as they did when outcomes were `Bool`s.
newtype Bit = Bit Bool
  deriving (Eq, Ord, Enum, Bounded, NFData)

instance Show Bit where
  show (Bit b) = if b then "1" else "0"

instance Num Bit where
  fromInteger 0 = Bit False
  fromInteger 1 = Bit True
  fromInteger n = error (show n ++ " is not a bit")
  Bit a + Bit b = Bit (a /= b)
  Bit a * Bit b = Bit (a && b)
  negate        = id
  abs           = id
  signum        = id
