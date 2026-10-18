-- | Computational-basis measurement on cyclotomic vectors. Collapse retains
-- the unnormalized projection, so no square root or field extension is needed.
module QPP.Exact.Internal.Measurement
  ( measure1, measure1WithLimits, probability, budget, correctBits
  ) where

import Control.Monad (unless, when)
import Data.Bits (testBit, setBit, xor)
import qualified Data.List as List
import qualified Data.Vector as V
import GHC.Stack (HasCallStack, withFrozenCallStack)
import QPP.Syntax (Nat, Outcomes, RNG)
import QPP.Exact.Types
import QPP.Exact.Internal.State (StateT(..))
import qualified QPP.Exact.Internal.State as St
import QPP.Exact.Internal.Buffer (Buffer(..))
import qualified QPP.Exact.Internal.Buffer as B
import QPP.Exact.Internal.Scalar (Field, Coeff, Cyclo(..))
import qualified QPP.Exact.Internal.Scalar as S
import qualified QPP.Exact.Internal.Order as Order

-- Input, chosen projection, branch weights, and comparison intervals. Scalar
-- polynomial workspace is additional, as in the unitary evaluators.
budget :: Limits -> Field -> Int -> Result ()
budget limits field n = do
  dim <- S.dimension limits n
  S.workspace limits field ((2 * toInteger dim + 8) * toInteger (S.degree field) + 64)

weights :: Limits -> StateT -> Nat -> Result (Coeff, Coeff)
weights limits state@(StateT n input@(Buffer field xs)) k = do
  St.validate limits state
  unless (k >= 0 && k < n) (Left (InvalidInput "measurement wire out of bounds"))
  budget limits field n
  (w0,w1) <- V.ifoldM' accumulate (S.zero field, S.zero field) xs
  when (S.isZeroC w0 && S.isZeroC w1)
    (Left (InvalidInput "cannot measure the zero vector"))
  pure (w0,w1)
  where
    accumulate (w0,w1) j a = do
      ac <- S.conjugateC limits (B.field input) a
      w <- S.mulC limits (B.field input) ac a
      if testBit j (n-1-k)
        then (w0,) <$> S.addC limits w1 w
        else (,w1) <$> S.addC limits w0 w

-- | Conditional probability for the requested bit, even on an unnormalized
-- vector. This division is for observation; sampling avoids the inverse.
probability :: Limits -> StateT -> Nat -> Bool -> Result Cyclo
probability limits state@(StateT _ (Buffer field _)) k outcome = do
  (w0,w1) <- weights limits state k
  total <- S.addC limits w0 w1
  inverse <- S.inverseC limits field total
  Cyclo field <$> S.mulC limits field (if outcome then w1 else w0) inverse

measure1 :: HasCallStack => (StateT, Outcomes, RNG) -> Nat -> (StateT, Outcomes, RNG)
measure1 input k = withFrozenCallStack $
  orFail "measure1" (measure1WithLimits defaultLimits input k)

-- | Consume one RNG value, prepend the outcome, and retain the projection.
-- A Double RNG value is compared as its exact dyadic rational. Equality at
-- the threshold selects outcome one, matching the numerical backend rule.
measure1WithLimits :: Limits -> (StateT, Outcomes, RNG) -> Nat
                   -> Result (StateT, Outcomes, RNG)
measure1WithLimits limits (state@(StateT n (Buffer field xs)), outcomes, rng) k = do
  (w0,w1) <- weights limits state k
  (r,remaining) <- case rng of
    [] -> Left (InvalidInput "measurement RNG is empty")
    r:rs -> do
      unless (not (isNaN r || isInfinite r) && r >= 0 && r < 1)
        (Left (InvalidInput "measurement RNG value must be finite and in [0,1)"))
      pure (r,rs)
  total <- S.addC limits w0 w1
  q <- S.checkRational limits (toRational r)
  threshold <- S.scaleC limits q total
  difference <- S.subC limits w0 threshold
  ordering <- Order.signReal limits field difference
  let outcome = ordering /= GT
      projected = V.imap (\j a -> if testBit j (n-1-k) == outcome then a else S.zero field) xs
  pure (StateT n (Buffer field projected), outcome:outcomes, remaining)

-- | Apply the tensor product of X gates selected by measurement/target XORs.
-- This is a basis-index permutation and introduces no scalar constants.
-- Wires are distinct and validated by checkStep before this internal call.
correctBits :: Limits -> [(Nat, Bool)] -> StateT -> Result StateT
correctBits limits corrections state@(StateT n (Buffer field xs)) = do
  let mask = List.foldl' (\m (k,flipBit) -> if flipBit then setBit m (n-1-k) else m) 0 corrections
  if mask == 0 then pure state else do
    S.reserve limits field S.LinearWork (2 * toInteger (V.length xs) * toInteger (S.degree field))
    pure (StateT n (Buffer field (V.generate (V.length xs) (\j -> xs V.! (j `xor` mask)))))
