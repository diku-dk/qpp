-- | Exact state shape wrapped around a field-owned coefficient buffer.
module QPP.Exact.Internal.State
  ( StateT(..), qubits, fieldIndex, amplitudes, amplitude, validate
  , fromAmplitudes, ket, scale, add, sub, tensor, inner, normSquared, equal
  ) where

import Control.Monad (unless)
import qualified Data.List as List
import qualified Data.Vector as V
import QPP.Exact.Types
import QPP.Syntax (HasQubits(..))
import QPP.Exact.Internal.Scalar (Cyclo(..))
import qualified QPP.Exact.Internal.Scalar as S
import QPP.Exact.Internal.Buffer (Buffer(..))
import qualified QPP.Exact.Internal.Buffer as B

data StateT = StateT !Int !Buffer deriving (Show)

instance HasQubits StateT where n_qubits = qubits

qubits :: StateT -> Int
qubits (StateT n _) = n

fieldIndex :: StateT -> Integer
fieldIndex (StateT _ b) = S.index (B.field b)

amplitudes :: StateT -> [Cyclo]
amplitudes (StateT _ b) = B.scalars b

amplitude :: StateT -> Int -> Result Cyclo
amplitude (StateT _ (Buffer f xs)) j
  | j >= 0 && j < V.length xs = pure (Cyclo f (xs V.! j))
  | otherwise = Left (InvalidInput "amplitude index out of bounds")

validate :: Limits -> StateT -> Result ()
validate l (StateT n b) = do
  d <- S.dimension l n
  unless (B.count b == toInteger d) (Left (InvalidInput "internal state dimension"))
  B.validate l b

fromAmplitudes :: Limits -> [Cyclo] -> Result StateT
fromAmplitudes l xs = do
  n <- S.widthOfLength (length xs)
  _ <- S.dimension l n
  StateT n <$> B.fromScalars l xs

ket :: Limits -> [Int] -> Result StateT
ket l bs = do
  unless (all (\b -> b == 0 || b == 1) bs)
    (Left (InvalidInput "ket expects bits equal to 0 or 1"))
  d <- S.dimension l (length bs)
  f <- S.makeField l 1
  S.storage l f (toInteger d)
  let j = List.foldl' (\acc b -> 2*acc + b) 0 bs
  pure (StateT (length bs) (Buffer f (V.generate d (\k -> if j == k then S.one f else S.zero f))))

sameWidth :: StateT -> StateT -> Result ()
sameWidth a b = unless (qubits a == qubits b) (Left (WidthMismatch (qubits a) (qubits b)))

scale :: Limits -> Cyclo -> StateT -> Result StateT
scale l c (StateT n b) = StateT n <$> B.scale l c b

add, sub :: Limits -> StateT -> StateT -> Result StateT
add l = combine (B.add l)
sub l = combine (B.sub l)

combine :: (Buffer -> Buffer -> Result Buffer) -> StateT -> StateT -> Result StateT
combine action a@(StateT n xs) b@(StateT _ ys) = do
  sameWidth a b
  StateT n <$> action xs ys

inner :: Limits -> StateT -> StateT -> Result Cyclo
inner l a@(StateT _ xs) b@(StateT _ ys) = sameWidth a b >> B.inner l xs ys

normSquared :: Limits -> StateT -> Result Cyclo
normSquared l v = inner l v v

equal :: Limits -> StateT -> StateT -> Result Bool
equal l (StateT _ a) (StateT _ b) = B.equal l a b

tensor :: Limits -> StateT -> StateT -> Result StateT
tensor l (StateT n a) (StateT m b) = do
  _ <- S.dimension l (n+m)
  StateT (n+m) <$> B.tensor l a b
