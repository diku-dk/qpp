-- | Small, self-contained helpers: bit twiddling, integer logarithms,
--   permutations, and a couple of fold/queue utilities.
--
--   This module deliberately depends on nothing else in QPP (and in
--   particular not on "QPP.Syntax"), so that syntax and semantics modules can
--   import it freely. QOp-level analysis lives in "QPP.Syntax", the shared
--   backend interpreter in "QPP.Semantics", and the hmatrix block-splitting
--   helpers in "QPP.MPS" (their only consumer).
module QPP.Util where

import Data.Bits(FiniteBits,finiteBitSize,countLeadingZeros,shiftL)
import Data.List (sort)
import Math.NumberTheory.Logarithms (integerLog2')
import qualified Data.PQueue.Prio.Min as PriorityQ
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as MV
import Control.Monad.ST (runST)


-- Bit manipulation

toBits :: (Integral a) => a -> [Int]
toBits 0 = []
toBits k = (toBits (k `div` 2)) ++ [fromIntegral (k `mod` 2)]

toBits' :: Integral t => Int -> t -> [Int]
toBits' n k = let
    bits = toBits k
    m    = length bits
  in
    (replicate (n-m) 0) ++ bits

fromBits :: [Int] -> Int
fromBits bs =
        foldl (\acc b -> (acc `shiftL` 1) + b) 0 bs  -- MSB-first decode

-- | Infinite list of powers of two, constructed with O(1) time per element
powersOfTwo :: [Int]
powersOfTwo = iterate (*2) 1

dotlists :: Num a => [a] -> [a] -> a
dotlists xs ys = sum $ zipWith (*) xs ys

bitIndexMSB, bitIndexLSB :: [Int] -> Int
bitIndexMSB ks = dotlists (reverse ks) powersOfTwo
bitIndexLSB ks = dotlists ks powersOfTwo


-- Integer logarithms and powers

-- | ilog2 m = floor (log2 m) for m >= 0
ilog2 :: (FiniteBits a, Integral a) => a -> Int
ilog2 m = finiteBitSize m - countLeadingZeros m - 1

integerlog2 :: Integer -> Int
integerlog2 = integerLog2'

pow2 :: Int -> Int
pow2 n = 1 `shiftL` n

-- | ceil_log2 m = ⌈log2 m⌉ for m >= 1; the smallest k with pow2 k >= m.
--   I.e. the number of bits needed to index m distinct items: 1→0, 2→1, 3→2, 4→2, 5→3, ...
--   Returns 0 for m <= 1 (defensive; log2 0 is undefined).
ceil_log2 :: (FiniteBits a, Integral a) => a -> Int
ceil_log2 m
  | m <= 1    = 0
  | otherwise = ilog2 (m - 1) + 1


-- List helpers

evenOdd :: [a] -> ([a],[a])
evenOdd [] = ([],[])
evenOdd [x] = ([x],[])
evenOdd (x:y:xs) = let (es,os) = evenOdd xs in (x:es,y:os)


-- | Working with permutations
permApply :: [Int] -> [a] -> [a]
permApply ks xs = [ xs !! k | k <- ks ]

permSupport :: [Int] -> [Int]
permSupport ks = [ i | (i,j) <- zip [0..] ks, i /= j ]

permInvert :: [Int] -> [Int]
permInvert ks = map snd $  -- For each index in the output, find its position in the input
    sort [ (k, i) | (i, k) <- zip [0..] ks ]


-- | Minimal adjacent-swap indices (0-based) sending permutation p to identity.
--   Swap index i means swapping positions i and i+1.
--
--   For v = 0..n-1, move value v left until it sits at position v.
--   Each step performs exactly the inversions involving v, so total swap count is minimal,
--   one of very few applications where bubble sort is optimal.
--
-- Uses locally mutable vectors to avoid the  O(n^2) worst case unless it is actually needed:
-- Θ(n+inv(p)) instead of Θ(n^2+inv(p)) = O(n^2).
permutationSwaps :: V.Vector Int -> [Int]
permutationSwaps p0 = runST $ do
  let n    = V.length p0
      -- inverse permutation: pos[v] = index where value v currently sits.
      -- V.indexed p0 = [(i, p0!i)]; we need pairs (p0!i, i) to update pos[v] = i.
      pos0 = V.update (V.replicate n 0) (V.imap (\i v -> (v, i)) p0)
  p   <- V.thaw p0
  pos <- V.thaw pos0

  let -- adjacent swap s_i: swap p[i],p[i+1] and update pos accordingly
      swapAt i = do -- Updates mutable vectors p and pos to swap values at positions i and i+1
        a <- MV.read p i; b <- MV.read p (i+1)
        MV.write p i b;   MV.write p (i+1) a
        MV.write pos a (i+1); MV.write pos b i

      -- bubble value v left by s_{k-1} ... s_v, where k = pos[v]
      bubble v k swaps
        | k <= v    = pure swaps
        | otherwise = let i = k-1 in swapAt i >> bubble v i (i:swaps)

      -- Bubble sort remaining values v..n-1, assuming positions 0..v-1 are already sorted:
      -- p[j]=j for all j < v,
      -- while producing the list of swaps performed.
      sortFrom v swaps
        | v >= n    = pure (reverse swaps)
        | otherwise = MV.read pos v
                  >>= \k -> bubble v k swaps
                  >>= sortFrom (v+1)

  sortFrom 0 []


-- Folds and vector helpers

-- | Fold a binary operator over a vector via a balanced binary-tree shape
--   rather than the linear left/right chain of `foldl1`/`foldr1`. Length must
--   be a positive power of 2. Useful when the operator's result-arity grows
--   with fold depth (e.g. `foldBalanced v DirectSum` has depth log₂(length v)
--   instead of length v − 1). No equivalent exists in the Haskell base
--   libraries; `mconcat`/`foldMap` and friends use a linear fold by default.
foldBalanced :: V.Vector t -> (t -> t -> t) -> t
foldBalanced v f = go v
  where
    go xs
      | V.length xs == 1 = V.head xs
      | otherwise        = go (V.generate (V.length xs `div` 2) $ \i ->
                                  f (xs V.! (2*i)) (xs V.! (2*i+1)))


-- Flyttes fra MatrixPreparation.hs ... SKAL TESTES
padToPowerOf2 :: Int -> a -> V.Vector a -> V.Vector a
padToPowerOf2 numQbits paddingObj vec
    | len == targetLen = vec -- Already 2^n
    | otherwise        = vec V.++ V.replicate (targetLen - len) paddingObj
  where
    len       = V.length vec
    targetLen = 2 ^ numQbits



-- HELPER DATA STRUCTURES
pushTopK :: Int
         -> Double
         -> a -> PriorityQ.MinPQueue Double a
         -> PriorityQ.MinPQueue Double a
pushTopK !k !key !val !q
  | k <= 0        = PriorityQ.empty
  | PriorityQ.size q < k = PriorityQ.insert key val q
  | otherwise     =
      let (!kmin, _) = PriorityQ.findMin q
      in  if key <= kmin then q else PriorityQ.insert key val (PriorityQ.deleteMin q)
