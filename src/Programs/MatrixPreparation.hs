module Programs.MatrixPreparation where
import QPP
import QPP.Util
import Data.Complex
import Data.Vector (Vector)
import qualified Data.Vector as V
import Data.Maybe

import Numeric.LinearAlgebra (conj, toRows, toList)
import QPP.Semantics.Matrix (CMat)

-- | The row-prep internals below consume Data.Vector, so convert from hmatrix's
--   storable rows once here. Remove when the row-prep is refactored too.
getRows :: CMat -> V.Vector (V.Vector ComplexT)
getRows = V.fromList . map (V.fromList . toList) . toRows

----------------------Matrixpreparation--------------------------

-- | Norm-encoding half of the U†V construction: prepares the row-norm vector
--   of M, tensored with an n-qubit identity for the system register.
unitaryV :: CMat -> QOp
unitaryV mat =
    let rowsVec = getRows mat
        rowLen v = sqrt . V.sum $ V.map (\x -> x * conjugate x) v
        allLens = V.map rowLen rowsVec
        retQOp = buildRowQOp allLens
        numQbits = ceil_log2 (V.length rowsVec)
    in
        retQOp ⊗ (Id numQbits)

-- | Row-encoding half: for each row of conj(M), build a row-prep unitary;
--   assemble them into a balanced DirectSum tree, then SWAP to put the
--   ancilla register first.
unitaryU :: CMat -> QOp
unitaryU mat =
    let
        -- Get the row encodings
        rowsVec = getRows (conj mat)
        uBlocks = (V.map buildRowQOp rowsVec)
        -- Pad with I_n to 2^n length
        numQbits = ceil_log2 (V.length uBlocks)
        uBlocksPadded = padToPowerOf2 numQbits (Id numQbits) uBlocks
        -- Build the direct sums
        dsQOp = foldBalanced uBlocksPadded DirectSum

        -- Building the SWAP
        swap = if numQbits == 1 then Permute [1,0] else Permute ([numQbits .. 2*(numQbits) - 1 ] ++ [0 .. (numQbits - 1)])
    in
        dsQOp ∘ swap

{- Moved to QPP.Util ... SKAL TESTES
padToPowerOf2 :: Int -> QOp -> V.Vector QOp -> V.Vector QOp
padToPowerOf2 numQbits paddingObj vec
    | len == (2^numQbits) = vec -- Already 2^n
    | otherwise   = vec V.++ V.replicate ((2^numQbits) - len) paddingObj
  where
    len  = V.length vec
-}


-- | Block-encode a complex matrix M as a 2n-qubit unitary U such that
--   the upper-left 2^n × 2^n block of U equals M / ‖M‖_F (Frobenius norm).
matrixPrep :: CMat -> QOp
matrixPrep mat = (Adjoint (unitaryU mat)) ∘ (unitaryV mat)

------------------------ Row preparation -------------------------

-- | Build a QOp that prepares the amplitude state v/‖v‖ from |0..0⟩.
--   The input vector is zero-padded to the next power of 2 implicitly via
--   `splitList`'s singleton chunks (each becomes an R Y 0 = I leaf).
buildRowQOp :: Vector ComplexT -> QOp
buildRowQOp vs = createRotations 0 (splitList vs)

createRotations :: Int -> [V.Vector ComplexT] -> QOp
createRotations level vs
    | length vs == 1 = (createQOp level (head vs) ⊗ (Id level))
    | otherwise =
        let -- One 1-qubit leaf gate per chunk. splitList can produce a non-power-of-2
            -- number of chunks (e.g. length-5 input → 3 chunks); pad the leaf array with
            -- Id 1 placeholders so foldBalanced can fold a balanced binary DirectSum
            -- tree. The placeholder branches correspond to indices with zero amplitude.
            leafOps    = V.fromList (map (createQOp level) vs)
            nQubits    = ceil_log2 (V.length leafOps)
            leafOpsPad = padToPowerOf2 nQubits I leafOps
            accQOp     = foldBalanced leafOpsPad DirectSum
            -- Inner-node norms for the next level of the rotation tree.
            pairLenLst = V.fromList (map pairNorm vs)
        in (accQOp ⊗ (Id level)) <> createRotations (level + 1) (splitList pairLenLst)
  where
    pairNorm pair =
        let val0 = pair V.! 0
            val1 = fromMaybe 0 (pair V.!? 1)
        in sqrt (val0 * conjugate val0 + val1 * conjugate val1)

createQOp :: Int -> V.Vector ComplexT -> QOp
createQOp level pairVector
    | level < 0  = error "Ups ... level er negativt"
    | level == 0 = if abs r1 < 1e-9 && abs r2 < 1e-9 then I else calculateComplexGate r1 phi1 r2 phi2
    | otherwise  = if abs r1 < 1e-9 && abs r2 < 1e-9 then I else calculateRealGate r1 r2
  where
    val0       = pairVector V.! 0
    val1       = fromMaybe 0 (pairVector V.!? 1)
    (r1, phi1) = polar val0
    (r2, phi2) = polar val1


calculateComplexGate :: Double -> Double -> Double -> Double -> QOp
calculateComplexGate r1 phi1 r2 phi2
    | abs r1 < 1e-9 && abs r2 < 1e-9 = I
    | otherwise =
        let phi    = toRational $   (phi2 - phi1)/pi
            lambda = toRational $  -(phi1 + phi2)/pi
            theta  = toRational $ if r1 == 0 then 1 else (2/pi) * acos (r1 / sqrt (r1**2 + r2**2))
        in
            -- Convention: R Y θ |0⟩ = cos(πθ/2)|0⟩ + sin(πθ/2)|1⟩. Positive θ encodes a
            -- positive amplitude on |1⟩. (The angles here used to be negated to compensate
            -- for an old sign bug in R that has since been fixed.)
            (R Z phi) ∘ (R Y theta) ∘ (R Z lambda)

calculateRealGate :: Double -> Double -> QOp
calculateRealGate r1 r2
    | abs r1 < 1e-9 && abs r2 < 1e-9 = I
    | otherwise =
        let theta = toRational $ if abs r1 < 1e-9 then 1 else (2/pi) * acos (r1 / sqrt (r1**2 + r2**2))
        in
            R Y theta

----------------- Helper funcs ----------

-- | Recursively split a vector into power-of-2 chunks of length ≤ 2.
splitList :: V.Vector ComplexT -> [V.Vector ComplexT]
splitList vs = go [vs]
  where
    go curr
        | all (\v -> V.length v <= 2) curr = curr
        | otherwise = go (concatMap splitStep curr)

splitStep :: V.Vector ComplexT -> [V.Vector ComplexT]
splitStep v
    | len <= 2 = [v]
    | otherwise =
        let n = floor (logBase 2 (fromIntegral len) :: Double) :: Integer
            -- If len is exactly a power of 2 (e.g., 4), split it in half (2, 2)
            -- Otherwise, split at the largest power of 2 (e.g., 6 becomes 4, 2)
            splitPoint = if 2 ^ n == len then 2 ^ (n - 1) else 2 ^ n
            left  = V.slice 0 splitPoint v
            right = V.slice splitPoint (len - splitPoint) v
        in [left, right]
  where len = V.length v
