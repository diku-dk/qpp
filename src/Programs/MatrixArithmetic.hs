module Programs.MatrixArithmetic where
import QPP
import QPP.Util
import Data.Vector (Vector)
import qualified Data.Vector as V

import Data.List (sortOn)
import Data.Ord (Down(..))
import qualified Data.Set as S


-- Kun Testet på reelle matricer
productEncoding :: Int -> QOp -> QOp -> QOp
productEncoding numQbits qOp1 qOp2 =
    let 
        swap = if numQbits == 1 then Permute [1,0] else Permute ([numQbits .. (2*numQbits - 1)] ++ [0 .. (numQbits - 1)])
    in
        ((Id numQbits ⊗ (qOp1 <> swap)) ∘ (qOp2 ⊗ Id numQbits) ∘ (Id numQbits ⊗ swap))

-- Kun Testet på reelle matricer
sumEncoding :: QOp -> QOp -> QOp -> QOp
sumEncoding normQOp qOp1 qOp2 =
    let 
        wMatrix = qOp1 ⊕ qOp2
        n = op_qubits wMatrix - op_qubits normQOp
    in
        (Adjoint normQOp ⊗ Id n) ∘ wMatrix ∘ (normQOp ⊗ Id n)


-- Dyadic non rotational LCU (linear combination of unitaries)
-- A variation of Gilyen et al Lemma 52
-- Uses NO rotations at the cost of instead extra qubits : 
-- A number of precision qubits and 1 ancilla
dyadicLCU :: Int -> V.Vector (Double, QOp) -> QOp
dyadicLCU precision terms =
    let 
        -- Split the single vector of pairs into two perfectly matched vectors
        (coeffs, qOps) = V.unzip terms

        -- Relevant Hilbert space qubits are found
        termQbits = ceil_log2 (V.length terms)
        sysQubits = if V.null qOps then 0 else op_qubits (V.head qOps)

        -- Sample coeffs with precision set to 2^precision 
        weights = cumulativeWeights precision coeffs 
        listOfSigns = signVector coeffs
        
        -- Padding so that vectors are balanced
        weightsPad = padToPowerOf2 termQbits (2^precision) weights
        listOfSignsPad = padToPowerOf2 termQbits (Phase 0) listOfSigns
        qOpsPad    = padToPowerOf2 termQbits I qOps

        -- Create a list of comparator operators based on precision and weights
        uCmps = comperatorQOps precision sysQubits weightsPad
 
        -- Create a list of controlled operators based on signs and qOps
        -- Controlled on a single ancilla
        -- Target is system qubits
        uCtrls = controlQOps precision (safeZip qOpsPad listOfSignsPad)

        -- The select oracle is finalized based on comparators and controls
        selectQOp =  createSelectQOp (safeZip uCmps uCtrls)

        -- The PREPARE unitary oracle is created
        -- It is just a diffusion operator on precision
        prepareQOp = diffusionOp precision ⊗ Id (1 + sysQubits) 
    in
        Adjoint prepareQOp ∘ selectQOp ∘ prepareQOp

-- DENNE ER DOBBELT DEFINERET findes også i SparseMatrixPreparation
diffusionOp :: Int -> QOp
diffusionOp n = foldr1 Tensor (replicate n H) 

-- Find the dyadic approximations
-- May have to optimize this later
-- But should be well suited for exponentially decaying coefficients
discretize :: Int -> Vector Double -> Vector Int
discretize m cs =
    let 
        n  = 2 ^ m
        nD = fromIntegral n

        -- ignore sign completely
        xs = V.map (\c -> abs c * nD) cs

        floors = V.map floor xs
        fracs  = V.zipWith (\x f -> x - fromIntegral f) xs floors

        r = n - V.sum floors

        -- indices of largest fractional parts
        idx =
            map fst $
            take r $
            sortOn (Down . snd) $
            zip [0..] (V.toList fracs)

        idxSet = S.fromList idx
    in 
        V.imap (\i f -> if S.member i idxSet then f + 1 else f) floors

-- Build
cumulativeWeights :: Int -> V.Vector Double -> V.Vector Int
cumulativeWeights precision coeffs =
    let
        discretizedCoeffs = discretize precision coeffs
    in
        -- Accumulate
        V.scanl1 (+) discretizedCoeffs

-- Extract the signs from the list of coefficients
signVector :: V.Vector Double -> V.Vector QOp
signVector = V.map (\x -> if x < 0 then Phase 1 else Phase 0)
    
-- Create a list of comparator operators based on the weights and signs
-- The vector weights is a look-up table
comperatorQOps :: Int -> Int -> V.Vector Int -> V.Vector QOp
comperatorQOps m sysQbits weights =
  let n = 2 ^ m
  in V.map (\k ->
        let prev = if k == 0 then 0 else weights V.! (k - 1)
            curr = weights V.! k

            cmpList =
                 V.replicate prev I
              V.++ V.replicate (curr - prev) X
              V.++ V.replicate (n - curr) I

        in foldBalanced cmpList DirectSum ⊗ Id sysQbits
     )
     (V.enumFromN 0 (V.length weights))

-- Create a list of controlled operators
-- Controlled on a single ancilla
-- Target is system qubits
controlQOps :: Int -> V.Vector (QOp,QOp) -> V.Vector QOp
controlQOps precision = V.map (\(qOp, signOp) ->  Id precision ⊗ C (signOp ⊗ qOp) ) 


-- The select oracle is finalized 
createSelectQOp :: V.Vector (QOp,QOp) -> QOp
createSelectQOp vPairs = V.foldl1 (∘) $ V.map (\(ck, uk) -> Adjoint ck ∘ uk ∘ ck) vPairs

----------- Helpers ----------

-- | Zip only if lengths match exactly
safeZip :: V.Vector a -> V.Vector b -> V.Vector (a, b)
safeZip v1 v2
  | V.length v1 == V.length v2 = V.zip v1 v2
  | otherwise                  = error "Cannot create operation: qOps and signs vector lengths mismatch!"