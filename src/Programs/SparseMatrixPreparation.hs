{-# LANGUAGE BangPatterns #-}
module Programs.SparseMatrixPreparation where

import QPP
import QPP.Util

import Programs.MatrixPreparation ( calculateComplexGate )

import Data.Complex (polar)
import Control.Monad.ST (runST)
import qualified Data.Vector as V          -- Boxed Vector
import qualified Data.Vector.Mutable as VM  -- Mutable Boxed Vector

----------------------------------------------------------------------------------------
-- Unitary encoding U_A
-- This is Theorem 4.1 Camps et al. arXiv:2203.10236v4
----------------------------------------------------------------------------------------

sparseUAEncoding :: Nat -> QOp -> QOp -> QOp
sparseUAEncoding sQubits valOracle posOracle =
    let noQubits = op_qubits valOracle
        n = noQubits - (sQubits + 1)
        dOp = I ⊗ ( diffusionOp sQubits ) ⊗ (Id n)
    in 
       dOp ∘ (I ⊗ posOracle) ∘ valOracle ∘ dOp
    
diffusionOp :: Int -> QOp
diffusionOp n = foldr1 Tensor (replicate n H) 

----------------------------------------------------------------------------------------
-- positionOracle :: Vector QOp -> QOp
-- This corresponds to the permutation operator O_c = ∑|j,c_jl> <l,j| 
-- (Theorem 4.1 Camps et al. arXiv:2203.10236v4)
-- Unitarity demands that the sparseColOracle c_jl is a permutation for each l
-- The implementation receives a vector of QOp permutations.
-- One permutation for each sparsity index l in [1,...,s = 2^n]
----------------------------------------------------------------------------------------
positionOracle :: V.Vector QOp -> QOp
positionOracle permutations = foldBalanced permutations DirectSum

----------------------------------------------------------------------------------------
-- valueOracle :: GenericSparseMatrix -> QOp
----------------------------------------------------------------------------------------

valueOracle :: V.Vector QOp -> QOp
valueOracle rotBlocks =
    let 
        -- Form the direct sum of rotation matrices
        dsQop = foldBalanced rotBlocks DirectSum

        -- Construct a basis change permutation 
        n = op_qubits dsQop
        pQOp = Permute ([1 .. (n-1)] ++ [0])
    in 
        (Adjoint pQOp) ∘ dsQop ∘ pQOp

----------------------------------------------------------------------------------------
-- Reading data from the [[ComplexT]] containing the sparse data
----------------------------------------------------------------------------------------

-- | Stepping through data rows first. 
-- Multiple repetitions for each row if repetitions set to > 1
buildRotVec :: Nat -> [[ComplexT]] -> V.Vector QOp
buildRotVec repetitions mat = runST $ do
    -- 1. Calculate total elements, factoring in the repetition count
    let !totalElements = sum (map length mat) * repetitions
    targetVec <- VM.new totalElements
    
    -- 2. Define the main recursive loop over rows
    let populate [] !idx = return idx
        populate (row:rows) !idx = do
            -- Inner helper to walk across columns of a single row
            let walkRow [] !i = return i
                walkRow (c:cs) !i = do
                    let !obj = buildRotQOp c
                    VM.write targetVec i obj
                    walkRow cs (i + 1)
            
            -- Loop helper to repeat the same row sequence multiple times
            let repeatRow 0 !i = return i
                repeatRow !remReps !i = do
                    !nextI <- walkRow row i
                    repeatRow (remReps - 1) nextI  -- Fixed typo here
                
            !nextIdx <- repeatRow repetitions idx
            populate rows nextIdx

    _ <- populate mat 0
    V.freeze targetVec

buildRotQOp :: ComplexT -> QOp
buildRotQOp val =
    let  
        (r1, phi1) = polar val
        r2 = sqrt(1 - r1*r1)
        phi2 = 0
    in
        calculateComplexGate r1 phi1 r2 phi2
