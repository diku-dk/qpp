module Programs.SVT where
import QPP
import Programs.Qubitization

import Data.Poly (VPoly)

import Data.Sequence (Seq(..), (><))
import qualified Data.Sequence as Seq

-- Calculates the SVT vector of phases
-- Input:
-- 1) p, q complex polynomials
-- 2) p(x)^2 + (1 - x^2)q(x)^2 = 1 for x in [-1,1] :
-- 3) (p even, q odd) or (p odd, q even)
-- 4) deg(p) = deg(q) + 1
-- Output:
-- Real k dimensional sequence where k = deg(p)
svtVectorPhi :: VPoly ComplexT -> VPoly ComplexT -> Seq Double
svtVectorPhi p q =
    let
        -- Generate a sequence
        qspVector = Seq.fromList (qspVectorPhi p q)

        -- 2. Modify to get the SVT vector
        (phi0, svtVector, phiLast) = case qspVector of
                    (p0 :<| (mid :|> pL))  -> (p0, mid, pL)
                    _                      -> error "Vector too short for SVT modification"

        d = (fromIntegral (Seq.length qspVector - 2) * (pi / 2) :: Double)
        modQspVector = fmap (subtract (pi / 2)) svtVector
    in
       Seq.fromList ([phi0 + phiLast + d ]) >< modQspVector

-- This function creates the Alternating Phase Modulation Operator 
-- It is the result of the SVT process.
-- The operator is unitary and has the transformed singular values
-- UA = ∑P(s_i)|v_i><v_i| when A = ∑s_i|v_i><v_i| is the input encoded matrix
-- Input:
-- matrixQubits n, col(A) = 2^n. So it is the log2 size of the original matrix
-- svtVec is a real vector of phases
-- blockEnc is the block encoding of A, that needs to have ||A|| <= 1 ... I think ...
-- Output:
-- 1) A unitary operator UA encoded in the same number of qubits as blockEnc
-- 2) The singular values should correspond to UA = ∑P(s_i)|v_i><v_i|
altPhaseMod :: Int -> Seq Double -> QOp -> QOp
altPhaseMod matrixQubits svtVec blockEnc =
    let
        -- Prøver med denne i stedet for encQubits
        encQubits = (op_qubits blockEnc) - matrixQubits

        -- Create the phase matrices ...
        phaseOperators = fmap (phaseOpBuilder matrixQubits encQubits) svtVec 
    in
        -- Create the phased alternating sequence
        phaseAltSeqBuilder (isEvenLength svtVec) phaseOperators blockEnc

phaseOpBuilder :: Int -> Int -> Double -> QOp
phaseOpBuilder matrixQubits encQubits phi = 
    let 
        phiPhaseOp = ((Id matrixQubits) ∘ Phase (toRational ( phi / pi))) ⊕ ((Id matrixQubits) ∘ Phase (toRational ( - phi / pi)))
    in
        phasePadder (matrixQubits + 1) (encQubits - 1) phi phiPhaseOp

phasePadder :: Int -> Int -> Double -> QOp -> QOp
phasePadder opQubits encQubits phi phiPhaseOp
    | encQubits == 0 = phiPhaseOp
    | otherwise =
        let 
            newphiPhaseOp = phiPhaseOp ⊕ (Id opQubits <> Phase (toRational ( - phi / pi)))
        in
            phasePadder (opQubits + 1) (encQubits - 1) phi newphiPhaseOp

phaseAltSeqBuilder :: Bool -> Seq QOp -> QOp -> QOp
phaseAltSeqBuilder isEven phaseOps blockEnc  
    | not isEven =  
        let
            (phaseOp, phaseOpsRest) = case phaseOps of
                (p0 :<| mid ) -> (p0, mid)
                _             -> error "Vector too short for SVT modification"

            part = phaseOp <> blockEnc
        in
            part <> phaseAltSeqBuilder True phaseOpsRest blockEnc

    | otherwise =
        let
            (phaseOp1, phaseOpsRest1) = case phaseOps of
                (p0 :<| rest ) -> (p0, rest)
                _              -> error "Vector too short for SVT modification"
            part1 = phaseOp1 <> (Adjoint blockEnc)

            (phaseOp2, phaseOpsRest2) = case phaseOpsRest1 of
                (p0 :<| rest ) -> (p0, rest)
                _              -> error "Vector too short for SVT modification"
            part2 = phaseOp2 <> blockEnc

            part = part1 <> part2
        
        in case phaseOpsRest2 of 
            Empty -> part 
            _  -> part <> phaseAltSeqBuilder True phaseOpsRest2 blockEnc
            

isEvenLength :: Seq a -> Bool
isEvenLength s = even (Seq.length s)