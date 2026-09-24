{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE FlexibleContexts #-}
module Main where 

import QPP
import QPP.Semantics.Statevector

import Programs.Qubitization
import Programs.SVT
import Programs.MatrixPreparation

import Numeric.LinearAlgebra (cmap,eigenvalues, singularValues, fromLists, norm_Frob)
import Data.Complex
import Data.Sequence (Seq(..))
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import Data.List (sortOn)
import Data.Ord (Down(..))
import Data.Poly (VPoly, eval, toPoly, unPoly)

main :: IO ()
main = do
    print (show (svtVectorPhiTest()))
    print ("evenPolySVTtest max error: " ++show (evenPolySVTtest()))
    print ("oddPolySVTtest max error: " ++show (oddPolySVTtest()))
    print ("oddHermitianSVTtest max error: " ++ show (oddHermitianSVTtest()))
    print ("evenHermitianSVTtest max error: " ++ show (evenHermitianSVTtest()))

svtVectorPhiTest :: () -> Seq Double
svtVectorPhiTest () =
    let 
        p = [0,0,0.5,0,0.5] -- p(x) = 1/2x^2+1/2x^4
        (pC,qC) = qspPolys p
    in
        svtVectorPhi pC qC

evenPolySVTtest :: () -> Double
evenPolySVTtest () =
    let 
        -- Et lige polynomium p(x) med p(x)^2 mindre end 1 for alle x i [-1,1] defineres
        p = [0,0,0.5,0,0.5] -- 0.5x^2 + 0.5x^4 
        (pC,qC) = qspPolys p
        svtVec = svtVectorPhi pC qC

        -- Define a matrix and set matrixQubits
        -- This is a non-diagonal matrix with singular values 0.8 and 0.6 and Frobeniusnorm 1
        matData =  [[((4*sqrt(6) - 3*sqrt(2))/20) :+ 0.0, ((-4*sqrt(6)-3*sqrt(2))/20) :+ 0.0], 
                    [((4*sqrt(2)+3*sqrt(6))/20) :+ 0.0, ((-4*sqrt(2)+3*sqrt(6))/20) :+ 0.0]]  
        matIn = fromLists matData
        matrixQubits = 1

        -- Lav block encoding
        blockEnc = matrixPrep matIn

        -- SVT udføres
        svtQOp = altPhaseMod matrixQubits svtVec blockEnc
        svtQOt = evalOp $ svtQOp

        -- SVT udføres med (-1 * svtVec) ... med det formål at ...
        svtVecMinus = fmap (* (-1)) svtVec
        svtConjQOp = altPhaseMod matrixQubits svtVecMinus blockEnc
        svtConjQOt = evalOp $ svtConjQOp

        -- ... tage gennemsnit af de tilstande man får ved anvendelse af 
        -- den alternerende fasemodulerende sekvens hørende til blockEnc med svtVec hhv (-svtVec) : 
        tempState1 = (0.5 :+ 0) .* ((apply svtQOt (ket [0,0])) .+ (apply svtConjQOt (ket [0,0])))
        tempState2 = (0.5 :+ 0) .* ((apply svtQOt (ket [0,1])) .+ (apply svtConjQOt (ket [0,1])))
        
        -- De singulære normeringskonstanter skal også kvadreres:
        normSV = 0.5*(norm_Frob matIn)^(2 :: Int) + 0.5*(norm_Frob matIn)^(4 :: Int)

        finalState1 = (normSV :+ 0) .* tempState1
        finalState2 = (normSV :+ 0) .* tempState2
     
        -- Der projiceres
        m11 = inner  (ket [0,0]) finalState1
        m12 = inner  (ket [0,0]) finalState2
        m21 = inner  (ket [0,1]) finalState1
        m22 = inner  (ket [0,1]) finalState2

        -- Den resulterende matrix ... har den mon de kvadrerede singulære værdier?
        matOutData = [[m11,m12], [m21,m22]]
        matOut = fromLists matOutData
        
        -- Calculating the transformed input SV and the output SV
        matInSV =  singularValues matIn
        matInSVT =  VS.map (polyEval p) matInSV
        matInSVTabs = sortDescendingReal $ VS.map (abs) matInSVT -- Singular values are positive!
        matOutSV = sortDescendingReal $ singularValues matOut
    in 
        -- What is the maximum deviation in the SVT?
        VS.maximum $ matInSVTabs - matOutSV

oddPolySVTtest :: () -> Double
oddPolySVTtest () =
    let 
        -- Et ulige polynomium p(x) med p(x)^2 mindre end 1 for alle x i [-1,1] defineres
        p = [0,(-1.5),0,0.5,0,2] -- 2x^5 + 0.5x^3 - 1.5x 
        (pC,qC) = qspPolys p
        svtVec = svtVectorPhi pC qC

        -- Definer en matrix og sæt matrixQubits
        -- Dette er en ikke-diagonal matrix med singulære værdier 0.8 og 0.6 and Frobeniusnorm 1
        matData =   [[((4*sqrt(6) - 3*sqrt(2))/20) :+ 0.0, ((-4*sqrt(6)-3*sqrt(2))/20) :+ 0.0], 
                        [((4*sqrt(2)+3*sqrt(6))/20) :+ 0.0, ((-4*sqrt(2)+3*sqrt(6))/20) :+ 0.0]]  

        matIn = fromLists matData
        matrixQubits = 1

        -- Lav block encoding
        blockEnc = matrixPrep matIn

        -- SVT udføres
        svtQOp = altPhaseMod matrixQubits svtVec blockEnc
        svtQOt = evalOp $ svtQOp

        -- SVT udføres med (-1 * svtVec) ... med det formål at ...
        svtVecMinus = fmap (* (-1)) svtVec
        svtConjQOp = altPhaseMod matrixQubits svtVecMinus blockEnc
        svtConjQOt = evalOp $ svtConjQOp

        -- ... tage gennemsnit af de tilstande man får ved anvendelse af 
        -- den alternerende fasemodulerende sekvens hørende til blockEnc med svtVec hhv (-svtVec) : 
        finalState1 = (0.5 :+ 0) .* ((apply svtQOt (ket [0,0])) .+ (apply svtConjQOt (ket [0,0])))
        finalState2 = (0.5 :+ 0) .* ((apply svtQOt (ket [0,1])) .+ (apply svtConjQOt (ket [0,1]))) 
     
        -- Der projiceres
        m11 = inner  (ket [0,0]) finalState1
        m12 = inner  (ket [0,0]) finalState2
        m21 = inner  (ket [0,1]) finalState1
        m22 = inner  (ket [0,1]) finalState2

        -- Den resulterende matrix
        matOutData = [[m11,m12], [m21,m22]]
        matOut = fromLists matOutData
        
        -- Calculating the transformed input SV and the output SV
        matInSV =  singularValues matIn
        matInSVT =  VS.map (polyEval p) matInSV
        matInSVTabs = sortDescendingReal $ VS.map (abs) matInSVT -- Singular values are positive!
        matOutSV = sortDescendingReal $ singularValues matOut
    in 
        -- What is the maximum deviation in the SVT?
        VS.maximum $ matInSVTabs - matOutSV

oddHermitianSVTtest :: () -> Double
oddHermitianSVTtest () =
    let 
        -- Et ulige polynomium p(x) med p(x)^2 mindre end 1 for alle x i [-1,1] defineres
        p = [0,(-1.5),0,0.5,0,2] -- 2x^5 + 0.5x^3 - 1.5x 
        (pC,qC) = qspPolys p
        svtVec = svtVectorPhi pC qC

        -- Definer en matrix og sæt matrixQubits
        matData =   [[0.0 :+ 0.0, sqrt(2)/4 :+ (-sqrt(2)/4)], 
                        [sqrt(2)/4 :+ sqrt(2)/4, (-1/sqrt(2)) :+ 0.0]]  

        matIn = fromLists matData
        matrixQubits = 1

        -- Lav block encoding
        blockEnc = matrixPrep matIn

        -- SVT udføres
        svtQOp = altPhaseMod matrixQubits svtVec blockEnc
        svtQOt = evalOp $ svtQOp

        -- SVT udføres med (-1 * svtVec) ... med det formål at ...
        svtVecMinus = fmap (* (-1)) svtVec
        svtConjQOp = altPhaseMod matrixQubits svtVecMinus blockEnc
        svtConjQOt = evalOp $ svtConjQOp

        -- ... tage gennemsnit af de tilstande man får ved anvendelse af 
        -- den alternerende fasemodulerende sekvens hørende til blockEnc med svtVec hhv (-svtVec) : 
        finalState1 = (0.5 :+ 0) .* ((apply svtQOt (ket [0,0])) .+ (apply svtConjQOt (ket [0,0])))
        finalState2 = (0.5 :+ 0) .* ((apply svtQOt (ket [0,1])) .+ (apply svtConjQOt (ket [0,1])))
     
        -- Der projiceres
        m11 = inner  (ket [0,0]) finalState1
        m12 = inner  (ket [0,0]) finalState2
        m21 = inner  (ket [0,1]) finalState1
        m22 = inner  (ket [0,1]) finalState2

        -- Den resulterende matrix 
        matOutData = [[m11,m12], [m21,m22]]
        matOut = fromLists matOutData

        -- Calculating the transformed input SV and the output SV
        -- In the case where the matrix is Hermitian we want to look at eigenvalues instead
        matInSV =  eigenvalues matIn
        matInSVT = sortDescendingComplex $ VS.map (polyEval p) matInSV
        matOutSV = sortDescendingComplex $ eigenvalues matOut
    in 
        -- What is the maximum deviation in the SVT (eigenvalues)?
        VS.maximum $ cmap magnitude (matInSVT - matOutSV)

evenHermitianSVTtest :: () -> Double
evenHermitianSVTtest () =
    let 
        -- Et lige polynomium p(x) med p(x)^2 mindre end 1 for alle x i [-1,1] defineres
        p = [0,0,0.5,0,0.5] :: VPoly Rational -- 0.5x^2 + 0.5x^4 

        (pC,qC) = qspPolys p
        svtVec = svtVectorPhi pC qC

        -- Definer en Hermitisk matrix og sæt matrixQubits
        matData =   [[0.0 :+ 0.0, sqrt(2)/4 :+ (-sqrt(2)/4)], 
                        [sqrt(2)/4 :+ sqrt(2)/4, (-1/sqrt(2)) :+ 0.0]]   

        matIn = fromLists matData
        matrixQubits = 1

        -- Lav block encoding
        blockEnc = matrixPrep matIn

        -- SVT udføres
        svtQOp = altPhaseMod matrixQubits svtVec blockEnc
        svtQOt = evalOp $ svtQOp

        -- SVT udføres med (-1 * svtVec) ... med det formål at ...
        svtVecMinus = fmap (* (-1)) svtVec
        svtConjQOp = altPhaseMod matrixQubits svtVecMinus blockEnc
        svtConjQOt = evalOp $ svtConjQOp

        -- ... tage gennemsnit af de tilstande man får ved anvendelse af 
        -- den alternerende fasemodulerende sekvens hørende til blockEnc med svtVec hhv (-svtVec) : 
        tempState1 = (0.5 :+ 0) .* ((apply svtQOt (ket [0,0])) .+ (apply svtConjQOt (ket [0,0])))
        tempState2 = (0.5 :+ 0) .* ((apply svtQOt (ket [0,1])) .+ (apply svtConjQOt (ket [0,1])))
        
        -- De singulære normeringskonstanter beregnes:
        normSV = 0.5*(norm_Frob matIn)^(2 :: Int) + 0.5*(norm_Frob matIn)^(4 :: Int)

        finalState1 = (normSV :+ 0) .* tempState1
        finalState2 = (normSV :+ 0) .* tempState2
     
        -- Der projiceres
        m11 = inner  (ket [0,0]) finalState1
        m12 = inner  (ket [0,0]) finalState2
        m21 = inner  (ket [0,1]) finalState1
        m22 = inner  (ket [0,1]) finalState2

        -- Den resulterende matrix ... har den mon de kvadrerede singulære værdier?
        matOutData = [[m11,m12], [m21,m22]]
        matOut = fromLists matOutData 

        -- Calculating the transformed input SV and the output SV
        -- For Hermitian matrices we want to look at eigenvalues instead
        matInSV =  eigenvalues matIn
        matInSVT = sortDescendingComplex $ VS.map (polyEval p) matInSV
        matOutSV = sortDescendingComplex $ eigenvalues matOut
    in 
        -- What is the maximum deviation in the SVT (eigenvalues)?
        VS.maximum $ cmap magnitude (matInSVT - matOutSV)


----------------------------- Helpers -----------------------------------------

polyEval :: (Real a, Eq b, Fractional b) => VPoly a -> b -> b
polyEval p v = eval (toPoly (V.map (fromRational . toRational) (unPoly p))) v

-- Use this helper when dealing with real-valued Double vectors
sortDescendingReal :: VS.Vector Double -> VS.Vector Double
sortDescendingReal vec = VS.fromList $ sortOn (Down . abs) (VS.toList vec)

-- Use this helper when dealing with complex-valued eigenvalues
sortDescendingComplex :: VS.Vector (Complex Double) -> VS.Vector (Complex Double)
sortDescendingComplex vec = VS.fromList $ sortOn (Down . magnitude) (VS.toList vec)
