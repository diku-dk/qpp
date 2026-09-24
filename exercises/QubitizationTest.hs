{-# LANGUAGE OverloadedLists #-}
module Main where 
import QPP
import Programs.Qubitization
import Data.Poly


main :: IO ()
main = do
    print (show (createRealBCTest())) 
    print (show (qspPolysTest()))
    print (show (qspVectorPhiTest()))

-----------------------------------------------------------------  
--
--   HUSK IMPLEMENTER ALLE BETINGELSER (Fra sætningerne)
--
-----------------------------------------------------------------

createRealBCTest :: () -> PolyPair
createRealBCTest () =
    let 
        -- Given a real even poly A(x) with deg(A) <= 2k and A(x)>= 0 on [-1,1] 
        -- return real polys B, C as in Lemma 6
        pp = [1,0,0,0,-0.25,0, -0.5,0, -0.25] :: VPoly Rational 
    in
        createRealPolyPair pp

qspVectorPhiTest :: () -> [Double]
qspVectorPhiTest () =
    let 
        p = [0,0,0.5,0,0.5] :: VPoly Rational
        (pC,qC) = qspPolys p
    in
        qspVectorPhi pC qC

qspPolysTest :: () -> (VPoly ComplexT,VPoly ComplexT)
qspPolysTest () = 
    let 
        p = [0,0,0.5,0,0.5] :: VPoly Rational
        (pC,qC) = qspPolys p
    in
        (pC,qC)