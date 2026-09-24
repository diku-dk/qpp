{-# LANGUAGE OverloadedLists #-}
module Main where 

import Programs.RootFinding

import Data.Complex
import Data.Poly
import qualified Data.Poly as P

----------------------------------------------------------------
-- Main
----------------------------------------------------------------

main :: IO ()
main = do
    putStrLn "Square-free factorisation:"
    print (show (yunTest()))

    putStrLn "Processing of yun factors:"
    print (show (processYunResultsTest()))

    putStrLn "Processing of Aberth factors:"
    print (show (runAberthTest()))

    putStrLn "Finding all roots and multiplicities:"
    print (show (getAllRootsMultTest()))

    putStrLn "Finding all SORTED roots and multiplicities:"
    print (show (getPartitionedFQRootsMultTest()))


yunTest :: () -> [(P.VPoly Rational, Int)]
yunTest () =
    let f :: VPoly Rational
        f = (X^2 + 1) * (X - 2)^2 * (X+1)^3
    in 
        yun f

processYunResultsTest :: () -> ([FinalRoot], [ReadyForAberth])
processYunResultsTest () =
    let f :: VPoly Rational
        f = (X^2 + 1) * (X - 2)^2 * (X+1)^3 * (X^3 - 2)^5
    in 
        processYunResults (yun f)

runAberthTest :: () -> [(Complex Double, Int)]
runAberthTest () =
    let f :: VPoly Rational
        f = (X^2 + 1) * (X - 2)^2 * (X+1)^3 * (X^3 - 2)^5
        (_,readyForAberth) = processYunResults (yun f)
    in 
        runAberth readyForAberth

getAllRootsMultTest :: () -> [(Complex Double, Int)]
getAllRootsMultTest () =
    let f :: VPoly Rational
        f = (X^2 + 1) * (X - 2)^2 * (X+1)^3 * (X^3 - 2)^5
    in 
        getAllRootsMult f


getPartitionedFQRootsMultTest :: () -> ([SubZeroRoot], [OtherFirstQuadRoot])
getPartitionedFQRootsMultTest () =
    let x = X :: VPoly Rational
        g = [-0.5, 1] :: VPoly Rational
        f = (x^2 + 1) * (x - 2)^2 * (x + 1)^3 * (x^3 - 2)^5 * x^2 * g^6 * (x-1)^7
    in 
        getPartitionedFQRootsMult f