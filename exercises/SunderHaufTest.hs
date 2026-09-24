{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE PatternSynonyms #-}

module Main where

import Programs.SunderhaufOptPoly

import Data.Poly (eval, unPoly)

-- ============================================================================
-- Example Harness
-- ============================================================================
main :: IO ()
main = do
    let kappa   = 5.0      -- Matrix condition number
        epsilon = 1e-3     -- Uniform error accuracy target
        poly    = sunderhaufOptPoly kappa epsilon
        coeffs  = generateChebyshevCoeffs kappa epsilon
    
    putStrLn "=== Sünderhauf Optimal Polynomial via Data.Poly ==="
    putStrLn $ "Condition Number (kappa): " ++ show kappa
    putStrLn $ "Error Tolerance (epsilon): " ++ show epsilon
    
    -- Print out the unboxed internal coefficient vector representation
    putStrLn $ "\nDense Vector Representation (Constant term up to X^d):"
    print (Data.Poly.unPoly poly)
    
    -- Verification test at a target singular value in the active domain
    let testX = 0.5 -- safely within [1/5, 1]
    putStrLn "\n--- Numerical Verification ---"
    putStrLn $ "Ideal Target Value (1/x) : " ++ show (1.0 / testX)
    putStrLn $ "Data.Poly Evaluated P(x) : " ++ show (eval poly testX)
    
    putStrLn "=== Sünderhauf Pure Chebyshev Coefficients ==="
    putStrLn $ "Condition Number (kappa): " ++ show kappa
    putStrLn $ "Error Tolerance (epsilon): " ++ show epsilon
    putStrLn $ "Total Odd Terms Generated: " ++ show (length coeffs)
    
    putStrLn "\nList of Coefficients (Stable and Ready for QSP Phase Solver):"
    mapM_ (\(deg, val) -> putStrLn $ "  T_" ++ show deg ++ " coefficient: " ++ show val) coeffs

