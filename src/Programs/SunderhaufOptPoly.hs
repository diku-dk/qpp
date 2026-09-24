{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE PatternSynonyms #-}

module Programs.SunderhaufOptPoly where

import Data.Poly (UPoly,pattern X, toPoly, unPoly,subst, monomial)
import qualified Data.Vector.Unboxed as V


chebyshevT :: Int -> UPoly Double
chebyshevT n = polys !! n
  where
    polys = 1 : X : zipWith (\t1 t2 -> 2 * X * t1 - t2) (tail polys) polys


-- | Computes the optimal polynomial for QSVT matrix inversion 
-- using Sünderhauf's exact analytical method.
-- From the article : 
-- Matrix inversion polynomials for the quantum singular value transformation
-- Sünderhauf et al.
sunderhaufOptPoly :: Double -> Double -> UPoly Double
sunderhaufOptPoly kappa epsilon =
    let a = 1.0 / kappa
        -- Optimal degree bound O(kappa * log(kappa / epsilon))
        d = oddDegreeCeil (kappa * log (kappa / epsilon))
        
        -- Analytical linear combination of mapped Chebyshev polynomials
        terms = [ scalePoly (coeffC k a d) (mappedChebyshev k a) | k <- [1, 3 .. d] ]
    in sum terms
  where
    oddDegreeCeil :: Double -> Int
    oddDegreeCeil x = let n = ceiling x in if odd n then n else n + 1

    -- Exact coefficient formula derived by Sünderhauf et al.
    coeffC :: Int -> Double -> Int -> Double
    coeffC k a d =
        let factor = 4.0 / (fromIntegral d * (1.0 - a*a))
            term1  = sin (fromIntegral k * pi / fromIntegral (2 * d))
        in factor * term1

    scalePoly :: Double -> UPoly Double -> UPoly Double
    scalePoly c p = toPoly (V.map (* c) (unPoly p))

    -- Maps interval x -> (2*x^2 - 1 - a^2) / (1 - a^2) inside the Chebyshev basis
    mappedChebyshev :: Int -> Double -> UPoly Double
    mappedChebyshev k a =
        let baseT = chebyshevT k
            -- Using monomial 2 1 generates exactly: 1.0 * x^2
            -- This explicitly avoids ambiguous literal inference errors for X^2
            xSquared = monomial 2 (1.0 :: Double)
            
            -- Construct the transformation poly matching the exact type constraint
            transformPoly = toPoly [-(1.0 + a*a) / (1.0 - a*a)] + scalePoly (2.0 / (1.0 - a*a)) xSquared
        in subst baseT transformPoly -- FIXED: Use subst instead of eval for polynomial composition


-- Returns pairs of (Degree, Coefficient) for odd terms only.
generateChebyshevCoeffs :: Double -> Double -> [(Int, Double)]
generateChebyshevCoeffs kappa epsilon =
    let a = 1.0 / kappa
        -- Optimal degree bound O(kappa * log(kappa / epsilon))
        d = oddDegreeCeil (kappa * log (kappa / epsilon))
    in [ (k, coeffC k a d) | k <- [1, 3 .. d] ]
  where
    oddDegreeCeil :: Double -> Int
    oddDegreeCeil x = let n = ceiling x in if odd n then n else n + 1

    -- Exact coefficient formula derived by Sünderhauf et al.
    coeffC :: Int -> Double -> Int -> Double
    coeffC k a d =
        let factor = 4.0 / (fromIntegral d * (1.0 - a*a))
            term1  = sin (fromIntegral k * pi / fromIntegral (2 * d))
        in factor * term1