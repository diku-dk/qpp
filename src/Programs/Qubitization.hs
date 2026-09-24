module Programs.Qubitization where

import QPP
import Programs.RootFinding

import Data.Complex (Complex((:+)), realPart, imagPart,magnitude,polar)
import qualified Data.Vector as V
import Data.List (foldl')

import Data.Poly (VPoly, toPoly, unPoly,leading, scale)



--------------------------------------------------------------------------------
------------- Quantum signal processing-----------------------------------
--------------------------------------------------------------------------------

-- qspVectorPhi calculates the quantum signal processing vector  
-- defined in Gilyen et al. Theorem 3.
-- Input:
-- 1) p, q complex polynomials
-- 2) p(x)^2 + (1 - x^2)q(x)^2 = 1 for x in [-1,1] :
-- 3) (p even, q odd) or (p odd, q even)
-- 4) deg(p) = deg(q) + 1
-- Output:
-- Real k+1 vector where k = deg(p)
qspVectorPhi :: VPoly ComplexT -> VPoly ComplexT -> [Double]
qspVectorPhi p q
    | Just (0, _) <- leading p =
        let 
            (_, phi) = polar (unPoly p V.! 0)
        in  
            [phi]           

    | otherwise =
        let
            -- Calculate phi value
            leadingOrderTermP = case leading p of
                                Just (_deg, c) -> c
                                Nothing        -> error "Polynomial is empty or zero"

            leadingOrderTermQ = case leading q of
                                Just (_deg, c) -> c
                                Nothing        -> error "Polynomial is empty or zero"

            lTermFrac = leadingOrderTermP / leadingOrderTermQ
            (_, theta) = polar lTermFrac
            phi = theta / 2

            xPoly =  toPoly $ V.fromList [0,1]
            x2Poly = toPoly $ V.fromList [1, 0, -1]

            pNextTerm1 = xPoly * p
            pNextTerm2 = scale 0 lTermFrac (x2Poly * q)
            pNext = scale 0 (exp (-(0.0 :+ phi))) (pNextTerm1 + pNextTerm2 )

            qNextTerm1 = scale 0 lTermFrac (xPoly * q)
            qNext = scale 0 (exp (-(0.0 :+ phi))) (qNextTerm1 - p)

            pNextFinal = removeSmallTrailing pNext
            qNextFinal = removeSmallTrailing qNext
        in
            (qspVectorPhi pNextFinal qNextFinal) ++ [phi] 

--------------------------------------------------------------------------------
------------- From real to complex polynomial-----------------------------------
--------------------------------------------------------------------------------

-- qspPolys calculates the pair of complex polynomials satisfying the conditions in Gilyen et al. (2018) Theorem 5.
-- The polynomials are used in the algorithm (qspVectorPhi) that finds the quantum signal processing vector.
-- Input conditions:
-- 1) P is a real polynomial
-- 2) P is an even or an odd pylynomial
-- 3) P(x)^2 <= 1 for x in [-1,1]
-- Output conditions
-- 1) P1, Q1 complex polynomials
-- 2) P1(x)^2 + (1 - x^2)Q1(x)^2 = 1 for x in [-1,1] :
-- 3) (P1 even, Q1 odd) or (P1 odd, Q1 even)
-- 4) deg(P) = deg(Q) + 1
qspPolys :: VPoly Rational -> (VPoly ComplexT, VPoly ComplexT)
qspPolys realParityPoly =
        let realPolyPair = createRealPolyPair (1 - realParityPoly * realParityPoly)
            realPoly1 = poly1 realPolyPair
            realPoly2 = poly2 realPolyPair

            complexEvenPosPoly = toPoly $ V.map (\r -> fromRational r :+ 0.0) (unPoly realParityPoly)
            complexPoly1       = toPoly $ V.map (:+ 0.0) (unPoly realPoly1)
            complexPoly2       = toPoly $ V.map (:+ 0.0) (unPoly realPoly2)

            pPoly = complexEvenPosPoly + scale 0 (0 :+ 1) complexPoly1
            qPoly = scale 0 (0 :+ 1) complexPoly2
        in
            (pPoly, qPoly)


-- Given a real even poly p(x) with deg(p) <= 2k and p(x)>= 0 on [-1,1] 
-- return real polys B, C as in Gilyen et al. Lemma 6.
-- The method is a step in the construction of a suitable complex polynomial
-- used in the qubitization algorithm
-- Input conditions:
-- 1) p is a real polynomial
-- 2) p is even
-- 3) p has degree 2k
-- 3) p(x)^2 >= 0 for x in [-1,1]
-- Output: a polypair (p1,q1) satisfying
-- 1) p1, q1 real polynomials
-- 2) (p1 even, q1 odd) or (p1 odd, q1 even)
-- 3) degree(p1) = k, degree(q1) = k-1
-- 4) p(x) = p1(x)^2 + (1-x^2)q1(x)^2
createRealPolyPair :: VPoly Rational -> PolyPair
createRealPolyPair p = 
    let 
        -- 1. Find the initialCoeff
        leadingCoeff = case leading p of
                Nothing         -> error "The polynomial is zero."
                Just (_, coeff) -> coeff

        val = sqrt (abs (fromRational leadingCoeff))
        initialCoeff = PolyPair (toPoly (V.singleton val), toPoly (V.singleton 0))

        -- 2. Find all roots
        (subOneRoots, otherFQroots) = getPartitionedFQRootsMult p

        -- 3. Process the subOneRoots
        subOnePairsProd = product $ map process subOneRoots
            where
                process (z, k)
                    | even k    = (rootToPolyPair z) ^ (k `div` 2)
                    | otherwise = error $ "subOnePairs: odd multiplicity " ++ show k ++ " for root " ++ show z

        -- 4. otherFQroots
        remainingPairsProd = product [ (rootToPolyPair z) ^ k | (z, k) <- otherFQroots ]

    in 
        --removeTrailingPPZeros (). Burde ikke være nødv med Data.Poly
        initialCoeff * subOnePairsProd * remainingPairsProd


-- Input: 
-- Complex number z i the 1.quadrant incl. the posive axis and origo
-- Output: A polypair (p,q)
-- 1) p, q real polynomials
-- 2) (p even, q odd) or (p odd, q even)
-- 3) degree(p) = k, degree(q) = k-1
-- 4) 0 <= p(x)^2 + (1-x^2)q(x)^2
rootToPolyPair :: ComplexT -> PolyPair
rootToPolyPair z
    -- Case: Re z > 0 and Im z > 0
    | re >= epsilon && im >= epsilon = 
        let c = re^(2 :: Int) + im^(2 :: Int) + sqrt( (2 * (re^(2 :: Int) + 1) * im^(2 :: Int)) + (re^(2 :: Int) - 1)^(2 :: Int) + im^(4 :: Int) )
            rePol = toPoly $ V.fromList [-(re^(2 :: Int) + im^(2 :: Int)), 0, c]
            imPol = toPoly $ V.fromList [0, sqrt(c^(2 :: Int) - 1)]
        in 
            PolyPair (rePol, imPol)

    -- Case: z = 0
    | abs re < epsilon && abs im < epsilon = 
        let
            rePol = toPoly $ V.fromList [0,1]
            imPol = toPoly $ V.fromList [0]
        in
            PolyPair (rePol, imPol)

    -- Case: z in (0, 1) (Real axis between 0 and 1)
    | abs im < epsilon && re >= epsilon && re <= 1 - epsilon = 
        let
            rePol = toPoly $ V.fromList [-re^(2 :: Int),0,1]
            imPol = toPoly $ V.fromList [0]
        in
            PolyPair (rePol, imPol)

    -- Case: z in [1, infinity) (Real axis >= 1)
    | abs im < epsilon && re > 1 - epsilon = 
        let
            rePol = toPoly $ V.fromList [0,sqrt((max re 1)^(2 :: Int) - 1)]
            imPol = toPoly $ V.fromList [re]
        in
            PolyPair (rePol, imPol)

    -- Case: Re z = 0 and Im z > 0 (Purely imaginary)
    | abs re < epsilon && im >= epsilon = 
        let
            rePol = toPoly $ V.fromList [0,sqrt(im^(2 :: Int) + 1)]
            imPol = toPoly $ V.fromList [im]
        in
            PolyPair (rePol, imPol)

    -- Fallback/Default
    | otherwise = error "Unspecified complex region"
  where
    re = realPart z
    im = imagPart z
    epsilon = 1e-9

-- Polypair is a pair of real polyomials with a distint multipication
-- used in Gilyen et al. Lemma 6.
newtype PolyPair = PolyPair (VPoly Double, VPoly Double)
    deriving (Show, Eq)

instance Num PolyPair where
    -- (p, q) + (s, t) = (p+s, q+t)
    (PolyPair (p, q)) + (PolyPair (s, t)) = 
        PolyPair ( p + s, q + t)

    -- (p, q) - (s, t) = (p-s, q-t)
    (PolyPair (p, q)) - (PolyPair (s, t)) = 
        PolyPair ( p - s,  q - t)

    -- (p, q) * (s, t) = (p*s - (1-x^2)*q*t, p*t + q*s)
    (PolyPair (p, q)) * (PolyPair (s, t)) = 
        let ps = p * s
            qt = q * t
            pt = p * t
            qs = q * s
            term2 =  oneMinusX2 * qt
        in PolyPair (ps - term2, pt + qs)

    -- Lifting an integer constant c to (c, 0)
    --fromInteger n = PolyPair ([fromInteger n], [0])
    fromInteger n = 
        let val = fromInteger n
        in PolyPair (toPoly (V.singleton val), toPoly (V.singleton 0))


    -- Abs and signum are typically required for Num but often left undefined 
    -- for complex-like polynomial structures.
    abs = error "abs not defined for PolyPair"
    signum = error "signum not defined for PolyPair"

-- Define the polynomial (1 - x^2) as a list of coefficients: 1 + 0x - 1x^2
oneMinusX2 :: VPoly Double
oneMinusX2 = toPoly $ V.fromList [1, 0, -1]

poly1 :: PolyPair -> VPoly Double
poly1 (PolyPair (p, _)) = p

poly2 :: PolyPair -> VPoly Double
poly2 (PolyPair (_, q)) = q

--------------------------------------------------------------------------------
------------- Helpers ----------------------------------------------------------
--------------------------------------------------------------------------------

-- isEven      
-- parityInt 

-- Typical degree calculation for dsp-style lists [a0, a1, ..., an]
class HasMagnitude a where
    getMag :: a -> Double -- or (RealFloat b => b)

instance HasMagnitude Double where
    getMag = abs

instance HasMagnitude (Complex Double) where
    getMag = magnitude

polyDegree :: (HasMagnitude a) => [a] -> Int
polyDegree as = foldl' findMaxIndex (-1) (zip [0..] as)
  where
    epsilon = 1e-3
    findMaxIndex acc (i, x) 

        | getMag x > epsilon = i
        | otherwise          = acc

getLeadingOrderTerm :: (Eq a, Num a) => [a] -> Maybe (a,Int)
getLeadingOrderTerm coeffs = 
    case dropWhile (== 0) (reverse coeffs) of
        []    -> Nothing                -- The polynomial is zero
        (x:_) -> Just (x, length (dropWhile (== 0) (reverse coeffs)) - 1)

removeSmallTrailing :: VPoly ComplexT -> VPoly ComplexT
removeSmallTrailing p = 
    let epsilon = 1e-9
        coeffs = unPoly p
        -- Find the index of the last element that is NOT small
        lastSignificant = V.findIndexR (\c -> magnitude c >= epsilon) coeffs
    in case lastSignificant of
        Nothing -> toPoly V.empty -- Everything was small
        Just i  -> toPoly (V.take (i + 1) coeffs)

