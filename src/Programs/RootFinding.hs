module Programs.RootFinding where

import Prelude hiding (gcd, quot)
import qualified Data.Poly as P

import Data.Euclidean (gcd, quot)
import Data.Complex
import qualified Data.Vector as V
import qualified Data.Vector.Generic as VG
import qualified Data.Vector.Unboxed as U


---------------------------------------------------------
-- Types used in RootFinding module
---------------------------------------------------------

-- A root with its original multiplicity
type FinalRoot = (Complex Double, Int)

-- A high-degree factor ready for Aberth-Ehrlich
type ReadyForAberth = (P.UPoly (Complex Double), Int)

-- A real root in [0,1[
type SubZeroRoot = (Complex Double, Int)
-- A complex root in the first quadrant, but not in [0,1[
type OtherFirstQuadRoot = (Complex Double, Int)



---------------------------------------------------------
-- Yuns algoritme til at finde square-free faktorer
-- Dvs. faktorerne består af irreducible delfaktorer der
-- højest optræder med multiplicitet 1
---------------------------------------------------------


-- This function will now accept both VPoly and UPoly
getDegree :: VG.Vector v a => P.Poly v a -> Int
getDegree p = case P.leading p of
    Nothing     -> -1
    Just (d, _) -> fromIntegral d

yun :: P.VPoly Rational -> [(P.VPoly Rational, Int)]
yun f

  | f == 0    = []
  | otherwise = go h g 1
  where
    f' = P.deriv f
    g  = gcd f f'
    h  = f `quot` g

    go currH currG i
      | getDegree currH <= 0 = [] 
      | otherwise =
            let y = gcd currH currG
                z = currH `quot` y

                monicZ = case P.leading z of
                         Nothing -> z
                         Just (_, c) -> z `quot` P.toPoly (V.singleton c)

                next = go y (currG `quot` y) (i + 1)
            in if getDegree z >= 1 
                then (monicZ, i) : next 
                else next

---------------------------------------------------------
-- Kode der processerer resulltater fra Yuns algo
-- Der deles i faktorer der nemt kan løses algebraisk
-- og faktorer der bedst løses numerisk
---------------------------------------------------------

-- | Helper: Convert VPoly Rational to UPoly (Complex Double)
toComplexPoly :: P.VPoly Rational -> P.UPoly (Complex Double)
toComplexPoly p = P.toPoly $ U.convert $ V.map conv (P.unPoly p)
  where
    conv r = (fromRational r :: Double) :+ 0.0

-- | Loop through Yun results and partition them
processYunResults :: [(P.VPoly Rational, Int)] -> ([FinalRoot], [ReadyForAberth])
processYunResults = foldr process ([], [])
  where
    process (poly, mult) (roots, pending)

      | deg == 1 || deg == 2 = (solveSmall poly deg mult ++ roots, pending)
      | deg > 2              = (roots, (toComplexPoly poly, mult) : pending)

      | otherwise            = (roots, pending) -- Skip constants
      where 
        deg = fromIntegral $ maybe 0 fst (P.leading poly)

-- | Exact/closed-form solver for degree 1 and 2
solveSmall :: P.VPoly Rational -> Int -> Int -> [FinalRoot]
solveSmall p deg mult

  | deg == 1 = 
      let a = getCoeff 1; b = getCoeff 0
      in [(fromRational (-b / a) :+ 0, mult)]
  | deg == 2 =
      let a = getCoeff 2; b = getCoeff 1; c = getCoeff 0
          disc = b*b - 4*a*c
          sqrtDisc = sqrt (fromRational disc :+ 0)
          root1 = (fromRational (-b) + sqrtDisc) / fromRational (2 * a)
          root2 = (fromRational (-b) - sqrtDisc) / fromRational (2 * a)
      in [(root1, mult), (root2, mult)]
  | otherwise = []
  where
    getCoeff i = if i < V.length (P.unPoly p) then (P.unPoly p) V.! i else 0

---------------------------------------------------------
-- Aberth's algoritme (Numerisk rodfinding)
-- Anvendes på de faktorer som ikke let kun løses
-- algebraisk i processeringen af Yun-faktorerne
---------------------------------------------------------

-- | 1. Calculate a bounding radius (Guggenheimer's Bound)
-- All roots of the polynomial are guaranteed to lie within this circle.
calcRadius :: P.UPoly (Complex Double) -> Double
calcRadius p = 1 + (U.maximum (U.map magnitude coeffs) / magnitude lead)
  where
    coeffs = P.unPoly p
    lead   = U.last coeffs -- Leading coefficient

-- | 2. Generate Initial Guesses
-- Places 'n' points evenly on the Aberth Circle. 
-- We add 0.4 to the phase to "break symmetry" and avoid axes.
initialGuesses :: Int -> Double -> [Complex Double]
initialGuesses n radius = 
    [ (radius :+ 0) * exp (0 :+ (2 * pi * fromIntegral i / fromIntegral n + 0.4)) 

    | i <- [0..n-1] ]

-- | 3. The Aberth Step
-- Moves all root guesses simultaneously using the repulsion formula.
aberthStep :: P.UPoly (Complex Double) -> P.UPoly (Complex Double) -> [Complex Double] -> [Complex Double]
aberthStep p p' zs = 
    [ z - (offset z / (1 - offset z * repulsion i z)) 

    | (i, z) <- zip [0 :: Int ..] zs ]
  where
    offset z = P.eval p z / P.eval p' z
    repulsion i z = sum [ 1 / (z - zj) | (j, zj) <- zip [0..] zs, i /= j ]

-- | 4. The Convergence Loop
-- Runs until the roots stop moving significantly or we hit a max iteration limit.
solveAberth :: P.UPoly (Complex Double) -> [Complex Double]
solveAberth p 

    | deg < 1   = []
    | otherwise = loop (initialGuesses deg radius) 0
  where
    deg    = fromIntegral (getDegree p)
    radius = calcRadius p
    p'     = P.deriv p
    eps    = 1e-12
    maxIt  = 100 :: Int

    loop zs count

        | count >= maxIt = zs
        | maxDiff < eps  = nextZs
        | otherwise      = loop nextZs (count + 1)
      where
        nextZs  = aberthStep p p' zs
        maxDiff = maximum $ zipWith (\a b -> magnitude (a - b)) zs nextZs

-- | 5. The Runner Function
-- Takes your list of (Polynomial, Multiplicity) and returns [(Root, Multiplicity)]
runAberth :: [(P.UPoly (Complex Double), Int)] -> [(Complex Double, Int)]
runAberth = concatMap (\(poly, mult) -> [(r, mult) | r <- solveAberth poly])

---------------------------------------------------------
-- Get all roots with multiplicities
---------------------------------------------------------

getAllRootsMult :: P.VPoly Rational -> [(Complex Double, Int)]
getAllRootsMult f =
    let (algRoots,readyForAberth) = processYunResults (yun f)
        numericalRoots = runAberth readyForAberth
    in 
        algRoots ++ numericalRoots

getPartitionedFQRootsMult :: P.VPoly Rational -> ([SubZeroRoot], [OtherFirstQuadRoot])
getPartitionedFQRootsMult f =
    let allRoots = getAllRootsMult f
    in
        sortRoots allRoots

-------------------------------------------------------
-- Kun rødder i første kvadrant tages med
-- Der rundes til 16.decimal ... da dette er 
-- præcisionen der regnes med i Complex Double
-------------------------------------------------------
sortRoots :: [(Complex Double, Int)] -> ([(Complex Double, Int)], [(Complex Double, Int)])
sortRoots = foldr classify ([], [])
  where
    eps = 1e-16
    roundTo x | abs x < eps = 0.0

              | otherwise   = x

    classify (z, mult) (subZero, otherFirst)
        -- Rebuild the item with rounded values
        | im == 0.0 && re >= 0.0 && re < 1.0 = (roundedItem : subZero, otherFirst)
        | re >= 0.0 && im >= 0.0             = (subZero, roundedItem : otherFirst)
        | otherwise                          = (subZero, otherFirst)
      where
        re = roundTo (realPart z)
        im = roundTo (imagPart z)
        roundedItem = (re :+ im, mult)

