module Main where

import QPP
import QPP.Semantics.MPS

import Data.Complex (Complex((:+)), magnitude)

{-| SOLUTIONS: the Quantum Fourier Transform.

This exercise accompanies the "Efficient Quantum Algorithms and the QFT"
problem sheet (Exercise 1, 3, 4).

Two QPP conventions do all the work here:

  * `Compose a b` (= `a ∘ b`) applies b FIRST. The left-to-right alias `>:`
    reads in circuit order, so `a >: b >: c` means "a, then b, then c". Every
    circuit below is written that way.

  * Qubit 0 is the most significant bit: `ket [1,0,0]` is |4>, and `C op`
    controls on qubit 0, the top of its window.
-}

-- | Compose a list of layers in circuit order: the first element acts first.
circuit :: [QOp] -> QOp
circuit = foldr1 (>:)

--------------------------------------------------------------------------
-- Exercise 1a: the controlled phase gate
--------------------------------------------------------------------------

{-| The 1-qubit phase gate diag(1, e^{i pi theta}).

`Phase theta` is scalar multiplication by e^{i pi theta} — a 0-qubit
operator — and `⊕` glues two 0-qubit operators into a 1-qubit diagonal
("if qubit 0 then the right one, else the left one"). So this literally
says: 1 on the |0> branch, e^{i pi theta} on the |1> branch.
-}
phaseGate :: Rational -> QOp
phaseGate theta = One ⊕ Phase theta

{-| Controlled phase rotation R_j = diag(1, 1, 1, e^{2 pi i / 2^j}).

This is the answer to the sheet's "adjust the angle formula if there's a
factor of 2 discrepancy" — the discrepancy is real, but it is not a factor
of 2, and it is not fixed by changing the angle at all:

    R Z theta  =  exp(-i pi theta Z / 2)  =  diag(e^(-i pi theta/2), e^(i pi theta/2))

which is the phase gate diag(1, e^(i pi theta)) times the global phase
e^(-i pi theta/2). A global phase on a 1-qubit gate is unobservable, but
`C` turns it into a *relative* phase between the |0..> and |1..> blocks:

    C (R Z theta)  =  diag(1, 1, e^(-i pi theta/2), e^(i pi theta/2))

which is not diag(1,1,1,e^{i pi theta}) for any theta. So we build the
phase gate directly instead of borrowing a Z-rotation. (`main` prints the
damage: with `C (R Z ...)` the QFT is exact on |00> and |10> and wrong on
|01> and |11>, which is easy to miss if you only ever test basis states.)
-}
cphase :: Int -> QOp
cphase j = cphaseAt j 1

{-| `cphaseAt j d` is R_j between qubit 0 and the qubit d places below it,
    leaving the d-1 qubits in between alone. Controlled phases are symmetric
    in control and target, so it doesn't matter which of the two we call the
    control — which is what lets a single `C` reach past idle qubits. -}
cphaseAt :: Int -> Int -> QOp
cphaseAt j d = C (Id (d - 1) ⊗ phaseGate (2 / 2 ^ j))

--------------------------------------------------------------------------
-- Exercise 1b: QFT on 2 and 3 qubits
--------------------------------------------------------------------------

-- | QFT on 2 qubits.
qft2 :: QOp
qft2 = H ⊗ I            -- H on qubit 0
    >: cphase 2         -- R_2 between qubits 0 and 1
    >: I ⊗ H            -- H on qubit 1
    >: Permute [1,0]    -- reverse the qubit order

-- | QFT on 3 qubits: the same pattern one row deeper.
qft3 :: QOp
qft3 = circuit
    [ H @> 2                -- H on qubit 0
    , cphaseAt 2 1 @> 1     -- R_2 between qubits 0 and 1
    , cphaseAt 3 2          -- R_3 between qubits 0 and 2 (qubit 1 idle)
    , 1 <@ H @> 1           -- H on qubit 1
    , 1 <@ cphaseAt 2 1     -- R_2 between qubits 1 and 2
    , 2 <@ H                -- H on qubit 2
    , Permute [2,1,0]       -- reverse the qubit order
    ]

-- | The circuit the sheet hands you as a warm-up, with the Z-rotation still in
--   it. Kept so `main` can show exactly how it fails.
qft2Naive :: QOp
qft2Naive = H ⊗ I >: C (R Z (1/2)) >: I ⊗ H >: Permute [1,0]

--------------------------------------------------------------------------
-- Checking a circuit against the closed form
--------------------------------------------------------------------------

-- | Big-endian bit list: `bitsOf 3 5 == [1,0,1]` (qubit 0 is most significant).
bitsOf :: Int -> Int -> [Int]
bitsOf n x = [ (x `div` 2 ^ i) `mod` 2 | i <- [n-1, n-2 .. 0] ]

-- | The amplitudes <y|psi> for y = 0 .. 2^n-1.
amplitudes :: Int -> StateT -> [ComplexT]
amplitudes n psi = [ inner (ket (bitsOf n y)) psi | y <- [0 .. 2 ^ n - 1] ]

-- | Column x of the DFT matrix: QFT|x> = (1/sqrt N) sum_y e^{2 pi i x y / N} |y>.
qftColumn :: Int -> Int -> [ComplexT]
qftColumn n x =
    let nn = 2 ^ n :: Int
        s  = 1 / sqrt (fromIntegral nn)
    in  [ (s :+ 0) * exp (0 :+ (2 * pi * fromIntegral (x * y) / fromIntegral nn))
        | y <- [0 .. nn - 1] ]

-- | Largest deviation between a circuit and the closed-form QFT, over all
--   2^n computational basis inputs. This is the honest test: agreeing on a
--   few basis states is not the same as being the same operator.
maxDeviation :: Int -> QOp -> Double
maxDeviation n op = maximum
    [ magnitude d
    | x <- [0 .. 2 ^ n - 1]
    , d <- zipWith (-) (amplitudes n (apply (evalOp op) (ket (bitsOf n x))))
                       (qftColumn n x) ]

showColumn :: [ComplexT] -> String
showColumn = unwords . map (\z -> pad (showRounded z))
  where
    pad s               = s ++ replicate (20 - length s) ' '
    showRounded (a :+ b) = show (r a) ++ (if r b < 0 then "-" else "+") ++ show (abs (r b)) ++ "i"
    r x                 = fromIntegral (round (x * 1000) :: Int) / 1000 :: Double

main :: IO ()
main = do
    putStrLn "-- Exercise 1b: QFT_2 as a circuit, checked against basis states --"
    let op2 = evalOp qft2
    mapM_ (\x -> putStrLn (" |" ++ concatMap show (bitsOf 2 x) ++ ">  |->  "
                            ++ showColumn (amplitudes 2 (apply op2 (ket (bitsOf 2 x))))))
          [0 .. 3]
    putStrLn $ "Closed form, column by column, for comparison:"
    mapM_ (\x -> putStrLn ("   x = " ++ show x ++ "        " ++ showColumn (qftColumn 2 x))) [0 .. 3]
    putStrLn $ "Max deviation from the closed form: " ++ show (maxDeviation 2 qft2)

    putStrLn "\n-- Why the Z-rotation version is wrong --"
    putStrLn $ "Max deviation for `C (R Z (1/2))`: " ++ show (maxDeviation 2 qft2Naive)
    mapM_ (\x -> putStrLn (" |" ++ concatMap show (bitsOf 2 x) ++ ">  |->  "
                            ++ showColumn (amplitudes 2 (apply (evalOp qft2Naive) (ket (bitsOf 2 x))))))
          [0 .. 3]
    putStrLn "The leftover global phase of R Z is a stray phase gate diag(1, e^(-i pi/4))"
    putStrLn "sitting on the control wire; the final swap moves it onto the low bit, so"
    putStrLn "every odd y picks up e^(-i pi/4). Every column is wrong -- but every"
    putStrLn "amplitude still has magnitude 1/2, so measuring the output of a single"
    putStrLn "basis state in the computational basis cannot see the error at all. Only"
    putStrLn "the phases give it away, which is exactly what the QFT is for."

    putStrLn "\n-- Exercise 3: QFT of a specific basis state --"
    putStrLn $ "QFT_2 |10> = " ++ showState (apply op2 (ket [1,0]))
    putStrLn "x = 2, N = 4: amplitudes (1/2) e^{2 pi i 2 y / 4} = (1/2)(1, -1, 1, -1). Matches."

    putStrLn "\n-- Exercise 4: toy period finding --"
    putStrLn $ "Max deviation of qft3 from the closed form: " ++ show (maxDeviation 3 qft3)
    let psi = (0.5 .* ket [0,0,0]) .+ (0.5 .* ket [0,1,0])
          .+ (0.5 .* ket [1,0,0]) .+ (0.5 .* ket [1,1,0])
    putStrLn $ "psi        = " ++ showState psi
    putStrLn $ "QFT_3 psi  = " ++ showState (apply (evalOp qft3) psi)
    putStrLn ""
    putStrLn "psi is the uniform superposition over x in {0,2,4,6}: period 2 in a"
    putStrLn "register of size 8. The QFT turns that into"
    putStrLn "    sum_{m=0..3} e^{2 pi i (2m) y / 8}  =  sum_m e^{2 pi i m y / 4}"
    putStrLn "which is 4 when y is a multiple of 4 and 0 otherwise -- so the output is"
    putStrLn "supported exactly on |000> and |100>, i.e. y in {0,4} = {0, N/period}."
    putStrLn "Reading the period off the *output* spacing rather than the input is the"
    putStrLn "whole trick behind Shor's algorithm."
