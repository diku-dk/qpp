-- Exercise scaffold: the placeholder bindings below are for the student to use.
{-# OPTIONS_GHC -Wno-unused-local-binds -Wno-unused-matches #-}
module Main where

import QPP
import QPP.Semantics.MPS

{-| Exercise: the Quantum Fourier Transform.

This exercise accompanies the "Efficient Quantum Algorithms and the QFT"
problem sheet (Exercise 1, 3, 4).

We build QFT_2 (on 2 qubits) fully worked out below as a warm-up (compare
it against your by-hand matrix from paper Exercise 1b using printS), and
leave QFT_3 (3 qubits) for you to complete.

Recall the standard circuit: for qubit k (1-indexed from the top), apply
H, then a controlled phase rotation R_j = diag(1, e^{2 pi i / 2^j}) from
every qubit below it, with j increasing by 1 each time; finally reverse
the qubit order with SWAPs.

Two QPP conventions you need here: rotation angles are `Rational` in
units of pi (so `R Z (1/2)` is a pi/2 rotation), and there is no SWAP
constructor -- a swap of two adjacent qubits is `Permute [1,0]`.
-}

-- | Controlled phase rotation R_j = diag(1, e^{2 pi i / 2^j}), as a QOp.
-- Note e^{i*phi} = cos(phi) + i sin(phi); since R Z is defined as a
-- rotation about the Z axis it will pick up a physically irrelevant
-- *global* phase relative to R_j -- for a *controlled* R Z this global
-- phase becomes a genuine relative phase on the |1>-controlled branch,
-- which is exactly the R_j behaviour we want. TODO: double check this
-- against the matrix you derived in the QFT paper exercise, and adjust
-- the angle formula below if there's a factor of 2 discrepancy (recall
-- from the Pauli/Bloch exercises how R Z theta relates to a Z-phase
-- gate diag(1, e^{i pi theta}) -- remember theta is in units of pi).
cphase :: Int -> QOp
cphase j = C (R Z (2 / 2 ^ j))  -- TODO: verify/fix the angle

-- | QFT on 2 qubits, fully worked out.
qft2 :: QOp
qft2 =
    Compose (Tensor H I)              -- H on qubit 1
    (Compose (cphase 2)               -- controlled-R_2 from qubit 2 onto qubit 1
             (Compose (Tensor I H)    -- H on qubit 2
                       (Permute [1,0]))) -- reverse qubit order
    -- TODO: check this Compose order matches the library's convention
    -- (does `Compose a b` mean "apply a then b", or "apply b then a"?)
    -- Test it against your hand-derived 4x4 matrix from the paper
    -- exercise by applying it to each of |00>,|01>,|10>,|11>.

-- | QFT on 3 qubits -- your turn!
qft3 :: QOp
qft3 = Id 3 -- TODO: generalize the qft2 pattern: H + controlled phases +
            -- final Permute to reverse all 3 qubits. You will need
            -- cphase 2 AND cphase 3 at different points in the circuit.


main :: IO ()
main = do
    putStrLn "-- Exercise 1b: QFT_2 as a circuit, checked against basis states --"
    let op2 = evalOp qft2
    mapM_ (\b -> putStrLn (showState b ++ "  |->  " ++ showState (apply op2 b)))
          [ket [0,0], ket [0,1], ket [1,0], ket [1,1]]
    putStrLn "Compare each output amplitude vector to your hand-derived QFT_4 matrix column."

    putStrLn "\n-- Exercise 3: QFT of a specific basis state --"
    putStrLn $ "QFT_2 |10> = " ++ showState (apply op2 (ket [1,0]))
    putStrLn "Compare to the closed-form formula you derived on paper for x=2."

    putStrLn "\n-- Exercise 4: toy period finding --"
    -- TODO: once qft3 works, build
    --   psi = 0.5 .* ket[0,0,0] .+ 0.5 .* ket[0,1,0] .+ 0.5 .* ket[1,0,0] .+ 0.5 .* ket[1,1,0]
    -- (i.e. equal superposition over 0,2,4,6 out of 8, matching the
    -- paper exercise) and apply qft3 to it, then inspect with printS.
    -- Confirm the result is supported only on |000> and |100> (i.e.
    -- decimal values 0 and 4).
    let psi = (0.5 .* ket [0,0,0]) .+ (0.5 .* ket [0,1,0]) .+ (0.5 .* ket [1,0,0]) .+ (0.5 .* ket [1,1,0])
    putStrLn $ "psi = " ++ showState psi
    let op3 = evalOp qft3
    putStrLn $ "QFT_3 psi = " ++ showState (apply op3 psi)
    putStrLn "Is this supported only on |000> and |100>, as predicted?"
