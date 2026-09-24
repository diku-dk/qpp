-- Exercise scaffold: the placeholder bindings below are for the student to use.
{-# OPTIONS_GHC -Wno-unused-local-binds -Wno-unused-matches #-}
module Main where

import QPP
import QPP.Semantics.MPS

{-| Exercise: State preparation and the Input Problem.

This exercise accompanies the "State Preparation and the Input Problem"
problem sheet (Exercise 1, 2, 5).

Rotation angles are `Rational` numbers in units of pi throughout, so
`R Y (1/2)` is a pi/2 rotation about the Y axis.

1. Complete `prepareBasisState`, which prepares a computational basis
   state |b1 b2 ... bn> from |00...0> using only X gates on the qubits
   where b_k = 1, built as a single n-qubit QOp via `Tensor` (paper
   Exercise 1).

2. Complete `controlledRotationPrep2`, a hand-built state-prep circuit
   for a *general real* 2-qubit state, using one Hadamard-like R_Y
   rotation on qubit 1 followed by a uniformly-controlled R_Y rotation
   on qubit 2 (i.e. two DIFFERENT rotation angles on qubit 2, selected
   by the value of qubit 1). This is the n=2 base case of the
   recursive construction from paper Exercise 2 -- work out what
   angles you need for a specific target state, and build the circuit
   with `R Y theta`, `C`, and `Compose`/`Tensor`.

3. Complete `gateCount`, a simple structural recursion over the `QOp`
   AST that counts how many "elementary" gates (leaves like X, Y, Z, H,
   SX, R _ _, Permute; a `C op` or `DirectSum op1 op2` counts its
   sub-circuit(s) plus itself) appear in a circuit. Use it to see how
   circuit size grows for the recursive state-prep construction
   (paper Exercise 2), by building it for n = 1, 2, 3 qubits.

4. (Discussion, paper Exercise 4 -- no code needed) Why would loading N
   arbitrary real numbers via a QRAM-style structure still cost you
   Omega(N) gates even though queries have O(log N) *depth*? Write your
   answer as a comment in `main`.
-}

-- | Prepare |b1..bn> from |00..0> using X on exactly the qubits with
-- b_k = 1. `bits` is e.g. [1,0,1] for |101>.
prepareBasisState :: [Int] -> QOp
prepareBasisState bits = I -- TODO: build via `Tensor`, folding X/I according to `bits`

-- | A hand-built 2-qubit real-amplitude state-prep circuit:
-- rotate qubit 1 by angle0, then apply a *different* R_Y rotation to
-- qubit 2 depending on qubit 1's value (a "uniformly controlled
-- rotation", the building block from paper Exercise 2).
controlledRotationPrep2 :: Rational -> Rational -> Rational -> QOp
controlledRotationPrep2 angle0 angle1if0 angle1if1 =
    Tensor (R Y angle0) I  -- TODO: this only handles the qubit-1 rotation.
    -- You still need to apply R Y angle1if0 to qubit 2 when qubit 1 = 0,
    -- and R Y angle1if1 when qubit 1 = 1. Hint: DirectSum applies its
    -- two arguments conditioned on a control qubit's value -- but check
    -- carefully which qubit plays the role of control here, and how to
    -- Compose this with the qubit-1 rotation above.

-- | Count "elementary" gates in a QOp circuit description.
gateCount :: QOp -> Int
gateCount I           = 0  -- identity doesn't really cost anything
gateCount op          = 1  -- TODO: this is wrong for everything except
                            -- single-qubit leaves! Pattern-match on X, Y,
                            -- Z, H, SX, R _ _, Permute, C, Tensor,
                            -- DirectSum, Compose separately and recurse
                            -- properly.


main :: IO ()
main = do
    putStrLn "-- Exercise 1: basis state preparation --"
    let target  = [1, 0, 1]
        prepOp  = evalOp (prepareBasisState target)
        psi000  = ket [0,0,0]
    putStrLn $ "Preparing |101>: " ++ showState (apply prepOp psi000)

    putStrLn "\n-- Exercise 2: a hand-built 2-qubit state-prep circuit --"
    -- TODO: pick a target 2-qubit real-amplitude state (e.g. from the
    -- tensor-products exercise sheet) and work out the three angles
    -- that `controlledRotationPrep2` needs to reach it starting from
    -- |00>. Then check your circuit with showState/printS.
    let circuit = evalOp (controlledRotationPrep2 0 0 0)  -- TODO: real angles
    putStrLn $ "Result: " ++ showState (apply circuit (ket [0,0]))

    putStrLn "\n-- Exercise 3: counting gates as n grows --"
    -- TODO: once you've generalized `controlledRotationPrep2` in your
    -- head to n qubits (you don't have to implement the full general
    -- version!), at least build the n=1 and n=2 cases here and print
    -- their gateCount, to sanity check against the T(n) formula you
    -- derived on paper.
    putStrLn $ "gateCount (n=2 example) = " ++ show (gateCount (controlledRotationPrep2 0 0 0))

    putStrLn "\n-- Exercise 4 (discussion) --"
    -- TODO: write your QRAM discussion answer here as a comment.
    return ()
