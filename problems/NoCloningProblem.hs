module Main where

import QPP
import QPP.Semantics.MPS

{-| Exercise: Demonstrate the no-cloning theorem by building a (non-physical) cloned qubit, and comparing to the quantum version of the classical XOR/CX based cloning.

1. write the function cx_cloned_state, which takes a 1Q-state |psi> and produces the 2Q-state CX|psi,0>

2. Inspect the result for |0>, |1>, H|0> = sqrt(1/2)*(|0> + |1>).
3. Inspect the result for |psi> = sqrt(1/3)|0> + sqrt(2/3)|1>
When does the "classical" cloning work to clone the quantum state?
4. Show that the differences in wave-function reflect different measurement probabilities (show first the probabilities for measuring qubit 1, then conditional probabilities for measuring qubit 0 given result 0 on Q1.
5. This proves the No-Cloning Theorem. Why? (Hint, you have made a linear function that matches the cloning function on the basis vectors |0> and |1>).

-}

{-|
  Clones the full state of a 1-qubit state |ψ> into a 2-qubit state |ψψ>
  Not a physically realizable operation due to the No-cloning theorem!
 -}
quantum_cloned_state :: StateT -> StateT
quantum_cloned_state ψ = ψ ⊗ ψ

cx_cloned_state :: StateT -> StateT
cx_cloned_state ψ = ψ ⊗ ket [0] -- Dummy placeholder. Write this function yourself!


main :: IO()
main = do
-- Show that CX correctly clones |0> and |1>
    let q00 = quantum_cloned_state $ ket [0]
        c00 = cx_cloned_state      $ ket [0]
    putStr $ "-- Can we clone |0> ? --\n" ++
             "Quantum clone |00> = " ++ (showState q00) ++ "\n" ++
             "XOR cloned    |00> = " ++ (showState c00) ++ "\n"

    --- et cetera: Fill in the rest.
