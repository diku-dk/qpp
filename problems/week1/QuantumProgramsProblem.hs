-- Exercise scaffold: the placeholder bindings below are for the student to use.
{-# OPTIONS_GHC -Wno-unused-local-binds -Wno-unused-matches #-}
module Main where

import QPP
import QPP.Semantics.MPS

{-| Exercise: Quantum Programs.

This exercise accompanies the "Quantum Programs" problem sheet
(Exercise 1, 2, 3, 4). It builds directly on the `QOp` grammar from the
lecture:

  data QOp = Id Nat | Phase Rational
           | X | Y | Z | H | SX | R QOp Rational
           | C QOp | Permute [Int]
           | Tensor QOp QOp | DirectSum QOp QOp | Compose QOp QOp
           | Adjoint QOp

with `I = Id 1` the 1-qubit identity, and rotation angles given as
`Rational` numbers in units of pi (so `R Z (1/2)` is a pi/2 rotation).

QPP does have a `Program` type for sequences of unitary and measurement
steps (you will meet it in the teleportation exercise). For this sheet
we don't need it: a "quantum program" is just a Haskell `do` block that
chains `apply` and `measureProjection` calls in order -- exactly like
the No-Cloning example you were given. Running `main` *is* running the
program.

1. Complete `bellProgram`, which runs a 2-step "program" (H on qubit 0,
   then CNOT) starting from |00>, using `Compose` to build a single QOp
   for the whole 2-qubit circuit before calling `evalOp`/`apply` once.

2. Complete `runWithMidCircuitMeasurement`, a 3-step "program": prepare
   the Bell state, measure qubit 0, then measure qubit 1, printing the
   state after each step (paper Exercise 3).

3. Complete `directSumCNOT`, building CNOT using `DirectSum I X` instead
   of `C X`, and check in `main` that it agrees with `C X` on all four
   computational basis states (paper Exercise 4b).

4. Complete `myUnitary`, an arbitrary single-qubit unitary built as a
   `Compose` chain of `H` and `R Z theta` QOps (recall H conjugates
   Z-rotations into X-rotations), and check what it does to |0> (paper
   Exercise 5a).
-}

bellProgram :: StateT
bellProgram =
    let psi0    = ket [0, 0]
        circuit = Compose (Tensor H I) (C X)   -- TODO: check the argument order/semantics of Compose
        op      = evalOp circuit
    in psi0 -- TODO: actually apply `op` to `psi0`

runWithMidCircuitMeasurement :: IO ()
runWithMidCircuitMeasurement = do
    let bell = bellProgram
        mP   = measureProjection
        p00  = mP 2 0 0
        p01  = mP 2 0 1
    putStrLn $ "Step 1 (after Bell prep): " ++ showState bell
    -- TODO: apply p00 to `bell`, normalize, print as "Step 2 (measure qubit 0 -> 0): ..."
    -- TODO: then measure qubit 1 of that result and print "Step 3: ..."
    return ()

directSumCNOT :: QOp
directSumCNOT = C X -- TODO: replace with `DirectSum I X`

-- | theta is a Rational in units of pi, e.g. 1/2 for a pi/2 rotation.
myUnitary :: Rational -> QOp
myUnitary theta = H -- TODO: replace with `Compose H (R Z theta)` (or similar) so it's a genuine 1-qubit unitary


main :: IO ()
main = do
    putStrLn "-- Exercise 1: Bell state as a single composed QOp --"
    putStrLn $ "bellProgram = " ++ showState bellProgram

    putStrLn "\n-- Exercise 3: program with mid-circuit measurement --"
    runWithMidCircuitMeasurement

    putStrLn "\n-- Exercise 4: DirectSum vs. C as ways of building CNOT --"
    let cnotOp = evalOp (C X)
        dsOp   = evalOp directSumCNOT
        basis  = [ket [0,0], ket [0,1], ket [1,0], ket [1,1]]

    mapM_ (\b -> putStrLn $
                    "input " ++ showState b ++
                    "  C X -> "        ++ showState (apply cnotOp b) ++
                    "  DirectSum -> "  ++ showState (apply dsOp b))
          basis

    putStrLn "\n-- Exercise 5: an arbitrary single-qubit unitary from Compose --"
    let u = evalOp (myUnitary (1/2))
    putStrLn $ "myUnitary(pi/2) |0> = " ++ showState (apply u (ket [0]))

    -- TODO (stretch, paper Exercise 5): using only H, R Z theta (for any
    -- theta) and C X, build (as a single QOp via Tensor/Compose/C) a
    -- 3-qubit circuit that prepares the GHZ state
    -- 1/sqrt(2) (|000> + |111>) from |000>, and check it with showState.
