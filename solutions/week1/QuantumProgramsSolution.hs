module Main where

import QPP
import QPP.Semantics.MPS

import Data.Complex (realPart)

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
        -- `∘` is right-to-left, so H ⊗ I acts first and C X second. The stub's
        -- `Compose (Tensor H I) (C X)` says the opposite and builds CX-then-H.
        circuit = (C X) ∘ (H ⊗ I)
        qop     = evalOp circuit
    in apply qop psi0


{-| Measure qubit k of a 2-qubit state in the computational basis, returning
    the probability of outcome v together with the renormalized
    post-measurement state. -}
measureQubit :: Int -> Int -> StateT -> (Double, StateT)
measureQubit k v psi =
    let branch = apply (measureProjection 2 k v) psi
    in  (realPart (inner psi branch), normalize branch)

{-| The three steps, with both branches shown at each measurement.

We enumerate the branches with explicit projectors rather than sampling one
of them, because the point here is to see the whole tree at once. If you
want an actual run with a coin flip instead, that is what `Program` /
`evalProg` and an `RNG` are for -- see the teleportation exercise. -}
runWithMidCircuitMeasurement :: IO ()
runWithMidCircuitMeasurement = do
    let bell = bellProgram
    putStrLn $ "Step 1 (after Bell prep): " ++ showState bell
    mapM_ (\v0 -> do
        let (p0, afterQ0) = measureQubit 0 v0 bell
        putStrLn $ "  Step 2 (qubit 0 -> " ++ show v0 ++ ", p = " ++ show p0 ++ "): "
                    ++ showState afterQ0
        mapM_ (\v1 -> do
            let (p1, afterQ1) = measureQubit 1 v1 afterQ0
            putStrLn $ "    Step 3 (qubit 1 -> " ++ show v1 ++ ", p = " ++ show p1 ++ "): "
                        ++ if p1 < 1e-12 then "impossible branch" else showState afterQ1)
          [0, 1])
      [0, 1]
    putStrLn "Each qubit-0 outcome has probability 1/2, but once it is known the"
    putStrLn "qubit-1 measurement is already decided: one branch has probability 1"
    putStrLn "and the other 0. That determinism *is* the entanglement."

directSumCNOT :: QOp
directSumCNOT = I ⊕ X -- = DirectSum I X: identity on the |0> block, X on the |1> block

-- | theta is a Rational in units of pi, e.g. 1/2 for a pi/2 rotation.
myUnitary :: Rational -> QOp
myUnitary theta = H ∘ R Z theta  -- i.e. `Compose H (R Z theta)`: a genuine 1-qubit unitary


-- | Exercise 5: GHZ = H on qubit 0, then CNOT 0->1, then CNOT 1->2.
--   `@>` and `<@` bind looser than `>:`, hence the parentheses.
ghz :: QOp
ghz = (H @> 2) >: (C X @> 1) >: (1 <@ C X)

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
    let u = evalOp (myUnitary (1 / 2))
    putStrLn $ "myUnitary(pi/2) |0> = " ++ showState (apply u (ket [0]))

    putStrLn "\n-- Exercise 5 (stretch): GHZ from H, R Z and C X --"
    putStrLn $ "ghz |000> = " ++ showState (apply (evalOp ghz) (ket [0,0,0]))
    putStrLn "The allowed set is {H, R Z theta, C X}, but no rotation is needed: H"
    putStrLn "and C X alone suffice. R Z earns its place when you want a GHZ-like"
    putStrLn "state with a relative phase, e.g. (|000> + e^(i pi theta)|111>)/sqrt 2,"
    putStrLn "which is `ghz` followed by R Z theta on any single qubit."
