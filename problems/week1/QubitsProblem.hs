-- Exercise scaffold: the placeholder bindings below are for the student to use.
{-# OPTIONS_GHC -Wno-unused-local-binds -Wno-unused-matches #-}
module Main where

import QPP
import QPP.Semantics.MPS

{-| Exercise: Single qubits -- normalization, change of basis, and global
phase invariance.

This exercise accompanies the "Qubits" problem sheet (Exercise 1, 2, 4).
Consult the document "QPP - What it does" if you get lost in code space :-)

1. Complete `normSquared`, which should compute <psi|psi> for a 1-qubit
   state. Use it to check which of the states in `main` are normalized.

2. Complete `basisChange`, which takes a 1-qubit state |psi> and returns
   its coordinates (b0, b1) in the {|+>, |-> } basis, i.e. such that
   |psi> = b0|+> + b1|-> . Build |+>, |-> using the Hadamard operator H
   applied to |0>, |1> (do not hard-code them by hand), then use `inner`
   to project.

3. Complete `applyGlobalPhase`, which multiplies a state by e^{i*alpha}
   for a real angle alpha, using (.*).

4. In `main`, use `applyGlobalPhase` to build psi' = e^{i*alpha} psi for
   some alpha of your choosing, and show numerically (via
   `measureProjection` and `inner`) that psi and psi' give identical
   measurement probabilities in the computational basis. This is the
   computational counterpart of paper Exercise 4 (global phase is
   unobservable).
-}

-- | a0|0> + a1|1>
mkState :: ComplexT -> ComplexT -> StateT
mkState a0 a1 = (a0 .* ket [0]) .+ (a1 .* ket [1])

-- | <psi|psi> -- should equal 1 for a properly normalized qubit state.
normSquared :: StateT -> ComplexT
normSquared psi = 0 -- TODO: replace with the correct computation using `inner` (See QPP/Syntax.hs for the definition)

-- | Coordinates (b0, b1) of |psi> in the {|+>, |->} basis.
basisChange :: StateT -> (ComplexT, ComplexT)
basisChange psi =
    let plusOp  = evalOp H
        plus    = apply plusOp (ket [0])   -- TODO: is this really |+>? check! (Use repl command to try)
        minus   = ket [1]                  -- TODO: replace with the real |->
        b0      = 0                        -- TODO: project psi onto |+>
        b1      = 0                        -- TODO: project psi onto |->
    in (b0, b1)

-- | Multiply a state by an overall phase e^{i*alpha}.
applyGlobalPhase :: Double -> StateT -> StateT
applyGlobalPhase alpha psi = psi -- TODO: fill in, using (.*) and a complex phase factor


main :: IO ()
main = do
    let psiA = ket [0]
    let psiB = mkState (sqrt (1/3)) (sqrt (2/3))
    let psiC = mkState (sqrt (1/5)) (sqrt (4/5))  -- from paper Exercise 1
    let bad  = mkState 1 1                        -- NOT normalized: sanity check

    putStrLn "-- Exercise 1: normalization --"
    mapM_ (\(name, psi) ->
              putStrLn (name ++ ": <psi|psi> = " ++ show (normSquared psi)))
          [("psiA", psiA), ("psiB", psiB), ("psiC", psiC), ("bad", bad)]

    putStrLn "\n-- Exercise 2: change of basis --"
    let (b0, b1) = basisChange psiB
    putStrLn $ "psiB in {|+>,|->} basis: b0 = " ++ show b0 ++ ", b1 = " ++ show b1
    putStrLn $ "Check: |b0|^2 + |b1|^2 should equal 1 -- fill in a check here"

    putStrLn "\n-- Exercise 4: global phase invariance --"
    let alpha   = pi / 3
        psiB'   = applyGlobalPhase alpha psiB
        mP      = measureProjection
        p0      = mP 1 0 0   -- projector: 1 qubit, qubit index 0, outcome 0
        p1      = mP 1 0 1

    putStrLn $ "psiB  = " ++ showState psiB
    putStrLn $ "psiB' = " ++ showState psiB'

    let prob0  = inner psiB  (apply p0 psiB)
        prob0' = inner psiB' (apply p0 psiB')
        prob1  = inner psiB  (apply p1 psiB)
        prob1' = inner psiB' (apply p1 psiB')

    putStrLn $ "P(0) for psiB  = " ++ show prob0  ++ "   P(0) for psiB' = " ++ show prob0'
    putStrLn $ "P(1) for psiB  = " ++ show prob1  ++ "   P(1) for psiB' = " ++ show prob1'
    putStrLn "These should match exactly, even though psiB /= psiB' as vectors."

    -- TODO (open-ended): repeat the phase-invariance check measuring in
    -- the {|+>,|->} basis instead, using your `basisChange` function.
