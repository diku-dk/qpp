module Main where

import QPP
import QPP.Semantics.MPS

import Data.Complex

{-| SOLUTIONS: Single qubits -- normalization, change of basis, and global
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
normSquared psi = inner psi psi -- returns the vector squared. SOLVED

-- | Coordinates (b0, b1) of |psi> in the {|+>, |->} basis.
basisChange :: StateT -> (ComplexT, ComplexT)
basisChange psi =
    let plusOp  = evalOp H
        plus    = apply plusOp (ket [0])   -- This is really |+>
        minus   = apply plusOp (ket [1])   -- replaced with the real |->
        b0      = inner plus psi           -- projected psi onto |+>
        b1      = inner minus psi          -- projected psi onto |->
    in (b0, b1)

-- | Multiply a state by an overall phase e^{i*alpha}.
applyGlobalPhase :: Double -> StateT -> StateT
applyGlobalPhase alpha psi = (exp (0 :+ alpha)) .* psi -- Filled in, using (.*) and a complex phase factor


main :: IO ()
main = do
    let psiA = ket [0]
    let psiB = mkState (sqrt (1/3)) (sqrt (2/3))
    let psiC = mkState (sqrt (1/5)) (sqrt (4/5))  -- from paper Exercise 1
    let bad  = mkState 1 1                        -- NOT normalized: sanity check

    putStrLn "-- Exercise 1: normalization --"
    mapM_ (\(name, psi) ->
              putStrLn ("<" ++ name ++ "|" ++ name ++ "> = " ++ show (normSquared psi)))
          [("psiA", psiA), ("psiB", psiB), ("psiC", psiC), ("bad", bad)]

    putStrLn "\n-- Exercise 2: change of basis --"
    let (b0, b1) = basisChange psiB
    putStrLn $ "psiB in {|+>,|->} basis: b0 = " ++ show b0 ++ ", b1 = " ++ show b1
    --putStrLn $ "Check: |b0|^2 + |b1|^2 should equal 1 -- filled in a check here"
    let lenB = conjugate b0 * b0 + conjugate b1 * b1
    putStrLn $ "Check: |b0|^2 + |b1|^2 = " ++ show lenB

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
    putStrLn "These should match exactly, even though psiB /= psiB' as vectors. ... And they do"

    putStrLn "\n-- Stretch: the same check in the {|+>,|->} basis --"
    let (c0,  c1)  = basisChange psiB
        (c0', c1') = basisChange psiB'
    putStrLn $ "psiB  : |b0|^2 = " ++ show (magnitude c0  ** 2)
                ++ ", |b1|^2 = " ++ show (magnitude c1  ** 2)
    putStrLn $ "psiB' : |b0|^2 = " ++ show (magnitude c0' ** 2)
                ++ ", |b1|^2 = " ++ show (magnitude c1' ** 2)
    putStrLn "Identical again. A global phase is invisible in EVERY basis, not just"
    putStrLn "the computational one: it factors out of <b|psi> as a scalar, and the"
    putStrLn "modulus throws it away. Only *relative* phases between branches survive."
