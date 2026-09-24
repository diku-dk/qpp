-- Exercise scaffold: the placeholder bindings below are for the student to use.
{-# OPTIONS_GHC -Wno-unused-local-binds -Wno-unused-matches #-}
module Main where

import QPP
import QPP.Semantics.MPS

{-| Exercise: the Output Problem and Amplitude Amplification.

This exercise accompanies the "The Output Problem and Amplitude
Amplification" problem sheet (Exercise 2, 3, 4). We build a full Grover
search for N=4 (2 qubits, 1 marked item |11>), and check the amplified
success probability against the sin^2((2k+1) theta) formula you derived
on paper.

Recipe for N=4 with the marked item |11>:
  * The "oracle" phase-flips |11> and leaves everything else alone. This
    is exactly `C Z` (control qubit 1, target qubit 2): when qubit 1 = 1
    it applies Z to qubit 2, which only flips the sign of the |1>
    component, i.e. exactly the |11> branch.
  * The "diffusion" operator reflects about the uniform superposition:
    D = H^{ox2} (2|00><00| - I) H^{ox2}, and (2|00><00| - I) itself can
    be built as (X^{ox2}) (C Z) (X^{ox2}) (flip everything, phase-flip
    |11>, flip back -- which phase-flips everything EXCEPT |00>, i.e.
    exactly 2|00><00| - I up to an overall sign you should track).

Note: `Id 2` below is the 2-qubit identity -- the placeholders have to
have the right arity or `apply` will reject them.
-}

oracle :: QOp
oracle = C Z   -- phase-flips |11>; nothing to fill in here, it's a freebie

diffusion :: QOp
diffusion = Id 2 -- TODO: build (H (x) H) . (X (x) X) . (C Z) . (X (x) X) . (H (x) H)
                 -- using Tensor and Compose. Careful with the Compose order
                 -- convention (see the QFT exercise sheet) and with an
                 -- overall sign -- an unobservable global phase is fine,
                 -- but make sure it's *global* and not relative.

groverIterate :: QOp
groverIterate = Compose oracle diffusion -- TODO: check this is the right order (oracle then diffusion, or vice versa?)

-- | Apply the Grover iterate k times to a starting state.
groverRun :: Int -> StateT -> StateT
groverRun 0 psi = psi
groverRun k psi = psi -- TODO: apply `evalOp groverIterate` k times, recursively


main :: IO ()
main = do
    putStrLn "-- Setting up: uniform superposition over 2 qubits --"
    let h2   = evalOp (Tensor H H)
        s0   = apply h2 (ket [0,0])
        mP   = measureProjection
        pMarked = mP 2 0 1  -- TODO: is this really the projector onto |*1> or |11>?
                            -- You may need to combine two single-qubit
                            -- projectors (one per qubit) to isolate |11>
                            -- specifically -- check measureProjection's
                            -- semantics and adjust if needed.

    putStrLn $ "Uniform superposition: " ++ showState s0
    putStrLn $ "P(marked) before any Grover iterations = " ++ show (inner s0 (apply pMarked s0))

    putStrLn "\n-- Exercise 3: sweep k = 0,1,2 and compare to theory --"
    let theta = pi / 6 :: Double  -- theta = arcsin(sqrt(1/4)) = pi/6, from paper Exercise 3a
    mapM_ (\k ->
              let psiK       = groverRun k s0
                  probK      = inner psiK (apply pMarked psiK)
                  theoryProb = (sin ((2 * fromIntegral k + 1) * theta)) ^ (2::Int)
              in putStrLn ("k = " ++ show k ++
                            "  simulated P(marked) = " ++ show probK ++
                            "  theory sin^2((2k+1)theta) = " ++ show theoryProb))
          [0, 1, 2 :: Int]

    putStrLn "\nDo the simulated and theoretical probabilities agree? At which k is P(marked) maximized?"

    -- TODO (stretch, paper Exercise 4): discuss in a comment why this
    -- sqrt(N) - style scaling (here: it took O(1) iterations for N=4)
    -- is a quadratic rather than exponential speedup over classical
    -- search, and why that's not enough to put NP-complete problems in
    -- BQP.
