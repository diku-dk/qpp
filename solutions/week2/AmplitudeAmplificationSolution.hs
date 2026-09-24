module Main where

import QPP
import QPP.Semantics.MPS

import Data.Complex (realPart)

{-| SOLUTIONS: the Output Problem and Amplitude Amplification.

This exercise accompanies "The Output Problem and Amplitude Amplification"
problem sheet (Exercise 2, 3, 4). Grover search for N = 4 (2 qubits, one
marked item |11>), checked against sin^2((2k+1) theta).

Remember that `>:` is left-to-right composition: `a >: b` means "a, then b".
The sheet's stub wrote `Compose oracle diffusion`, which is `oracle ∘
diffusion` and therefore applies the *diffusion* first -- one of the two
things you were asked to check.
-}

--------------------------------------------------------------------------
-- Exercise 2: the two reflections
--------------------------------------------------------------------------

-- | Phase-flips the marked item |11> and nothing else: diag(1,1,1,-1) = C Z.
oracle :: QOp
oracle = C Z

{-| Reflection about |00>, up to an overall sign.

`C Z` phase-flips |11>; conjugating by X on both wires moves that flip onto
|00>, giving diag(-1,1,1,1) = -(2|00><00| - I). The sign is a global phase
(it multiplies the whole operator, not one branch), so it is harmless -- but
you should know it is there rather than assume it away.
-}
reflectAbout00 :: QOp
reflectAbout00 = X ⊗ X >: C Z >: X ⊗ X

-- | The diffusion operator D = H^{ox2} (2|00><00| - I) H^{ox2}, inheriting the
--   same overall minus sign from `reflectAbout00`.
diffusion :: QOp
diffusion = H ⊗ H >: reflectAbout00 >: H ⊗ H

-- | One Grover iteration: oracle first, then diffusion.
groverIterate :: QOp
groverIterate = oracle >: diffusion

-- | Apply the Grover iterate k times. `evalOp` runs once, outside the loop.
groverRun :: Int -> StateT -> StateT
groverRun k psi = iterate (apply (evalOp groverIterate)) psi !! k

--------------------------------------------------------------------------
-- Measuring "did we find it?"
--------------------------------------------------------------------------

{-| Projector onto the marked item |11>.

`measureProjection n k v` only ever constrains a *single* qubit -- `mP 2 0 1`
is the projector onto |1*>, which covers |10> as well as |11>. To isolate
|11> you need both single-qubit projectors. They commute (they act on
different qubits), so the order does not matter.
-}
projectMarked :: StateT -> StateT
projectMarked = apply (measureProjection 2 0 1) . apply (measureProjection 2 1 1)

-- | P(marked) = <psi|P|psi>, real because P is a projector.
probMarked :: StateT -> Double
probMarked psi = realPart (inner psi (projectMarked psi))

main :: IO ()
main = do
    putStrLn "-- Setting up: uniform superposition over 2 qubits --"
    let s0 = apply (evalOp (H ⊗ H)) (ket [0,0])
    putStrLn $ "Uniform superposition: " ++ showState s0
    putStrLn $ "P(marked) before any Grover iterations = " ++ show (probMarked s0)
    putStrLn "(1/4, as it must be: one marked item out of N = 4.)"

    putStrLn "\n-- Exercise 3: sweep k = 0,1,2 and compare to theory --"
    let theta = pi / 6 :: Double   -- theta = arcsin(sqrt(1/4)) = pi/6
    mapM_ (\k ->
              let simulated = probMarked (groverRun k s0)
                  theory    = sin ((2 * fromIntegral k + 1) * theta) ** 2
              in putStrLn ("k = " ++ show k
                            ++ "   simulated P(marked) = " ++ show simulated
                            ++ "   sin^2((2k+1)theta) = " ++ show theory))
          [0, 1, 2 :: Int]

    putStrLn ""
    putStrLn "Simulation and theory agree. P(marked) is maximised at k = 1, where it"
    putStrLn "reaches 1 exactly: for N = 4 a single Grover iteration finds the marked"
    putStrLn "item with certainty. Going on to k = 2 *overshoots* -- the amplitude"
    putStrLn "rotates past the target and P(marked) falls back to 1/4. Grover is a"
    putStrLn "rotation, not a monotone improvement, so running it too long is as bad"
    putStrLn "as not running it long enough."

    putStrLn "\n-- Exercise 4: why quadratic is not enough --"
    putStrLn "Each iteration rotates the state by 2 theta in the 2-dimensional span of"
    putStrLn "|marked> and |unmarked>, with sin theta = sqrt(M/N). We need (2k+1) theta"
    putStrLn "= pi/2, so k = O(sqrt(N/M)) -- a quadratic saving over the O(N) classical"
    putStrLn "scan, and it is optimal: no quantum algorithm does better with a black-box"
    putStrLn "oracle (the BBBV bound)."
    putStrLn ""
    putStrLn "For an NP-complete problem, N = 2^n over n-bit certificates, so sqrt(N) ="
    putStrLn "2^(n/2): still exponential. Quadratic speedup halves the exponent, it does"
    putStrLn "not remove it -- it buys you a longer key length, not membership in BQP."
    putStrLn "And the speedup only applies to *unstructured* search: a real SAT solver"
    putStrLn "exploits structure the oracle model deliberately hides, so the sqrt(N)"
    putStrLn "bound is not even the right comparison against a good classical solver."
