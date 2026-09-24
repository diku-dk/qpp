-- Exercise scaffold: the placeholder bindings below are for the student to use.
{-# OPTIONS_GHC -Wno-unused-local-binds -Wno-unused-matches #-}
module Main where

import QPP
import QPP.Semantics.MPS

{-| Exercise: Classical Simulability and the Limits of Quantum Advantage.

This exercise accompanies the "Classical Simulability and the Limits of
Quantum Advantage" problem sheet (Exercise 1, 2, 3).

Key trick used throughout: instead of reaching for the `dagger` function,
note that H, X, Y, Z, and CNOT (= C X) are all self-inverse, and
`R op theta` is inverted by `R op (-theta)`. So wherever you'd
conceptually write U^dagger, you can write it out explicitly using this
fact. Remember that theta is a `Rational` in units of pi, so the T gate
is `R Z (1/4)` and the S (phase) gate is `R Z (1/2)`.

1. Complete `conjugate`, which computes (informally) "U P U^dagger
   applied to psi" by applying U^dagger, then P, then U to a test state
   -- and use it to verify a couple of the Clifford conjugation
   identities from paper Exercise 1 by comparing against the *claimed*
   conjugated Pauli, using your `sameState` helper from the Pauli/Bloch
   exercise sheet (feel free to copy it in here).

2. Complete `bellViaConjugation`, which tracks a Bell-state preparation
   by only ever applying single-qubit conjugation rules on paper/in
   comments (paper Exercise 2) -- then, separately, `bellViaStatevector`
   actually runs the circuit on the full state vector so you can sanity
   check your by-hand stabilizer bookkeeping against the ground truth.

3. Complete `checkTBreaksPauli`, which computes T X T^dagger (using the
   T = R Z (1/4) trick above) applied to test states, and checks
   whether the result matches ANY of I, X, Y, Z (up to phase) applied to
   those same test states -- demonstrating that T's conjugation action
   is not a Pauli operator (paper Exercise 3).

4. Complete `sizeComparison`, a plain (non-quantum) Haskell function
   printing, for a range of n, the size 2^n of a full amplitude vector
   versus the size n*n of an n-qubit stabilizer tableau -- a simple,
   very concrete way to see why Clifford circuits are cheap to simulate
   classically (paper Exercise 2c / 4a).
-}

sameState :: StateT -> StateT -> Bool
sameState phi psi = False -- TODO: as in the Pauli/Bloch exercise sheet:
                           -- compare |<phi|psi>| to 1

-- | Apply a QOp (thought of as a "Pauli" P) sandwiched between U^dagger
-- and U, to a test state: computes U P U^dagger |psi>.
-- `uInverse` should be the QOp you've worked out is U's inverse.
conjugate :: QOp -> QOp -> QOp -> StateT -> StateT
conjugate u uInverse p psi = psi -- TODO: apply uInverse, then p, then u, in that order

checkConjugationRule :: String -> QOp -> QOp -> QOp -> QOp -> [StateT] -> IO ()
checkConjugationRule name u uInverse p claimedResult testStates =
    mapM_ (\psi ->
              let lhs = conjugate u uInverse p psi
                  rhs = apply (evalOp claimedResult) psi
              in putStrLn (name ++ " on " ++ showState psi ++
                           ": matches claim? " ++ show (sameState lhs rhs)))
          testStates

-- | Track the Bell-state circuit by direct application to |00>.
bellViaStatevector :: StateT
bellViaStatevector =
    let h0 = evalOp (Tensor H I)
        cx = evalOp (C X)
    in apply cx (apply h0 (ket [0,0]))

checkTBreaksPauli :: [StateT] -> IO ()
checkTBreaksPauli testStates = do
    let t    = R Z (1/4)
        tinv = R Z (negate (1/4))
        candidates = [("I", I), ("X", X), ("Y", Y), ("Z", Z)]
    mapM_ (\psi -> do
              let lhs = conjugate t tinv X psi   -- T X T^dagger |psi>
              putStrLn $ "T X T^dagger on " ++ showState psi ++ " = " ++ showState lhs
              mapM_ (\(name, p) ->
                        let rhs = apply (evalOp p) psi
                        in putStrLn ("  matches " ++ name ++ "? " ++ show (sameState lhs rhs)))
                    candidates)
          testStates

sizeComparison :: [Int] -> IO ()
sizeComparison ns =
    mapM_ (\n -> putStrLn $
                    "n = " ++ show n ++
                    "   full statevector size 2^n = " ++ show ((2 :: Integer) ^ n) ++
                    "   stabilizer tableau size ~ n^2 = " ++ show (n * n))
          ns


main :: IO ()
main = do
    putStrLn "-- Exercise 1: checking Clifford conjugation rules --"
    let testStates = [ket [0], ket [1], apply (evalOp H) (ket [0])]
    -- TODO: fill in `conjugate`'s three applications above, then:
    checkConjugationRule "H X H =?= Z" H H X Z testStates
    checkConjugationRule "H Z H =?= X" H H Z X testStates
    -- TODO: also check S X S^dagger =?= Y, using S = R Z (1/2), with
    -- the appropriate inverse QOp.

    putStrLn "\n-- Exercise 2: Bell state, statevector ground truth --"
    putStrLn $ "bellViaStatevector = " ++ showState bellViaStatevector
    putStrLn "Compare this to the {X (x) X, Z (x) Z}-stabilized state you derived by hand."

    putStrLn "\n-- Exercise 3: T breaks the Pauli tableau --"
    checkTBreaksPauli testStates

    putStrLn "\n-- Exercise 4: statevector vs. stabilizer-tableau size --"
    sizeComparison [1, 2, 4, 8, 16, 32, 64]
