module Main where

import QPP
import QPP.Semantics.MPS

import Data.Complex (magnitude)

{-| SOLUTIONS: Classical Simulability and the Limits of Quantum Advantage.

This exercise accompanies the "Classical Simulability and the Limits of
Quantum Advantage" problem sheet (Exercise 1, 2, 3).

H, X, Y, Z and CNOT are all self-inverse, and `R op theta` is inverted by
`R op (-theta)`, so every U^dagger below is written out explicitly. Angles
are `Rational` in units of pi, so S = R Z (1/2) and T = R Z (1/4).
-}

-- | Comparison threshold. The library's shared `tol` is 1e-12; we use a
--   looser one here because these are accumulated floating-point results,
--   and because `QPP.Semantics.MPS` exports an `EvalCfg` field also called
--   `tol`, which would make the name ambiguous.
eps :: Double
eps = 1e-9

-- | Equal up to an overall (unobservable) phase: |<phi|psi>| = 1.
--   Assumes both states are normalized.
sameState :: StateT -> StateT -> Bool
sameState phi psi = abs (1 - magnitude (inner phi psi)) < eps

-- | Equal on the nose, phase included. Needed for stabilizers, where the
--   whole point is the difference between eigenvalue +1 and eigenvalue -1.
equalState :: StateT -> StateT -> Bool
equalState phi psi = norm (phi .- psi) < eps

--------------------------------------------------------------------------
-- Exercise 1: Clifford conjugation rules
--------------------------------------------------------------------------

{-| U P U^dagger applied to a state, with `uInverse` the QOp you have worked
    out is U^dagger. Read the applications right to left, as in the formula:
    U^dagger acts first. -}
conjugate :: QOp -> QOp -> QOp -> StateT -> StateT
conjugate u uInverse p psi =
    apply (evalOp u) (apply (evalOp p) (apply (evalOp uInverse) psi))

checkConjugationRule :: String -> QOp -> QOp -> QOp -> QOp -> [StateT] -> IO ()
checkConjugationRule name u uInverse p claimedResult testStates =
    mapM_ (\psi ->
              let lhs = conjugate u uInverse p psi
                  rhs = apply (evalOp claimedResult) psi
              in putStrLn ("  " ++ name ++ " on " ++ showState psi ++
                           ": matches claim? " ++ show (sameState lhs rhs)))
          testStates

--------------------------------------------------------------------------
-- Exercise 2: the Bell state, tracked two ways
--------------------------------------------------------------------------

-- | The Bell circuit run forwards on the full state vector: the ground truth.
bellViaStatevector :: StateT
bellViaStatevector =
    apply (evalOp (C X)) (apply (evalOp (H ⊗ I)) (ket [0,0]))

{-| The same circuit in the Heisenberg picture. Instead of carrying 4
    amplitudes we carry 2 stabilizer generators and push them through the
    circuit with the conjugation rules from Exercise 1:

      |00> is stabilized by  Z(x)I  and  I(x)Z
      H on qubit 0:          H Z H = X, and I(x)Z is untouched
      CNOT (0 -> 1):         CNOT (X(x)I) CNOT = X(x)X
                             CNOT (I(x)Z) CNOT = Z(x)Z

    Two 2-qubit Pauli strings instead of four amplitudes -- and for n qubits,
    n strings of 2n+1 bits instead of 2^n amplitudes. That is Gottesman-Knill.
-}
bellViaConjugation :: [(String, [QOp])]
bellViaConjugation =
    [ ("|00>",                [Z ⊗ I, I ⊗ Z])
    , ("after H on qubit 0",  [X ⊗ I, I ⊗ Z])
    , ("after CNOT (0 -> 1)", [X ⊗ X, Z ⊗ Z])
    ]

-- | Does S fix psi (eigenvalue +1, not just +-1)?
stabilizes :: QOp -> StateT -> Bool
stabilizes s psi = equalState (apply (evalOp s) psi) psi

--------------------------------------------------------------------------
-- Exercise 3: T is not a Clifford gate
--------------------------------------------------------------------------

checkTBreaksPauli :: [StateT] -> IO ()
checkTBreaksPauli testStates = do
    let t    = R Z (1/4)
        tinv = R Z (negate (1/4))
        candidates = [("I", I), ("X", X), ("Y", Y), ("Z", Z)]
    mapM_ (\psi -> do
              let lhs = conjugate t tinv X psi   -- T X T^dagger |psi>
              putStrLn $ "  T X T^dagger on " ++ showState psi ++ " = " ++ showState lhs
              putStrLn $ "    matches: " ++ unwords
                  [ name ++ "=" ++ show (sameState lhs (apply (evalOp p) psi))
                  | (name, p) <- candidates ])
          testStates

--------------------------------------------------------------------------
-- Exercise 4: how big is the classical description?
--------------------------------------------------------------------------

sizeComparison :: [Int] -> IO ()
sizeComparison ns =
    mapM_ (\n -> putStrLn $
                    "  n = " ++ pad 3 (show n) ++
                    "  statevector 2^n = " ++ pad 22 (show ((2 :: Integer) ^ n)) ++
                    "  tableau ~ n^2 = " ++ show (n * n))
          ns
  where pad w s = s ++ replicate (w - length s) ' '

main :: IO ()
main = do
    let zero  = ket [0]
        one   = ket [1]
        plus  = apply (evalOp H) (ket [0])
        testStates = [zero, one, plus]

    putStrLn "-- Exercise 1: checking Clifford conjugation rules --"
    checkConjugationRule "H X H =?= Z" H H X Z testStates
    checkConjugationRule "H Z H =?= X" H H Z X testStates
    -- S = R Z (1/2) is its own story: R Z is a rotation, so it carries a
    -- global phase relative to the textbook diag(1, i). Conjugation cancels
    -- that phase (it appears once as U and once as U^dagger), so the rule
    -- comes out clean anyway.
    checkConjugationRule "S X S^dagger =?= Y" (R Z (1/2)) (R Z (-1/2)) X Y testStates
    checkConjugationRule "S Z S^dagger =?= Z" (R Z (1/2)) (R Z (-1/2)) Z Z testStates

    putStrLn "\n-- Exercise 2: Bell state, both pictures --"
    putStrLn $ "  statevector: " ++ showState bellViaStatevector
    mapM_ (\(step, gens) -> putStrLn ("  " ++ step ++ ": stabilized by "
                                       ++ unwords (map showOp gens)))
          bellViaConjugation
    let finalGens = snd (last bellViaConjugation)
    putStrLn $ "  do those generators actually fix the state? "
                ++ show (map (`stabilizes` bellViaStatevector) finalGens)
    putStrLn "  Note we check with `equalState`, not `sameState`: a stabilizer must"
    putStrLn "  have eigenvalue +1, and -X(x)X would pass an up-to-phase test."

    putStrLn "\n-- Exercise 3: T breaks the Pauli tableau --"
    checkTBreaksPauli testStates
    putStrLn ""
    putStrLn "  T X T^dagger = (X + Y)/sqrt 2, which is not a Pauli. But look at the"
    putStrLn "  rows: on |0> and |1> it reports a match with both X and Y. That is not"
    putStrLn "  a bug -- X|0>, Y|0> and (X+Y)/sqrt2 |0> are all |1> up to a phase, so a"
    putStrLn "  basis state simply cannot tell Paulis apart. |+> is what settles it:"
    putStrLn "  there the answer is a genuine superposition of |+> and |->, matching"
    putStrLn "  none of I, X, Y, Z. Choose test states that can see the difference."

    putStrLn "\n-- Exercise 4: statevector vs. stabilizer-tableau size --"
    sizeComparison [1, 2, 4, 8, 16, 32, 64]
    putStrLn ""
    putStrLn "  A Clifford circuit on 64 qubits needs a few thousand bits of tableau;"
    putStrLn "  the amplitude vector would need 2^64 complex numbers. Clifford circuits"
    putStrLn "  can entangle heavily and still be simulated in polynomial time, so"
    putStrLn "  entanglement alone is not what makes quantum computing hard to fake."
    putStrLn "  The T gate is what leaves the tableau: it is the expensive resource in"
    putStrLn "  fault-tolerant architectures for exactly this reason."
