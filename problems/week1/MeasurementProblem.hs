-- Exercise scaffold: the placeholder bindings below are for the student to use.
{-# OPTIONS_GHC -Wno-unused-local-binds -Wno-unused-matches #-}
module Main where

import QPP
import QPP.Semantics.MPS

{-| Exercise: Measurement.

This exercise accompanies the "Measurement" problem sheet
(Exercise 1, 2, 4).

Note on types: `measureProjection` hands you an already-evaluated
backend operator (`OpT`), not a `QOp` syntax tree, so you feed it
straight to `apply` -- no `evalOp` in between.

1. Complete `prob`, computing the probability of a given measurement
   outcome as <psi| P |psi> for a (normalized) state psi and a projector
   P built with `measureProjection`.

2. In `main`, reproduce paper Exercise 1+2: measure
   psi = sqrt(1/5)|0> + sqrt(4/5)|1> both in the computational basis and
   in the {|+>,|->} basis (apply H first, then measure computationally),
   and print all four probabilities.

3. Complete `postMeasurementState`, which applies a projector and then
   *normalizes* the result -- i.e. the actual physical post-measurement
   state (not the unnormalized "branch amplitude" used for probability
   computations).

4. In `main`, reproduce the entangled-state partial-measurement analysis
   from paper Exercise 4, but for the state
     phi = 1/2 (|00> + |01> + |10> - |11>)
   Measure qubit 0 in the computational basis; for each outcome, print
   the (normalized) post-measurement 2-qubit state, and comment (as a
   putStrLn) on whether qubit 1 ends up in a definite classical value or
   in a superposition.

5. Complete `biasedCoin`, preparing
     psi(p) = sqrt(p)|0> + sqrt(1-p)|1>
   for p in [0,1], and check numerically that measuring it in the
   computational basis gives P(0) = p (paper Exercise 5).
-}

-- | Probability of a measurement outcome, given the corresponding
-- projector and a NORMALIZED state.
prob :: OpT -> StateT -> ComplexT
prob projOp psi = 0 -- TODO: use `apply` and `inner`

-- | The actual, renormalized state of the system after a projective
-- measurement (as opposed to the unnormalized "branch amplitude").
postMeasurementState :: OpT -> StateT -> StateT
postMeasurementState projOp psi = psi -- TODO: apply the projector, then `normalize`

biasedCoin :: Double -> StateT
biasedCoin p = ket [0] -- TODO: build sqrt(p)|0> + sqrt(1-p)|1>


main :: IO ()
main = do
    putStrLn "-- Exercise 1+2: basis-dependent measurement statistics --"
    let psi = (sqrt (1/5) .* ket [0]) .+ (sqrt (4/5) .* ket [1])
        mP  = measureProjection
        p0  = mP 1 0 0
        p1  = mP 1 0 1

    putStrLn $ "psi = " ++ showState psi
    putStrLn $ "P(0) computational = " ++ show (prob p0 psi)
    putStrLn $ "P(1) computational = " ++ show (prob p1 psi)

    let psiPM = apply (evalOp H) psi   -- change of basis trick: measure H|psi> computationally
    putStrLn $ "P(+) i.e. computational-P(0) of H|psi> = " ++ show (prob p0 psiPM)
    putStrLn $ "P(-) i.e. computational-P(1) of H|psi> = " ++ show (prob p1 psiPM)

    putStrLn "\n-- Exercise 4: partial measurement of an entangled state --"
    let phi = (0.5 .* ket [0,0]) .+ (0.5 .* ket [0,1]) .+ (0.5 .* ket [1,0]) .+ ((-0.5) .* ket [1,1])
        q00 = mP 2 0 0
        q01 = mP 2 0 1

    putStrLn $ "phi = " ++ showState phi
    putStrLn $ "P(qubit0 = 0) = " ++ show (prob q00 phi)
    putStrLn $ "P(qubit0 = 1) = " ++ show (prob q01 phi)

    let phi0 = postMeasurementState q00 phi
        phi1 = postMeasurementState q01 phi
    putStrLn $ "post-measurement state (outcome 0) = " ++ showState phi0
    putStrLn $ "post-measurement state (outcome 1) = " ++ showState phi1
    putStrLn "TODO: is qubit 1 left in a classical bit value or a superposition, in each branch?"

    putStrLn "\n-- Exercise 5: biased coin --"
    mapM_ (\p ->
              let psiP = biasedCoin p
              in putStrLn ("p = " ++ show p ++ " -> P(0) = " ++ show (prob p0 psiP)))
          [0.0, 0.25, 0.5, 0.75, 1.0]
