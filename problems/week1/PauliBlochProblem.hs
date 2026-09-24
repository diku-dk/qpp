-- Exercise scaffold: the placeholder bindings below are for the student to use.
{-# OPTIONS_GHC -Wno-unused-local-binds -Wno-unused-matches #-}
module Main where

import QPP
import QPP.Semantics.MPS

{-| Exercise: Pauli matrices and the Bloch sphere.

This exercise accompanies the "Pauli Matrices and the Bloch Sphere"
problem sheet (Exercise 1, 2, 5, 6).

Note on angles: in QPP, `R axis theta` takes theta as a `Rational` in
units of pi, so `R Z (1/2)` is a rotation by pi/2 about the Z axis. Keep
rotation angles rational -- don't route them through `Double`.

1. Complete `sameState`, a helper that checks whether two 1-qubit states
   are equal up to global phase (needed because several checks below are
   only true "up to phase"). Hint: two states are equal up to phase iff
   |<phi|psi>| = 1.

2. Complete `checkPauliSquares`, which applies X, Y, Z twice in a row to
   a test state and checks the result equals the original state (up to
   phase). This checks X^2 = Y^2 = Z^2 = I numerically (paper Ex. 1a).

3. Complete `blochVector`, which computes the Bloch vector
   (<psi|X|psi>, <psi|Y|psi>, <psi|Z|psi>) of a 1-qubit state, using
   `evalOp`, `apply`, and `inner`.

4. In `main`, compute the Bloch vector of |0>, |1>, |+>, |-> and check
   they come out as (0,0,1), (0,0,-1), (1,0,0), (-1,0,0) -- confirming
   that orthogonal states land on *antipodal* Bloch-sphere points
   (paper Exercise 6).

5. Use the `R` operator (rotation about a Pauli axis) to rotate |+>
   about the Z-axis by a range of angles theta, and print the resulting
   Bloch vector for each. Confirm it traces a circle around the z-axis
   (paper Exercise 5b).
-}

-- | Are two 1-qubit states equal up to an overall (unobservable) phase?
sameState :: StateT -> StateT -> Bool
sameState phi psi = False -- TODO: compare |<phi|psi>| to 1 (within a small tolerance)

-- | Apply the same 1-qubit QOp twice and report whether we're back to the
-- original state (up to phase).
checkPauliSquares :: StateT -> IO ()
checkPauliSquares psi = do
    let opX  = evalOp X
        opY  = evalOp Y
        opZ  = evalOp Z
        xx   = apply opX (apply opX psi)  -- TODO: is this correct? (order matters!)
        yy   = psi                        -- TODO: replace with Y applied twice
        zz   = psi                        -- TODO: replace with Z applied twice
    putStrLn $ "X^2 |psi> == |psi>? " ++ show (sameState xx psi)
    putStrLn $ "Y^2 |psi> == |psi>? " ++ show (sameState yy psi)
    putStrLn $ "Z^2 |psi> == |psi>? " ++ show (sameState zz psi)

-- | (<psi|X|psi>, <psi|Y|psi>, <psi|Z|psi>)
blochVector :: StateT -> (ComplexT, ComplexT, ComplexT)
blochVector psi =
    let opX = evalOp X
        opY = evalOp Y
        opZ = evalOp Z
        rx  = 0 -- TODO: <psi|X|psi>
        ry  = 0 -- TODO: <psi|Y|psi>
        rz  = 0 -- TODO: <psi|Z|psi>
    in (rx, ry, rz)


main :: IO ()
main = do
    putStrLn "-- Exercise 1: Pauli matrices square to the identity --"
    checkPauliSquares (ket [0])
    checkPauliSquares (apply (evalOp H) (ket [0]))   -- try it on |+> too

    putStrLn "\n-- Exercise 3+6: Bloch vectors of the standard states --"
    let zero  = ket [0]
        one   = ket [1]
        plus  = apply (evalOp H) (ket [0])
        minus = apply (evalOp H) (ket [1])   -- TODO: check this really gives |->!

    mapM_ (\(name, psi) -> putStrLn (name ++ ": Bloch vector = " ++ show (blochVector psi)))
          [("|0>", zero), ("|1>", one), ("|+>", plus), ("|->", minus)]

    putStrLn "\nExpected (up to numerical rounding): (0,0,1), (0,0,-1), (1,0,0), (-1,0,0)"
    putStrLn "Note |0> and |1> are orthogonal but their Bloch vectors are ANTIPODAL, not perpendicular."

    putStrLn "\n-- Exercise 5: rotating the Bloch vector with R --"
    -- Angles are rationals in units of pi: 0, pi/4, pi/2, 3pi/4, pi.
    let angles = [0, 1/4, 1/2, 3/4, 1] :: [Rational]
    mapM_ (\theta ->
              let rz    = evalOp (R Z theta)   -- TODO: confirm this is the right QOp constructor/argument order
                  psi   = apply rz plus
                  bloch = blochVector psi
              in putStrLn ("theta = " ++ show theta ++ " pi -> Bloch vector = " ++ show bloch))
          angles

    -- TODO (stretch, paper Exercise 7): find delta, theta, and an axis n
    -- such that R n theta (up to a global phase e^{-i delta/2}) equals
    -- the Hadamard gate H, by comparing apply (evalOp (R ...)) (ket [0])
    -- and apply (evalOp (R ...)) (ket [1]) against apply (evalOp H) ...
    -- using `sameState`.
