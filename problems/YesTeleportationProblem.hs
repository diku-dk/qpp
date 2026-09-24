-- Exercise scaffold: `teleport` and `rng0` stay unused until the commented-out
-- driver at the bottom of `main` is filled in.
{-# OPTIONS_GHC -Wno-unused-top-binds -Wno-unused-local-binds #-}
module Main where

import QPP
import QPP.Semantics.MPS
import System.Random(mkStdGen, randoms)
-- import Data.List(intercalate)

{-| Exercise: Implement the quantum teleportation protocol from N&C Section 1.3.7 (Fig. 1.11), and demonstrate that teleportation of the full quantum state, contrary to cloning, *is* possible.

Use the operators of the QOp language in src/QPP/Syntax.hs.

- How do you implement X^M2, Z^M1 without dynamic control? (Hint: A^0 = I for any operator A).

- Our simple language doesn't have qubit numbering. How can you implement the flow from measurement M1 to qubit 2? (Hint: op can be multi-qubit in our controlled operation 'C op')
-}
teleprog :: Program
teleprog = let
        op1 = Unitary $
            Id 3 -- Replace by your own program

        op2 = Unitary $
            Id 3 -- Replace by your own program
    in
        [op1, Measure [0,1], op2]

{-| 'teleport psi rng' runs the teleprog program on the input 1Q-state psi. rng is a
list of numbers [0..1] for random measurements. See below for how to make an infinite list of random numbers, or supply a fixed list for debugging.

returns: (final_state, outcomes, rng') where outcomes are the measurement results
(most recent first) and rng' is the remaining random number stream.
-}
teleport :: StateT -> RNG -> (StateT, Outcomes, RNG)
teleport psi rng = let
        bell     = ket [0,0] -- replace by actual definition, N&C Section 1.3.6
        psi_bell = psi ⊗ bell
    in
        evalProg teleprog psi_bell rng


main :: IO()
main = do
    let rng0 = randoms (mkStdGen 42) :: [Double]
    --let rng0 = [0,0,0] -- (always measure 0)
    -- let rng0 = [0,1,0] -- (First meausure 0, then measure 1)

    putStr $ "\nTeleport program:\n" ++ showProgram teleprog ++ "\n\n"

--    putStr $ "|ψ>   = "++(showState psi) ++ "\n" -- ". Running teleport program!\n"
    --putStr $ "|ψbb> = "++(showState psi_bell)++".\nRunning teleport program!\n"

    -- let (end_state,outcomes,_) = teleport psi rng0

    -- putStr $ "Measurement outcomes: " ++ (show outcomes) ++ "\n"
    -- putStr $ "Final 3-qubit state:\n" ++ (showState end_state) ++ "\n\n"
