-- Exercise scaffold: the placeholder bindings below are for the student to use.
{-# OPTIONS_GHC -Wno-unused-local-binds -Wno-unused-matches #-}
module Main where

import QPP
import QPP.Semantics.MPS

{-| Exercise: Tensor products of states and operators.

This exercise accompanies the "Tensor Products" problem sheet
(Exercise 1, 2, 4, 5).

1. Complete `productState`, building |psi1> ⊗ |psi2> using (⊗).

2. Complete `bellPhiPlus`, preparing the Bell state
   1/sqrt(2) (|00> + |11>) starting from |00>, using H and C X  and  the
   evalOp as well as the apply command.

3. Complete `applyOnQubit0` / `applyOnQubit1`, which apply a 1-qubit QOp
   to only the first / second qubit of a 2-qubit state, using the
   `Tensor` QOp constructor (paper Exercise 4a).

4. In `main`, verify that (X ⊗ I) and (I ⊗ X) act differently on the
   Bell state, and check whether the results are still entangled
   (informally: is the state still equal, up to relabelling, to a Bell
   state? use `showState`/`printS` to inspect).

5. In `main`, use an entanglement witness: measure qubit 0 of the Bell
   state in the computational basis, and separately in the {|+>,|->}
   basis (apply H first). Compare the conditional post-measurement state
   of qubit 1 in each case with what you would see for the *product*
   state |+> ⊗ |0> treated the same way. This is a hands-on version of
   paper Exercise 2 (Bell states are not product states).
-}

productState :: StateT -> StateT -> StateT
productState psi1 psi2 = psi1 -- TODO: fix, should combine both states

bellPhiPlus :: StateT
bellPhiPlus =
    let psi0 = ket [0, 0]        -- TODO: is |00> built like this, or via (⊗)? check both work
        h0   = evalOp H          -- TODO: this needs to act on qubit 0 only of a 2-qubit state -- see applyOnQubit0
        cx   = evalOp (C X)
    in psi0 -- TODO: replace with H on qubit 0, then CX

applyOnQubit0 :: QOp -> StateT -> StateT
applyOnQubit0 op psi = psi -- TODO: use `Tensor op I`, evalOp, apply

applyOnQubit1 :: QOp -> StateT -> StateT
applyOnQubit1 op psi = psi -- TODO: use `Tensor I op`, evalOp, apply


main :: IO ()
main = do
    putStrLn "-- Exercise 1: simple product states --"
    let plus  = apply (evalOp H) (ket [0])
        pp    = productState plus plus
    putStrLn $ "|+> ⊗ |+> = " ++ showState pp

    putStrLn "\n-- Exercise 2+5: Bell state vs. product state --"
    let bell = bellPhiPlus
    putStrLn $ "|Phi+> = " ++ showState bell

    let prod = productState plus (ket [0])   -- |+> (x) |0>, NOT entangled
    putStrLn $ "|+>|0> = " ++ showState prod

    let mP  = measureProjection
        p00 = mP 2 0 0   -- project qubit 0 to |0>
        p01 = mP 2 0 1   -- project qubit 0 to |1>

    putStrLn "\nMeasuring qubit 0 of |Phi+> in computational basis:"
    let bell0 = apply p00 bell
        bell1 = apply p01 bell
    putStrLn $ "  outcome 0 -> (unnormalized) " ++ showState bell0
    putStrLn $ "  outcome 1 -> (unnormalized) " ++ showState bell1
    putStrLn "  TODO: is qubit 1 now in a *definite* classical state in each branch?"

    putStrLn "\nMeasuring qubit 0 of |+>|0> in computational basis:"
    let prod0 = apply p00 prod
        prod1 = apply p01 prod
    putStrLn $ "  outcome 0 -> (unnormalized) " ++ showState prod0
    putStrLn $ "  outcome 1 -> (unnormalized) " ++ showState prod1
    putStrLn "  TODO: compare to the entangled case above -- what's the key difference?"

    putStrLn "\n-- Exercise 4: (X ⊗ I) vs (I ⊗ X) on the Bell state --"
    let bellX0 = applyOnQubit0 X bell
        bellX1 = applyOnQubit1 X bell
    putStrLn $ "(X ⊗ I)|Phi+> = " ++ showState bellX0
    putStrLn $ "(I ⊗ X)|Phi+> = " ++ showState bellX1
    putStrLn "TODO: are these the same state? Are they still entangled?"

    -- TODO (stretch, paper Exercise 6): build a 3-qubit GHZ state
    -- (1/sqrt 2)(|000> + |111>) using Tensor/Compose/C on 3-qubit QOps.
    -- Hint: you'll need `H ⊗ Id 2` (or `Tensor H (Id 2)`) and two CNOTs.
