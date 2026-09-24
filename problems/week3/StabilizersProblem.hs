-- Exercise scaffold: the placeholder bindings below are for the student to use.
{-# OPTIONS_GHC -Wno-unused-local-binds -Wno-unused-matches #-}
module Main where

import QPP
import QPP.Semantics.MPS

{-| Exercise: Stabilizers, check matrices, and the Clifford group.

This exercise accompanies the "Stabilizers and the Clifford Group" problem
sheet. It has two independent halves:

  PART A (Exercises 6-9 on paper): pure classical bit-vector arithmetic on
  stabilizer tableaux. This does NOT use QPP -- it's plain Haskell, since
  the check-matrix formalism is itself a classical data structure. Fill in
  `pauliMultiply`, `symplecticProduct`, `commutes`, and `rankF2`.

  PART B (Exercises 10-11 on paper): verifying Clifford conjugation
  identities (H X H^dagger = Z, etc.) numerically against the QPP
  operators, the same way QubitsProblem.hs checked global-phase invariance
  by comparing measurement statistics rather than raw matrices.

Consult "QPP - What it does" if you get lost in code space :-)

Two things to know about the QPP operator language before you start:
there are no separate `S` and `T` constructors -- they are the rotations
`R Z (1/2)` and `R Z (1/4)`, since rotation angles are `Rational` in
units of pi -- and `C op` controls a whole (possibly multi-qubit) `op`
on the qubit above it, so there are no explicit qubit indices to pass.
-}

--------------------------------------------------------------------------
-- PART A: stabilizer tableaux as bit vectors (no QPP needed)
--------------------------------------------------------------------------

-- | A Pauli group element on n qubits, encoded as (s, c, xs, zs):
--   the operator is (-1)^s * i^c * X^(xs!!0) Z^(zs!!0) (x) ... (x) X^(xs!!(n-1)) Z^(zs!!(n-1))
--   Bits are represented as Bool (False = 0, True = 1).
data PauliCode = PauliCode
  { signBit  :: Bool   -- ^ s
  , phaseBit :: Bool   -- ^ c
  , xBits    :: [Bool] -- ^ x_0 .. x_{n-1}
  , zBits    :: [Bool] -- ^ z_0 .. z_{n-1}
  } deriving (Eq, Show)

xor :: Bool -> Bool -> Bool
xor = (/=)

-- | Exercise 6: encode a single-qubit Pauli symbol into its (x,z) pair.
--   [I] = (0,0), [Z] = (0,1), [X] = (1,0), [Y] = (1,1).
encodeSymbol :: Char -> (Bool, Bool)
encodeSymbol 'I' = (False, False)
encodeSymbol 'X' = (True,  False)
encodeSymbol 'Z' = (False, True)
encodeSymbol 'Y' = (True,  True)   -- TODO: double check against Y = iXZ convention used in lecture
encodeSymbol c   = error ("encodeSymbol: not a Pauli symbol: " ++ [c])

-- | Exercise 6: pretty-print a PauliCode back as a signed tensor-product string,
--   e.g. "-i * Y (x) I (x) Z". Used to check your by-hand decoding.
decodePauli :: PauliCode -> String
decodePauli (PauliCode s c xs zs) =
  signStr ++ phaseStr ++ concatWith " (x) " (map symbolOf (zip xs zs))
  where
    signStr  = if s then "-" else ""
    phaseStr = if c then "i * " else ""
    symbolOf (False, False) = "I"
    symbolOf (True,  False) = "X"
    symbolOf (False, True)  = "Z"
    symbolOf (True,  True)  = "Y"
    concatWith sep = foldr1 (\a b -> a ++ sep ++ b)

-- | Exercise 7: multiply two same-length PauliCodes using the XOR rule
--     [P1 P2] = [P1] xor [P2] xor (sign correction),
--     sign correction bit = (z1 . x2) mod 2   (dot product over F2)
--   TODO: fill in. Don't forget the phase bits (c1,c2) also add mod 4 in
--   general -- for this exercise you may assume c1 = c2 = 0 (i.e. inputs
--   are +-1 times a tensor product of X's and Z's with no lone i factor)
--   and only worry about the sign bit s, as in the lecture's worked example.
pauliMultiply :: PauliCode -> PauliCode -> PauliCode
pauliMultiply (PauliCode s1 c1 x1 z1) (PauliCode s2 c2 x2 z2) =
  PauliCode s1 c1 x1 z1  -- TODO: replace with the real product; this is a placeholder

-- | Exercise 8: the symplectic product of two same-length bit-vector pairs,
--   ignoring sign/phase: <(x1,z1),(x2,z2)> = x1.z2 + z1.x2  (mod 2).
--   Returns True if the product is 1 (i.e. the operators ANTI-commute).
symplecticProduct :: [Bool] -> [Bool] -> [Bool] -> [Bool] -> Bool
symplecticProduct x1 z1 x2 z2 = False
  -- TODO: replace with the real computation, e.g.
  --   foldr xor False (zipWith (&&) x1 z2) `xor` foldr xor False (zipWith (&&) z1 x2)

-- | Exercise 8: do two Pauli strings commute? (Ignores overall phase, which
--   never affects commutativity.)
commutes :: PauliCode -> PauliCode -> Bool
commutes p q = not (symplecticProduct (xBits p) (zBits p) (xBits q) (zBits q))

-- | Exercise 9: rank of a list of (x|z) bit-vectors over F2, via Gaussian
--   elimination mod 2. A set of n generators on n qubits is independent
--   (and hence pins down a genuine stabilizer state, cf. Theorem 1 in the
--   lecture) iff this rank equals n.
--   TODO: implement Gaussian elimination over F2. A reasonable approach:
--     1. Represent each row as [Bool] (concatenation of xBits ++ zBits).
--     2. Repeatedly pick a pivot column, XOR it out of all other rows that
--        have a 1 there, and recurse on the remaining rows/columns.
rankF2 :: [[Bool]] -> Int
rankF2 rows = length rows  -- TODO: WRONG on purpose -- replace with real rank computation

--------------------------------------------------------------------------
-- PART B: verifying Clifford conjugation with QPP
--------------------------------------------------------------------------

-- | Apply U, then P, then U^dagger to a state, i.e. compute (U P U^dagger) |phi>.
--   `dagger` (from QPP.Syntax) pushes the adjoint down into the operator tree.
--   TODO: does `Compose` associate left-to-right or right-to-left? Double
--   check against a known identity (e.g. H ∘ H = I) in the REPL before
--   trusting the order below.
conjugateApply :: QOp -> QOp -> StateT -> StateT
conjugateApply u p phi =
  apply (evalOp u) (apply (evalOp p) (apply (evalOp (dagger u)) phi))
  -- TODO: this composes U . P . U^dagger acting *rightmost-first* on phi,
  -- i.e. computes U (P (U^dagger |phi>)) -- confirm this is what you want
  -- (it should equal (U P U^dagger) |phi>, which is NOT the same operator
  -- as U^dagger P U -- watch the order!).

-- | Exercise 10/11: numerically check H X H^dagger = Z by comparing
--   (H X H^dagger)|phi> against Z|phi> for a few test states |phi>, the
--   same trick QubitsProblem.hs used to check global-phase invariance via
--   measurement statistics rather than raw matrix equality.
checkHXH :: [StateT] -> [(StateT, StateT)]
checkHXH testStates =
  [ (conjugateApply H X phi, apply (evalOp Z) phi) | phi <- testStates ]
  -- TODO: after computing this list, use `inner` to check each pair is
  -- equal up to global phase (recall Exercise 4 of the Qubits sheet!).

main :: IO ()
main = do
    putStrLn "-- Part A: tableau arithmetic (fill in pauliMultiply, symplecticProduct, rankF2) --"

    let p1 = PauliCode False False [True, False, True] [False, True, True]   -- (110|011) placeholder, adjust to match paper Ex.7
    let p2 = PauliCode False False [True, False, True] [True, True, False]

    putStrLn $ "P1 = " ++ decodePauli p1
    putStrLn $ "P2 = " ++ decodePauli p2
    putStrLn $ "P1 P2 (via XOR rule) = " ++ decodePauli (pauliMultiply p1 p2)
    putStrLn $ "Do P1, P2 commute?    = " ++ show (commutes p1 p2)

    putStrLn "\n-- Exercise 9: independence check for the paper's P1,P2,P3 --"
    let rows =
          [ [True,True,False, False,False,True]   -- P1 = (110|001)
          , [False,True,True, True,False,False]   -- P2 = (011|100)
          , [True,False,True, True,False,True]    -- P3 = (101|101)
          ]
    putStrLn $ "rank = " ++ show (rankF2 rows) ++ " (should be 3 for independence, on 3 qubits)"

    putStrLn "\n-- Part B: Clifford conjugation checks --"
    let testStates = [ket [0], ket [1], evalOp H `apply` ket [0]]
    mapM_ (\(lhs, rhs) ->
              putStrLn $ "  H X H^dagger |phi> = " ++ showState lhs
                      ++ "   vs.   Z |phi> = " ++ showState rhs)
          (checkHXH testStates)

    -- TODO (open-ended, Exercise 12): implement a small "push T through
    -- Clifford" routine that takes a list of gates like [H, T, H, T, H]
    -- and returns the rewritten form N1' N2' C described on the paper
    -- exercise, by repeatedly conjugating trailing non-Clifford gates
    -- through the Clifford gates that follow them.
