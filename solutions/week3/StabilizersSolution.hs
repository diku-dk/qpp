module Main where

import QPP
import QPP.Semantics.MPS

{-| SOLUTIONS: Stabilizers, check matrices, and the Clifford group.

This exercise accompanies the "Stabilizers and the Clifford Group" problem
sheet. Two independent halves:

  PART A (Exercises 6-9): bit-vector arithmetic on stabilizer tableaux. No
  QPP needed -- the check-matrix formalism is a classical data structure.

  PART B (Exercises 10-11): verifying Clifford conjugation identities
  numerically against the QPP operators.

The bridge between them is `sameOp` at the end of Part A: it checks the
XOR sign rule against the actual matrix product, so the bit twiddling is
not taken on faith.
-}

-- | Comparison threshold for the numerical checks.
eps :: Double
eps = 1e-9

--------------------------------------------------------------------------
-- PART A: stabilizer tableaux as bit vectors (no QPP needed)
--------------------------------------------------------------------------

{-| A Pauli group element on n qubits, encoded as (s, c, xs, zs). The
    operator it denotes is

      (-1)^s i^c  (X^x_0 Z^z_0) (x) ... (x) (X^x_(n-1) Z^z_(n-1))

    Bits are Bool (False = 0, True = 1). Note the X-before-Z ordering: it
    is what makes the multiplication rule below a pure XOR plus one sign
    bit, and it is why the phase bits s and c are needed at all.
-}
data PauliCode = PauliCode
  { signBit  :: Bool   -- ^ s
  , phaseBit :: Bool   -- ^ c
  , xBits    :: [Bool] -- ^ x_0 .. x_{n-1}
  , zBits    :: [Bool] -- ^ z_0 .. z_{n-1}
  } deriving (Eq, Show)

xor :: Bool -> Bool -> Bool
xor = (/=)

-- | F2 inner product: the parity of the bitwise AND.
dotF2 :: [Bool] -> [Bool] -> Bool
dotF2 u v = foldr xor False (zipWith (&&) u v)

{-| Exercise 6: encode a single-qubit Pauli symbol into its (x,z) pair.

Careful with Y. In this encoding (1,1) means X Z, and Y = i X Z, so
X Z = -i Y: the symbol table below is exact only if the missing factor is
tracked in the phase bits. Nothing in Exercises 7-9 depends on it --
commutation and rank ignore phases entirely -- but `pauliToQOp` at the end
of Part A does have to get it right, and it does (see `symbolOp`).
-}
encodeSymbol :: Char -> (Bool, Bool)
encodeSymbol 'I' = (False, False)
encodeSymbol 'X' = (True,  False)
encodeSymbol 'Z' = (False, True)
encodeSymbol 'Y' = (True,  True)
encodeSymbol c   = error ("encodeSymbol: not a Pauli symbol: " ++ [c])

-- | Exercise 6: pretty-print a PauliCode as a signed tensor-product string.
--   (Prints (1,1) as "Y", so read it as "X Z up to the phase bits".)
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

{-| Exercise 7: multiply two PauliCodes.

Per qubit, pushing the second X past the first Z is the only place a sign
can appear:

    (X^x1 Z^z1)(X^x2 Z^z2) = (-1)^(z1 x2) X^(x1+x2) Z^(z1+z2)

because Z X = -X Z. Summing the exponent over all qubits, the total sign
correction is the F2 inner product z1 . x2, and the x and z rows are just
XORed. Assuming c1 = c2 = 0, as the sheet allows, the product needs no
lone i factor either, so c stays 0.
-}
pauliMultiply :: PauliCode -> PauliCode -> PauliCode
pauliMultiply (PauliCode s1 c1 x1 z1) (PauliCode s2 c2 x2 z2) =
  PauliCode (s1 `xor` s2 `xor` dotF2 z1 x2)
            (c1 `xor` c2)
            (zipWith xor x1 x2)
            (zipWith xor z1 z2)

{-| Exercise 8: the symplectic product, ignoring sign and phase:

    <(x1,z1),(x2,z2)> = x1 . z2 + z1 . x2  (mod 2)

This is exactly the sign that survives when you compare P1 P2 against
P2 P1 -- True means they ANTI-commute. -}
symplecticProduct :: [Bool] -> [Bool] -> [Bool] -> [Bool] -> Bool
symplecticProduct x1 z1 x2 z2 = dotF2 x1 z2 `xor` dotF2 z1 x2

-- | Exercise 8: do two Pauli strings commute?
commutes :: PauliCode -> PauliCode -> Bool
commutes p q = not (symplecticProduct (xBits p) (zBits p) (xBits q) (zBits q))

{-| Exercise 9: rank over F2 by Gaussian elimination.

Take the leftmost column, use any row with a 1 there as pivot, XOR it into
every row below that also has a 1, drop the column, recurse. A set of n
generators on n qubits is independent -- and so pins down a stabilizer
state -- exactly when the rank is n.
-}
rankF2 :: [[Bool]] -> Int
rankF2 rows
  | null rows        = 0
  | null (head rows) = 0
  | otherwise = case span (not . head) rows of
      (_, [])              -> rankF2 (map tail rows)   -- column is all zero
      (zeroes, pivot:more) ->
          1 + rankF2 (map tail (zeroes ++ map (eliminate pivot) more))
  where
    eliminate p r = if head r then zipWith xor p r else r

-- | The (x|z) row of a PauliCode, as `rankF2` wants it.
checkRow :: PauliCode -> [Bool]
checkRow p = xBits p ++ zBits p

--------------------------------------------------------------------------
-- The bridge: does the bit arithmetic agree with the matrices?
--------------------------------------------------------------------------

-- | The QOp denoted by a PauliCode. This is where the X Z = -i Y convention
--   has to be spelled out honestly: symbol (1,1) becomes `X ∘ Z`, not `Y`.
pauliToQOp :: PauliCode -> QOp
pauliToQOp (PauliCode s c xs zs)
  | null xs   = Phase (phaseOf s c)
  | otherwise = Phase (phaseOf s c) ⊗ foldr1 (⊗) (zipWith symbolOp xs zs)
  where
    -- (-1)^s i^c = e^{i pi (s + c/2)}
    phaseOf a b = (if a then 1 else 0) + (if b then 1/2 else 0)
    symbolOp x z = case (x, z) of
        (False, False) -> I
        (True,  False) -> X
        (False, True ) -> Z
        (True,  True ) -> X ∘ Z

-- | Do two n-qubit QOps act identically? Checked on every basis state.
sameOp :: Int -> QOp -> QOp -> Bool
sameOp n a b = all agree [0 .. 2 ^ n - 1]
  where
    agree x = let psi = ket (bitsOf n x)
              in norm (apply (evalOp a) psi .- apply (evalOp b) psi) < eps

-- | Big-endian bit list: `bitsOf 3 5 == [1,0,1]`.
bitsOf :: Int -> Int -> [Int]
bitsOf n x = [ (x `div` 2 ^ i) `mod` 2 | i <- [n-1, n-2 .. 0] ]

--------------------------------------------------------------------------
-- PART B: verifying Clifford conjugation with QPP
--------------------------------------------------------------------------

{-| (U P U^dagger) |phi>, read right to left: U^dagger acts on phi first.

`dagger` (from QPP.Syntax) pushes the adjoint structurally down the tree,
so we do not have to know U's inverse by hand. And `Compose`/`∘` is
right-to-left, so `a ∘ b` applies b first -- consistent with the formula
as written, which is why the nesting below reads the way it does.
-}
conjugateApply :: QOp -> QOp -> StateT -> StateT
conjugateApply u p phi =
  apply (evalOp u) (apply (evalOp p) (apply (evalOp (dagger u)) phi))

-- | Exercise 10/11: check U P U^dagger = Q by comparing both sides on a few
--   test states. Equality is exact here, not up to phase: conjugation by a
--   Clifford maps Paulis to Paulis on the nose, signs included.
checkConjugation :: String -> QOp -> QOp -> QOp -> [StateT] -> IO ()
checkConjugation name u p claimed testStates =
    putStrLn ("  " ++ name ++ ": " ++ show (all agree testStates))
  where
    agree phi = norm (conjugateApply u p phi .- apply (evalOp claimed) phi) < eps

main :: IO ()
main = do
    putStrLn "-- Exercise 7: multiplying Pauli strings with the XOR rule --"
    -- (101|011) and (101|110): X (x) Z (x) XZ  and  XZ (x) Z (x) X
    let p1 = PauliCode False False [True, False, True] [False, True, True]
        p2 = PauliCode False False [True, False, True] [True,  True, False]

    putStrLn $ "  P1       = " ++ decodePauli p1
    putStrLn $ "  P2       = " ++ decodePauli p2
    putStrLn $ "  P1 P2    = " ++ decodePauli (pauliMultiply p1 p2)
    putStrLn $ "  agrees with the matrix product? "
                ++ show (sameOp 3 (pauliToQOp (pauliMultiply p1 p2))
                                  (pauliToQOp p1 ∘ pauliToQOp p2))
    putStrLn "  The minus sign is real: z1 . x2 = 1, from the third qubit where"
    putStrLn "  X Z meets X and Z X = -X Z."

    putStrLn "\n-- Exercise 8: commutation --"
    putStrLn $ "  Do P1, P2 commute? " ++ show (commutes p1 p2)
    putStrLn $ "  cross-check on the matrices: "
                ++ show (sameOp 3 (pauliToQOp p1 ∘ pauliToQOp p2)
                                  (pauliToQOp p2 ∘ pauliToQOp p1))
    putStrLn "  x1.z2 = 1 and z1.x2 = 1, so the symplectic product is 0: they commute."

    putStrLn "\n-- Exercise 9: independence --"
    let rows =
          [ [True,  True,  False, False, False, True ]   -- P1 = (110|001)
          , [False, True,  True,  True,  False, False]   -- P2 = (011|100)
          , [True,  False, True,  True,  False, True ]   -- P3 = (101|101)
          ]
    putStrLn $ "  rank = " ++ show (rankF2 rows) ++ " out of 3 rows on 3 qubits"
    putStrLn "  The rank is 2, not 3: row1 XOR row2 = row3, i.e. P3 = P1 P2 up to"
    putStrLn "  phase. These three generators are NOT independent, so they do not"
    putStrLn "  pin down a single stabilizer state -- they fix a 2-dimensional"
    putStrLn "  subspace (a 1-qubit code space), not a point."

    putStrLn "\n-- Exercises 10/11: Clifford conjugation, checked numerically --"
    let testStates = [ket [0], ket [1], apply (evalOp H) (ket [0])]
    checkConjugation "H X H^dagger = Z" H X Z testStates
    checkConjugation "H Z H^dagger = X" H Z X testStates
    checkConjugation "H Y H^dagger = -Y" H Y (Phase 1 ⊗ Y) testStates
    checkConjugation "S X S^dagger = Y" (R Z (1/2)) X Y testStates
    checkConjugation "S Z S^dagger = Z" (R Z (1/2)) Z Z testStates
    putStrLn "  H Y H^dagger needs the explicit -1 (written `Phase 1 ⊗ Y`): the sign"
    putStrLn "  is part of the answer, and this is exactly the bookkeeping the sign"
    putStrLn "  bit s in PauliCode exists to carry."

    putStrLn "\n-- Exercise 12 (open-ended): pushing T through Cliffords --"
    putStrLn "  The rewrite is C T = T' C with T' = C T C^dagger. For a Clifford C"
    putStrLn "  and T = R Z (1/4), T' is a pi/4 rotation about the *conjugated* axis"
    putStrLn "  C Z C^dagger, which is a Pauli -- so a T gate stays a pi/4 rotation"
    putStrLn "  about some Pauli axis as it moves right, and only its axis changes."
    putStrLn "  Sweeping right to left through [H,T,H,T,H] therefore collects every"
    putStrLn "  non-Clifford factor at the front and leaves one Clifford at the back:"
    putStrLn "  the N1' N2' C form. The count of T gates never drops, which is the"
    putStrLn "  point -- T-count is invariant under this rewriting, and it is the"
    putStrLn "  honest measure of how far a circuit sits outside Gottesman-Knill."
