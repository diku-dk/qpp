module Main where

import QPP
import QPP.Semantics.Matrix

import Data.Complex (Complex((:+)), magnitude)
import Data.Ratio (approxRational)

{-| Exercise: State preparation and the Input Problem.

This exercise accompanies the "State Preparation and the Input Problem"
problem sheet (Exercise 1, 2, 5).

1. Complete `prepareBasisState`, which prepares a computational basis
   state |b1 b2 ... bn> from |00...0> using only X gates on the qubits
   where b_k = 1, built as a single n-qubit QOp via `Tensor` (paper
   Exercise 1).

2. Complete `controlledRotationPrep2`, a hand-built state-prep circuit
   for a *general real* 2-qubit state, using one Hadamard-like R_Y
   rotation on qubit 1 followed by a uniformly-controlled R_Y rotation
   on qubit 2 (i.e. two DIFFERENT rotation angles on qubit 2, selected
   by the value of qubit 1). This is the n=2 base case of the
   recursive construction from paper Exercise 2 -- work out what
   angles you need for a specific target state, and build the circuit
   with `R Y theta`, `C`, and `Compose`/`Tensor`.

3. Complete `gateCount`, a simple structural recursion over the `QOp`
   AST that counts how many "elementary" gates (leaves like X, Y, Z, H,
   SX, R _ _, Permute; a `C op` or `DirectSum op1 op2` counts its
   sub-circuit(s) plus itself) appear in a circuit. Use it to see how
   circuit size grows for the recursive state-prep construction
   (paper Exercise 2), by building it for n = 1, 2, 3 qubits.

   Hmmm: Dette er en sjov ... men ikke helt triviel opgave .... Hvad gør vi ved den? 


4. (Discussion, paper Exercise 4 -- no code needed) Why would loading N
   arbitrary real numbers via a QRAM-style structure still cost you
   Omega(N) gates even though queries have O(log N) *depth*? Write your
   answer as a comment in `main`.
-}

-- | Prepare |b1..bn> from |00..0> using X on exactly the qubits with
-- b_k = 1. `bits` is e.g. [1,0,1] for |101>.
prepareBasisState :: [Int] -> QOp
prepareBasisState [] = error "Cannot fold an empty list of integers!"
prepareBasisState bits = foldr1 Tensor (map chooseOp bits)
  where
    chooseOp n = if n `mod` 2 == 0 then I else X -- Handles operation modulo 2 ... other options can be chosen


-- | A hand-built 2-qubit real-amplitude state-prep circuit (angles are
-- `Rational`, in units of pi):
-- rotate qubit 1 by angle0, then apply a *different* R_Y rotation to
-- qubit 2 depending on qubit 1's value (a "uniformly controlled
-- rotation", the building block from paper Exercise 2).
controlledRotationPrep2 :: Rational -> Rational -> Rational -> QOp
controlledRotationPrep2 angle0 angle1if0 angle1if1 =
   let 
      firstRot = (R Y angle0) ⊗ I  -- This only handles the qubit-1 rotation.
    
      -- You still need to apply R Y angle1if0 to qubit 2 when qubit 1 = 0,
      -- and R Y angle1if1 when qubit 1 = 1. Hint: DirectSum applies its
      -- two arguments conditioned on a control qubit's value -- but check
      -- carefully which qubit plays the role of control here, and how to
      -- Compose this with the qubit-1 rotation above.
      
      -- |0><0|⊗(R Y angle1if0) + |1><1|⊗(R Y angle1if1)
      -- Hvilket er det samme som (R Y angle1if0) ⊕ (R Y angle1if1) ... pga Kronecker-produktet.
      secondRot = (R Y angle1if0) ⊕ (R Y angle1if1)
   in
      secondRot ∘ firstRot
      
{-| Count "elementary" gates in a QOp circuit description.

Note `Id _` rather than the `I` the stub started from: `I` is a pattern
synonym for `Id 1`, so matching on it would have missed `Id 2`, `Id 3`, ...
and charged a gate for every idle wire.
-}
gateCount :: QOp -> Int
gateCount (Id _)          = 0   -- identities are wires, not gates
gateCount (Phase _)       = 0   -- an unobservable scalar costs nothing
gateCount X               = 1
gateCount Y               = 1
gateCount Z               = 1
gateCount H               = 1
gateCount SX              = 1
gateCount (R _ _)         = 1
gateCount (Permute _)     = 1
gateCount (C op)          = 1 + gateCount op
gateCount (Tensor a b)    = gateCount a + gateCount b
gateCount (DirectSum a b) = 1 + gateCount a + gateCount b
gateCount (Compose a b)   = gateCount a + gateCount b
gateCount (Adjoint op)    = gateCount op

--------------------------------------------------------------------------
-- Exercise 2: which angles reach a given target state?
--------------------------------------------------------------------------

{-| The three R_Y angles (in units of pi) that prepare a real, normalized
2-qubit amplitude vector [a00, a01, a10, a11] from |00>.

Run the circuit symbolically. With alpha = pi*angle0/2 and
beta_i = pi*angle1if_i/2, and R Y theta |0> = cos(pi theta/2)|0> + sin(pi theta/2)|1>:

    after R_Y(angle0) on qubit 0:  cos(alpha)|0> + sin(alpha)|1>
    after the controlled R_Y:      cos(alpha)cos(beta0)|00> + cos(alpha)sin(beta0)|01>
                                 + sin(alpha)cos(beta1)|10> + sin(alpha)sin(beta1)|11>

So alpha splits the weight between the two halves of the register, and
beta0/beta1 split it again inside each half. Inverting that is just three
arctangents. This is the n = 2 case of the general recursion: at each level
you spend one uniformly controlled rotation to split every branch you have
so far.

A caveat worth noticing: QPP angles are `Rational` (in units of pi), but the
angles a general target needs are irrational, so `approxRational` is doing
real work here. Exact state preparation of an arbitrary vector is not
available at any finite gate count -- you always approximate.
-}
prepAngles :: [Double] -> (Rational, Rational, Rational)
prepAngles [a00, a01, a10, a11] =
    ( angle (sqrt (a10*a10 + a11*a11)) (sqrt (a00*a00 + a01*a01))
    , angle a01 a00
    , angle a11 a10 )
  where
    -- theta with (cos(pi theta/2), sin(pi theta/2)) parallel to (c, s)
    angle s c = approxRational (2 / pi * atan2 s c) 1e-15
prepAngles _ = error "prepAngles: expected exactly four amplitudes"

-- | The 2-qubit state with the given real amplitudes.
realState :: [Double] -> StateT
realState as = foldr1 (.+)
    [ (a :+ 0) .* ket b | (a, b) <- zip as [[0,0],[0,1],[1,0],[1,1]] ]

-- | Amplitudes <b|psi> for b = |00>, |01>, |10>, |11>.
amplitudes :: StateT -> [ComplexT]
amplitudes psi = [ inner (ket b) psi | b <- [[0,0],[0,1],[1,0],[1,1]] ]



main :: IO ()
main = do
    putStrLn "-- Exercise 1: basis state preparation --"
    let target  = [1, 0, 1]
        prepOp  = evalOp (prepareBasisState target)
        psi000  = ket [0,0,0]
    putStrLn $ "Preparing |101>: " ++ showState (apply prepOp psi000)

    putStrLn "\n-- Exercise 2: a hand-built 2-qubit state-prep circuit --"
    -- Target: amplitudes proportional to (1,2,3,4), normalized.
    let want            = map (/ sqrt 30) [1, 2, 3, 4]
        (a0, b0, b1)    = prepAngles want
        prep2           = controlledRotationPrep2 a0 b0 b1
        got             = apply (evalOp prep2) (ket [0,0])
    putStrLn $ "Target amplitudes:  " ++ show want
    putStrLn $ "Angles (units of pi): " ++ show (fromRational a0 :: Double)
                                        ++ ", " ++ show (fromRational b0 :: Double)
                                        ++ ", " ++ show (fromRational b1 :: Double)
    putStrLn $ "Prepared state:     " ++ showState got
    putStrLn $ "Fidelity |<target|prepared>| = "
                ++ show (magnitude (inner (realState want) got))

    putStrLn "\n-- Exercise 3: counting gates as n grows --"
    -- n = 1 is a bare rotation; n = 2 adds one uniformly controlled rotation.
    let prep1 = R Y a0
    putStrLn $ "gateCount (n = 1) = " ++ show (gateCount prep1)
    putStrLn $ "gateCount (n = 2) = " ++ show (gateCount prep2)
    putStrLn ""
    putStrLn "Level k of the recursion has to split every branch produced so far, so"
    putStrLn "it costs 2^(k-1) rotations: T(n) = T(n-1) + 2^(n-1), T(1) = 1, giving"
    putStrLn "T(n) = 2^n - 1 rotations. That is exactly one angle per amplitude of the"
    putStrLn "target (minus the one fixed by normalization) -- the circuit cannot be"
    putStrLn "smaller, because it has to carry that much information."
    putStrLn ""
    putStrLn "(n = 2 prints 4, not 3: the sheet's definition of gateCount also charges"
    putStrLn "one for each DirectSum or C node, so this is 3 rotations plus the one"
    putStrLn "multiplexer that selects between them. The 2^n growth is in the rotations.)"

    putStrLn "\n-- Exercise 4 (discussion): why QRAM does not rescue you --"
    putStrLn "A QRAM query has O(log N) *depth*: the address tree is only log N deep,"
    putStrLn "so a single lookup is fast. But depth is not gate count. The tree has"
    putStrLn "Omega(N) nodes, every one of which needs hardware holding one of the N"
    putStrLn "data values, and building or loading that structure touches all of them."
    putStrLn ""
    putStrLn "The information-theoretic version of the same point: N arbitrary reals"
    putStrLn "carry Omega(N) bits, and a circuit of g gates over a fixed gate set"
    putStrLn "carries O(g log g) bits of description. No circuit with g = o(N) can"
    putStrLn "encode an arbitrary N-entry input, whatever gates you allow it."
    putStrLn ""
    putStrLn "So an algorithm advertising an exponential speedup in N is only honest"
    putStrLn "if its input is *computed* rather than *loaded* -- given by a formula or"
    putStrLn "a sparse oracle, not read off a list. This is the Input Problem, and it"
    putStrLn "is where a lot of claimed quantum advantage quietly goes."
