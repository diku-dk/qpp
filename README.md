# QPP — Quantum Programming Playground

QPP is a Haskell framework in which a quantum program is an expression in an
algebra of unitary operators. It was written for the DIKU MSc course *Hybrid
Quantum Programming*, whose aim is to train people who can build quantum
software — algorithms, simulators, compilers, analyses — and not only use it.
Contrary to large industrial platforms such as Qiskit or Pytket, QPP's core
fits on a page: a program is a value of one small data type, a compiler pass
is a function on that type, and a simulation backend is a module that
interprets it. d, simplified, converted to a

## Operators

An operator is a term of type `QOp`, defined in `QPP.Syntax`:

```haskell
data QOp
  = Id Nat              -- n-qubit identity; Id 0 = One is the scalar unit, Id 1 = I
  | Phase Rational      -- global phase e^{iπθ}
  | X | Y | Z | H | SX  -- single-qubit gates
  | R QOp Rational      -- rotation by θπ about an axis given by a Pauli string
  | C QOp               -- controlled operator, Id n ⊕ op
  | Permute [Int]       -- permutation of qubits
  | Tensor QOp QOp      -- a ⊗ b
  | DirectSum QOp QOp   -- a ⊕ b: "if the first qubit is 0 then a, else b"
  | Compose QOp QOp     -- a ∘ b: b first, as in mathematics
  | Adjoint QOp         -- adj a
```

Angles are rational multiples of π, so `R Z (1/2)` is exactly S and
`R Z (1/4)` exactly T. There is no SWAP constructor; the swap of two adjacent
qubits is `Permute [1,0]`.

The indices in `Permute` refer to the operator's own qubits, not to positions
in the register it will end up in, so a permutation can be tensored,
controlled or direct-summed without renumbering. To place an operator in a
larger register, `k <@ op @> l` abbreviates `Id k ⊗ op ⊗ Id l`: `op` at qubit
`k`, with `l` idle qubits after it.

The combinators `⊗`, `⊕`, `∘` and `adj` (in ASCII, `<.>`, `<+>` and `<>`) are
methods of the `Operator` class, implemented by `QOp` and by the dense
backend's matrix type `CMat`. A function written against the class builds a
term or a matrix depending on the type it is used at.

## Programs

A program is a list of steps:

```haskell
data Step
  = Unitary QOp
  | Initialize [Nat] [Bool]   -- set the given qubits to classical values
  | Measure    [Nat]

type Program  = [Step]
type Outcomes = [Bool]        -- most recent first
type RNG      = [Double]      -- an infinite stream in [0,1)
```

Measurement draws from a random stream that the caller supplies, so a run can
be repeated exactly.

## Example: the quantum Fourier transform

The QFT on *n* qubits is the QFT on the first *n−1* qubits, after a layer
that puts a Hadamard on the last qubit and then a phase on it controlled by
each of the earlier qubits. QPP lets you write exactly that:

```haskell
import QPP
import QPP.Semantics.Matrix        -- exactly one backend; swap this line to switch

-- | Phase gate diag(1, e^{2πi/2^k}): a Z-rotation times a global phase.
p :: Int -> QOp
p k = Phase (1/2^k) ⊗ R Z (2/2^k)

-- | Hadamard on the last of n qubits, then a phase on it controlled by each earlier qubit.
layer :: Int -> QOp
layer n = foldr (∘) ((n-1) <@ H) [ k <@ C ((n-k-2) <@ p (n-k)) | k <- [0..n-2] ]

-- | QFT on n qubits, output bits in reverse order.
qftrev :: Int -> QOp
qftrev 0 = One
qftrev n = (qftrev (n-1) ⊗ I) ∘ layer n

-- | QFT on n qubits.
qft :: Int -> QOp
qft n = Permute [n-1, n-2 .. 0] ∘ qftrev n

main :: IO ()
main = do
  putStrLn (showOp (cleanop (qft 2)))          -- the simplified term
  printM (evalOp (qft 2))                      -- its matrix
  let (_, outcomes, _) = evalProg [Unitary (qft 3), Measure [0,1,2]] (ket [0,0,1]) (cycle [0.3, 0.7])
  print outcomes                               -- a run, with an explicit random stream
```

```
(Permute [1,0] ∘ ((H ⊗ I) ∘ (C ((Phase (1 % 4) ⊗ R Z (1 % 2))) ∘ (I ⊗ H))))
0.5 *
 [[1, 1, 1, 1],
  [1, 1,-1,-1],
  [1,-1, 0, 0],
  [1,-1, 0, 0]]

0.5 *
i *
 [[0,0, 0, 0],
  [0,0, 0, 0],
  [0,0, 1,-1],
  [0,0,-1, 1]]

[False,True,False]
```

The same `qft` is in `Programs.QFT`, where `Programs.QuantumAdder` builds
addition from it.

## Backends

A backend is a module, not a class instance. Every backend defines the same
two types and six functions:

```haskell
StateT, OpT                                          -- states and operators
ket      :: [Int] -> StateT                          -- a computational-basis state
evalOp   :: QOp -> OpT                               -- interpret a term
apply    :: OpT -> StateT -> StateT
measure1 :: (StateT, Outcomes, RNG) -> Nat -> (StateT, Outcomes, RNG)
evalStep :: (StateT, Outcomes, RNG) -> Step -> (StateT, Outcomes, RNG)
evalProg :: Program -> StateT -> RNG -> (StateT, Outcomes, RNG)
```

| Module | Representation |
|---|---|
| `QPP.Semantics.Matrix` | dense `hmatrix` matrices; the reference semantics |
| `QPP.Semantics.Statevector` | `massiv` state vectors |
| `QPP.Semantics.MPS` | matrix product states with truncation |
| `QPP.Semantics.Stabilizer` | empty; writing the stabilizer-tableau backend is a course exercise |

A program imports `QPP` and exactly one backend. `QPP` re-exports the syntax,
rewriting and pretty-printing layers but no backend, so the choice of backend
is always visible in the import list.

A new backend defines the two types, implements `ket`, `evalOp`, `apply` and
`measure1`, and takes the two interpreters from `QPP.Semantics`:

```haskell
evalStep = evalStepWith (apply . evalOp) measure1
evalProg = evalProgWith (apply . evalOp) measure1
```

To have it validated, give `OpT` a `Convertible OpT CMat` instance and add
the module to `qpp/test/BackendValidation.hs`, which compares each backend it
imports against the dense reference on random operators. The convention is
spelled out in the header comment of `QPP.Semantics`.

## Compiler passes and analyses

A compiler pass is a function `QOp -> QOp`. `QPP.Rewrite` contains the
semantics-preserving ones that come with the library — merging identities,
pushing adjoints inward, flattening compositions — and `simplifyFixpoint`, which
applies a list of passes until nothing changes; `cleanop` runs all of the
built-in ones. Two of them, `pushComposes` and `liftComposes`, bring a term
into the normal forms "tensor products of compositions" and "compositions of
tensor products", the natural starting point for a translation to a
stabilizer or ZX representation. An analysis is a function that reads a
term: `QPP.Syntax` has `op_qubits`, `op_support` and `dagger`, and
`QPP.PrettyPrint` has the printers used above.

## Building and running

QPP needs GHC 9.6 or later, cabal 3.10 or later, and BLAS/LAPACK for
`hmatrix`.

```sh
git clone git@github.com:diku-dk/QPP.git
cd QPP
cabal build lib:qpp          # the library
cabal build all              # library, exercises, problem sheets, test suites
cabal test all               # the six test suites
cabal run NoCloning          # a worked exercise
cabal run QubitsProblem      # a weekly problem sheet
cabal repl lib:qpp           # the library in GHCi
```

## Layout

```
qpp/
  src/QPP/           the library: Syntax, Rewrite, Semantics, the backends, PrettyPrint, Util
  src/Programs/      example programs: QFT, quantum adder, repeater protocol,
                     block encodings, QSVT phase sequences
  exercises/         worked course demos, one executable each
  problems/          weekly problem sheets: they compile and run, with the
                     interesting functions left as stubs (see problems/README.md)
  solutions/         answer keys; not part of the cabal build
  test/              tasty/QuickCheck suites (see test/README.md)
  limbo/             parked designs and legacy code; not built
  CHANGELOG.md       HQPlayground → QPP module and symbol map
docs/                notes on block encodings, QSVT and matrix inversion; design notes
```

## Coming from HQPlayground

QPP was previously named HQPlayground. The module names changed with the
rename; the mapping, including moved symbols, is in
[`qpp/CHANGELOG.md`](qpp/CHANGELOG.md). The rationale for `Id n`, `Phase`,
rational angles and relative `Permute` indices is in
[`docs/v2-to-v3-changes.md`](docs/v2-to-v3-changes.md).

## Authors and license

QPP is developed by James Avery, Thomas Mork, Michael Kirkedal Thomasen, and Fritz Henglein
at the Department of Computer Science, University of Copenhagen (DIKU).
It is released under the BSD 3-Clause license.
