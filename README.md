# QPP — Quantum Programming Playground

QPP (Quantum Programming Playground) is an open platform for teaching and research, by which
students and researchers can build new quantum software (algorithms, compilers and analyses)
from a small set of building blocks that combine according to a few clear algebraic rules.

## Operators

An operator is a term of type `QOp`, defined in `QPP.Syntax`:

```haskell
data QOp
  = Id Nat              -- n-qubit identity; Id 0 = One is the scalar unit, Id 1 = I
  | Phase Rational      -- global phase e^{iπθ}
  | X | Y | Z | H | SX  -- single-qubit gates
  | R QOp Rational      -- rotation by θπ; the axis must be Hermitian unless θ = 0
  | C QOp               -- controlled operator, Id n ⊕ op
  | Permute [Nat]       -- permutation of qubits
  | Tensor QOp QOp      -- a ⊗ b
  | DirectSum QOp QOp   -- a ⊕ b: "if the first qubit is 0 then a, else b"
  | Compose QOp QOp     -- a ∘ b: b first, as in mathematics
  | Adjoint QOp         -- adj a
```

Angles are rational multiples of π. `R Z (1/2)` and `R Z (1/4)` equal S and T
up to global phases; `C (Phase (1/2))` and `C (Phase (1/4))` give S and T
exactly. There is no SWAP constructor; the swap of two adjacent
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

### General finite-dimensional unitaries

`QPP.UOp` provides the smaller `UOp` syntax for unitaries on arbitrary
finite-dimensional spaces, including dimension zero. Products and direct
sums take lists; `Compose d []` is the identity of matrix dimension `d`.
Use a qualified import alongside `QPP`, since constructor names overlap:

```haskell
import qualified QPP.UOp as U

x = U.PermuteSum [1,1] [1,0]
z = U.directSum [U.identity 1, U.phase 1]
u = U.compose 2 [U.R x (-1/2), z]
-- U.showUOp u:
-- "(R(PermuteSum [1,1] [1,0], -1/2) ∘ (I[1] ⊕ R(I[1], -2)))"
```

`R a t` means `exp(-i*pi*t*A/2)`, with the same sign and half angle as QOp,
and a Hermitian axis unless `t` is zero. Scalar phase `U.phase t` still means
`exp(i*pi*t)` and is represented by `R (identity 1) (-2*t)`.
`U.dim` counts basis vectors,
not qubits. The shared operator notation also works for `UOp`; its builders
flatten nested lists without performing algebraic simplification.
`U.fromQOp :: QOp -> Either TypeError U.UOp` checks and lowers a QOp,
preserving its global phase. Apply `U.normalize` to the result to flatten
associations, remove identities, collect scalar phases, simplify rotations,
and fuse or cancel adjacent composition factors. It preserves relative phases
inside direct sums. Non-scalar rotation angles use `(-2,2]` modulo four;
scalar phase parameters use `(-1,1]` modulo two before encoding as `R`.
Interchange lifts compositions out of tensors and direct sums, padding shorter
operand lists with identities. This produces a sequence of parallel layers.
It is a structural normal form, not a complete test for operator equality;
tensor products of direct sums stay factored.

Normalization assumes well-formed input; it does not validate raw UOp syntax.
General UOp validation and semantics are not yet implemented. Run
`cabal run QFTSyntax -- 4` to write the QOp pretty print, raw and pretty UOp,
and `uop.normalized.raw.txt` / `uop.normalized.pretty.txt` into `output/qft4/`.
The normalized pretty print puts each outer composition factor on its own line.
Use `16` for the corresponding sixteen-qubit example.

### TikZ operator diagrams

`showTikZQOp` and `showTikZUOp` produce LaTeX fragments for both syntaxes:

```haskell
import qualified QPP as Q
import qualified QPP.UOp as U

writeFile "qop.tex" (Q.showTikZQOp (Q.Compose (Q.C Q.X) (Q.Tensor Q.H Q.I)))
writeFile "uop.tex" (U.showTikZUOp (U.DirectSum [U.identity 1, U.phase (1/4)]))
```

Copy `latex/qpp-tikz.sty` beside the document, load `\usepackage{qpp-tikz}`,
then use `\input{qop.tex}`. Each operator is a node in a syntax tree, with its
operands as children from left to right. Outlined nodes use pale constructor
colors, joined by compact horizontal branches with vertical drops. The rightmost
child of a composition acts first. Printing preserves the input syntax,
including arbitrary UOp dimensions and zero-dimensional spaces, and does not
type-check or normalize it.

Every node has a gray qubit count outside its outline, such as `[⊕]: 4`.
Counts assume well-formed input and are computed from the children without
building matrices.
For UOp, dimension `2^n` gives count `n`; other dimensions are explicitly labelled
`d=...`. Thus dimension one has zero qubits, while dimension zero is labelled
`d=0`. Use `showTikZQOpWithoutQubits` or `showTikZUOpWithoutQubits` to omit
the outside annotations.

Restyle diagrams with ordinary TikZ keys, without changing generated files:

```latex
\tikzset{
  qpp/gate/.append style={fill=yellow!20},
  qpp/node/.append style={font=\footnotesize},
  qpp/edge/.append style={thick, draw=blue!60!black},
  qpp/annotation/.append style={text=black!70}
}
\forestset{qpp tree/.append style={for tree={l sep=5mm, s sep=7mm, fork sep=2.5mm}}}
```

Shared TikZ styles are `qpp/node`, `qpp/edge`, and `qpp/annotation`;
the Forest style `qpp tree` controls layout: `l sep` sets the vertical gap,
`s sep` the horizontal spacing,
and `fork sep` the drop from the parent to the horizontal branch. Constructor
styles are `gate`, `identity`, `phase`, `permutation`,
`tensor`, `sum`, `compose`, `rotation`, `control`, and `adjoint`,
each under `qpp/`. Both syntaxes use `R(t)` for the same rotation convention;
identity subscripts and UOp composition subscripts are matrix dimensions.

The implementation in `QPP.PrettyPrint.TikZ` translates each syntax into
`Data.Tree`, then recursively emits [Forest](https://ctan.org/pkg/forest)'s bracket
notation. Forest handles spacing and edges; the `.sty` file holds the styles
and measures annotation widths to keep neighboring nodes clear.
The Haskell code never estimates text widths or calculates coordinates.

Run `cabal run TikZOperators` from the repository root to generate examples
and QFT fragments in `output/tikz/`. Run `make -C latex` to generate them and
compile the gallery to `output/pdf/operator-diagrams.pdf` (requires `latexmk`
and a TeX installation with Forest). The gallery source is
`latex/operator-diagrams.tex`. Large terms retain their full structure;
they may need a larger page or explicit decomposition for presentation.

### Type-checking

`typeCheck`, exported by `QPP` and `QPP.TypeCheck`, checks a `QOp` and returns
its qubit count or the first `TypeError`:

```haskell
typeCheck (C X)                    -- Right 2
typeCheck (R (X <.> Z) (1/3))      -- Right 2
typeCheck (X <+> Id 2)             -- Left (ArityMismatch (DirectSum X (Id 2)) 1 2)
typeCheck (R SX (1/3))             -- Left (UnprovenHermitian SX)
typeCheck (R SX 0)                 -- Right 1: zero-angle rotations are identity
```

Identity widths must be nonnegative; composition and direct sum require
equal operand widths. `Permute ks` must contain each index in
`[0 .. length ks - 1]` exactly once, so `Permute []` is the zero-qubit
identity. Tensor, control, and direct sum widths are also checked for
`Int` overflow.

A rotation axis must be provably Hermitian when the angle is nonzero.
At exactly zero, the rotation is identity for any well-formed axis. The
axis is still checked recursively for valid widths, permutations, and
nested rotations; a zero angle cannot hide an invalid subterm.
The checker propagates the more general property
`A† = exp(i*pi*p) A`, using exact rational phases modulo 2. The exported
`adjointPhase :: QOp -> Either TypeError (Maybe Rational)` exposes this proof:
`Just 0` proves Hermiticity, another `Just p` proves non-Hermiticity, and
`Nothing` means unknown. An invalid subterm still returns `Left`.

A **phase-labelled Pauli expression** has the form
`exp(i*pi*q) (P₀ ⊗ … ⊗ Pₙ₋₁)`, where each `Pᵢ` is `I`, `X`, `Y`, or `Z`.
For example, `X <> Y = iZ`, while `(X ⊗ X) <> (Z ⊗ Z) = -(Y ⊗ Y)`.
Products can therefore introduce phases even without an explicit `Phase`.
The checker normalizes these expressions to a rational phase and a sparse
map of nonidentity Pauli factors. Every bare Pauli string is Hermitian,
so the expression is Hermitian **exactly when `q` is an integer**.

This gives a complete Hermiticity decision for the well-formed sublanguage
generated by `Id`, `Phase`, `X`, `Y`, `Z`, `Tensor`, `Compose`, and `Adjoint`,
with arbitrary nesting. It does not claim to recognize every expression
outside that grammar which happens to denote a Pauli operator.

For other operators, the sufficient rules include:

| Construction | Adjoint-phase rule (all phases modulo 2) |
|---|---|
| `a ⊗ b` | Add the operand phases. |
| `a ⊕ b` | Both blocks must have the same phase; retain it. |
| `adj a` | Negate the phase. |
| `C a` | The operand's phase must be 0 because the other block is identity. |
| `a <> b` | Add the phases and the commutation phase, when commutation or anticommutation is proved. |
| `u <> a <> adj u` | Retain `a`'s phase, for either association. |
| `a <> R k t`, or `R k t <> a` | Retain `a`'s phase if `a` anticommutes with `k`. |

Commutation is proved by exact Pauli arithmetic, equal operands, scalar
operands, or disjoint support. The checker also recognizes inverse pairs,
`SX <> SX`, involutive permutations, integral-angle rotations, and adjacent
rotations about the same syntactic axis whose angles sum to an integer.
Some controls, direct sums, and rotations also normalize to Pauli forms.

```haskell
adjointPhase (X <> Y)                        -- Right (Just 1): anti-Hermitian
typeCheck (R ((X <.> X) <> (Z <.> Z)) (1/3)) -- Right 2
typeCheck (R (X <> R Z (1/3)) (1/5))         -- Right 1
```

Soundness follows by induction: each stored Pauli form denotes the original
operator, and each phase rule preserves `A† = exp(i*pi*p) A`. For example,
if `b a = exp(i*pi*w) a b`, then
`(a b)† = exp(i*pi*(p_a+p_b+w)) a b`. Nonzero-angle rotations are admitted
only after their axes have phase 0; zero-angle rotations are identity.
Thus all checked operators remain unitary.

Outside the complete fragment, `UnprovenHermitian` can still mean a Hermitian
axis needs additional algebra. For example, `(H ⊗ H) <> (Y ⊗ Y)` is Hermitian
but is not currently proved. Analysis is polynomial in the expression size
and rational bit lengths; no matrices or expansions into sums of Pauli
strings are constructed. Even large `Id n` nodes use constant-size Pauli forms.

Call the checker explicitly before evaluation or rewriting. It checks
`QOp` well-formedness; individual backends can have additional restrictions
on supported operators or rotation axes.

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

### Exact cyclotomic evaluation

`QPP.Semantics.Exact.Matrix` constructs exact dense matrices;
`QPP.Semantics.Exact.Statevector` applies the same `QOp` terms directly to
amplitudes. Both use `QPP.Exact.Scalar` and the shared `QPP.Exact.State`
representation. `QOp` and its `Rational` angles are unchanged.

```haskell
import QPP.Syntax (QOp(..))
import QPP.Exact.Types
import qualified QPP.Exact.Scalar as S
import qualified QPP.Semantics.Exact.Statevector as E

output = E.apply (E.evalOp (R H (1/3))) (E.ket [0])

-- Observation and linear algebra retain their checked scalar API.
example = do
  squaredNorm <- E.normSquared defaultLimits output
  pure (E.amplitudes output, S.asRational squaredNorm)  -- norm squared = Just 1
```

Switch the qualified backend import to `QPP.Semantics.Exact.Matrix` to
run the same example through dense evaluation. The ordinary `ket`, `evalOp`,
`apply`, and `evalProg` signatures match the numerical backends and use
`defaultLimits`. `ket` accepts `[Int]` containing only zeroes and ones;
qubit 0 is the most significant bit. Failures raise errors with a call stack.
Use `ketWithLimits`, `evalOpWithLimits`, `applyWithLimits`, and
`evalProgWithLimits` for explicit limits and recoverable failures:
these return `Result a = Either ExactError a`.

`evalProg program input rng` preflights the whole program and supports
`Unitary`, `Measure`, and `Initialize`, returning `(output, outcomes, remainingRng)`.
The checked form is `evalProgWithLimits limits program input rng`.
`evalStep` and `measure1` also have the conventional backend signatures and
checked `WithLimits` variants. Step errors have zero-based indices; standalone
`evalStepWithLimits` uses index zero. An empty program preserves its input and RNG.

Measurement retains the **unnormalized projected vector**. For nonzero input
`v`, outcome `b` has probability `normSquared (P_b v) / normSquared v`, and the
stored output is `P_b v`. Normalizing this vector could leave the cyclotomic
field. `amplitudes`, `inner`, `normSquared`, and `equal` still refer to the
stored vector, not a normalized vector or an equivalence class of states.
Use `E.probability limits state wire outcome` for a conditional probability;
squared amplitudes alone are probabilities only when the stored norm is one.
For example, measuring the Bell state can leave `[h,0,0,0]`, with `h=1/sqrt(2)`
and squared norm `1/2`; measuring the same bit again returns zero with probability one.

Each measurement consumes one finite RNG value in `[0,1)`, even for a certain
outcome. The sampler treats that `Double` as its exact rational value and
compares it with the exact probability using certified interval refinement;
it never rounds a probability to `Double`. Threshold equality selects one.
This retains the other backends' finite-precision RNG convention. `Measure ks`
processes `reverse ks`, prepending outcomes so that the new block is in `ks`
order. `Initialize ks values` measures those wires, applies X where the sample
differs from the target, and hides its measurement outcomes while consuming RNG.
Initialization wires must be distinct and match the number of target values.
Empty wire lists are no-ops. A zero vector cannot be measured. Reset returns one
sampled branch; averaging over its hidden outcomes can produce a mixed state.
Hidden samples affect the stored squared norm; it need not equal the probability
of the recorded outcomes.

Statevector `evalOp` compiles gate constants once; `apply` reuses that plan
or embeds its constants into a larger common field. Both backends share
coefficient-buffer operations while retaining their own quantum shapes.
The scalar-parameterized plan derives folds and traversals for checking and
embedding its constants.

To read the direct evaluator, start with `runPlan` in
[`Exact/Statevector.hs`](src/QPP/Semantics/Exact/Statevector.hs): its cases
give the action of each prepared operator on a vector.
[`Internal/Vector.hs`](src/QPP/Exact/Internal/Vector.hs) implements the array
operations, including the row/column algorithm for tensor application.
[`Internal/Plan.hs`](src/QPP/Exact/Internal/Plan.hs) prepares exact constants,
pushes adjoints to the leaves, and accounts for execution storage.

The symbolic type checker validates each operator, including composite
rotation axes. After full validation, zero-angle rotations become identities
before field selection: unused axes cannot cause field growth. They accept
any well-formed axis, as in `QPP.TypeCheck`. Every global phase is retained,
including zero-qubit
`Phase` terms and identity-axis rotations under control.

The linear API provides `scale`, `add`, `sub`, `tensor`, `inner`,
`normSquared`, and `equal`, along with `fromAmplitudes`, `amplitudes`,
and checked `amplitude` access. Vectors need not be normalized. Matrix
also provides `entry`, `entries`, `fromRows`, `addOp`, `scaleOp`,
`composeOp`, `adjointOp`, `innerOp`, `normSquaredOp`, and `equalOp`.
Its operator inner product is Frobenius; arbitrary matrix linear
combinations need not be unitary. The existing numerical `HilbertSpace`
class is unchanged: these checked operations do not require `Floating`,
`Complex RealT`, square roots, or normalization.

Scalars have private, reduced rational coordinates modulo the cyclotomic
polynomial, with a distinguished complex embedding. Use `S.rational`,
`S.rootOfUnity`, `S.phase`, or `S.fromCoordinates` to construct values;
`S.add`, `S.mul`, `S.inverse`, `S.conjugate`, and `S.equal` perform checked
exact arithmetic. Equality across fields promotes to a common field;
resource exhaustion is an error, not inequality. `S.compareReal` orders real
cyclotomic values with certified bounds and rejects nonreal inputs; exhausted
comparison precision is also a resource error. There is no `Eq`, `Ord`, or
`Floating` instance. `S.approximate` is a lossy observation function and
is never used by evaluation, comparison, or zero tests.
Inversion uses a normalized polynomial Euclidean algorithm: each nonzero
remainder and its Bézout coefficient are divided by the remainder's leading
coefficient together. This curbs avoidable intermediate coefficient growth.

This prototype uses a dense power basis and one field per state/matrix.
Field degree can grow independently of qubit count or entanglement:
rotation denominators 101 and 103 give a sufficient field of degree
40,800. `defaultLimits` bounds root indices at 4096, degree at 256,
conservative coefficient-workspace estimates at 2,000,000 rational slots,
and integer coefficients/intermediates at 4096 bits. These are explicit
representation limits, not wall-time or heap-byte guarantees. Validation
checks retained storage; operations separately reserve their temporary
buffers. Empty programs and same-field equality need no arithmetic workspace.
Direct evaluation counts retained plan constants separately from live vectors:
sequential composition takes the maximum of the child workspace bounds, while
rotations and tensor operations account for saved vectors. Subtraction combines
coordinates directly, without allocating a negated input vector.
Existing field contexts are reused, including when constructed under larger
limits; their retained values must still fit the current limits. Cancellation
may make an estimate larger than the ultimately needed result. Raise
limits deliberately, or use smaller angle denominators. Sparse storage,
cached embeddings, local fields, and field descent remain student-project
directions; the prototype never substitutes a numerical approximation.

Run `cabal run ExactEvaluation` for the phased-Bell example from the paper,
including exact matrix/direct agreement, probabilities `3/4` and `1/4`
after interference, and measurement followed by reset. Run `cabal test ExactTest` for the exact validation suite,
including polynomial identities, larger-field inverses, embedding composition,
and resource-boundary regressions. Random generators shrink failing examples
while preserving field coordinates and quantum shapes.

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
