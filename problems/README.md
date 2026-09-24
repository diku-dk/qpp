# problems/

The weekly problem sheets, as runnable Haskell. Each file is a `main`-carrying
scaffold that compiles and runs *as handed out* — but the interesting functions
are TODO stubs, so the printed answers are wrong until the student fills them
in. That's the point: you run it first, see the nonsense, then fix it.

Shared deps live in the `problems-common` cabal stanza; each sheet is its own
`executable` stanza in `qpp.cabal`, named after the file. Run one with
`cabal run QubitsProblem`, or explore it interactively with
`cabal repl QubitsProblem`.

These sit alongside two related directories: `exercises/` holds the *worked*
demos shown in class (and a few larger test drivers), and `solutions/` holds
reference answers. `problems/` is the unsolved tier.

Every sheet has an answer key. The ten weekly ones live in `solutions/week1/`
through `solutions/week3/`; `NoCloningProblem` and `YesTeleportationProblem` are
answered by `solutions/NoCloning.hs` and `solutions/YesTeleportation.hs`. The
keys check themselves where they can — the QFT key prints its deviation from the
closed-form DFT, the stabilizer key validates its XOR sign rule against the
actual matrix product — so a convention drifting in the library shows up as a
number, not as silent rot.

Keys are deliberately outside the cabal build. Typecheck one with
`cabal exec -- ghc -fno-code -Wall qpp/solutions/week1/QubitsSolution.hs`, run it
with `cabal exec -- runghc <file>`.

## Sheets

| Executable | Sheet | What it exercises |
|---|---|---|
| `QubitsProblem` | week 1, Qubits | Normalization via `inner`, change of basis into {\|+>,\|->}, global-phase invariance of measurement statistics. |
| `PauliBlochProblem` | week 1, Pauli Matrices and the Bloch Sphere | X²=Y²=Z²=I numerically, equality up to global phase, Bloch vectors of the standard states, `R Z θ` tracing a circle. |
| `TensorProductsProblem` | week 1, Tensor Products | Product states via `⊗`, Bell-state preparation, applying a 1-qubit op to one qubit of a pair, entanglement witnessed by partial measurement. |
| `MeasurementProblem` | week 1, Measurement | Born-rule probabilities from `measureProjection`, normalized post-measurement states, partial measurement of an entangled pair, biased coin. |
| `QuantumProgramsProblem` | week 1, Quantum Programs | Composing a whole circuit as one `QOp`, mid-circuit measurement, `DirectSum I X` vs `C X`, building a 1-qubit unitary from `H` and `R Z θ`. |
| `StatePreparationProblem` | week 2, State Preparation and the Input Problem | Basis-state prep from X gates, uniformly controlled R_Y rotations, structural gate counting over the `QOp` AST, the QRAM Ω(N) discussion. |
| `QFTProblem` | week 2, Efficient Quantum Algorithms and the QFT | Controlled phase rotations, QFT₂ given as a warm-up (with a deliberate order/angle question), QFT₃ to build, toy period finding. |
| `AmplitudeAmplificationProblem` | week 2, The Output Problem and Amplitude Amplification | Grover for N=4: oracle, diffusion operator, iterating, checking P(marked) against sin²((2k+1)θ). |
| `GottesmanKnillProblem` | week 2, Classical Simulability and the Limits of Quantum Advantage | Clifford conjugation rules checked numerically, Bell-state stabilizer bookkeeping, T falling outside the Pauli group, statevector vs tableau size. |
| `StabilizersProblem` | week 3, Stabilizers and the Clifford Group | Part A: Pauli codes, the XOR product rule, symplectic product, F₂ rank — plain Haskell, no QPP. Part B: verifying Clifford conjugation against the QPP operators. |
| `NoCloningProblem` | no-cloning | The unsolved form of `exercises/NoCloning.hs`: build CX\|ψ,0> and compare it to \|ψ>⊗\|ψ>. |
| `YesTeleportationProblem` | teleportation | The unsolved form of `exercises/YesTeleportation.hs`: write the two `Unitary` steps and the Bell pair of the teleportation `Program`. |

## Conventions to remind students of

The three that trip people up, all documented in the top-of-file comments:

- **Angles are `Rational` in units of π.** `R Z (1/2)` is a π/2 rotation. Don't
  route angles through `Double`.
- **`Id n` is the n-qubit identity**, with `I = Id 1` and `One = Id 0`. Placeholder
  stubs have to have the right arity or `apply` rejects them — that's why some
  stubs say `Id 2` or `Id 3` rather than `I`.
- **There is no `SWAP`, `S` or `T` constructor.** A swap of adjacent qubits is
  `Permute [1,0]`; S and T are `R Z (1/2)` and `R Z (1/4)`.

Every sheet imports `QPP` (syntax, rewriting, pretty-printing) plus exactly one
backend, `QPP.Semantics.MPS`. Swapping that single import line to
`QPP.Semantics.Matrix` or `QPP.Semantics.Statevector` switches backends.
