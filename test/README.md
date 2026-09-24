# test/

Validation test suites for QPP. Built with `tasty` + `tasty-hunit` +
`tasty-quickcheck`. Shared deps live in the `test-common` cabal stanza; each
suite is a `test-suite` stanza in `qpp.cabal`.

Run everything with `cabal test all`. Run one suite with `cabal test AlgebraTest`,
or with tasty options via `cabal run AlgebraTest -- -p '/pushComposes/'`.

Backends are abbreviated MS (`QPP.Semantics.Matrix`, the dense reference),
MPS (`QPP.Semantics.MPS`) and SV (`QPP.Semantics.Statevector`), matching the
qualified import aliases used in the test files.

## Test suites

| Suite | What it covers |
|---|---|
| `BackendValidation` | Cross-backend matrix equivalence (MS ↔ MPS ↔ SV), unitarity, adjoint involution, rotation algebra, permutation inverse. Has caught six MPS/SV bugs to date. |
| `MPSTest` | MPS-internal invariants the cross-backend suites can't see: canonical-center validity (isometry check per site) after random ops/measurements, *binding* truncation (`approxCfg` with maxBond below the exact bond dimension), profiling stats, deterministic program-level oracle vs MS on a shared RNG, `R` with `Phase`/`Adjoint`-carrying axes. |
| `MatrixPrepTest` | Pure helpers in `Programs.MatrixPreparation` (`toBitString`, `splitList`, `padToPowerOf2`); `buildRowQOp` row preparation; `matrixPrep` block encoding correctness and unitarity; three-backend agreement on prepared circuits. |
| `SamplingTest` | `sampleAll` / `measure1` / diagonal-MPS sampling: determinism, index mapping, Born-rule agreement, cross-method and cross-backend histograms, truncation safety. |
| `AlgebraTest` | Simplifier preservation, adjoint distribution, and operator-algebra laws, each checked under all three backends; also checks that the shared shrinkers terminate. |
| `QFTTest` | QFT-specific identities and three-backend agreement on `qft n`. |

## Generators (in `Generators.hs`)

| Generator | Produces |
|---|---|
| `RandomQOp` | An `n`-qubit `QOp` with recursion depth ≤ `d` (both ∈ [1,4]). Arity-preserving shrink. |
| `genQOpAt n d` | Same, parameterised — used inside `RandomQOp`'s `arbitrary`. |
| `genQOpBase n` | A leaf op of arity `n`: primitive (X, Y, Z, H, SX, Id), `Phase q`, `Permute π`, or `R ax θ`. |
| `genPauli n` | An `n`-qubit Pauli string (tensor of `Id 1`/`X`/`Y`/`Z`). Used as `R` axis. |
| `genPermutation n` | A permutation of `[0..n-1]`, via QuickCheck's `shuffle`. |
| `genRat` | A rational in `[-7,7]/[1,8]`. |
| `RandomVec`, `RandomMatrix` | Random complex vectors / square matrices (small integer entries). In `MatrixPrepTest`. |

QC budget: `QuickCheckTests 50`, `QuickCheckMaxSize 5` (4 in MatrixPrepTest).
The shrinkers are arity-preserving so the shrunk counterexample is still a
valid input for the property. Every shrink candidate must also be strictly
smaller than its input. Otherwise a failing property never finishes
shrinking and the suite looks hung. AlgebraTest's "Shrinkers terminate"
group checks this against a size measure.

## Properties asserted today

### `BackendValidation`

| Property | Idea |
|---|---|
| `MPS matches MS`, `SV matches MS` | Full-matrix equivalence via `to (MPS.evalOp op) :: CMat` vs `MS.evalOp op` (and likewise SV). |
| `MPS output is unitary` | `U · U† = I` for the materialized MPS matrix. |
| `op ∘ adj op = I (MPS)` | End-to-end inverse. |
| `adj² (op) = op (MPS)` | Adjoint involution. |
| `R ax 0 = Id` | Rotation zero point. |
| `R ax 2θ = (R ax θ)²` | Angle additivity. |
| `R ax θ ∘ R ax (-θ) = I` | Rotation inverse. |
| `Permute π ∘ π⁻¹ = Id` | Permutation inverse. |

### `MatrixPrepTest`

| Property | Idea |
|---|---|
| Block extraction | Upper-left `2^n × 2^n` of `evalOp (matrixPrep M)` equals `M / ‖M‖_F`. |
| Unitarity | `matrixPrep` always produces a unitary. |
| Three-backend match | MS = MPS = SV for `matrixPrep`. |
| Row prep correctness | `buildRowQOp v |0..0⟩ = v / ‖v‖`. |

### `AlgebraTest`

Every property compares the two sides under MS, MPS and SV (`eqMat`), so a
failure names the backend that disagrees.

| Group | Properties |
|---|---|
| Simplifier preservation | `cleanOnes`, `cleanAdjoints`, `liftComposes`, `doComposes`, `pushComposes` each preserve `evalOp`; `pushComposes` leaves `(a⊗b)∘(c⊗d)` alone when the tensor splits differ. |
| Shrinkers terminate | Every `shrinkRat` and `shrinkQOp` candidate is strictly smaller under a size measure and keeps the arity. |
| Adjoint distribution | `adj (a ∘ b) = adj b ∘ adj a`; `adj` distributes over `⊗`, `⊕`, `C`; `adj (R ax θ) = R ax (-θ)` for Pauli axes; `adj (Permute π) = Permute (permInvert π)`; `adj (adj a) = a`. |
| Algebraic structure | `Compose`/`Tensor` associativity; identity laws for `Id n` and `One`; bifunctoriality `(a⊗b)∘(c⊗d) = (a∘c)⊗(b∘d)`; `(a⊕b)∘(c⊕d) = (a∘c)⊕(b∘d)`; `C a = Id n_a ⊕ a`. |
| Pauli / Clifford identities | `X² = Y² = Z² = H² = I`; `(C X)² = I`; `H X H = Z`; `H Z H = X`; `(I⊗H) CX (I⊗H) = CZ`; `R Z θ ∘ X = X ∘ R Z (-θ)`. |
| Permutation algebra | `Permute [0..n-1] = Id n`; `Permute π ∘ Permute σ = Permute (σ∘π)`; conjugating `X` on qubit 0 by `Permute π` moves it to qubit `π⁻¹[0]`. |
| Rotation algebra | `R ax 4 = Id n`; `R ax (θ₁+θ₂) = R ax θ₁ ∘ R ax θ₂`; `Tr(R ax θ) = 2 cos(πθ/2)` for 1-qubit axes. |
| Phase identities | `Phase a ∘ Phase b = Phase (a+b)`; `Phase 0 = Phase 2 = Id 0`. |
| Bifunctoriality bisect | HUnit cases pinning down the minimal `C(P)∘P` trigger of a former MPS/MS disagreement. |

### `QFTTest`

| Property | Idea |
|---|---|
| `qft 1 = H`, `qft 2 = DFT_4`, `qft 3 = DFT_8` | Fixed cases against explicit matrices. |
| `qft n = analytic DFT` | `ω^{jk}/√N` entrywise. |
| `qft n |0⟩ = uniform` | Uniform superposition from the zero ket. |
| `qft n is unitary`, `qft n ∘ adj (qft n) = I` | Unitarity two ways. |
| `(qft n)⁴ = I`, `(qft n)² = index reversal` | Order-4 and the bit-reversal permutation. |
| `MS / MPS / SV agree on qft` | Three-backend agreement. |

`MPSTest` and `SamplingTest` document their own properties in their module
headers.

## Bugs caught so far

The current property suite has surfaced (and led to fixes for):

1. `permutationSwaps`'s `pos` array initialised as `p0` instead of `p0⁻¹`
   (only manifested on non-self-inverse permutations like 3-cycles).
2. `R axis θ` single-qubit case applied the bare Pauli, ignoring θ.
3. `R` multi-qubit sign error (`+i·s·P` instead of `−i·s·P`).
4. `SX` matrix `(p,p,m,p)` — non-unitary typo.
5. `dagger SX = Adjoint SX` infinite recursion.
6. `supportInterval` empty-support fallback out-of-bounds.
7. `MatrixPreparation.calculate{Real,Complex}Gate` sign convention.
8. `createRotationsDS` left-leaning DirectSum tree (dim-mismatch on
   non-power-of-2 chunk counts).
9. `op_support` for `Compose (Permute ks) b` returning ∅ when `b` had
   empty support.
10. `svd_compact` returned the *untruncated* singular values next to
    truncated U/V — any binding `maxBond` crashed on a dimension mismatch
    (the truncating branch had never been exercised).
11. Stale `center_site` bookkeeping (`addLocal` windows not covering the
    inputs' centers; `swapSites` and `tensorMPS` not maintaining the
    center) — wrong *measurement probabilities* while the state vector
    itself stayed correct, so `BackendValidation` couldn't see it.
12. `R` with a phase-carrying axis (`Phase q ⊗ P`, `q` non-integer) used
    the two-term formula `c·I − i·s·φP`, which equals the matrix
    exponential only when `φ² = 1`; MS and MPS disagreed.
13. 1-qubit `R` rotations went through the LCU-plus-SVD path, breaking
    closure of the diagonal-core MPS class under 1-qubit gates
    (SamplingTest's "Closure under 1q + Permute + Phase").
14. Truncation statistics stored inside the state's `cfg` dropped one
    branch's counters at every `addLocal` merge — `C a` reported *fewer*
    SVDs than evaluating `a` alone. Fixed by moving stats to a `Profile`
    monoid emitted per-SVD and threaded through the evaluation
    (`applyProfiled` / `evalProgProfiled`; see
    docs/MPS-profiling-redesign.md).
15. `pushComposes` fused `(a⊗b)∘(c⊗d)` into `(a∘c)⊗(b∘d)` even when `a`
    and `c` act on different numbers of qubits, producing an ill-formed
    term. It showed up as a hang rather than a failure because the
    shrinkers never terminated.

## Backlog of invariants to add

Each item is a one-liner property that can be encoded as either an HUnit
case or a QuickCheck property over `RandomQOp`. The earlier backlog items
(simplifier preservation, adjoint distribution, algebraic structure, gate,
permutation, rotation and phase identities) now live in `AlgebraTest`.

**Simplifier fixpoints** — the individual rules are covered; the composites
are not.

- `evalOp (cleanop op) ≈ evalOp op` (fixpoint composition). Written but
  disabled in `AlgebraTest`: it fails because `cleanOnes` rewrites
  `R (Id n) θ` to `Id n` and drops the phase.
- `evalOp (simplifyFixpoint n rules op) ≈ evalOp op` for any rule list.

**Norm and trace preservation**

- `inner (apply U ψ) (apply U ψ) = inner ψ ψ` for unitary U
- `Tr(U M U†) = Tr(M)`
- `|det U| = 1` for unitary U

**Under `Truncate` cfg in MPS** — now covered by `MPSTest` (bond cap,
profile population, loose-truncation agreement with the dense reference).

**Apply-path (state) properties** — useful once we want to test n > 4
qubits where full-matrix materialisation is expensive.

- `apply (Tensor a b) (ket bs) = (apply a (ket b_a)) ⊗ (apply b (ket b_b))`
- `apply (Compose a b) ψ = apply a (apply b ψ)`
- Probabilities sum to 1 after a Born-rule measurement.

## Adding a property

1. Decide where it lives: cross-backend / general → `BackendValidation`;
   algebraic identities and rewrite rules → `AlgebraTest`; matrix-prep
   specific → `MatrixPrepTest`; MPS internals → `MPSTest`; measurement and
   sampling → `SamplingTest`; QFT → `QFTTest`.
2. Write the property as `forAll genQOpAt …` or against `RandomQOp` /
   `RandomMatrix` / `RandomVec` from the generators above.
3. Match matrices via the `(~~)` operator (L1 diff < `tol`), or `eqMat` in
   `AlgebraTest` to check under all three backends at once.
4. Register in the suite's `defaultMain` group.
