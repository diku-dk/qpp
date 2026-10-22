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
| `AlgebraTest` | Simplifier preservation, adjoint distribution, and operator-algebra laws, each checked under all three backends; the angle reduction of `cleanOnes`; the gate identities of `fuseOps` (adjoint pairs, `SX ∘ SX = X` in both spellings, same-axis rotations, permutations), with a soundness property under all three backends; also checks that the shared shrinkers terminate. |
| `TypeCheckTest` | QOp widths, malformed subterms, permutations, exact Pauli-fragment decisions, and conservative adjoint-phase proofs; checks phases and accepted rotation axes against dense matrix semantics, including anticommutation and conjugation. |
| `UOpTest` | General unitary syntax, arbitrary and zero dimensions, list construction, permutation adjoints with unequal block sizes, pretty printing, and symbolic QOp lowering/normalization through `qft 16`. Normalization checks dimensions, idempotence, and phase-preserving vector action. Raw and normalized lowering are also compared with QOp matrices, including rotation signs, half angles, period four, and scalar phases. |
| `TikZTest` | Syntax-tree structure and operand order, exact angles, protected Forest labels, empty and singleton lists, symbolic dimensions, optional qubit counts, and separation of syntax from styling. `make -C latex` additionally compiles the tree gallery. |
| `ExactTest` | Cyclotomic polynomial identities, embedding composition, normalized inversion in larger fields, and conjugation; independent exact gate and DFT values; dense/direct agreement on unnormalized inputs and compiled-plan promotion; composite axes, zero-qubit phases, linear algebra, rewrite phase preservation, program preflight, exact real ordering and threshold ties, unnormalized measurement, repeated measurements and subset resets with conventional outcome/RNG ordering, long-composition workspace, zero-rotation field pruning, and resource boundaries. Generators shrink while preserving valid shapes. No floating-point tolerances are used. |
| `QFTTest` | QFT-specific identities and three-backend agreement on `qft n`. |
| `QCircTest` | Flat circuits (`QPP.QCirc`): fixed gate extractions (time order, wire windows, phases on no wires, `qft 3`), `toQOp` placement by identity padding and by permutation conjugation, round-trip and idempotence properties over `RandomQOp`, agreement of the circuit `Operator` instances with QOp's, and hand-built circuits on shuffled wires checked for well-formedness, unitarity, and an independent CNOT bit-rule specification. Rewrite passes (`flatten`, `compact`, `absorbPermutations`, `fuse`, `collectPhases`, `simplify`): fixed cases, matrix preservation on extracted, hand-built and permutation-heavy circuits, structural post-conditions (elementary gates, full support, no permutation but a trailing one, a single leading phase), the `qft 4` gate count, and idempotence of `simplify`. Diagonal gates: canonical trees from sums, compositions and controlled phases, idle-wire detection, the structural fallback above `maxDiagonalWires`, and fusion across spellings (`Z` with `C (Phase 1)`, a rotation with a phase gate, controlled phases in either wire order) and past other diagonal gates; wire ordering before fusion (controls, tensor factors, rotation axes, permutations on shuffled wires); composite diagonals on three and four shuffled wires (`DiagCirc`): matrix preservation by every pass, no two fused gates on one wire set, idempotence of `simplify`. Layers: ASAP layering keeps every gate and the order on each wire, layers are wire-disjoint (span layers row-disjoint), depth lies between the busiest wire and the gate count, fixed cases including `qft 3`. Diagrams: pinned text and quantikz drawings of controls across idle wires, swaps, a Toffoli with its target between its controls, phase gates, numbered boxes on non-consecutive wires, the one-box fallback, scalars and the empty circuit; shape properties (a row per wire and gap of one width; a header, a layer row and a bare-wire row per span layer in the vertical drawing; a comment per scalar, a row per wire, one quantikz cell per span layer). |

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
| `RandomQCirc` | A circuit on 1–4 wires with up to four gates, each on up to two distinct shuffled wires and a depth-1 payload of matching arity. Shrinks by dropping a gate or shrinking a payload. In `QCircTest`. |
| `PermCirc` | A circuit on 2–4 wires mixing `Permute` gates on shuffled wire subsets with random one-qubit gates. Shrinks by dropping a gate. In `QCircTest`. |
| `DiagCirc` | A circuit on 2–4 wires of composite diagonal gates (`genDiagonal`: phases, `Z`, rotations about `Z`/`Id 1` tensor axes, under controls, sums, products, compositions and adjoints) on two to four shuffled wires each. Shrinks by dropping a gate. In `QCircTest`. |

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

Most properties compare both sides under MS, MPS and SV (`eqMat`), so a
failure names the backend that disagrees. Identity-axis phase regressions
compare the rewritten term in all three backends with the original dense
semantics; this also covers zero-qubit rotations that MPS cannot evaluate
before rewriting.

| Group | Properties |
|---|---|
| Simplifier preservation | `cleanOnes`, `cleanAdjoints`, `liftComposes`, `doComposes`, `pushComposes`, and the combined `cleanop` preserve `evalOp`. Identity-axis rotations retain their phase under tensor, control, and direct sum, including after identity factors merge; `pushComposes` leaves `(a⊗b)∘(c⊗d)` alone when the tensor splits differ. |
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
16. `cleanOnes` discarded the phase of `R (Id n) q`, changing controlled
    and direct-sum operators. It now retains `Phase (-q/2) ⊗ Id n`;
    the previously disabled `cleanop` preservation property is enabled.

## Backlog of invariants to add

Each item is a one-liner property that can be encoded as either an HUnit
case or a QuickCheck property over `RandomQOp`. The earlier backlog items
(simplifier preservation, adjoint distribution, algebraic structure, gate,
permutation, rotation and phase identities) now live in `AlgebraTest`.

**Simplifier fixpoints** — individual rules and `cleanop` are covered;
arbitrary combinations remain.

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
