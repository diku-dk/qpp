# Changelog

## 0.1.0.0 — QPP (2026-09)

HQPlayground (the ATPL 2025/26 "HQP" framework, v3) becomes QPP, the Quantum
Programming Playground. Same language and semantics; new package name, module
layout, and cabal setup. If you have code written against HQP v3, this is the
map.

### Package and directories

| Before | After |
|---|---|
| `HQPlayground/`, package `hqp`, `library core` | `qpp/`, package `qpp`, default library (`build-depends: qpp`) |
| `HQPlayground/exe/` | `qpp/exercises/` |
| test executables (`cabal run X`) | `test-suite` stanzas: `cabal test all` runs all six; `cabal run X -- <tasty opts>` still works |
| cabal commands from inside the package dir only | root `cabal.project`; commands work from the repo root or `qpp/` |

### Modules

| Before | After |
|---|---|
| `HQP` | `QPP` — re-exports `QPP.Syntax`, `QPP.Rewrite`, `QPP.PrettyPrint` only (no `Util`, no backend) |
| `HQP.QOp` (re-exported Syntax, Simplify, HelperFunctions) | dissolved — `import QPP` covers Syntax and Rewrite; the helpers need an explicit `import QPP.Util` |
| `HQP.QOp.Syntax` | `QPP.Syntax` |
| `HQP.QOp.Simplify` | `QPP.Rewrite` |
| `HQP.QOp.HelperFunctions` | split: see "Symbols" below |
| — | `QPP.Semantics` (new) — the backend *convention* + shared `evalStepWith`/`evalProgWith`. There is no backend class; the `SemanticsBackend` sketch was never in the build and lives in `limbo/Semantics.hs` |
| `HQP.QOp.MatrixSemantics` | `QPP.Semantics.Matrix` |
| `HQP.QOp.StatevectorSemantics` | `QPP.Semantics.Statevector` |
| `HQP.QOp.MPSSemantics` | `QPP.Semantics.MPS` |
| `HQP.QOp.MPS` | `QPP.MPS` |
| `HQP.QOp.DiagMPS` | `QPP.MPS.Diagonal` |
| `HQP.QOp.StabilizerSemantics` | `QPP.Semantics.Stabilizer` (exercise hole) |
| `HQP.PrettyPrint`, `.PrettyOp`, `.PrettyMatrix` | `QPP.PrettyPrint`, `.Op`, `.Matrix` |
| `Programs.*` | unchanged |

### Symbols

| Symbol | Was in | Now in |
|---|---|---|
| `op_dimension`, `op_support`, `step_qubits`, `prog_qubits`, `dagger` | `HelperFunctions` | `QPP.Syntax` |
| `evalStepWith`, `evalProgWith`, `evalStepWithW`, `evalProgWithW` | `HelperFunctions` | `QPP.Semantics` |
| `Convertible`, `SparseMat` | `Syntax` | `QPP.Semantics` |
| `HasWork` | `StatevectorSemantics` | `QPP.Semantics`, generalised with an associated type `Work t` |
| `CMatable` | `MatrixSemantics` | removed — use `Convertible _ CMat` (`to`) |
| `invertPerm` | `HelperFunctions` | `QPP.Util.permInvert` |
| bits / permutations / folds (`toBits`, `permApply`, `foldBalanced`, …) | `HelperFunctions` | `QPP.Util` (hmatrix-free; not re-exported by `QPP`, import it explicitly) |
| `split2x2`, `split2x1`, `split1x2`, `firstBelow` | `HelperFunctions` | `QPP.MPS` |

### Typical import block, before and after

```haskell
-- HQP v3
import HQP
import HQP.QOp.MPSSemantics

-- QPP
import QPP
import QPP.Semantics.MPS
import QPP.Util (toBits)        -- only if you use the helpers
```
