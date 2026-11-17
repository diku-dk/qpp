# Changelog

## Unreleased

### The state-vector backend is replaced

`QPP.Semantics.Statevector` is a new module with the same backend names
(`StateT`, `OpT`, `ket`, `evalOp`, `apply`, `measure1`, `evalStep`,
`evalProg`). It defines the semantics by one rule per constructor of `QOp`,
on `hmatrix`, in 181 lines against the 621 of the `massiv` module. Programs
that only use the backend names need no change.

| What | Before | Now |
|---|---|---|
| `StateT` | a `massiv` array | `CMat`, the state type of `QPP.Semantics.Matrix` |
| Rotation axes | Pauli strings (with `Phase`, `Permute`, `∘`, `Adjoint`) | any Hermitian operator |
| `OpT` | a closure with `n_qubits` only | also an `Operator` and a `HilbertSpace`: `<>`, `∘`, `⊗`, `⊕`, `adj`, `.*`, `.+`, `.-`, `inner`, `norm`. `adj`, `inner` and `norm` form the operator's matrix |
| The `massiv` module | `src/QPP/Semantics/Statevector.hs` | `attic/StatevectorMassiv.hs`; not built |
| `HasWork`, `toWork`, `fromWork` | in `QPP.Semantics`; used only by the `massiv` module | removed |
| `massiv` | a dependency of the library | not needed |
| Tests | cross-backend suites | those, and the new suite `StatevectorTest` |

`applyOp e ψ` takes a state ψ of exactly the qubits of `e`, like `apply`;
another size is an error. Until 2026-11-17 it also accepted a state of more
qubits, acted on the first ones and carried the rest along. Write `e ⊗ Id p`
for that now. The columns of ψ may still be several states side by side, and
`applyOp e ψ` is then the matrix product `M(e) ψ`.

### Rotations have one formula in every backend

`R a θ` evaluates to `cos t·I − i sin t·A`, with `t = πθ/2` and `A` the
matrix of the axis `a`. For a Hermitian axis, such as a Pauli string or `H`,
`A² = I` and this is `exp(−it·A)`, as before.

| Backend | Before | Now |
|---|---|---|
| `QPP.Semantics.Matrix` | the matrix exponential `exp(−it·A)` | `cos t·I − i sin t·A` |
| `QPP.Semantics.MPS` | for an axis `φ·P`, with `P` a Pauli string: `cos(tφ)·I − i sin(tφ)·P`, which is `exp(−it·A)` | `cos t·I − i sin t·φP` |

The results differ only for an axis that is not Hermitian, such as
`Phase q ⊗ X` with `q` not an integer. On such an axis the matrix and MPS
backends gave `exp(−it·A)`, and the state-vector backend gave the formula,
which is then a different matrix. All three now give the formula.

### Classical bits have a type of their own

`Bit` is the type of a classical bit. It is defined in the new module
`QPP.Bit` and comes with `import QPP`. It replaces `Int` in `ket` and `Bool`
in `Initialize` and in measurement outcomes. A bit is written `0` or `1` and
prints that way. `+` and `*` on bits are arithmetic modulo 2, so `+` is
exclusive or.

| What | Before | Now |
|---|---|---|
| `ket` | `[Int] -> StateT`; `ket [0,2]` gave a wrong vector | `[Bit] -> StateT`; `ket [0,0,1]` is unchanged, `ket [0,2]` is the error "2 is not a bit" |
| `Initialize` | `Initialize [0,2] [False,True]` | `Initialize [0,2] [0,1]` |
| `Outcomes` | `[Bool]`, printed as `[False,True]` | `[Bit]`, printed as `[0,1]` |
| A test on an outcome | `if b then …`, `bits == [True]` | `if b == 1 then …`, `bits == [1]` |
| Exclusive or of two bits | `xor a b`, `a /= b` | `a + b` |
| `measureProjection n k v` | `v :: Int` | `v :: Bit`; `measureProjection n k 0` is unchanged |
| `toBits`, `toBits'`, `fromBits` in `QPP.Util` | lists of `Int` | lists of `Bit` |
| Tests | | the new suite `BitTest` |

To convert: `Bit (x > y)` makes a bit from a `Bool`, and the pattern `Bit b`
gives the `Bool` back. `fromEnum b` is the `Int` 0 or 1, and `fromIntegral n`
makes a bit from an integer that is 0 or 1. `Bit` is an instance of `NFData`,
so `force` and `deepseq` work on outcomes as before; the library now depends
on `deepseq` for that.

The sum of a list of bits is its parity. To count the ones, write
`length (filter (== 1) bits)`. GHC does not know that 0 and 1 are all the
bits, so a function defined by the two cases `f 0` and `f 1` is reported as
incomplete: write the second case as `f _`.

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
