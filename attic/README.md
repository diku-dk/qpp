# attic/

Retired code, kept for reference. Nothing here is built, and nothing in the
library depends on it.

## `StatevectorMassiv.hs`

`QPP.Semantics.Statevector` as it was until 2026-11-12: the state-vector
backend on `massiv` arrays, delayed while computing and unboxed at rest. The
file is unchanged, its module name included. The module of the same name in
`src/QPP/Semantics/` replaced it.

| | `StatevectorMassiv.hs` | its replacement |
|---|---|---|
| Size | 621 lines | 176 lines |
| Arrays | `massiv` | `hmatrix`; `StateT = CMat`, as in the matrix backend |
| Rotation axes | Pauli strings (with `Phase`, `Permute`, `∘`, `Adjoint`); others raise an error | any Hermitian operator |
| Tests | the cross-backend suites | those, and `StatevectorTest` |

### To build it again

1. Copy it to `src/QPP/Semantics/StatevectorMassiv.hs`, change its module
   name to `QPP.Semantics.StatevectorMassiv`, and list it in `qpp.cabal`.
2. Add `massiv` to the library's `build-depends`.
3. Restore the class it used. It was removed from `QPP.Semantics` together
   with this backend, its only user:

   ```haskell
   class HasWork t where
     type Work t
     toWork   :: t -> Work t
     fromWork :: Work t -> t
   ```

## `StatevectorBench.hs`

The comparison that was run before the retirement. Its header says how to
build it, and `statevector-bench-2026-11-12.log` is its raw output. Where
both backends finished, their results agreed: the amplitudes up to rounding
(compared up to 14 qubits) and the outcomes of the programs.

### Results, 2026-11-12

Seconds to apply a circuit to |0…0⟩ and force every amplitude. One core
(`+RTS -N1`), Apple M3 Max, GHC 9.6.7, `-O1`. Times under 2 s are the best
of three runs. "> 60" is a run stopped after 60 s.

| Circuit | Qubits | massiv | new |
|---|---:|---:|---:|
| QFT | 12 | 0.078 | 0.0044 |
| QFT | 16 | 2.1 | 0.18 |
| QFT | 20 | 53 | 5.9 |
| QFT | 22 | > 60 | 47 |
| GHZ preparation | 16 | 0.17 | 0.017 |
| GHZ preparation | 20 | 3.3 | 0.47 |
| GHZ preparation | 22 | 14 | 3.5 |
| 3 Trotter steps of an Ising chain | 4 | 9.0 | 0.0001 |
| 3 Trotter steps of an Ising chain | 6 | > 60 | 0.0002 |
| 3 Trotter steps of an Ising chain | 20 | > 60 | 4.5 |
| `H ⊗ … ⊗ H` | 12 | 1.3 | 0.0009 |
| `H ⊗ … ⊗ H` | 14 | 21 | 0.0036 |
| `H ⊗ … ⊗ H` | 22 | > 60 | 2.3 |
| one `H` per qubit, composed | 12 | 3.0 | 0.0010 |
| one `H` per qubit, composed | 14 | 48 | 0.0042 |
| one `H` per qubit, composed | 22 | > 60 | 3.4 |
| one `R X` per qubit, composed | 12 | 3.0 | 0.0011 |
| one `R X` per qubit, composed | 14 | 49 | 0.0049 |
| one `R X` per qubit, composed | 22 | > 60 | 4.1 |
| 18 rotations `R Z` in a row on one qubit | 10 | 28 | 0.0005 |
| bit reversal, one `Permute` | 20 | 0.60 | 0.22 |
| bit reversal, one `Permute` | 24 | 16 | 4.3 |
| `X` under n−1 controls, after `H` on every qubit | 12 | 1.3 | 0.0009 |
| `X` under n−1 controls, after `H` on every qubit | 20 | > 60 | 0.40 |
| program: GHZ, `Measure` all, `Initialize` all | 16 | 0.33 | 0.14 |
| program: GHZ, `Measure` all, `Initialize` all | 20 | 6.3 | 3.6 |

Many small operators, as the tests and exercises use them. Total seconds for
random operators of depth at most 4, each applied to |0…0⟩:

| Qubits | Operators | massiv | new | matrix backend |
|---:|---:|---:|---:|---:|
| 1 | 2000 | 0.038 | 0.0064 | 0.013 |
| 3 | 2000 | 0.045 | 0.020 | 0.060 |
| 5 | 1000 | 0.031 | 0.019 | 0.15 |
| 8 | 200 | 0.030 | 0.0098 | 1.6 |

Converting such operators to their matrices with `to (evalOp op)`, as the
test suites do, took 0.30 s against 0.032 s for 100 operators on 6 qubits.

The massiv module works on delayed arrays, which are recomputed at every
use. Its time grows by a factor of 16 for every two more qubits in a layer of
one-qubit gates, and doubles with every further rotation in a row. The new
module was faster on every case measured.

More cores did not help the massiv module. With `+RTS -N14`, QFT on 18
qubits took 14 s, and QFT on 20 qubits did not finish in 60 s.

Peak memory of one run, as "total memory in use" from `+RTS -s`:

| Circuit | Qubits | massiv | new |
|---|---:|---:|---:|
| QFT | 16 | 230 MiB | 21 MiB |
| QFT | 20 | 3192 MiB | 109 MiB |
| GHZ preparation | 22 | 1464 MiB | 462 MiB |
| bit reversal, one `Permute` | 20 | 42 MiB | 66 MiB |

The bit reversal is the one case in which the new module used more.
