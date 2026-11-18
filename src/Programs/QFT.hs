module Programs.QFT where
import QPP

-- | Quantum Fourier Transform on n qubits with MSB input and MSB output:
--   |j⟩ ↦ 1/√(2^n) Σ_k e^{2πi jk/2^n} |k⟩
qft :: Nat -> QOp
qft n = Permute [n-1, n-2 .. 0] ∘ qftrev n

-- | QFT with MSB input and LSB output (reverse order):
--   |j_1 … j_n⟩ ↦ ⊗_{l=1..n} (|0⟩ + e^{2πi 0.j_l … j_n} |1⟩)/√2
qftrev :: Nat -> QOp
qftrev 0 = One
qftrev n = (I ⊗ qftrev (n-1)) ∘ C (fraction n) ∘ (H ⊗ Id (n-1))

-- | |j_2 … j_n⟩ ↦ e^{2πi 0.0j_2 … j_n} |j_2 … j_n⟩
fraction :: Nat -> QOp
fraction n = foldr (⊗) One [ p k | k <- [2..n] ]

-- | P-gate (relative z-phase change) by angle 2π/2^k
p :: Nat -> QOp
p k = C (Phase (2/2^k))
