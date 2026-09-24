module Programs.QFT where
import QPP

-- | Quantum Fourier Transform on n qubits with LSB input and LSB output.
qft :: Int -> QOp
qft n = Permute [n-1, n-2 .. 0] ∘ qftrev n 

-- | QFT with LSB input and MSB output (reverse order)
qftrev :: Int -> QOp
qftrev 0 = One
qftrev 1 = H
qftrev n = (qftrev (n-1) ⊗ I) ∘ layer n


-- | QFT layer acting on qubit n
layer :: Int -> QOp  -- acts on n qubits; targets last qubit
layer n = let
        phases_at_n = foldr (∘) (Id n) [ cP k n | k <- [1..n-1] ] 
        h_at_n    = (Id (n-1)) ⊗ H
    in 
        phases_at_n ∘ h_at_n

-- | Bare Z-rotation by angle 2π/2^k
rz :: Int -> QOp
rz k = R Z (2 / (2^k))

-- | P-gate (relative z-phase change) by angle 2π/2^k
p :: Int -> QOp
p k = Phase (1/2^k) ⊗ rz k

-- | controlled Rz(2π/2^(n-k+1)) with control qubit k and target qubit n. 
cP :: Int -> Int -> QOp
cP k n = (Id (k-1)) ⊗ C ((Id (n-k-1)) ⊗ p (n-k+1)) 

