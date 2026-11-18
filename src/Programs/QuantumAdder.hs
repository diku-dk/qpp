module Programs.QuantumAdder where
import Programs.QFT
import QPP

-- | Quantum addition of a classical constant (mod 2^n) using QFT: |x⟩ ↦ |x + a mod 2^n⟩
addConstant :: Nat -> Integer -> QOp
addConstant n a = adj (qft n) ∘ phaseAdd [1..n] a ∘ qft n

-- | Quantum addition of a register (mod 2^nx) using QFT: |a⟩ ⊗ |x⟩ ↦ |a⟩ ⊗ |x + a mod 2^nx⟩
addRegister :: Nat -> Nat -> QOp
addRegister na nx =
    (Id na ⊗ adj (qft nx)) ∘ phaseAddRegister [1..nx] na ∘ (Id na ⊗ qft nx)

-- | Quantum addition of a classical constant (mod 2^n) using QFT without bit reversal.
addConstant' :: Nat -> Integer -> QOp
addConstant' n a = adj (qftrev n) ∘ phaseAdd [n, n-1 .. 1] a ∘ qftrev n

-- | Quantum addition of a register (mod 2^nx) using QFT without bit reversal.
addRegister' :: Nat -> Nat -> QOp
addRegister' na nx =
    (Id na ⊗ adj (qftrev nx)) ∘ phaseAddRegister [nx, nx-1 .. 1] na ∘ (Id na ⊗ qftrev nx)

-- | Addition of a in the Fourier basis, |y⟩ ↦ e^{2πi ay/2^n} |y⟩: one P-gate per qubit,
--   by angle 2πa/2^k for the qubit listed as k, which holds the 2^(n-k) bit of y.
phaseAdd :: [Nat] -> Integer -> QOp
phaseAdd ks a = foldr (⊗) One [ C (Phase (2 * fromIntegral a / 2^k)) | k <- ks ]

-- | Addition of the register a_1 … a_na in the Fourier basis: qubit a_j controls adding 2^(na-j).
phaseAddRegister :: [Nat] -> Nat -> QOp
phaseAddRegister ls na = foldr (∘) (Id (na + length ls))
    [ Id (j-1) ⊗ C (Id (na-j) ⊗ phaseAdd ls (2^(na-j))) | j <- [1..na] ]
