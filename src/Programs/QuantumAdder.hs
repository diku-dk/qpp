module Programs.QuantumAdder where
import Programs.QFT
import QPP
import Data.Bits (shiftR)

-- | Quantum addition of a classical constant (mod 2^n) using QFT
addConstant :: Int -> Integer -> QOp
addConstant n a =
    let 
        -- |0> + e^{2πi a_k} <1|, where a_k is the k'th binary digit of a
        rzs    = [ R Z (fromIntegral $ a `shiftR` k) | k <- [0..n-1] ] 
        d_a = foldr (⊗) One rzs
    in
        (adj $ qft n) ∘ d_a ∘ (qft n)


addRegister :: Int -> Int -> QOp
addRegister na nx =
    let 
        qftx = qft nx
        -- Control on qubit a_j, rotation targets on all x_t for t in [0..nx-j-1]
        -- -> prefix I^(j-1), and inside control a prefix I^(na-j-1). 
        rotations j      = foldr (⊗) One [ rz (j+t) | t <- [1..nx-j] ] 
        rotation_layer j = (Id $ j-1) ⊗    -- Id prefix before control a_j
                           C((Id (na-j-1)) -- Id prefix inside control (rest of a)
                              ⊗ rotations j) ⊗ -- Rz θ_{j,0} ⊗ Rz θ_{j,1} ⊗ ... ⊗ Rz θ_{j,nx-j-1} 
                           (Id $ j-1)      -- Id suffix after controlled rotations
                                           -- Total length should be na + nx
        -- | Π_{j=0}^{na-1} C( Π_{t=0}^{nx-j-1} (R^{x_t} Z (θ_{j,t})  )
        d_a = foldr (∘) One [ rotation_layer j | j <- [0..na-1] ]
    in 
        (adj qftx) ∘ d_a ∘ qftx

-- | Quantum addition of a classical constant (mod 2^n) using QFT without bit reversal.
addConstant' :: Int -> Integer -> QOp
addConstant' n a =
    let 
        rsz = [ R Z (fromIntegral $ a `shiftR` (n-k)) | k <- [0..n-1] ]
        d_a = foldr (⊗) One rsz
    in
        (adj $ qftrev n) ∘ d_a ∘ (qftrev n)


addRegister' :: Int -> Int -> QOp
addRegister' na nx =
    let 
        qftx = qftrev nx
        -- Control on qubit a_j, rotation targets on all x_t for t in [0..nx-j-1]
        -- -> prefix I^(j-1), and inside control a prefix I^(na-j-1). 
        rotations j      = foldr (⊗) One [ rz (nx-t-j) | t <- [0..nx-j-1] ] 
        rotation_layer j = (Id $ j-1) ⊗    -- Id prefix before control a_j
                           C((Id (na-j-1)) -- Id prefix inside control (rest of a)
                              ⊗ rotations j) ⊗ -- Rz θ_{j,0} ⊗ Rz θ_{j,1} ⊗ ... ⊗ Rz θ_{j,nx-j-1} 
                           (Id $ j-1)      -- Id suffix after controlled rotations
                                           -- Total length should be na + nx
        -- | Π_{j=0}^{na-1} C( Π_{t=0}^{nx-j-1} (R^{x_t} Z (θ_{j,t})  )
        d_a = foldr (∘) One [ rotation_layer j | j <- [0..na-1] ]
    in 
        (adj qftx) ∘ d_a ∘ qftx




