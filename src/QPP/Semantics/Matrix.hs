{-# LANGUAGE TypeFamilies #-}

module QPP.Semantics.Matrix where

import QPP.Syntax
import QPP.PrettyPrint.Op
import QPP.Util
import QPP.Semantics
import Numeric.LinearAlgebra hiding(normalize,step,(<>)) -- the hmatrix library
import Data.Bits(shiftL)

type CMat = Matrix ComplexT
type CVec = Vector ComplexT
type RMat = Matrix RealT
type StateT = CMat
type OpT    = CMat

{-| 
 evalOp :: Op -> CMat
 (evalOp op) evaluates the pure quantum operator op (defined in Circuit.hs) in the matrix semantics, and produces a complex matrix of dimension 2^n x 2^n (where n is the qubit-length of op)
-}

apply :: OpT -> StateT -> StateT
apply op state = op <> state

measure1 :: (CMat, Outcomes, RNG) -> Nat -> (CMat, Outcomes, RNG)
measure1 (state, outcomes, (r:rng)) k = let
      n = ilog2 (rows state)
      proj0 = measureProjection n k 0
      proj1 = measureProjection n k 1

      s0 = proj0 <> state
      s1 = proj1 <> state

      prob0 = realPart $ inner state s0
      prob1 = realPart $ inner state s1

      outcome = if (r < prob0) then False else True
      collapsed_state = normalize $ if(outcome) then s1 else s0

      in
          if (abs(prob1+prob0-1)>tol) then
              error $ "Probabilities don't sum to 1: " ++ (show (prob0,prob1))
          else
              (collapsed_state, outcome:outcomes, rng)

measure1 (_,_,[]) _ = error "No more random numbers. This never happens."

  {-| bra [v0,...,v_{n-1}] (where v_j is 0 or 1) is the adjoint state <v0 v1 ... v_{n-1}|. In the matrix representation, this is a row-vector in C^{2^n}, i.e. of dimension (1 >< 2^n). -}
ket :: [Int] -> CMat
ket []     = (1><1) [1]
ket (v':vs) = let 
      v = fromIntegral v' :+ 0
  in 
      (2><1) [1-v,v] ⊗ (ket vs)

evalOp :: QOp -> CMat
evalOp op = case op of
  Id n -> ident (2^n)

  Phase q -> let theta = (realToFrac q * pi) :+ 0 
             in  scalar (exp (ii*theta))

  X -> (2 >< 2) [0,1,
                 1,0]

  Y -> (2 >< 2) [0, -ii,
                 ii, 0]

  Z -> (2 >< 2) [1,0,
                 0,-1]

  H -> let s = 1/sqrt 2
       in  s * ((2><2) [1, 1,
                       1,-1])

  SX -> let p = 0.5:+0.5
            m = 0.5:+(-0.5)
        in (2><2) [p, m, m, p]   -- IBM-convention SX; matches MPS/SV. SX² = X.

  R axis q -> 
      let mat  = evalOp axis
          theta = (realToFrac q * pi) :+ 0 -- Haskell won't multiply real and complex numbers
      in
          matFunc exp ( (-ii*theta/2) .* mat )

  Permute ks -> let
        n    = length ks
        dims = 1 `shiftL` n

        outIndex i = fromBits $ map (toBits' n i !!) ks   -- output bits in order [b_{ks0},..,b_{ks(n-1)}]

        indices = [ (outIndex i, i) | i <- [0..dims-1] ]
    in assoc (dims,dims) (0 :+ 0) (zip indices (replicate dims (1 :+ 0)))
                                          
  C op1          ->  let 
                          mop = evalOp op1
                          mI  = ident (rows mop)
                      in
                          mI <+> mop -- |0><0| ⊗ I^n + |1><1| ⊗ op1
  
  Tensor    op1 op2 -> (evalOp op1)  ⊗  (evalOp op2)
  Compose   op1 op2 | (op_qubits op1 == op_qubits op2) -> (evalOp op1)  ∘  (evalOp op2)  
                    | otherwise -> error $ 
                     "\n\nDim-mismatch: " ++ showOp op1 ++ " ∘ " ++ showOp op2 ++ "\n"
  DirectSum op1 op2 | (op_qubits op1 == op_qubits op2) -> (evalOp op1)  <+> (evalOp op2)  
                    | otherwise -> error $ 
                     "\n\nDim-mismatch: " ++ showOp op1 ++ "<+>" ++ showOp op2 ++ "\n"
  Adjoint op1       -> adj $ evalOp op1


-- | Unitary application with a friendly arity check (the raw matrix product
--   would fail with an opaque hmatrix dimension error).
applyUnitary :: QOp -> StateT -> StateT
applyUnitary op st
  | n_qubits op == n_qubits st = apply (evalOp op) st
  | otherwise = error $ "Dim-mismatch between " ++ showOp op ++ " and n=" ++ show (n_qubits st)

evalStep :: (StateT, Outcomes, RNG) -> Step -> (StateT, Outcomes, RNG)
evalStep = evalStepWith applyUnitary measure1

  -- | evalProg simply executes a program step by step
evalProg :: Program -> StateT -> RNG -> (StateT, Outcomes, RNG)
evalProg = evalProgWith applyUnitary measure1



{-| measure1 measures a single qubit, projecting the state on the subspace corresponding to the 
    outcome, and prepends the outcome to the list of outcomes -}
 

{-|
 - The tensor product of two operators in matrix representation a :: (m><n) and b :: (p><q) is the kronecker product c :: (m*p >< n*q) with entries c_{i*p+k, j*q+l} = a_{ij}*b_{kl} 

 - The direct sum of two operators in matrix representation a :: (m><n) and b :: (p><q) is the block-diagonal matrix c :: ((m+p) >< (n+q)) with a and b as blocks. Note that the additive dimension necessitates that a and b have the same dimensions for it to be realizable on qubits. 

 - The operator adjoint in matrix representation is the conjugate transpose.
-}    

instance HasTensorProduct CMat where (⊗) = kronecker
instance HasAdjoint CMat where adj = tr -- HMatrix confusingly defines conjugate transpose as 'tr' 
instance HasQubits CMat where  n_qubits op = ilog2 (rows op)
instance HasDirectSum CMat where
    (⊕) a b = fromBlocks [[a, zeros (rows a) (cols b)],
                          [zeros (rows b) (cols a), b]]

instance Operator CMat

{-| We implement ket-states as (n >< 1) matrices (= column vectors), and bra-states as (1 >< n) matrices (= row vectors). Thus StateT inherits the operator operations (recall that states
are also linear operators by Riesz representation), but we can also treat them like vectors,
so we implement the VectorSpace operations. -}
instance HilbertSpace CMat where
    type Scalar  CMat = ComplexT 
    type Realnum CMat = RealT
    
    (.*) = scale
    (.+) a b = a+b
    (.-) a b = a-b

    -- | inner a b is the usual dot product, with the adjoint coefficients complex conjugated
    inner a b = let 
            (m,n,p,q) = (rows a, cols a, rows b, cols b) 
        in
            if (m,n) == (p,q) then 
                sumElements (conj a * b)        
            else 
                error $ "inner: incompatible shapes " ++ (show (m,n)) ++ " != " ++ (show (p,q))

    normalize :: CMat -> CMat
    normalize s = let 
          n = (norm_2 (flatten s))
        in 
          if (n < tol) then s 
          else              s / scalar (n :+ 0)


{-|
 measureProjection n k result 
 
Builds the projection operator onto the subspace of C^{2^n} in which qubit k has value 'result'.
 
Construction of the operator:
I ⨷ ... ⨷ I ⨷ P ⨷ I ⨷ ... ⨷ I 
\-- k-1 --/                   /
 \------------ n ------------/
-}
measureProjection :: Int -> Int -> Int -> CMat
measureProjection n k v' = let     
        v = fromIntegral v' :+ 0
        p = (2><2) [1-v,0,
                    0,  v]
    in
        evalOp(Id k) ⊗ p ⊗ evalOp (Id (n-k-1))

-- Auxiliary definitions -- move to internal module?
-- (The numerical tolerance `tol` is shared by all backends and lives in QPP.Semantics.)
ii, one :: ComplexT
ii  = 0 :+ 1 
one = 1     

zeros :: Int -> Int -> CMat
zeros m n = konst 0 (m,n)

        


-- apply :: Op -> Vector -> Vector -- How do define abstract state instead, for easy backend switching?

realM, imagM :: CMat -> RMat
realM m = cmap realPart m
imagM m = cmap imagPart m

realMI, imagMI :: CMat -> [[Int]]
realMI = map (map (floor.realPart)) . toLists 
imagMI = map (map (floor.imagPart)) . toLists 

{-| sparseMat takes an (m >< n) matrix and returns ((m,n), nonzeros) where nonzeros is a list of every nonzero index paired with the corresponding value. -}
sparseMat :: CMat -> SparseMat
sparseMat mat = 
    let 
        (m,n) = (rows mat, cols mat) 
        full_list = case (m,n) of
            (1,1) -> []
            (_,1) -> [((fromIntegral i,0), mat `atIndex` (i,0)) | i <- [0..m-1]]
            (1,_) -> [((0,fromIntegral j), mat `atIndex` (0,j)) | j <- [0..n-1]]
            _     -> error $ show "Use sparseOp for operators"
    in
        SparseMat ((fromIntegral m, fromIntegral n), filter (\(_,v) -> (magnitude v > tol)) full_list)


instance Convertible CMat SparseMat where
    to mat = sparseMat mat

    from _ = error "CMat from SparseMat: not implemented yet."

-- | Identity Convertible instance. Exists so backend-agnostic user code can write
--   `to (evalOp op) :: CMat` uniformly across backends — in this backend OpT is already
--   CMat, so this is a no-op; in lambda backends (MPS, Statevector) the same call dispatches
--   to a canonical-basis materialization. `Convertible t CMat` is also what the dense
--   pretty-printers (`showCMat`, `printM`) require of a state.
instance Convertible CMat CMat where
    to   = id
    from = id
        