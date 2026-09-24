-- Exercise scaffold: Data.BitVector is unused until task 4 is done.
{-# OPTIONS_GHC -Wno-unused-imports #-}
module Main where
import QPP
import QPP.Semantics (tol)
import QPP.Semantics.Matrix
import Numeric.LinearAlgebra(find,scale,magnitude,atIndex,Complex((:+)))
import Data.BitVector

{-|
 1. Write a function to check if a given QOp is in the Pauli group (start with 1-qubit operators, then write recursive n-qubit test).

    You can use the following facts about the Pauli group P_n to help you:

    The tensor product         A ⊗ B ∈ P_{m+n} iff A ∈ P_m and B ∈ P_n. 
    The sequential composition A ∘ B ∈ P_n iff A,B ∈ P_n.
    The direct sum           A <+> B ∈ P_{n+1} iff A ∈ P_n and B==+/-A 
                                              (A<+>A == I⊗A, A<+>(-A) == Z ⊗ A).
                            Special case: C A ∈ P_{n+1} iff A==+/-I_n
    P_n is closed under adjoint.
    
    2. Write a function to check if a given QOp is a Clifford operator (start with 1Q operators, then write recursive n-qubit test).

    You can use the following facts about the Clifford group C_n to help you:

    The tensor product         A ⊗ B ∈ C_{m+n} iff A ∈ C_m and B ∈ C_n.
    The sequential composition A ∘ B ∈ C_n iff A,B ∈ C_n.
    The direct sum           A <+> B ∈ C_{n+1} iff A,B ∈ C_n 
                                              and A ∘ adj B ∈ P_n
                            Special case: (C A) = I<+>A ∈ C_{n+1} iff A ∈ P_n
    C_n is closed under adjoint.
    
    3. Write a function that takes a 1Q Clifford operator G and computes the effect of conjugation on the Paulis: GXG^\dagger, GZG^\dagger. Do you need GYG^\dagger as well? (Why/why not?)
    4. Generalize to n-qubit Clifford operators. Use BitVector (from Data.BitVector) to represent the result. 

    5. Use this to implement a stabilizer simulator that tracks the evolution of stabilizers under Clifford operations.
-}

{-| Divide out global phase from a matrix such that the first non-zero entry is real and positive.
    This allows us to compare matrices up to global phase. -}
factor_phase :: CMat -> CMat
factor_phase m = case first_nonzero m of
    Nothing -> m
    Just x  -> let
        phase = x / (magnitude x :+ 0)
      in
        scale (1 / phase) m

{-| Find the first non-zero entry in a matrix. We use tol to allow for floating point inaccuracy. -}
first_nonzero :: CMat -> Maybe ComplexT
first_nonzero m = 
    case find (\x -> magnitude x > tol) m of
        []     -> Nothing
        (ij:_) -> Just (m `atIndex` ij)

main :: IO ()
main = putStrLn "Stabilizers exercise: implement the functions described in the comment block at the top of this file."
