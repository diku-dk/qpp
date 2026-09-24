module PI where

import QPP hiding (fixpoint)
import Programs.QFT
import qualified QPP.Semantics.Matrix as MatSem
import Numeric.LinearAlgebra(rows,norm_2)
import Data.List (sort)


{-| 
   In src/Programs/QFT.hs we have implemented the Quantum Fourier Transform (QFT) as a quantum program. First have a look at the QFT implementation there.

   In this exercise, we will get comfortable working with programs as data structures that we can manipulate symbolically, performing program transformations that preserve their semantics.

   1. First generate the QFT program on 1,2, and 3 qubits, and inspect it using both show and the pretty printer. 
   -}
[qft1, qft2, qft3] = [qft n | n <- [1,2,3]]
   
pp    = putStrLn . showOp   
 {-| Look at the ASTs interactively with cabal repl ProgramInversion. -}



{-|
   2. Evaluate the QFT program on 1, 2, 4, 10 qubits using the matrix semantics from MatrixSemantics.hs. Verify that the resulting matrices are unitary. What would happen if you tried this with 20 qubits? Use 'cabal repl ProgramInversion' to explore this interactively.
   -}
[mqft1,mqft2,mqft3, mqft4, mqft10] = map (MatSem.evalOp . qft) [1,2,3,4,10]

show2 :: IO ()
show2 = do
  printM (mqft1 <> adj mqft1) -- should be identity
  printM (mqft2 <> adj mqft2) -- should be identity
  printM (mqft4 <> adj mqft4) -- should be identity
  putStrLn $ ("QFT10 matrix is " ++ show (rows mqft10)) ++ " x " ++ show (rows mqft10)

{-| Evaluating the QFT on 20 qubits results in a 1e6 x 1e6 matrix, which would require a 8 petabytes of memory. -}

{-|
   Look in src/QPP/Rewrite.hs for an example of a program transformation that simplifies QOp syntax trees by removing One operators. You can use the cleanOnes function from there to simplify your QFT programs before printing or evaluating them. Check that the simplification preserves the semantics of the program.

   3. Now write a function (QOp -> QOp) that performs the inversion of a quantum program. I.e., given Adjoint (qft n), this rewrite function should actually perform the Adjoint operation recursively on all sub-operators, resulting in a new QOp that is the inverse of the original QFT program. 

   To achieve this: 1) Find out what the adjoints of primitive operators I, X, Y, X, Rz(Phi) are.  Can they be expressed as operators in our operator language?  2) What are the rules for what the adjoints of compositions, tensor products, direct sums and permutations are?  
   
   Every operator term in our operator language has an equivalent operator term without any occurrences of the  Adjoint constructor.
-}
permInvert :: [Int] -> [Int]
permInvert ks = map snd $  -- For each index in the output, find its position in the input
    sort [ (k, i) | (i, k) <- zip [0..] ks ]

simplifyAdjoints :: QOp -> QOp
simplifyAdjoints op = case op of
    Adjoint One                -> One
    Adjoint I                 -> I
    Adjoint X                 -> X
    Adjoint Y                 -> Y
    Adjoint Z                 -> Z
    Adjoint H                 -> H
    Adjoint (R a theta)       -> R (simplifyAdjoints a) (-theta)
    Adjoint (C a)             -> C (simplifyAdjoints (Adjoint a))
    Adjoint (Adjoint a)       -> simplifyAdjoints a
    Adjoint (Permute ks)    -> Permute (permInvert ks)
    --(AB)^-1 = B^-1 A^-1
    Adjoint (Compose a b)     -> Compose (simplifyAdjoints (Adjoint b)) (simplifyAdjoints (Adjoint a))
    Adjoint (Tensor a b)      -> Tensor  (simplifyAdjoints (Adjoint a)) (simplifyAdjoints (Adjoint b))
    Adjoint (DirectSum a b)   -> DirectSum (simplifyAdjoints (Adjoint a)) (simplifyAdjoints (Adjoint b))
      
      -- No more rewrite rules. Recurse over non-leaf constructors.
    Tensor  a b               -> Tensor    (simplifyAdjoints a) (simplifyAdjoints b)
    DirectSum a b             -> DirectSum (simplifyAdjoints a) (simplifyAdjoints b)
    Compose a b               -> Compose   (simplifyAdjoints a) (simplifyAdjoints b)
    C a                       -> C         (simplifyAdjoints a)
    R a phi                   -> R         (simplifyAdjoints a) phi
    _                         -> op    

{-|    
    4. Use this function to construct the inverse QFT program by applying it to the QFT program on n qubits. Verify that the inverse QFT program is indeed the inverse by evaluating the composition of QFT and its inverse.
-}
invqft n = cleanOnes $ Adjoint (qft n)
[siqft1, siqft2, siqft3]    =  [simplifyAdjoints (invqft n) | n <- [1,2,3]]
[mi1, mi2, mi3]             =  [MatSem.evalOp    (invqft n) | n <- [1,2,3]]
[msi1, msi2, msi3]          =  [MatSem.evalOp (simplifyAdjoints (invqft n)) | n <- [1,2,3]]

show4 :: IO ()
show4 = do
   putStrLn "Check that siqft3 has no Adjoint constructors:"
   pp siqft3

   putStrLn "Is the matrix representation of Adjoint(qft 3) equal to its simplified version?"
   putStrLn $ "||mi3 - msi3|| = " ++ show (norm_2 (mi3 - msi3))  -- should be 0
   
   putStrLn "Is the composition of qft 3 and its inverse equal to identity?"
   printM $ msi3 <> mqft3  -- should be identity
   
      
{-|
    5. Write a function that takes a list of rewrite rules (e.g. [cleanOnes,cleanAdjoint]) and applies them all to a QOp. Then write a function that applies such a list of rewrite rules repeatedly until a fixpoint is reached (i.e., applying the rules does not change the QOp anymore). You can use this to combine multiple optimization passes in your future work.
 -}

-- | Apply a list of rewrite rules once.
simplifyPass :: [o -> o] -> o -> o
simplifyPass rewriterules op = foldr (\f acc -> f acc) op rewriterules

-- | Fixpoint with a guard to ensure termination after at most n iterations.
fixpoint :: Eq a => Int -> (a -> a) -> a -> a 
fixpoint 0 _ x = x
fixpoint n f x = let 
      fx = f x
   in 
      if fx == x then x else fixpoint (n-1) f fx

-- | Apply a list of rewrite rules repeatedly until a fixpoint is reached, or at most n iterations.
simplifyFixpoint :: Eq o => Int -> [o -> o] -> o -> o
simplifyFixpoint n rewriterules op = fixpoint n (simplifyPass rewriterules) op

