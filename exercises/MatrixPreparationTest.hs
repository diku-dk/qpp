module Main where

import QPP
import QPP.Semantics.Statevector ( apply, evalOp, ket, StateT )
import QPP.Semantics.Matrix (CMat)

import Programs.MatrixPreparation
import Programs.MatrixArithmetic

import Numeric.LinearAlgebra (fromLists, norm_Frob,dispcf)

import Control.Monad (forM_)
 
import qualified Data.Vector as V
import Data.Complex


printbuildRowQOpTester :: IO ()
printbuildRowQOpTester = do
    let results = buildRowQOpTester ()
    forM_ results $ \res -> do
        printS res

main :: IO ()
main = do
    -- Number of decimals set in call to dispcf
    let decimals = 6
    putStrLn $ "matrixPrepTester2: " ++ dispcf decimals (matrixPrepTester2 ())
    putStrLn $ "matrixPrepTester3: " ++ dispcf decimals (matrixPrepTester3 ()) 
    putStrLn $ "matrixPrepTester4: " ++ dispcf decimals (matrixPrepTester4 ())
    putStrLn $ "matrixPrepTester5: " ++ dispcf decimals (matrixPrepTester5 ())
    putStrLn $ "productTester: " ++ dispcf decimals (productTester ())
    putStrLn $ "sumPrepTester: " ++ dispcf decimals (sumTester ())
    putStrLn $ "dyadicLCUTester: " ++ dispcf decimals (dyadicLCUTester ())
    putStrLn $ "discretizeTest: " ++ show (discretizeTest ())
    --print $ show (dyadicLCUTester ())
    --printbuildRowQOpTester

discretizeTest :: () -> V.Vector Int
discretizeTest () = discretize 3 (V.fromList [9/16,-4/16,2/16,-1/16,0])






sumTester :: () -> CMat
sumTester () =
    let
        -- Matrix data
        matData1 =   [[1,2,3], [4,5,6], [7,8,9]]
        matData2 =   [[10,11,12], [13,14,15],[16,17,18]]
        matIn1 = fromLists matData1
        matIn2 = fromLists matData2
        
        -- Norms
        norm1 = (norm_Frob matIn1 :+ 0)
        norm2 = (norm_Frob matIn2 :+ 0)

        -- Norm QOp
        normQOp = buildRowQOp (V.fromList [sqrt(norm1), sqrt(norm2)])  

        matData1norm = map (map (/ norm1)) [[1,2,3], [4,5,6], [7,8,9]]
        matData2norm = map (map (/ norm2)) [[10,11,12], [13,14,15],[16,17,18]]
        
        -- Normerede data matricer
        matIn1norm = fromLists matData1norm
        matIn2norm = fromLists matData2norm

        -- Sum encoding
        mat1QOp = matrixPrep matIn1norm
        mat2QOp = matrixPrep matIn2norm
        sumEncQOp = sumEncoding normQOp mat1QOp mat2QOp
        matQOt = evalOp $ sumEncQOp

        normFactor = (norm_Frob matIn1) + (norm_Frob matIn2)

        finalState1 = (normFactor :+ 0) .* (apply matQOt (ket [0,0,0,0,0]))
        finalState2 = (normFactor :+ 0) .* (apply matQOt (ket [0,0,0,0,1]))
        finalState3 = (normFactor :+ 0) .* (apply matQOt (ket [0,0,0,1,0]))
     
        m11 = inner  (ket [0,0,0,0,0]) finalState1
        m12 = inner  (ket [0,0,0,0,0]) finalState2
        m13 = inner  (ket [0,0,0,0,0]) finalState3

        m21 = inner  (ket [0,0,0,0,1]) finalState1
        m22 = inner  (ket [0,0,0,0,1]) finalState2
        m23 = inner  (ket [0,0,0,0,1]) finalState3
 
        m31 = inner  (ket [0,0,0,1,0]) finalState1
        m32 = inner  (ket [0,0,0,1,0]) finalState2
        m33 = inner  (ket [0,0,0,1,0]) finalState3
 
        matOutData = [[m11,m12,m13], [m21,m22,m23], [m31,m32,m33]]
        matOut = fromLists matOutData
    in 
        matOut - (matIn1 + matIn2)

productTester :: () -> CMat
productTester () =
    let
        -- Matrix data
        matData1 =   [[1,2], [3,4]]
        matData2 =   [[0,1], [1,0]]
        matIn1 = fromLists matData1
        matIn2 = fromLists matData2
        
        -- Product encoding
        mat1QOp = matrixPrep matIn1
        mat2QOp = matrixPrep matIn2
        prodEncQOp = productEncoding 1 mat1QOp mat2QOp
        matQOt = evalOp $ prodEncQOp

        normFactor = (norm_Frob matIn1 :+ 0) * (norm_Frob matIn2 :+ 0)

        finalState1 = normFactor .* (apply matQOt (ket [0,0,0]))
        finalState2 = normFactor .* (apply matQOt (ket [0,0,1]))

        m11 = inner (ket [0,0,0]) finalState1
        m12 = inner (ket [0,0,0]) finalState2
        
        m21 = inner (ket [0,0,1]) finalState1  
        m22 = inner (ket [0,0,1]) finalState2                                               

        matOutData =   [[m11,m12], [m21,m22]]
        matOut = fromLists matOutData 
        matInProduct = fromLists ([[2,1], [4,3]])
    in 
        matOut - matInProduct

buildRowQOpTester :: () -> [StateT] 
buildRowQOpTester =
    let 
        -- 2 dim
        resQOp2 = buildRowQOp (V.fromList [1:+1,2])      
        rowQOt2 = evalOp resQOp2
        -- Calculate the normalization factor and apply the operator
        finalState2 = sqrt(1**2 + 1**2 + 2**2) .* (apply rowQOt2 (ket ([0])))

        -- 3 dim
        resQOp3 = buildRowQOp (V.fromList [1,2,3])      
        rowQOt3 = evalOp resQOp3
        -- Calculate the normalization factor and apply the operator
        finalState3 = sqrt(1**2 + 2**2 + 3**2) .* (apply rowQOt3 (ket ([0,0])))

        -- 4 dim
        resQOp4 = buildRowQOp (V.fromList [1,2,3,4])      
        rowQOt4 = evalOp resQOp4
        -- Calculate the normalization factor and apply the operator
        finalState4 = sqrt(1**2 + 2**2 + 3**2 + 4**2) .* (apply rowQOt4 (ket ([0,0])))

        -- 5 dim
        resQOp5 = buildRowQOp (V.fromList [1,2,3,4,5])      
        rowQOt5 = evalOp resQOp5
        -- Calculate the normalization factor and apply the operator
        finalState5 = sqrt(1**2 + 2**2 + 3**2 + 4**2 + 5**2) .* (apply rowQOt5 (ket ([0,0,0])))

    in
        return [finalState2,finalState3,finalState4,finalState5]


matrixPrepTester2 :: () -> CMat 
matrixPrepTester2 () =
     let
        matIn = fromLists [[ 0.6 :+ 0.0, 0.0 :+ 0.0 ], 
                           [ 0.0 :+ 0.0, 0.3 :+ 0.0 ]]      

        -- Unitary encoding of the matrix
        matQOt = evalOp $ matrixPrep matIn

        -- Final states on input 00 and 01
        finalState1 =  (apply matQOt (ket [0,0]))
        finalState2 =  (apply matQOt (ket [0,1]))

        -- The encoding is normalized
        normMatIn = norm_Frob matIn

        -- Checking if the matrix has been upper left corner encoded
        -- Making sure the matrix is de-normalized
        m11 = (normMatIn :+ 0.0) * inner (ket [0,0]) finalState1
        m12 = (normMatIn :+ 0.0) * inner (ket [0,0]) finalState2
        m21 = (normMatIn :+ 0.0) * inner (ket [0,1]) finalState1
        m22 = (normMatIn :+ 0.0) * inner (ket [0,1]) finalState2
        
        -- Conversion to MatClassic and rounding (16 decimals)
        matOutData = [[m11,m12], [m21,m22]]
        matOut = fromLists matOutData
    in 
        matOut - matIn   

matrixPrepTester3 :: () -> CMat 
matrixPrepTester3 () =
     let
        matIn = fromLists  [[ 1.0 :+ 2.0, 3.0 :+ 4.0 , 5.0 :+ 0.0  ], 
                            [ 3.0 :+ 0.0, 4.0 :+ 0.0 , 5.0 :+ 6.0  ],
                            [ 0.0 :+ 0.0, 0.0 :+ 0.0 , 9.0 :+ 0.0  ]]

        matQOt = evalOp $ matrixPrep matIn

        finalState1 = (norm_Frob matIn :+ 0) .* (apply matQOt (ket [0,0,0,0]))
        finalState2 = (norm_Frob matIn :+ 0) .* (apply matQOt (ket [0,0,0,1]))
        finalState3 = (norm_Frob matIn :+ 0) .* (apply matQOt (ket [0,0,1,0]))
     
        m11 = inner  (ket [0,0,0,0]) finalState1
        m12 = inner  (ket [0,0,0,0]) finalState2
        m13 = inner  (ket [0,0,0,0]) finalState3

        m21 = inner  (ket [0,0,0,1]) finalState1
        m22 = inner  (ket [0,0,0,1]) finalState2
        m23 = inner  (ket [0,0,0,1]) finalState3
 
        m31 = inner  (ket [0,0,1,0]) finalState1
        m32 = inner  (ket [0,0,1,0]) finalState2
        m33 = inner  (ket [0,0,1,0]) finalState3
 
        matOutData = [[m11,m12,m13], [m21,m22,m23], [m31,m32,m33]]
        matOut = fromLists matOutData
    in 
        matOut - matIn 

matrixPrepTester4 :: () -> CMat
matrixPrepTester4 () =
     let
        matIn = fromLists  [[ 1.0 :+ 2.0, 3.0 :+ 4.0 , 5.0 :+ 0.0 ,  6.0 :+ 0.0 ], 
                            [ 3.0 :+ 0.0, 4.0 :+ 0.0 , 5.0 :+ 6.0 ,  7.0 :+ 0.0 ],
                            [ 0.0 :+ 0.0, 0.0 :+ 0.0 , 9.0 :+ 0.0 ,  0.0 :+ 0.0 ],
                            [ 0.0 :+ 0.0, 0.0 :+ 0.0 , 0.0 :+ 0.0 , 10.0 :+ 0.0 ]]

        matQOt = evalOp $ matrixPrep matIn

        finalState1 = (norm_Frob matIn :+ 0) .* (apply matQOt (ket [0,0,0,0]))
        finalState2 = (norm_Frob matIn :+ 0) .* (apply matQOt (ket [0,0,0,1]))
        finalState3 = (norm_Frob matIn :+ 0) .* (apply matQOt (ket [0,0,1,0]))
        finalState4 = (norm_Frob matIn :+ 0) .* (apply matQOt (ket [0,0,1,1]))
     
        m11 = inner  (ket [0,0,0,0]) finalState1
        m12 = inner  (ket [0,0,0,0]) finalState2
        m13 = inner  (ket [0,0,0,0]) finalState3
        m14 = inner  (ket [0,0,0,0]) finalState4

        m21 = inner  (ket [0,0,0,1]) finalState1
        m22 = inner  (ket [0,0,0,1]) finalState2
        m23 = inner  (ket [0,0,0,1]) finalState3
        m24 = inner  (ket [0,0,0,1]) finalState4

        m31 = inner  (ket [0,0,1,0]) finalState1
        m32 = inner  (ket [0,0,1,0]) finalState2
        m33 = inner  (ket [0,0,1,0]) finalState3
        m34 = inner  (ket [0,0,1,0]) finalState4
     
        m41 = inner  (ket [0,0,1,1]) finalState1
        m42 = inner  (ket [0,0,1,1]) finalState2
        m43 = inner  (ket [0,0,1,1]) finalState3
        m44 = inner  (ket [0,0,1,1]) finalState4

        matOutData = [[m11,m12,m13,m14], [m21,m22,m23,m24], [m31,m32,m33,m34],[m41,m42,m43,m44]]
        matOut = fromLists matOutData
    in 
        matOut - matIn  

matrixPrepTester5 :: () -> CMat
matrixPrepTester5 () = 
    let
        matIn = fromLists  [[ 1.0 :+ 2.0, 3.0 :+ 4.0 , 5.0 :+ 0.0 ,  6.0 :+ 0.0 ,    0.0 :+ 7.0 ], 
                            [ 3.0 :+ 0.0, 4.0 :+ 0.0 , 5.0 :+ 6.0 ,  7.0 :+ 0.0 , (-8.0) :+ 0.0 ],
                            [ 0.0 :+ 0.0, 0.0 :+ 0.0 , 9.0 :+ 0.0 ,  0.0 :+ 0.0 ,    0.0 :+ 0.0 ],
                            [ 0.0 :+ 0.0, 0.0 :+ 0.0 , 0.0 :+ 0.0 , 10.0 :+ 0.0 ,    0.0 :+ 0.0 ],
                            [ 1.0 :+ 2.0, 3.0 :+ 4.0 , 5.0 :+ 0.0 ,  0.0 :+ 6.0 ,   11.0 :+ 0.0 ]]

        matQOt = evalOp $ matrixPrep matIn

        finalState1 = (norm_Frob matIn :+ 0) .* apply matQOt (ket [0,0,0,0,0,0])
        finalState2 = (norm_Frob matIn :+ 0) .* apply matQOt (ket [0,0,0,0,0,1])
        finalState3 = (norm_Frob matIn :+ 0) .* apply matQOt (ket [0,0,0,0,1,0])
        finalState4 = (norm_Frob matIn :+ 0) .* apply matQOt (ket [0,0,0,0,1,1])
        finalState5 = (norm_Frob matIn :+ 0) .* apply matQOt (ket [0,0,0,1,0,0])

        m11 = inner  (ket [0,0,0,0,0,0]) finalState1
        m12 = inner  (ket [0,0,0,0,0,0]) finalState2
        m13 = inner  (ket [0,0,0,0,0,0]) finalState3
        m14 = inner  (ket [0,0,0,0,0,0]) finalState4
        m15 = inner  (ket [0,0,0,0,0,0]) finalState5

        m21 = inner  (ket [0,0,0,0,0,1]) finalState1
        m22 = inner  (ket [0,0,0,0,0,1]) finalState2
        m23 = inner  (ket [0,0,0,0,0,1]) finalState3
        m24 = inner  (ket [0,0,0,0,0,1]) finalState4
        m25 = inner  (ket [0,0,0,0,0,1]) finalState5

        m31 = inner  (ket [0,0,0,0,1,0]) finalState1
        m32 = inner  (ket [0,0,0,0,1,0]) finalState2
        m33 = inner  (ket [0,0,0,0,1,0]) finalState3
        m34 = inner  (ket [0,0,0,0,1,0]) finalState4
        m35 = inner  (ket [0,0,0,0,1,0]) finalState5
        
        m41 = inner  (ket [0,0,0,0,1,1]) finalState1
        m42 = inner  (ket [0,0,0,0,1,1]) finalState2
        m43 = inner  (ket [0,0,0,0,1,1]) finalState3
        m44 = inner  (ket [0,0,0,0,1,1]) finalState4
        m45 = inner  (ket [0,0,0,0,1,1]) finalState5

        m51 = inner  (ket [0,0,0,1,0,0]) finalState1
        m52 = inner  (ket [0,0,0,1,0,0]) finalState2
        m53 = inner  (ket [0,0,0,1,0,0]) finalState3
        m54 = inner  (ket [0,0,0,1,0,0]) finalState4
        m55 = inner  (ket [0,0,0,1,0,0]) finalState5

        matOutData = [[m11,m12,m13,m14,m15], [m21,m22,m23,m24,m25], [m31,m32,m33,m34,m35],[m41,m42,m43,m44,m45],[m51,m52,m53,m54,m55]]
        matOut = fromLists matOutData
    in 
        matOut - matIn 

 
dyadicLCUTester :: () -> CMat 
dyadicLCUTester () =
     let
        matIn1 = fromLists  [[ 1.0 :+ 0.0, 0.0 :+ 0.0 ], 
                             [ 0.0 :+ 0.0, 1.0 :+ 0.0 ]]
        
        matIn2 = fromLists  [[ 0.0 :+ 0.0, 1.0 :+ 0.0 ], 
                             [ 1.0 :+ 0.0, 0.0 :+ 0.0 ]]
        
        -- Unitary encoding of the matrices
        -- and call to nonRotLCU to 
        precision = 3
        lcuQOp = dyadicLCU precision (V.fromList [(0.75, matrixPrep matIn1), (0.25, matrixPrep matIn2)])
        lcuQOt = evalOp $ lcuQOp

        -- Final states on input 00 and 01
        finalState1 =  apply lcuQOt (ket [0,0,0,0,0,0])
        finalState2 =  apply lcuQOt (ket [0,0,0,0,0,1])

        -- The encoding is normalized
        -- normMatIn = norm_Frob matIn

        -- Checking if the matrix has been upper left corner encoded
        -- Making sure the matrix is de-normalized
        m11 = inner (ket [0,0,0,0,0,0]) finalState1 --(normMatIn :+ 0.0) *
        m12 = inner (ket [0,0,0,0,0,0]) finalState2 --(normMatIn :+ 0.0) *
        m21 = inner (ket [0,0,0,0,0,1]) finalState1 --(normMatIn :+ 0.0) *
        m22 = inner (ket [0,0,0,0,0,1]) finalState2 --(normMatIn :+ 0.0) *
        
        -- Conversion to MatClassic and rounding (16 decimals)
        matOutData = [[m11,m12], [m21,m22]]
        matOut = fromLists matOutData
    in 
        matOut --     - (2 * matIn1 - matIn2)
