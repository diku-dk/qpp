module QPP.PrettyPrint.Op where
import QPP.Syntax
import Data.List (intercalate)

showOp :: QOp -> String
showOp op = case op of
    C a             -> "C ("++ showOp a ++ ")"
    a `Tensor`    b -> "(" ++ showOp a ++ " ⊗ " ++ showOp b ++ ")"
    a `DirectSum` b -> "(" ++ showOp a ++ " ⊕ "  ++ showOp b ++ ")"
    a `Compose`   b -> "(" ++ showOp a ++ " ∘ " ++ showOp b ++ ")"
    Adjoint a       -> "(adj " ++ showOp a ++ ")"
    Id 0     -> "One"
    Id 1     -> "I"
    _        -> show op
    
showStep :: Step -> String
showStep (Unitary op) = "Unitary $ " ++ showOp op
showStep step = show step

showProgram :: Program -> String
showProgram steps = intercalate "\n" [ "step" ++ show i ++ " = " ++ (showStep step)
                     | (i :: Int,step) <- zip [1..] steps
                    ]

printOp :: QOp -> IO ()
printOp = putStrLn . showOp