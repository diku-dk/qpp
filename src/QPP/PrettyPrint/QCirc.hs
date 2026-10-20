-- | Gate listings of flat circuits: one gate per line, in time order.
module QPP.PrettyPrint.QCirc (showQCirc, printQCirc) where

import Data.List (intercalate)
import QPP.PrettyPrint.Op (showOp)
import QPP.QCirc.Syntax

-- | @qubits: n@ followed by one @op \@ wires@ line per gate, first gate first.
showQCirc :: QCirc -> String
showQCirc (QCirc n gs) = intercalate "\n" (("qubits: " ++ show n) : map showGate gs)
  where
    showGate (Gate g ws) = showOp g ++ " @ " ++ show ws

printQCirc :: QCirc -> IO ()
printQCirc = putStrLn . showQCirc
