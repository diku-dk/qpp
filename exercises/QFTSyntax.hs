-- | Inspect the QFT as QOp, UOp and flat-circuit syntax, without evaluating its matrix.
-- Run with `cabal run QFTSyntax -- 16` (16 is also the default).
module Main where

import Data.List (intercalate)
import Programs.QFT (qft)
import qualified QPP as Q
import qualified QPP.QCirc as QC
import qualified QPP.UOp as U
import System.Directory (createDirectoryIfMissing)
import System.Environment (getArgs)
import System.FilePath ((</>))
import Text.Read (readMaybe)

main :: IO ()
main = do
  args <- getArgs
  n <- case args of
    [] -> pure 16
    [arg] | Just k <- readMaybe arg, k >= 0 -> pure k
    _ -> fail "Usage: QFTSyntax [nonnegative qubit count]"
  let q = qft n
  u <- either (fail . show) pure (U.fromQOp q)
  let normalized = U.normalize u
      circuit = QC.fromQOp q
      directory = "output" </> ("qft" ++ show n)
      outputs =
        [ ("qop.pretty.txt", Q.showOp q)
        , ("uop.raw.txt", show u)
        , ("uop.pretty.txt", U.showUOp u)
        , ("uop.normalized.raw.txt", show normalized)
        , ("uop.normalized.pretty.txt", showLayers normalized)
        , ("qcirc.txt", QC.showQCirc circuit)
        ]
      summary = unlines
        [ "qft " ++ show n ++ ": " ++ show n ++ " qubits; dimension " ++ show (U.dim u)
        , "QOp: " ++ show (nodes qChildren q) ++ " syntax nodes; depth " ++ show (depth qChildren q)
        , "UOp: " ++ show (nodes uChildren u) ++ " syntax nodes; depth " ++ show (depth uChildren u)
        , "Normalized UOp: " ++ show (nodes uChildren normalized) ++ " syntax nodes; depth " ++ show (depth uChildren normalized)
        , "Normalized top-level composition: " ++ show (compositionLength normalized) ++ " factors"
        , "QCirc: " ++ show (length (QC.circ_gates circuit)) ++ " gates on " ++ show (QC.circ_qubits circuit) ++ " wires"
        , "Counts expand every tree occurrence, including shared subexpressions."
        , "No operator matrix was constructed."
        ]
  createDirectoryIfMissing True directory
  writeFile (directory </> "summary.txt") summary
  putStr summary
  mapM_ (writeOutput directory) outputs

writeOutput :: FilePath -> (FilePath, String) -> IO ()
writeOutput directory (name, contents) = do
  let path = directory </> name
  writeFile path (contents ++ "\n")
  putStrLn (path ++ " (" ++ show (length contents) ++ " characters)")

nodes :: (a -> [a]) -> a -> Int
nodes children a = 1 + sum (map (nodes children) (children a))

depth :: (a -> [a]) -> a -> Int
depth children a = 1 + maximum (0 : map (depth children) (children a))

compositionLength :: U.UOp -> Int
compositionLength (U.Compose _ xs) = length xs
compositionLength _ = 1

-- Display the outer composition one factor per line, in mathematical order
-- (the last factor acts first). Each factor retains the ordinary UOp notation.
showLayers :: U.UOp -> String
showLayers (U.Compose _ xs@(_:_:_)) =
  "(\n  " ++ intercalate "\n  ∘ " (map U.showUOp xs) ++ "\n)"
showLayers op = U.showUOp op

qChildren :: Q.QOp -> [Q.QOp]
qChildren q = case q of
  Q.R a _ -> [a]
  Q.C a -> [a]
  Q.Tensor a b -> [a,b]
  Q.DirectSum a b -> [a,b]
  Q.Compose a b -> [a,b]
  Q.Adjoint a -> [a]
  _ -> []

uChildren :: U.UOp -> [U.UOp]
uChildren u = case u of
  U.R a _ -> [a]
  U.Tensor xs -> xs
  U.DirectSum xs -> xs
  U.Compose _ xs -> xs
  _ -> []
