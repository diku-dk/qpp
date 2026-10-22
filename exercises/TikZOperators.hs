-- | Generate reusable diagram fragments and a small LaTeX gallery.
-- Run from the repository root with `cabal run TikZOperators`.
module Main where

import Programs.QFT (qft)
import qualified QPP as Q
import qualified QPP.QCirc as QC
import qualified QPP.UOp as U
import System.Directory (copyFile, createDirectoryIfMissing)
import System.FilePath ((</>))

main :: IO ()
main = do
  translatedExample <- lower qopExample
  qft2 <- lower (qft 2)
  qft4 <- lower (qft 4)
  let directory = "output" </> "tikz"
      fragments =
        [ ("qop-example", Q.showTikZQOp qopExample)
        , ("uop-example", U.showTikZUOp translatedExample)
        , ("qop-annotated", Q.showTikZQOp qopExample)
        , ("qop-sum-annotated", Q.showTikZQOp (Q.DirectSum (Q.Id 3) (Q.Id 3)))
        , ("uop-annotated", U.showTikZUOp arbitraryDimensionExample)
        , ("uop-empty-annotated", U.showTikZUOp (U.DirectSum [U.Tensor [], U.DirectSum []]))
        , ("qop-gates", Q.showTikZQOp (foldr1 Q.Compose [Q.X,Q.Y,Q.Z,Q.H,Q.SX]))
        , ("qop-scalars", Q.showTikZQOp (Q.Compose (Q.Id 0) (Q.Phase (-1/3))))
        , ("qop-adjoint", Q.showTikZQOp (Q.Adjoint (Q.DirectSum Q.X Q.Y)))
        , ("qop-permute", Q.showTikZQOp (Q.Permute [2,0,1]))
        , ("uop-permute", U.showTikZUOp (U.PermuteTensor [2,3,4] [2,0,1]))
        , ("uop-empty-tensor", U.showTikZUOp (U.Tensor []))
        , ("uop-empty-sum", U.showTikZUOp (U.DirectSum []))
        , ("uop-identity", U.showTikZUOp (U.identity 0))
        , ("qft2-qop", Q.showTikZQOp (qft 2))
        , ("qft2-qop-annotated", Q.showTikZQOp (qft 2))
        , ("qft2-uop", U.showTikZUOp (U.normalize qft2))
        , ("qft2-uop-without-qubits", U.showTikZUOpWithoutQubits (U.normalize qft2))
        , ("qft4-qop", Q.showTikZQOp (qft 4))
        , ("qft4-uop", U.showTikZUOp qft4)
        , ("qft4-uop-normalized", U.showTikZUOp (U.normalize qft4))
        , ("circuit-bell", QC.showQuantikz (QC.QCirc 2 [QC.Gate Q.H [0], QC.Gate (Q.C Q.X) [0,1]]))
        , ("circuit-elements", QC.showQuantikz circuitElements)
        , ("qft3-circuit", QC.showQuantikz (QC.fromQOp (qft 3)))
        , ("qft4-circuit-simplified", QC.showQuantikz (QC.simplify (QC.fromQOp (qft 4))))
        ]
  createDirectoryIfMissing True directory
  copyFile "latex/qpp-tikz.sty" (directory </> "qpp-tikz.sty")
  mapM_ (\(name, body) -> writeFile (directory </> name ++ ".tex") body) fragments
  copyFile "latex/operator-diagrams.tex" (directory </> "operator-diagrams.tex")
  putStrLn ("Wrote TikZ fragments and operator-diagrams.tex to " ++ directory)
  where
    lower = either (fail . show) pure . U.fromQOp

qopExample :: Q.QOp
qopExample = Q.Compose (Q.Tensor Q.H (Q.Id 1)) (Q.C (Q.R Q.Z (1/4)))

arbitraryDimensionExample :: U.UOp
arbitraryDimensionExample = U.Compose 3
  [ U.PermuteSum [1,2] [1,0]
  , U.DirectSum [U.phase (1/3), U.R (U.PermuteSum [1,1] [1,0]) (1/2)]
  ]

-- | One of each element the circuit printers draw: a global phase, a swap, a
--   controlled swap, a Toffoli with its target between its controls, a
--   controlled phase with its control below, a phase gate, a box on
--   non-consecutive wires in reversed order, and an adjoint.
circuitElements :: QC.QCirc
circuitElements = QC.QCirc 3
  [ QC.Gate (Q.Phase (3/4)) []
  , QC.Gate (Q.Permute [1,0]) [0,1]
  , QC.Gate (Q.C (Q.Permute [1,0])) [0,1,2]
  , QC.Gate (Q.C (Q.C Q.X)) [0,2,1]
  , QC.Gate (Q.C (Q.C (Q.Phase (1/4)))) [2,0]
  , QC.Gate (Q.C (Q.Phase (1/2))) [1]
  , QC.Gate (Q.R (Q.Tensor Q.X Q.Y) (1/3)) [2,0]
  , QC.Gate (Q.Adjoint Q.SX) [1]
  ]
