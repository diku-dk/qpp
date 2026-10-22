-- | Gate listings and wire diagrams of flat circuits.
--
--   'showQCirc' lists the gates one per line. 'showQCircDiagram' (wires as
--   rows), 'showQCircDiagramVertical' (wires as columns, time running down)
--   and 'showQuantikz' draw the circuit on its wires, as text and as a quantikz
--   LaTeX fragment. The drawings share one model: gates on no wires are
--   set aside as scalars, the other gates are placed in columns as soon as
--   possible without two gates covering a common row ('layersBy' 'gateSpan'),
--   and each gate is drawn as controls around one element, the ⊕ of a
--   controlled NOT, the crosses of a swap, the dot of a phase gate, or a box.
module QPP.PrettyPrint.QCirc
  ( showQCirc, printQCirc
  , showQCircDiagram, printQCircDiagram
  , showQCircDiagramVertical, printQCircDiagramVertical
  , showQuantikz
  ) where

import Data.Char (isSpace)
import Data.List (dropWhileEnd, elemIndex, intercalate, partition, sort)
import qualified Data.Map.Strict as M
import Data.Maybe (mapMaybe)
import Data.Ratio (denominator, numerator)
import qualified Data.Set as S
import QPP.PrettyPrint.Op (showOp)
import QPP.QCirc.Syntax
import QPP.Syntax

--------------------------------------------------------------------------------
-- Listing
--------------------------------------------------------------------------------

-- | @qubits: n@ followed by one @op \@ wires@ line per gate, first gate first.
showQCirc :: QCirc -> String
showQCirc (QCirc n gs) = intercalate "\n" (("qubits: " ++ show n) : map showGate gs)
  where
    showGate (Gate g ws) = showOp g ++ " @ " ++ show ws

printQCirc :: QCirc -> IO ()
printQCirc = putStrLn . showQCirc

--------------------------------------------------------------------------------
-- Drawing model shared by the diagrams
--------------------------------------------------------------------------------

-- | What a gate puts on its wires. A gate has one element besides its
--   controls, and its controls lie outside the rows of that element.
data Element
  = Control Wire
  | Target Wire             -- ^ The ⊕ of a controlled NOT.
  | Swap Wire Wire          -- ^ The crossed wires of a swap, upper row first.
  | PhaseDot Wire Rational  -- ^ The phase gate diag(1, e^{iπθ}): the last wire of @C^k (Phase θ)@.
  | Box Wire [Wire] QOp     -- ^ A labelled box: the wire of local qubit 0, the other wires in local order, the operator.

-- | A gate as controls around one element. The controls of @C^k a@ are peeled
--   off when they all lie above or below the rows of @a@; otherwise the whole
--   gate is one box, so that nothing is ever drawn inside a box. A controlled
--   @X@ gets a ⊕, a swap of two wires gets crosses, @C^k (Phase θ)@ gets a
--   phase dot on its last wire, and anything else a box.
elements :: Gate -> [Element]
elements (Gate g ws) = case (core, targets) of
  (Phase t, [])            | c : cs <- reverse controls -> dots (reverse cs) ++ [PhaseDot c t]
  (_, [])                  | c : cs <- reverse controls -> dots (reverse cs) ++ [Box c [] (C core)]
  (X, [t])                 | k > 0, clear [t]           -> dots controls ++ [Target t]
  (Permute [1, 0], [a, b]) | clear [a, b]               -> dots controls ++ [Swap (min a b) (max a b)]
  (_, t : ts)              | clear targets              -> dots controls ++ [Box t ts core]
  _ -> case ws of
    w : more -> [Box w more g]
    []       -> []
  where
    (k, core) = peelControls g
    (controls, targets) = splitAt k ws
    dots = map Control
    clear ts = all (\c -> c < minimum ts || c > maximum ts) controls

-- | The scalars (gates on no wires) and the columns of the other gates: span
--   layers, each gate with its elements.
layout :: QCirc -> ([QOp], [[(Gate, [Element])]])
layout (QCirc _ gs) =
  ([ g | Gate g [] <- scalars ], [ [ (g, elements g) | g <- column ] | column <- layersBy gateSpan wired ])
  where
    (scalars, wired) = partition (null . gateWires) gs

-- | The rows a gate covers, from its lowest wire to its highest.
rowsOf :: Gate -> (Wire, Wire)
rowsOf (Gate _ ws) = (minimum ws, maximum ws)

--------------------------------------------------------------------------------
-- Labels
--------------------------------------------------------------------------------

-- | How a diagram spells the operators in its boxes and scalar lines.
data Notation = Notation
  { notePi       :: String
  , noteTensor   :: String
  , noteSum      :: String
  , noteCompose  :: String
  , noteSqrtX    :: String
  , noteAdjoint  :: String -> String
  , noteIdentity :: Nat -> String
  , notePhase    :: String -> String            -- ^ From the angle.
  , noteRotation :: String -> String -> String  -- ^ From the axis and the angle.
  , notePermute  :: [Nat] -> String
  }

textNotation :: Notation
textNotation = Notation
  { notePi = "π", noteTensor = " ⊗ ", noteSum = " ⊕ ", noteCompose = " ∘ ", noteSqrtX = "√X"
  , noteAdjoint = (++ "†")
  , noteIdentity = \n -> case n of { 0 -> "1"; 1 -> "I"; _ -> "Id " ++ show n }
  , notePhase = \t -> "Phase(" ++ t ++ ")"
  , noteRotation = \a t -> "R(" ++ a ++ ", " ++ t ++ ")"
  , notePermute = \ks -> "Perm" ++ show ks
  }

latexNotation :: Notation
latexNotation = Notation
  { notePi = "\\pi", noteTensor = " \\otimes ", noteSum = " \\oplus ", noteCompose = " \\circ ", noteSqrtX = "\\sqrt{X}"
  , noteAdjoint = (++ "^{\\dagger}")
  , noteIdentity = \n -> case n of { 0 -> "1"; 1 -> "I"; _ -> "I_{" ++ show n ++ "}" }
  , notePhase = \t -> "e^{i\\," ++ t ++ "}"
  , noteRotation = \a t -> "R_{" ++ a ++ "}(" ++ t ++ ")"
  , notePermute = \ks -> "\\mathrm{Perm}" ++ show ks
  }

-- | An angle in units of π: @0@, @π@, @-π/2@, @3π/4@.
angle :: Notation -> Rational -> String
angle nt t
  | t == 0    = "0"
  | otherwise = coefficient ++ notePi nt ++ (if d == 1 then "" else "/" ++ show d)
  where
    n = numerator t
    d = denominator t
    coefficient
      | n == 1    = ""
      | n == -1   = "-"
      | otherwise = show n

-- | An operator as a box label, with binary operands parenthesized.
label :: Notation -> QOp -> String
label nt = go
  where
    go op = case op of
      Id n          -> noteIdentity nt n
      X             -> "X"
      Y             -> "Y"
      Z             -> "Z"
      H             -> "H"
      SX            -> noteSqrtX nt
      Phase t       -> notePhase nt (angle nt t)
      R a t         -> noteRotation nt (go a) (angle nt t)
      C a           -> "C(" ++ go a ++ ")"
      Permute ks    -> notePermute nt ks
      Adjoint a     -> noteAdjoint nt (operand a)
      Tensor a b    -> operand a ++ noteTensor nt ++ operand b
      DirectSum a b -> operand a ++ noteSum nt ++ operand b
      Compose a b   -> operand a ++ noteCompose nt ++ operand b
    operand a
      | binary a  = "(" ++ go a ++ ")"
      | otherwise = go a
    binary op = case op of
      Tensor _ _    -> True
      DirectSum _ _ -> True
      Compose _ _   -> True
      _             -> False

--------------------------------------------------------------------------------
-- Text diagram
--------------------------------------------------------------------------------

-- | The circuit drawn on its wires in line-drawing characters: a row per wire
--   and a connector row between neighbouring wires. Controls are @●@, the
--   target of a controlled NOT @⊕@, swapped wires @×@, a phase gate @P(θ)@,
--   and a wire crossed by a gate's vertical line @┼@. A box on several wires
--   has its label on the wire of its local qubit 0 and @#i@ on the wire of
--   local qubit @i@. Gates on no wires come first, as @scalar:@ lines.
showQCircDiagram :: QCirc -> String
showQCircDiagram circuit@(QCirc n _) =
  intercalate "\n" (map scalarLine scalars ++ map trimEnd (interleave wireRows gapRows))
  where
    (scalars, cols) = layout circuit
    columns = map textColumn cols
    scalarLine g = "scalar: " ++ label textNotation g
    name w = "q" ++ show w ++ ":"
    nameWidth = length (name (n - 1))
    wireRows =
      [ padRight nameWidth (name w) ++ " ─"
          ++ concat [ "─" ++ centre '─' width (M.findWithDefault "" w cells) ++ "─" | (width, cells, _) <- columns ]
          ++ "─"
      | w <- [0 .. n - 1] ]
    gapRows =
      [ replicate nameWidth ' ' ++ "  "
          ++ concat [ " " ++ centre ' ' width (if S.member w bars then "│" else "") ++ " " | (width, _, bars) <- columns ]
      | w <- [0 .. n - 2] ]

printQCircDiagram :: QCirc -> IO ()
printQCircDiagram = putStrLn . showQCircDiagram

-- | The same drawing with the wires as columns and time running down the
--   page, for long circuits: a header of wire names, then a row per column of
--   'showQCircDiagram' with a row of bare wires between. Each wire's column is
--   as wide as its widest cell, and a gate's connector runs horizontally.
showQCircDiagramVertical :: QCirc -> String
showQCircDiagramVertical circuit@(QCirc n _) =
  intercalate "\n" (map scalarLine scalars ++ map trimEnd (header : wires : concatMap (\c -> [layerRow c, wires]) columns))
  where
    (scalars, cols) = layout circuit
    columns = map textColumn cols
    scalarLine g = "scalar: " ++ label textNotation g
    name w = "q" ++ show w
    widths = M.fromList [ (w, maximum (length (name w) : [ length s | (_, cells, _) <- columns, Just s <- [M.lookup w cells] ]))
                        | w <- [0 .. n - 1] ]
    width w = M.findWithDefault 1 w widths
    row cell gap = concat [ cell w ++ (if w < n - 1 then gap w else "") | w <- [0 .. n - 1] ]
    header = row (\w -> centre ' ' (width w) (name w)) (const "  ")
    wires = row (\w -> centre ' ' (width w) "│") (const "  ")
    layerRow (_, cells, bars) = row cell gap
      where
        cell w = centreWith (fill (w - 1)) (fill w) (width w) (M.findWithDefault "│" w cells)
        gap w = if S.member w bars then "──" else "  "
        fill w = if S.member w bars then '─' else ' '

printQCircDiagramVertical :: QCirc -> IO ()
printQCircDiagramVertical = putStrLn . showQCircDiagramVertical

-- | A column's width, its cell texts by wire, and the wires with a vertical
--   line between them and the next wire.
textColumn :: [(Gate, [Element])] -> (Int, M.Map Wire String, S.Set Wire)
textColumn column = (maximum (1 : map length (M.elems cells)), cells, bars)
  where
    cells = M.unions [ M.fromList (concatMap textCells es) `M.union` M.fromList [ (w, "┼") | w <- gateSpan g ]
                     | (g, es) <- column ]
    bars = S.fromList [ w | (g, _) <- column, let (top, bottom) = rowsOf g, w <- [top .. bottom - 1] ]

textCells :: Element -> [(Wire, String)]
textCells e = case e of
  Control w    -> [(w, "●")]
  Target w     -> [(w, "⊕")]
  Swap a b     -> [(a, "×"), (b, "×")]
  PhaseDot w t -> [(w, "P(" ++ angle textNotation t ++ ")")]
  Box w ws g   -> (w, label textNotation g) : [ (v, '#' : show i) | (i, v) <- zip [1 :: Int ..] ws ]

centre :: Char -> Int -> String -> String
centre fill = centreWith fill fill

-- | Pad a text to a width, with separate fillers on each side.
centreWith :: Char -> Char -> Int -> String -> String
centreWith leftFill rightFill width s = replicate left leftFill ++ s ++ replicate (missing - left) rightFill
  where
    missing = width - length s
    left = missing `div` 2

padRight :: Int -> String -> String
padRight width s = s ++ replicate (width - length s) ' '

trimEnd :: String -> String
trimEnd = dropWhileEnd isSpace

interleave :: [a] -> [a] -> [a]
interleave (x : xs) (y : ys) = x : y : interleave xs ys
interleave xs []             = xs
interleave [] ys             = ys

--------------------------------------------------------------------------------
-- quantikz
--------------------------------------------------------------------------------

-- | The circuit as a @quantikz@ environment for a document loading the
--   @quantikz@ package: one row per wire, led by @\\lstick{$q_i$}@ and padded
--   with @\\qw@. Controls are @\\ctrl@, each linked to the next control
--   towards its target; the target of a controlled NOT is @\\targ@, a swap
--   @\\swap@ with @\\targX@, a phase gate @\\phase@, and a box @\\gate@. A box
--   on several wires spans them with @\\gate[wires=k]@, puts @\\linethrough@
--   on wires it does not act on, and numbers its inputs with @\\gateinput@
--   unless its wires are consecutive from the top. Gates on no wires become
--   @%@ comment lines before the environment.
showQuantikz :: QCirc -> String
showQuantikz circuit@(QCirc n _) =
  intercalate "\n" (map scalarComment scalars ++ ["\\begin{quantikz}"] ++ rows ++ ["\\end{quantikz}"])
  where
    (scalars, cols) = layout circuit
    columns = map (M.unions . map latexCells) cols
    scalarComment g = "% scalar: " ++ label latexNotation g
    rows =
      [ intercalate " & " (("\\lstick{$q_{" ++ show w ++ "}$}") : [ M.findWithDefault "\\qw" w cells | cells <- columns ] ++ ["\\qw"])
          ++ (if w < n - 1 then " \\\\" else "")
      | w <- [0 .. n - 1] ]

-- | A gate's cells by wire.
latexCells :: (Gate, [Element]) -> M.Map Wire String
latexCells (_, es) = M.fromList (concatMap cell es)
  where
    controls = sort [ c | Control c <- es ]
    (top, bottom) = case mapMaybe target es of
      r : _ -> r
      []    -> (0, 0)
    target e = case e of
      Control _    -> Nothing
      Target t     -> Just (t, t)
      Swap a b     -> Just (a, b)
      PhaseDot w _ -> Just (w, w)
      Box w ws _   -> Just (minimum (w : ws), maximum (w : ws))
    -- A control links to the next control towards the target, or to the target.
    link c
      | c < top   = case [ d | d <- controls, c < d, d < top ] of
          d : _ -> d - c
          []    -> top - c
      | otherwise = case reverse [ d | d <- controls, bottom < d, d < c ] of
          d : _ -> d - c
          []    -> bottom - c
    cell e = case e of
      Control c    -> [(c, "\\ctrl{" ++ show (link c) ++ "}")]
      Target t     -> [(t, "\\targ{}")]
      Swap a b     -> [(a, "\\swap{" ++ show (b - a) ++ "}"), (b, "\\targX{}")]
      PhaseDot w t -> [(w, "\\phase{" ++ angle latexNotation t ++ "}")]
      Box w ws g   -> box (w : ws) g
    box wires g = [ (r, boxCell r) | r <- [lo .. hi] ]
      where
        lo = minimum wires
        hi = maximum wires
        numbered = wires /= [lo .. hi]
        input r = case elemIndex r wires of
          Just i | numbered -> "\\gateinput{" ++ show i ++ "}"
          _                 -> ""
        boxCell r
          | r == lo        = gate ++ input r
          | r `elem` wires = input r
          | otherwise      = "\\linethrough"
        gate
          | lo == hi  = "\\gate{" ++ label latexNotation g ++ "}"
          | otherwise = "\\gate[wires=" ++ show (hi - lo + 1) ++ "]{" ++ label latexNotation g ++ "}"
