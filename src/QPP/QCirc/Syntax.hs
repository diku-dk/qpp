-- | Flat quantum circuits on a fixed register: a wire count and a time-ordered
--   list of gates, each a 'QOp' placed on explicitly listed wires.
--
--   This is the gate-level representation, for listing, counting and emitting.
--   It is bridged to the compositional 'QOp' syntax by 'toQOp' here and by
--   'QPP.QCirc.FromQOp.fromQOp'. Functions assume the invariants documented on
--   'Gate' and 'QCirc', just as QOp functions assume 'typeCheck'; 'checkQCirc'
--   validates them.
module QPP.QCirc.Syntax
  ( Wire, Gate(..), QCirc(..), QCircError(..)
  , checkQCirc
  , toQOp
  , peelControls
  , gateWires, gateSpan, layersBy, layers, depth
  ) where

import Control.Monad (unless)
import Data.List ((\\), nub)
import qualified Data.Map.Strict as M
import QPP.Syntax
import QPP.TypeCheck (TypeError, typeCheck)
import QPP.Util (permInvert)

--------------------------------------------------------------------------------
-- Core syntax
--------------------------------------------------------------------------------

-- | Absolute wire index in a circuit's register.
-- ANCHOR: qcirc-syntax
type Wire = Nat

-- | @Gate op ws@ applies @op@ with its local qubit @i@ on wire @ws !! i@.
--   Invariants: @op_qubits op == length ws@ and the wires are distinct. The list
--   may be empty (a global 'Phase'), non-consecutive (an idle wire between
--   control and target), or decreasing (control below target).
data Gate = Gate QOp [Wire]
  deriving (Show, Eq)

-- | A circuit on @circ_qubits@ wires. Gates are in time order: the head acts
--   first. Every gate's wires lie in @[0, circ_qubits)@.
data QCirc = QCirc { circ_qubits :: Nat, circ_gates :: [Gate] }
  deriving (Show, Eq)
-- ANCHOR_END: qcirc-syntax

--------------------------------------------------------------------------------
-- Shared operator notation
--------------------------------------------------------------------------------

-- Composition is in mathematical order, as for QOp: @a <> b@ runs @b@ first.
instance Semigroup QCirc where
  a <> b = QCirc (circ_qubits a) (circ_gates b ++ circ_gates a)

-- @a ⊗ b@ puts @a@ on the first wires and @b@ on the following ones.
instance HasTensorProduct QCirc where
  a ⊗ b = QCirc (na + circ_qubits b) (circ_gates a ++ map shift (circ_gates b))
    where
      na = circ_qubits a
      shift (Gate g ws) = Gate g (map (+ na) ws)

-- @a ⊕ b@ conditions on a new wire 0: it equals @(a ⊕ I) ∘ (I ⊕ b)@, and
-- @I ⊕ b = C b@ distributes over b's gates, @C (g_k ∘ … ∘ g_1) = C g_k ∘ … ∘ C g_1@,
-- while @a ⊕ I@ is @C a@ conjugated by X on the control wire.
instance HasDirectSum QCirc where
  a ⊕ b = QCirc (circ_qubits a + 1) (anti (circ_gates a) ++ map control (circ_gates b))
    where
      control (Gate g ws) = Gate (C g) (0 : map (+ 1) ws)
      anti [] = []
      anti gs = Gate X [0] : map control gs ++ [Gate X [0]]

-- The structural 'dagger' keeps every gate free of a leading 'Adjoint'.
instance HasAdjoint QCirc where
  adj (QCirc n gs) = QCirc n (reverse [ Gate (dagger g) ws | Gate g ws <- gs ])

instance HasQubits QCirc where n_qubits = circ_qubits
instance Operator QCirc

--------------------------------------------------------------------------------
-- Validation
--------------------------------------------------------------------------------

data QCircError
  = NegativeWidth Nat
  | WireOutOfRange Gate
  | DuplicateWire Gate
  | GateArity Gate Nat Nat   -- ^ The gate, its operator's qubit count, and its wire count.
  | GateType Gate TypeError
  deriving (Show, Eq)

-- | Check the invariants of 'Gate' and 'QCirc', including 'typeCheck' of every
--   gate's operator. Reports the first violation in time order.
checkQCirc :: QCirc -> Either QCircError ()
checkQCirc (QCirc n gs)
  | n < 0     = Left (NegativeWidth n)
  | otherwise = mapM_ checkGate gs
  where
    checkGate gate@(Gate g ws) = do
      m <- either (Left . GateType gate) Right (typeCheck g)
      unless (m == length ws) $ Left (GateArity gate m (length ws))
      unless (all (\w -> 0 <= w && w < n) ws) $ Left (WireOutOfRange gate)
      unless (nub ws == ws) $ Left (DuplicateWire gate)

--------------------------------------------------------------------------------
-- Controls
--------------------------------------------------------------------------------

-- | The leading controls of an operator, counted, and the operator they
--   control: @peelControls (C (C a)) = (2, a)@ when @a@ is not itself a 'C'.
peelControls :: QOp -> (Int, QOp)
peelControls (C a) = let (k, b) = peelControls a in (k + 1, b)
peelControls a     = (0, a)

--------------------------------------------------------------------------------
-- Layers
--------------------------------------------------------------------------------

-- | The wires a gate acts on, in its local qubit order.
gateWires :: Gate -> [Wire]
gateWires (Gate _ ws) = ws

-- | Every wire from the lowest to the highest one a gate touches: the rows a
--   gate covers when drawn as one column of a diagram.
gateSpan :: Gate -> [Wire]
gateSpan (Gate _ []) = []
gateSpan (Gate _ ws) = [minimum ws .. maximum ws]

-- | As-soon-as-possible layering by a footprint: a gate joins the layer after
--   the latest layer of an earlier gate whose footprint meets its own, and a
--   gate with an empty footprint joins the first layer. Layers keep time order,
--   and their concatenation reorders the gates only across disjoint footprints.
layersBy :: (Gate -> [Wire]) -> [Gate] -> [[Gate]]
layersBy footprint gs =
  -- Each gate is consed onto its layer, so the layers come out reversed.
  map reverse (M.elems (M.fromListWith (++) [ (k, [g]) | (k, g) <- schedule M.empty gs ]))
  where
    schedule _ [] = []
    schedule used (g : more) = (k, g) : schedule used' more
      where
        ws = footprint g
        k = maximum (0 : [ i + 1 | w <- ws, Just i <- [M.lookup w used] ]) :: Int
        used' = foldr (`M.insert` k) used ws

-- | The gates in layers of pairwise disjoint wires, as soon as possible.
layers :: QCirc -> [[Gate]]
layers = layersBy gateWires . circ_gates

-- | The number of layers that hold a gate on at least one wire.
depth :: QCirc -> Int
depth = length . filter (not . all (null . gateWires)) . layers

--------------------------------------------------------------------------------
-- Embedding into QOp
--------------------------------------------------------------------------------

-- | The circuit as a 'QOp' on @circ_qubits@ wires, right-nested in mathematical
--   order (the last gate is the leftmost factor). A gate on consecutive
--   increasing wires is placed with identity padding, @k <@ g @> l@ without
--   zero-width identities. Any other wire list is conjugated by the permutation
--   that brings its wires to the front in the listed order.
toQOp :: QCirc -> QOp
toQOp (QCirc n gs) = case map (place n) gs of
  []  -> Id n
  ops -> foldr1 (∘) (reverse ops)

place :: Nat -> Gate -> QOp
place n (Gate g ws)
  | ws == [k .. k + m - 1] = pad k g (n - k - m)
  | otherwise = Permute (permInvert p) ∘ pad 0 g (n - m) ∘ Permute p
  where
    m = length ws
    k = case ws of
      (w : _) -> w
      []      -> 0
    p = ws ++ ([0 .. n - 1] \\ ws)  -- Permute p: position i receives wire (ws !! i).

-- | @Id k ⊗ g ⊗ Id l@, leaving out zero-width identities.
pad :: Nat -> QOp -> Nat -> QOp
pad k g l = wrap k (if l > 0 then Tensor g (Id l) else g)
  where
    wrap 0 op = op
    wrap i op = Tensor (Id i) op
