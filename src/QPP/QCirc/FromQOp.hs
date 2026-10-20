-- | Extract a flat circuit from a QOp by reading gates off its layered form.
module QPP.QCirc.FromQOp (fromQOp) where

import QPP.QCirc.Syntax
import QPP.Rewrite
  ( assocComposes, assocTensors, cleanAdjoints, cleanOnes, liftComposes, simplifyFixpoint )
import QPP.Syntax

-- | Flat circuit of a well-formed QOp. Compositions are lifted out of tensor
--   products and direct sums, adjoints are pushed to the leaves, and trivial
--   identities are removed. Each factor of the resulting outer composition is
--   one parallel layer, and every non-identity tensor factor in it becomes a
--   gate on its consecutive window of wires. Gates are neither fused, split,
--   nor compacted to their support: @C (Id 2 ⊗ X)@ stays one four-wire gate,
--   and @H ∘ H@ gives two gates. A zero-qubit 'Phase' factor becomes a gate on
--   no wires. Assumes well-formed input, like 'op_qubits'; run 'typeCheck'
--   first to validate.
fromQOp :: QOp -> QCirc
fromQOp op =
  QCirc (op_qubits op) (concatMap layerGates (reverse (assocComposes (layerForm op))))

-- | Compositions of tensor layers. 'pushComposes' is left out because it undoes
--   the lifting; 'doComposes' is left out to keep the gate list faithful.
layerForm :: QOp -> QOp
layerForm = simplifyFixpoint 10000 [cleanOnes, cleanAdjoints, liftComposes]

-- | Gates of one parallel layer, with a running wire offset.
layerGates :: QOp -> [Gate]
layerGates = go 0 . assocTensors
  where
    go _ [] = []
    go k (Id m : rest) = go (k + m) rest
    go k (g : rest) = Gate g [k .. k + m - 1] : go (k + m) rest
      where m = op_qubits g
