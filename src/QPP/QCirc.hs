-- | Flat quantum circuits: syntax, extraction from QOp, embedding into QOp,
--   and printing.
--
--   Import qualified alongside "QPP", as with "QPP.UOp": 'fromQOp' and 'toQOp'
--   are generic names, and the circuit type uses the same operator typeclasses.
module QPP.QCirc
  ( module QPP.QCirc.Syntax
  , module QPP.QCirc.FromQOp
  , module QPP.PrettyPrint.QCirc
  ) where

import QPP.QCirc.Syntax
import QPP.QCirc.FromQOp
import QPP.PrettyPrint.QCirc
