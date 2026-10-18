-- | General finite-dimensional unitary syntax, lowering, normalization, and printing.
--
-- Import qualified alongside "QPP", since the two syntaxes share constructor
-- names such as Tensor and Compose. They use the same operator typeclasses.
-- 'fromQOp' checks its input before lowering it. General UOp type-checking
-- and semantics are separate work.
module QPP.UOp
  ( module QPP.UOp.Syntax
  , module QPP.UOp.FromQOp
  , module QPP.UOp.Normalize
  , module QPP.PrettyPrint.UOp
  , showTikZUOp, showTikZUOpWithoutQubits
  ) where

import QPP.UOp.Syntax
import QPP.UOp.FromQOp
import QPP.UOp.Normalize
import QPP.PrettyPrint.UOp
import QPP.PrettyPrint.TikZ (showTikZUOp, showTikZUOpWithoutQubits)
