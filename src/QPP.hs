-- | Umbrella module for the symbolic layer: syntax, rewriting and
--   pretty-printing.
--
--   No semantics is re-exported here on purpose — a program runs against
--   exactly one backend, and the choice must be visible in the import list.
--   Import one `QPP.Semantics.*` module (Matrix, Statevector, MPS,
--   Stabilizer) alongside this one; see "QPP.Semantics" for the convention
--   those modules follow.
module QPP(
    module QPP.Syntax,
    module QPP.Rewrite,
    module QPP.PrettyPrint
  ) where

import QPP.Syntax
import QPP.Rewrite
import QPP.PrettyPrint
