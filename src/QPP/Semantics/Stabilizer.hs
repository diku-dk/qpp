-- Exercise scaffold: the imports below are unused until the hole is filled.
{-# OPTIONS_GHC -Wno-unused-imports #-}
-- | Stabilizer semantics — an exercise hole. Fill in the backend interface
--   described in "QPP.Semantics": `StateT`, `OpT`, `ket`, `evalOp`, `apply`,
--   `measure1`, and `evalStep`/`evalProg` (the latter two via
--   `evalStepWith`/`evalProgWith`).
module QPP.Semantics.Stabilizer where

import QPP.Syntax
import QPP.Semantics
