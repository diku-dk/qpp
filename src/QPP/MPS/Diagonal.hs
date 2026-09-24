-- | Diagonal-core MPS: a specialization for states that are sums of χ
--   separable terms,
--       |ψ⟩ = scalar · Σ_α ⊗_p (f_p(0,α) |0⟩ + f_p(1,α) |1⟩)
--   All sites share one bond dimension χ; per site we store the
--   coefficients f_p(0,·), f_p(1,·) ∈ C^χ directly (no χ² storage).
--   See `perfect-sampling.tex` § "Diagonal-MPS specialization".
--
--   Layered API:
--     * `DiagMPS` is the standalone diagonal-core data structure.
--       `sampleDiagMPS` / `measureDiagMPS` operate on it directly.
--     * `sampleAllDiag` / `measureAllDiag` are MPS wrappers: assert
--       diagonality via `isDiagonalMPS`, extract via `toDiagMPS`, run
--       the raw entry, remap outcomes to logical order via `phys2log`.
module QPP.MPS.Diagonal where

import Data.Complex (Complex(..), conjugate, magnitude, realPart)
import qualified Data.Vector as V
import Data.Vector ((!), (//))
import Data.List (foldl')
import QPP.Syntax
import QPP.MPS
import QPP.Semantics.Matrix (CMat, CVec)

import Numeric.LinearAlgebra ((#>), atIndex, rows, cols, dot, cmap)
import qualified Numeric.LinearAlgebra as H
import GHC.Stack (HasCallStack)

-- | Per-site diagonal coefficients (f_p(0,·), f_p(1,·)).
data DiagSite = DiagSite { d0 :: !CVec, d1 :: !CVec } deriving (Show, Eq)

-- | A diagonal-core MPS. Invariant (enforced by `mkDiagMPS`): every
--   `DiagSite` has CVecs of length `dBondDim`.
data DiagMPS = DiagMPS
  { dScalar  :: !ComplexT
  , dBondDim :: !Int
  , dSites   :: !(V.Vector DiagSite)
  } deriving (Show, Eq)

dNSites :: DiagMPS -> Int
dNSites = V.length . dSites

instance HasQubits DiagMPS where n_qubits = dNSites

-- | Validating constructor: every site's CVecs must have equal length.
mkDiagMPS :: HasCallStack => ComplexT -> V.Vector DiagSite -> DiagMPS
mkDiagMPS s ss
  | V.null ss = DiagMPS s 0 ss
  | otherwise =
      let chi = H.size (d0 (V.head ss))
          ok  = V.all (\(DiagSite v0 v1) -> H.size v0 == chi && H.size v1 == chi) ss
      in if ok then DiagMPS s chi ss
         else error "mkDiagMPS: inconsistent bond dimension across sites"

-- | Extract diagonal coefficients from an MPS Site. Boundary sites are 1×χ
--   or χ×1; internal sites are taken to be χ×χ diagonal (caller is
--   responsible for `isDiagonalSite`-checking first).
siteToDiag :: Site -> DiagSite
siteToDiag (Site x0 x1)
  | rows x0 == 1 = DiagSite (H.flatten x0)  (H.flatten x1)
  | cols x0 == 1 = DiagSite (H.flatten x0)  (H.flatten x1)
  | otherwise    = DiagSite (H.takeDiag x0) (H.takeDiag x1)

-- | A site is diagonal-or-boundary if it's 1×χ, χ×1, or square with all
--   off-diagonals below @eps@ in magnitude.
isDiagonalSite :: Double -> Site -> Bool
isDiagonalSite eps (Site x0 x1) =
  let isDiag m =
        let r = rows m; c = cols m
            offDiagSmall =
              and [ magnitude (m `atIndex` (i,j)) < eps
                  | i <- [0 .. r-1], j <- [0 .. c-1], i /= j ]
        in (r == 1) || (c == 1) || (r == c && offDiagSmall)
  in isDiag x0 && isDiag x1

-- | Whole-MPS predicate. Boundary sites must be 1×χ / χ×1; internal sites
--   must be square and (numerically) diagonal under @tol (cfg psi)@.
isDiagonalMPS :: MPS -> Bool
isDiagonalMPS psi =
  let n      = nSites psi
      sV     = sites psi
      eps    = tol (cfg psi)
      sLeft  = sV ! 0
      sRight = sV ! (n-1)
  in case n of
       0 -> True
       1 -> True   -- 1×1
       _ -> rows (a0 sLeft) == 1
         && cols (a0 sRight) == 1
         && all (\p -> isDiagonalSite eps (sV ! p)) [1 .. n-2]

-- | Convert an MPS to a DiagMPS. Asserts diagonality. The DiagMPS sites
--   are in *physical* order; the MPS's `phys2log` mapping (if non-identity)
--   is not preserved here — `sampleAllDiag` / `measureAllDiag` apply the
--   remap when wrapping back to the MPS API.
toDiagMPS :: HasCallStack => MPS -> DiagMPS
toDiagMPS psi
  | not (isDiagonalMPS psi) = error "toDiagMPS: state is not diagonal"
  | otherwise =
      let n   = nSites psi
          chi = if n == 0 then 0 else cols (a0 (sites psi ! 0))
      in DiagMPS { dScalar  = scalar psi
                 , dBondDim = chi
                 , dSites   = V.map siteToDiag (sites psi)
                 }

-- | Embed a DiagMPS as a generic MPS in standard chain form: 1×χ at the
--   left boundary, χ×χ diagonal internally, χ×1 at the right boundary
--   (1×1 for n=1). `log2phys` is the identity.
--
--   NOTE: the result is generally *not* in canonical form (diagonal cores are
--   not isometries), so `center_site` is nominal. It is safe for
--   center-agnostic consumers (`mpsToDenseVec`, `sampleAllDiag`,
--   `measureAllDiag`, applying unitaries); the center-based probability paths
--   (`measure1`, `sampleAll`) require genuinely canonical input.
fromDiagMPS :: DiagMPS -> MPS
fromDiagMPS dm =
  let n   = dNSites dm
      mk p (DiagSite v0 v1)
        | p == 0     = Site (H.asRow v0)    (H.asRow v1)
        | p == n-1   = Site (H.asColumn v0) (H.asColumn v1)
        | otherwise  = Site (H.diag v0)     (H.diag v1)
      sV  = V.imap mk (dSites dm)
      idm = V.generate n id
  in MPS { scalar      = dScalar dm
         , sites       = sV
         , center_site = 0
         , log2phys    = idm
         , phys2log    = idm
         , cfg         = defaultCfg
         }

-- | G_p[α,β] = Σ_s f_p(s,α)* · f_p(s,β) — Hermitian, rank ≤ 2.
--   (Conjugation on the *first* index: this is what falls out of |⟨s|ψ⟩|²
--   when integrating out one site.)
gMatrix :: DiagSite -> CMat
gMatrix (DiagSite v0 v1) =
  let outerC u = H.outer (cmap conjugate u) u
  in outerC v0 + outerC v1

-- | Right-chain @R[p][α,β] = ∏_{q > p} G_q[α,β]@, computed as elementwise
--   products (a right scan). @R[n-1]@ is the all-ones χ×χ matrix.
buildRChain :: Int -> V.Vector DiagSite -> V.Vector CMat
buildRChain chi sV =
  let jOnes = H.konst (1:+0) (chi, chi)
  in V.scanr' (\site acc -> gMatrix site * acc) jOnes (V.drop 1 sV)

-- | Below this total weight the sampling recurrences give up: the state (or
--   the conditioned branch) is numerically zero. Deliberately far below any
--   truncation tolerance — a legitimate small-probability branch should
--   sample, not error.
probFloor :: Double
probFloor = 1e-300

-- | Sample all qubits from a DiagMPS. Outcomes returned head-most-recent,
--   matching @Measure [0..n-1]@. Cost: O(n · χ²); working memory O(n · χ²)
--   for the right chain.
sampleDiagMPS :: HasCallStack => DiagMPS -> RNG -> (Outcomes, RNG)
sampleDiagMPS dm rng0
  | n == 0    = ([], rng0)
  | otherwise =
      let chi    = dBondDim dm
          rChain = buildRChain chi sV
          lInit  = H.konst (1:+0) chi
          (_, rng', bitsRev) =
            foldl' (step rChain) (lInit, rng0, []) [0 .. n-1]
      in (reverse bitsRev, rng')
  where
    sV = dSites dm
    n  = V.length sV
    step :: V.Vector CMat
         -> (CVec, RNG, [Bool])
         -> Int
         -> (CVec, RNG, [Bool])
    step _ (_, [], _) _ = error "sampleDiagMPS: empty RNG"
    step rChain (l, r:rs, acc) p =
      let DiagSite f0 f1 = sV ! p
          h0  = l * f0
          h1  = l * f1
          r_p = rChain ! p
          w0  = realPart (dot h0 (r_p #> h0))
          w1  = realPart (dot h1 (r_p #> h1))
          tot = w0 + w1
      in if tot < probFloor
           then error "sampleDiagMPS: prob ~ 0"
           else let bit = r * tot >= w0
                    l'  = if bit then h1 else h0
                in (l', rs, bit : acc)

-- | Projective measure-all on a DiagMPS. Returns a new DiagMPS whose
--   `dScalar` is set to 1/|⟨s|ψ_internal⟩|, so the result has unit norm
--   regardless of the input scalar (phase is unobservable post-collapse).
--   The unchosen branch is zeroed at every site.
measureDiagMPS :: HasCallStack => DiagMPS -> RNG -> (DiagMPS, Outcomes, RNG)
measureDiagMPS dm rng0
  | n == 0    = (dm, [], rng0)
  | otherwise =
      let chi    = dBondDim dm
          rChain = buildRChain chi sV
          zerof  = H.konst (0:+0) chi
          lInit  = H.konst (1:+0) chi
          (lFinal, rng', bitsRev, updates) =
            foldl' (step rChain zerof) (lInit, rng0, [], []) [0 .. n-1]
          amp     = H.sumElements lFinal
          renorm  = (1 / magnitude amp) :+ 0
          dm'     = DiagMPS { dScalar = renorm, dBondDim = chi, dSites = sV V.// updates }
      in (dm', reverse bitsRev, rng')
  where
    sV = dSites dm
    n  = V.length sV
    step :: V.Vector CMat
         -> CVec
         -> (CVec, RNG, [Bool], [(Int, DiagSite)])
         -> Int
         -> (CVec, RNG, [Bool], [(Int, DiagSite)])
    step _ _ (_, [], _, _) _ = error "measureDiagMPS: empty RNG"
    step rChain zerof (l, r:rs, acc, upd) p =
      let DiagSite f0 f1 = sV ! p
          h0  = l * f0
          h1  = l * f1
          r_p = rChain ! p
          w0  = realPart (dot h0 (r_p #> h0))
          w1  = realPart (dot h1 (r_p #> h1))
          tot = w0 + w1
      in if tot < probFloor
           then error "measureDiagMPS: prob ~ 0"
           else let bit = r * tot >= w0
                    l'  = if bit then h1 else h0
                    s'  = if bit then DiagSite zerof f1 else DiagSite f0 zerof
                in (l', rs, bit : acc, (p, s') : upd)

-- | Sample-only all-qubit measurement on a diagonal MPS state. Asserts
--   diagonality; outcomes returned in logical-qubit order.
sampleAllDiag :: HasCallStack => MPS -> RNG -> (Outcomes, RNG)
sampleAllDiag psi rng =
  let dm               = toDiagMPS psi
      (physOuts, rng') = sampleDiagMPS dm rng
      n                = nSites psi
      bits = V.replicate n False //
               [ (phys2log psi ! p, b) | (p, b) <- zip [0..] physOuts ]
      outs = [ bits ! k | k <- [0 .. n-1] ]
  in (outs, rng')

-- | Projective all-qubit measurement on a diagonal MPS. The collapsed state
--   is the computational-basis product state given by the outcomes (the
--   post-collapse phase is unobservable), so the post-state is rebuilt as a
--   bond-dimension-1 basis ket: unit norm and canonical by construction.
measureAllDiag :: HasCallStack => MPS -> RNG -> (MPS, Outcomes, RNG)
measureAllDiag psi rng =
  let dm                 = toDiagMPS psi
      (_, physOuts, rng') = measureDiagMPS dm rng
      n    = nSites psi
      bits = V.replicate n False //
               [ (phys2log psi ! p, b) | (p, b) <- zip [0..] physOuts ]
      outs = [ bits ! k | k <- [0 .. n-1] ]
      psi' = (ket (map fromEnum outs)) { cfg = cfg psi }
  in (psi', outs, rng')
