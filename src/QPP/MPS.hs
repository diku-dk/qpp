{-# LANGUAGE BangPatterns #-}

-- | Matrix-product-state representation: the data type, its canonical-form
--   machinery (center moves, compression, site swaps), state arithmetic
--   (tensor, add, inner), measurement/sampling, and dense/sparse conversions.
--
--   The `QOp` evaluator built on top of this lives in "QPP.Semantics.MPS";
--   the diagonal-core specialization in "QPP.MPS.Diagonal". This module is the
--   implementation layer and exports everything.
module QPP.MPS where

import Data.Complex (Complex(..), conjugate, magnitude, realPart)
import qualified Data.Vector as V
import Data.Vector ((!), Vector, (//))
import qualified Data.Vector.Generic as G
import Data.Bits (shiftL, testBit, (.|.))
import Data.List (foldl')
import QPP.Syntax
import QPP.Util
import QPP.Semantics hiding (tol)   -- this module's `tol` is the EvalCfg field
import QPP.Semantics.Matrix (CMat, CVec)
import qualified QPP.Semantics.Matrix as MS
import qualified Data.PQueue.Prio.Min as PriorityQ

import Numeric.LinearAlgebra (
  (><), (<#), (#>), atIndex, tr, rows, cols, size, diagBlock, compactSVD, dot,
  Matrix, Element, subMatrix
  )
import qualified Numeric.LinearAlgebra as H
import Numeric.IEEE(epsilon)
import GHC.Stack (HasCallStack)

-- MPS representation
data Interval = Ival !Int !Int deriving (Show,Eq)
hull :: Interval -> Interval -> Interval
hull (Ival l r) (Ival l' r') = Ival (min l l') (max r r')

singletonIval :: Int -> Interval
singletonIval p = Ival p p

-- A : Dl x 2 x Dr as (A0,A1), each Dl x Dr.
data Site = Site { a0 :: !CMat, a1 :: !CMat } deriving (Show,Eq)

-- | MPS in physical chain order; `Permute` updates only the wire maps.
--
--   Center invariant: `center_site` is a valid orthogonality center — every
--   site left of it is a left isometry (A₀†A₀ + A₁†A₁ = I) and every site
--   right of it a right isometry (A₀A₀† + A₁A₁† = I). All probability
--   computations (`measure1`, `sampleAll`) rely on this; every operation
--   exported from this module preserves it (`QPP.MPS.Diagonal.fromDiagMPS` is
--   the documented exception). The state's norm and phase live in the center
--   site and `scalar`.
data MPS = MPS
  { scalar     :: !ComplexT
  , sites      :: !(Vector Site)
  , center_site :: !Int
  , log2phys   :: !(Vector Int)
  , phys2log   :: !(Vector Int)
  , cfg        :: !EvalCfg
  } deriving (Show,Eq)

nSites :: MPS -> Int
nSites = V.length . sites

instance HasQubits MPS where n_qubits = nSites

invertVec :: Vector Int -> Vector Int
invertVec v =
  let n = V.length v
  in  (V.replicate n 0) // [ (p,q) | (q,p) <- zip [0..] (V.toList v) ]

-- | The class method discards the profile (fixed signature); use `tensorMPS`
--   directly when the statistics matter.
instance HasTensorProduct MPS where (⊗) a b = fst (tensorMPS a b)

{-| tensor product a ⊗ b for MPS a,b -}
tensorMPS :: MPS -> MPS -> (MPS, Profile)
tensorMPS a b
  | nSites a == 0 = (scaleMPS (scalar a) b { cfg = cfg' }, mempty)
  | nSites b == 0 = (scaleMPS (scalar b) a { cfg = cfg' }, mempty)
  | otherwise =
      let (na,nb) = (nSites a, nSites b)
          (a', pa) = moveCenterToPhys (na-1) a
          (b', pb) = moveCenterToPhys 0 b
          l2p = V.generate (na+nb) $ \q -> if q < na then log2phys a ! q else na + log2phys b ! (q-na)
          merged = MPS { scalar      = scalar a * scalar b,
                         sites       = sites a' V.++ sites b',
                         center_site = na-1,
                         log2phys = l2p,
                         phys2log = invertVec l2p,
                         cfg = cfg'
                         }
          -- a' (center at na−1) and b' (center at seam site na) each contribute a
          -- non-isometric site; unless b happens to be normalized, the seam site
          -- breaks the canonical form. One moveRight across the seam absorbs b's
          -- center into a's, leaving a single valid center at na.
          (out, pSeam) = moveRight (na-1) merged
      in (out, pa <> pb <> pSeam)
  where
    cfg' = truncMax (cfg a) (cfg b) -- Worst accuracy guarantee determines overall accuracy

{-| Combined truncation level: If we combine two MPS, use the least accurate truncation level -}
truncMax :: EvalCfg -> EvalCfg -> EvalCfg
truncMax c1 c2 = EvalCfg
  { trunc = case (trunc c1, trunc c2) of
      (Exact, t) -> t
      (t, Exact) -> t
      (Truncate m1 e1, Truncate m2 e2) -> Truncate (min m1 m2) (max e1 e2)
  , tol = max (tol c1) (tol c2)
  }

-- | Bring two states into the same wire frame (normalizing both site orders
--   if they disagree), reporting the SVD work done.
withSameFrame :: String -> MPS -> MPS -> ((MPS, MPS), Profile)
withSameFrame ctx x y
  | nSites x /= nSites y      = error (ctx ++ ": arity mismatch")
  | log2phys x == log2phys y  = ((x, y), mempty)
  | otherwise                 = let (x', px) = normalizeSiteOrder x
                                    (y', py) = normalizeSiteOrder y
                                in ((x', y'), px <> py)

scaleMPS :: ComplexT -> MPS -> MPS
scaleMPS c m = m { scalar = c * scalar m }

absorbScalarAt :: Int -> MPS -> MPS
absorbScalarAt p m
  | scalar m == (1:+0) = m
  | p < 0 || p >= nSites m =
      error ("absorbScalarAt: site " ++ show p ++ " of " ++ show (nSites m))
  | otherwise =
      let c = scalar m
          s = sites m ! p
          s' = Site (c .* (a0 s)) (c .* (a1 s))
      in m { scalar = 1:+0, sites = sites m // [(p,s')] }

(.*.) :: CMat -> CMat -> CMat
(.*.) = (H.<>) -- HMatrix <> conflicts with Prelude.<>

-- | Inner product. Any site-order normalization work is discarded (the
--   result is a scalar, and `inner`'s class signature has nowhere to put a
--   profile).
innerMPS :: HasCallStack => MPS -> MPS -> ComplexT
innerMPS x0 y0 =
  let ((x, y), _) = withSameFrame "innerMPS" x0 y0
      n  = nSites x
      (sx, sy) = (scalar x, scalar y)
      e0 = foldl' step ((1><1) [1:+0]) [n-1, n-2 .. 0] -- Full sweep for contraction
      -- e is the right environment ("environment" in MPS lingo): a map from
      -- x's right-bond space to y's right-bond space. Site matrices map
      -- left bond -> right bond, so contracting site p sandwiches e as
      -- y_s · e · x_s†, yielding the environment one bond to the left.
      step e p =
        let Site xa0 xa1 = sites x ! p
            Site ya0 ya1 = sites y ! p
            term xs ys = ys .*. e .*. tr xs
        in term xa0 ya0 + term xa1 ya1 -- Contract over physical index s ∈ {0,1}
  in conjugate sx * sy * (e0 `atIndex` (0,0))

-- | The class methods have fixed signatures (`v -> v -> v`), so state
--   arithmetic done through this instance discards its (small) profile;
--   call `addMPS` / `tensorMPS` directly when the statistics matter.
instance HilbertSpace MPS where
  type Realnum MPS = Double
  type Scalar  MPS = ComplexT

  (.*) c ψ = scaleMPS c ψ
  (.+) ψ φ = fst (addMPS ψ φ)
  (.-) ψ φ = fst (addMPS ψ (((-1):+0) .* φ))

  inner = innerMPS
  normalize ψ = let
     nrm = norm ψ
    in
      if nrm < tol (cfg ψ) then ψ else ((1/nrm):+0) .* ψ

data Trunc = Exact | Truncate { maxBond :: !Int, svd_r :: !Double }
        deriving (Show,Eq)

data EvalCfg    = EvalCfg { trunc :: !Trunc, tol :: !Double } deriving (Show,Eq)

defaultCfg :: EvalCfg
defaultCfg = EvalCfg { trunc = Exact, tol = 1e-12 }

approxCfg :: Int -> Double -> EvalCfg
approxCfg maxbond svd_r = EvalCfg
  { trunc = Truncate { maxBond = maxbond, svd_r = svd_r }
  , tol = 1e-12
  }

-- | Statistics of an evaluation *run* — not of a state. Emitted once per SVD
--   by `svd_compact` and combined monoidally along the actual execution
--   order, so nothing is dropped at branch merges and shared history is never
--   double-counted. (See docs/MPS-profiling-redesign.md.)
data Profile = Profile
  { nSVDs           :: !Int
  , maxBondDim      :: !Int    -- ^ largest bond dimension seen at any SVD
  , discardedWeight :: !Double -- ^ Σ δᵢ, δᵢ = Σ_{k>χ} σ_k² per truncation
  , normError       :: !Double -- ^ Σ √δᵢ: 2-norm error bound across unitary
                               --   segments and adds; heuristic across
                               --   measurement renormalization
  } deriving (Show, Eq)

instance Semigroup Profile where
  Profile n b w e <> Profile n' b' w' e' =
    Profile (n + n') (max b b') (w + w') (e + e')

instance Monoid Profile where
  mempty = Profile 0 0 0 0

-- | Basis ket (logical MSB-first). All sites are 1×1 and norm 1, so any
--   center position is valid.
ket :: [Int] -> MPS
ket bs =
  let n = length bs
      mk b = let v0 = if b==0 then 1:+0 else 0:+0
             in Site ((1><1) [v0]) ((1><1) [1-v0])
      ss  = V.fromList (map mk bs)
      idm = V.generate n id
  in MPS { scalar = 1:+0, sites = ss, center_site = 0, log2phys = idm, phys2log = idm, cfg = defaultCfg }

-- wire permutation: logical relabeling only
applyPermute :: HasCallStack => Int -> [Int] -> MPS -> MPS
applyPermute base π m =
  let n    = length π
      l2p  = log2phys m
      sl   = V.slice base n l2p
      sl'  = V.fromList [ sl ! k | k <- π ]
      l2p' = l2p // [ (base+i, sl' ! i) | i <- [0..n-1] ]
      p2l' = invertVec l2p'
  in m { log2phys = l2p',
         phys2log = p2l' }

-- small matrix helpers
hcat, vcat :: CMat -> CMat -> CMat
hcat = (H.|||)
vcat = (H.===)

-- | Block splits of a dense matrix. Used to cut the two-site tensor Θ back
--   into a pair of `Site`s after an SVD (`split2x1` on U, `split1x2` on SV†)
--   and to swap Θ's off-diagonal blocks in `swapSites`.
split2x2 :: Element e
         => Int      -- nr row split
         -> Int      -- nc column split
         -> Matrix e -- m  input matrix
         -> (Matrix e, Matrix e, Matrix e, Matrix e)
split2x2 nr nc m = let
    (r,c) = (rows m, cols m)
  in
    ( subMatrix (0,0)   (nr,  nc)   m,
      subMatrix (0,nc)  (nr,  c-nc) m,
      subMatrix (nr,0)  (r-nr,nc)   m,
      subMatrix (nr,nc) (r-nr,c-nc) m  )

split2x1 :: Element e
         => Int      -- nr row split
         -> Matrix e -- m  input matrix
         -> (Matrix e, Matrix e)
split2x1 nr m = let
    (r,c) = (rows m, cols m)
  in
    ( subMatrix (0,0)   (nr,  c) m,
      subMatrix (nr,0)  (r-nr,c) m  )

split1x2 :: Element e
         => Int      -- nc column split
         -> Matrix e -- m  input matrix
         -> (Matrix e, Matrix e)
split1x2 nc m = let
    (r,c) = (rows m, cols m)
  in
    ( subMatrix (0,0)   (r,  nc)   m,
      subMatrix (0,nc)  (r,  c-nc) m  )

-- | s descending. Return first i with s!i < cutoff, or n if none.
--   Locates the truncation rank in a singular-value spectrum.
firstBelow :: (G.Vector v a, Ord a) => a -> v a -> Int
firstBelow cutoff s = go 0 (G.length s)
  where
    go !lo !hi
      | lo >= hi           = lo
      | s G.! mid < cutoff = go lo mid
      | otherwise          = go (mid+1) hi
      where
        !mid = (lo+hi) `div` 2

{-| Dense matrix representation of two adjacent sites
               [ A0 ]
     Θ2(A,B) = [ A1 ] [B0 B1] -}
theta2 :: Site -> Site -> CMat
theta2 a b = vcat (a0 a)  (a1 a) .*. hcat (a0 b) (a1 b)

diagMulLeft :: H.Vector Double -> CMat -> CMat
diagMulLeft v m = (H.complex . H.asColumn $ v) * m

diagMulRight :: CMat -> H.Vector Double -> CMat
diagMulRight m v = m * (H.complex . H.asRow $ v)

-- | Truncating SVD. Emits its `Profile` delta — the only place profiles are
--   created; everything else merely accumulates them (`<>`) in execution
--   order. Always decomposes via `compactSVD`: hmatrix computes the full thin
--   SVD in every variant (`compactSVDTol` is a column post-filter over
--   `thinSVD`), so the discarded weight is available at zero extra cost and
--   error tracking is not a mode.
svd_compact :: EvalCfg -> CMat -> ((CMat, H.Vector Double, CMat), Profile)
svd_compact cfg m = case trunc cfg of
  Exact -> let (u,s,v) = compactSVD m
           in ((u, s, v), Profile 1 (size s) 0 0)

  Truncate{maxBond,svd_r} ->
    let
      (u,s,v) = compactSVD m
      -- Truncate internal bond dimension to χ. All three factors must be cut
      -- consistently: callers contract u/s/v against each other.
      u' = H.takeColumns chi u
      v' = H.takeColumns chi v
      s' = H.subVector 0 chi s
      δ  = H.sumElements (H.subVector chi (size s - chi) (s*s))

      -- Same relative threshold hmatrix's compactSVDTol would use.
      svd_tol = svd_r*g*epsilon*k where g = H.norm_Inf s
                                        k = fromIntegral (max (rows m) (cols m))
      chi     = min maxBond (firstBelow svd_tol s)
    in ((u', s', v'), Profile 1 chi δ (sqrt δ))

moveRight :: Int -> MPS -> (MPS, Profile)
moveRight j m
  | j < 0 || j+1 >= nSites m =
      error ("moveRight: site " ++ show j ++ " of " ++ show (nSites m))
  | otherwise =
      let
          (a, b)        = (sites m ! j, sites m ! (j+1))
          (dl, dr)      = (rows (a0 a), cols (a0 b)) -- External bond dimensions
          ((u,s,v), p)  = svd_compact (cfg m) (theta2 a b)

          -- Absorb singular values into B'. Θ = U S V† = A' B'  =>  B' = S V†
          sv   = diagMulLeft s (tr v)
          a'   = uncurry Site (split2x1 dl u) -- Left Isometry A'
          b'   = uncurry Site (split1x2 dr sv)

      in (m { sites = sites m // [(j,a'),(j+1,b')], center_site = j+1 }, p)

moveLeft :: Int -> MPS -> (MPS, Profile)
moveLeft j m
  | j <= 0 || j >= nSites m =
      error ("moveLeft: site " ++ show j ++ " of " ++ show (nSites m))
  | otherwise =
      let i = j-1
          (a,b)        = (sites m ! i, sites m ! j)
          (dl, dr)     = (rows (a0 a), cols (a0 b))
          ((u,s,v), p) = svd_compact (cfg m) (theta2 a b)

          -- Absorb singular values into A'. Θ = U S V† = A' B'  =>  A' = U S
          us  = diagMulRight u s
          a'  = uncurry Site (split2x1 dl us)
          b'  = uncurry Site (split1x2 dr (tr v))
      in (m { sites = sites m // [(i,a'),(j,b')], center_site = j-1 }, p)


-- TODO: This can be done cheaper by QR-decomposition.
--       If inv(p) is quadratic, we can reduce to O(n^2) QR's and O(n) SVD's.
-- | Swap physical sites j and j+1 (contents *and* wire maps). The two-site
--   SVD is only optimal — and the result's canonical form only valid — when
--   the orthogonality center lies in the swapped pair, so the center is moved
--   to j first; it ends at j+1 (B' = SV† holds the singular values).
swapSites :: MPS -> Int -> (MPS, Profile)
swapSites psi0 j
  | j < 0 || j+1 >= nSites psi0 =
      error ("swapSites: index " ++ show j ++ " of " ++ show (nSites psi0) ++ " sites")
  | otherwise =
      let (psi, pMove) = moveCenterToPhys j psi0
          (a,  b)  = (sites psi ! j, sites psi ! (j+1))
          (dl, dr) = (rows (a0 a), cols (a0 b))

          (a0b0,a0b1,
           a1b0,a1b1) = split2x2 dl dr $ theta2 a b

          theta' = H.fromBlocks [[a0b0,a1b0],
                                 [a0b1,a1b1]] -- Swap off-diagonal blocks

          ((u,s,v), pSvd) = svd_compact (cfg psi) theta'

          a' = uncurry Site (split2x1 dl u)
          sv = diagMulLeft s (tr v)
          b' = uncurry Site (split1x2 dr sv)

          (qa, qb) = (phys2log psi ! j, phys2log psi ! (j+1))
          l2p'     = log2phys psi // [(qa, j+1), (qb, j)]
          p2l'     = phys2log psi // [(j, qb), (j+1, qa)]
      in ( psi { sites = sites psi // [(j,a'),(j+1,b')]
               , center_site = j+1
               , log2phys = l2p', phys2log = p2l' }
         , pMove <> pSvd )

permutePhysicalSwaps :: MPS -> [Int] -> (MPS, Profile)
permutePhysicalSwaps = -- inv(pi) swaps (w/ SVD), so O(n^2 χ^3) worst case.
  foldl' (\(m, p) j -> let (m', p') = swapSites m j in (m', p <> p'))
    . (, mempty)

-- | Reorder the "physical" qubit sites to match the logical qubit order. This is needed before
--   adding two MPS together.
normalizeSiteOrder :: MPS -> (MPS, Profile)
normalizeSiteOrder psi =
  let swaps      = permutationSwaps (phys2log psi)
      (psi', p)  = permutePhysicalSwaps psi swaps
      ident      = V.generate (nSites psi) id
  in (psi' { log2phys = ident, phys2log = ident }, p)

moveCenterToPhys :: Int -> MPS -> (MPS, Profile)
moveCenterToPhys p0 = go mempty where
  go !acc !m
    | center_site m < p0 = step acc (moveRight (center_site m) m)
    | center_site m > p0 = step acc (moveLeft  (center_site m) m)
    | otherwise          = (m, acc)
  step acc (m', p) = go (acc <> p) m'

compressRange :: Interval -> MPS -> (MPS, Profile)
compressRange (Ival l r) m
  | r <= l    = (m, mempty)
  | otherwise = foldl' (\(acc, p) j -> let (acc', p') = moveRight j acc
                                       in (acc', p <> p'))
                       (m, mempty) [l..r-1]

-- local 1-qubit gate on a physical site: A'_s = sum_t U_{s,t} A_t
type Gate1 = (ComplexT,ComplexT,ComplexT,ComplexT) -- (u00,u01,u10,u11)

-- | A unitary Gate1 preserves left/right isometry of the site it acts on
--   (conjugation by U on the physical index), so the canonical form and
--   center survive unchanged.
apply1Phys :: Gate1 -> Int -> MPS -> MPS
apply1Phys (u00,u01,u10,u11) p m =
  let s = sites m ! p
      a0' = u00 .* (a0 s) + u01 .* (a1 s)
      a1' = u10 .* (a0 s) + u11 .* (a1 s)
  in m { sites = sites m // [(p, Site a0' a1')] }

apply1Logical :: Int -> Int -> Gate1 -> MPS -> MPS
apply1Logical base k u m = apply1Phys u (log2phys m ! (base+k)) m

cisPi :: Rational -> ComplexT
cisPi q = let t = pi * fromRational q in cos t :+ sin t

-- | Physical hull of positions where `phys2log[p] /= p`: exactly the sites
--   `normalizeSiteOrder` can mutate via `swapSites`, since the bubble sort
--   only sweeps between the lowest and highest misplaced positions. Reflects
--   ALL non-canonical state accumulated so far (e.g. log2phys inherited from
--   an outer `Permute`), so it can extend beyond any single op's
--   `op_support`. Returns `Nothing` when the site order is already identity.
normalizeTouchInterval :: MPS -> Maybe Interval
normalizeTouchInterval m =
  let p2l = phys2log m
      n   = nSites m
      bad = [ p | p <- [0 .. n - 1], p2l ! p /= p ]
  in case bad of
       [] -> Nothing
       _  -> Just (Ival (minimum bad) (maximum bad))

-- local addition with branch-index only on a hull, followed by local compression
addLocal :: Interval -> MPS -> MPS -> (MPS, Profile)
addLocal (Ival l0 r0) ψ0 φ0
    -- `scalar` is frame-independent, so the zero shortcuts fire before any
    -- normalization work is done (or counted).
  | scalar ψ0 == 0 = (φ0, mempty)
  | scalar φ0 == 0 = (ψ0, mempty)
  | otherwise =
    -- The caller's `[l0, r0]` is derived from `op_support` — sites the current
    -- op touches. It must be widened before merging:
    --
    --  * When ψ0 and φ0 disagree on `log2phys`, `withSameFrame` invokes
    --    `normalizeSiteOrder`, whose `swapSites` bubble sort mutates the full
    --    range of phys2log misplacement (its `normalizeTouchInterval`) plus
    --    the path of the center moves it performs — including misplacement
    --    inherited from ancestor ops, beyond any single op's support. After
    --    normalize, ψ and φ may disagree at sites the caller never named, so
    --    the merge's "outside [l, r] = identical" promise would break.
    --
    --  * Both inputs' orthogonality centers must lie inside the merged-and-
    --    recompressed window: sites outside it are inherited verbatim, and
    --    they are only left/right isometries — as the final `center_site = r`
    --    claims — if each input's own center was inside [l, r]. (A stale
    --    center gives *wrong measurement probabilities* later, even though
    --    the state vector itself stays correct.)
    --
    -- Zero-cost when both inputs are canonical with centers in [l0, r0].
    let Ival l r =
          foldr hull (Ival l0 r0) $
            [ singletonIval (center_site ψ0), singletonIval (center_site φ0) ]
            ++ maybe [] (:[]) (normalizeTouchInterval ψ0)
            ++ maybe [] (:[]) (normalizeTouchInterval φ0)

        ((ψ, φ), pFrame) = withSameFrame "addLocal" ψ0 φ0

        (ψ1,φ1) = (absorbScalarAt l ψ, absorbScalarAt l φ)
        (sψ, sφ) = (sites ψ1, sites φ1)

        -- Inside [l, r], combine ψ.site and φ.site by a position-dependent
        -- combinator: singleton support adds; left edge concatenates
        -- horizontally (entering Dl×2Dr); right edge concatenates vertically
        -- (exiting 2Dl×Dr); internal sites form a 2×2 block diagonal.
        -- Outside [l, r] the caller promised ψ.site p == φ.site p, so we
        -- take ψ.
        mk p
          | p < l || p > r = sψ ! p
          | otherwise =
              let a = sψ ! p; b = sφ ! p
                  combine = case (compare p l, compare p r) of
                    (EQ, EQ) -> (+)                      -- singleton support
                    (EQ, _ ) -> hcat                     -- left edge
                    (_,  EQ) -> vcat                     -- right edge
                    _        -> \x y -> diagBlock [x, y] -- internal
              in Site (combine (a0 a) (a0 b)) (combine (a1 a) (a1 b))

        out = ψ1 { scalar = 1:+0,
                   sites  = V.generate (nSites ψ1) mk }

        (out', pCompress) = compressRange (Ival l r) out
  in (out', pFrame <> pCompress)

-- | General state addition: the inputs may differ anywhere, so merge over the
--   full chain. (The evaluator's internal adds go through `addLocal` with the
--   op-support interval instead.)
addMPS :: MPS -> MPS -> (MPS, Profile)
addMPS ψ φ = addLocal (Ival 0 (max 0 (nSites ψ - 1))) ψ φ

-- beam-search sparse extraction (fast for low-entanglement states)
toSparseMat :: Double -> Int -> MPS -> SparseMat
toSparseMat eps maxTerms t0 =
  let mps  = fst (moveCenterToPhys 0 t0) -- ensures partial path amplitudes are strict bounds (yielding exact largest amplitudes)
      n    = nSites mps
      dim  = 2^(fromIntegral n :: Integer)
      eps2 = eps * eps

      step :: [(Integer, CVec)] -> Int -> [(Integer, CVec)]
      step states p =
        let Site a0 a1 = sites mps ! p
            q          = phys2log mps ! p
            bitpos     = (n - 1) - q
            branches   = [(0, a0), (1, a1)] :: [(Integer,CMat)]

            extend1
              :: (Integer, CMat)
              -> (Integer, CVec)
              -> PriorityQ.MinPQueue Double (Integer, CVec)
              -> PriorityQ.MinPQueue Double (Integer, CVec)
            extend1 (s,a) (!idx,!v) beam0 =
              let !v' = v <# a
                  !w2 = realPart $ dot v' v'
              in  if w2 <= eps2
                    then beam0 -- Continuing on this path cannot yield amplitude larger than eps2
                    else let !idx' = idx .|. (s `shiftL` bitpos)   -- Feasible candidate:
                         in  pushTopK maxTerms w2 (idx', v') beam0 -- push to top-k priority queue

            extendState
              :: PriorityQ.MinPQueue Double (Integer, CVec)
              -> (Integer, CVec)
              -> PriorityQ.MinPQueue Double (Integer, CVec)
            extendState beam0 st = foldl' (\b br -> extend1 br st b) beam0 branches

            beam = foldl' extendState PriorityQ.empty states
        in  map snd (PriorityQ.toDescList beam)

      finals = -- Branch and bound on path through sites from left to right
        foldl' step [(0 :: Integer, H.fromList [scalar mps])] [0 .. n-1]

      nz = [ ((i,0), v `atIndex` 0) | (i,v) <- finals ]
  in SparseMat ((dim,1), nz)


-- | Beam width for the default sparse extraction: keep at most this many
--   basis states. Amplitudes beyond the beam are dropped silently, so `to` is
--   only faithful for states with at most this many significant amplitudes.
defaultSparseTerms :: Int
defaultSparseTerms = 100

instance Convertible MPS SparseMat where
  to   mps = toSparseMat (tol $ cfg mps) defaultSparseTerms mps
  from = \sm ->
    let (SparseMat ((m,_), nonzeros)) = sm
        n = integerlog2 m
        kets = [ a .* ket (toBits' n k) | ((k,_),a) <- nonzeros ] :: [MPS]
    in case kets of
         [] -> error "fromSparseMat: empty"
         (x:xs) -> foldl' (.+) x xs


instance Convertible MPS CMat where
  to   = mpsToDenseVec
  from psimat =
    let psi_sparse = MS.sparseMat psimat
    in from psi_sparse :: MPS

-- | Contract an MPS into a dense (2^n × 1) column vector, indexed by logical bits MSB-first
--   (matching QPP.Semantics.Matrix's ket convention). Faithfully reflects the MPS's stored
--   amplitudes — does not modify the MPS or add any truncation beyond what is already in it.
--   O(2^n · χ²); intended for small n.
mpsToDenseVec :: MPS -> CMat
mpsToDenseVec mps =
  let n = nSites mps
      step !t p =
        let Site x0 x1 = sites mps ! p
            t0 = t .*. H.tr' x0
            t1 = t .*. H.tr' x1
        in t0 H.=== t1
      vPhys = H.scale (scalar mps) $
                foldl' step ((1><1) [1:+0]) [n-1, n-2 .. 0]
      l2p   = log2phys mps
      logToPhys !iLog = foldl' (.|.) 0
        [ if testBit iLog (n-1-q) then 1 `shiftL` (n-1-(l2p ! q)) else 0
        | q <- [0..n-1] ]
      dim   = pow2 n
  in H.asColumn $ H.fromList
       [ vPhys `atIndex` (logToPhys iLog, 0) | iLog <- [0 .. dim-1] ]

-- measurement
frob2 :: CMat -> Double
frob2 x =
  let v = H.flatten x
  in realPart $ dot v v -- Hmatrix conjugates left argument to dot

-- | Move the center to physical site p and project it there.
projectCtrl :: Int -> Bool -> MPS -> (MPS, Profile)
projectCtrl p one m =
  let (m', prof) = moveCenterToPhys p m
  in (projectCenter p one m', prof)

-- | Project the *center* site to |b> (b=False => |0>, True => |1>).
--   Produces an *unnormalized* post-measurement state. Zeroing one physical
--   branch of the center site leaves all other sites' isometries untouched,
--   so the center stays valid.
projectCenter :: Int -> Bool -> MPS -> MPS
projectCenter p b st =
  let Site x0 x1 = sites st ! p
      (y0,y1)    = if b then (zeros_like x0, x1) else (x0, zeros_like x1)
  in st { sites = sites st // [(p, Site y0 y1)] }

zeros_like :: CMat -> CMat
zeros_like x = H.konst (0:+0) (rows x, cols x)

measure1 :: HasCallStack => (MPS, Outcomes, RNG) -> Int -> (MPS, Outcomes, RNG)
measure1 acc k = fst (measure1Profiled acc k)

-- | `measure1` with the SVD work of its center move reported (under
--   `Truncate` those SVDs can truncate, so they belong in the error bound).
measure1Profiled :: HasCallStack
                 => (MPS, Outcomes, RNG) -> Int -> ((MPS, Outcomes, RNG), Profile)
measure1Profiled (st, outs, u:us) k =
  let
      p   = log2phys st ! k
      (st1, prof) = moveCenterToPhys p st
      Site x0 x1 = sites st1 ! p
      s2  = let a = magnitude (scalar st1) in a*a
      p0  = s2 * frob2 x0
      p1  = s2 * frob2 x1
      tot = p0 + p1
  in if tot < (tol $ cfg st) then error ("measure1: prob~0 measuring qubit " ++ show k) else
     let b    = (u*tot >= p0)
         pb   = if b then p1 else p0
         inv  = (1 / sqrt pb) :+ 0
         y0   = if b then zeros_like x0 else inv .* x0
         y1   = if b then inv .* x1 else zeros_like x1
         st2  = st1 { sites = sites st1 // [(p, Site y0 y1)] }
     in ((st2, b:outs, us), prof)
measure1Profiled (_,_,[]) _ = error "measure1: empty RNG"

-- | Sample a computational-basis outcome for *every* qubit without producing a
--   post-measurement state. Uses the Ferris--Vidal perfect-sampling recurrence:
--   given a mixed-canonical MPS at center @c@, sweep right from @c@ to @n-1@
--   maintaining a running matrix product @B@ as the left boundary, then sweep
--   left from @c-1@ to @0@ where the boundary has collapsed to a vector. No
--   SVDs; matmul for the right sweep, matvec for the left sweep.
--   Outcomes are returned head-most-recent, matching @Measure [0..n-1]@:
--   @outs = [b_0, b_1, ..., b_{n-1}]@ with @b_0@ at the head (last bit
--   processed by the @reverse ks@ fold).
sampleAll :: HasCallStack => MPS -> RNG -> (Outcomes, RNG)
sampleAll psi rng0
  | n == 0    = ([], rng0)
  | otherwise =
      let chiC  = rows (a0 (sites psi ! c))
          bInit = H.complex (H.ident chiC :: H.Matrix Double) :: CMat
          (bEnd, rng1, rightPairs) =
            foldl' rightStep (bInit, rng0, []) [c .. n-1]
          vInit = H.flatten bEnd
          (_v,   rng2, leftPairs)  =
            foldl' leftStep  (vInit, rng1, []) [c-1, c-2 .. 0]
          bits = V.replicate n False // (rightPairs ++ leftPairs)
          outs = [ bits ! k | k <- [0 .. n-1] ]
      in (outs, rng2)
  where
    n     = nSites psi
    c     = center_site psi
    tolP  = tol (cfg psi)

    rightStep :: (CMat, RNG, [(Int,Bool)]) -> Int
              -> (CMat, RNG, [(Int,Bool)])
    rightStep (_, [], _)     _ = error "sampleAll: empty RNG"
    rightStep (b, r:rs, acc) p =
      let Site x0 x1 = sites psi ! p
          m0  = b .*. x0
          m1  = b .*. x1
          w0  = frob2 m0
          w1  = frob2 m1
          tot = w0 + w1
      in if tot < tolP
           then error "sampleAll: prob ~ 0"
           else let bit = r * tot >= w0
                    b'  = if bit then m1 else m0
                    q   = phys2log psi ! p
                in (b', rs, (q, bit) : acc)

    leftStep :: (CVec, RNG, [(Int,Bool)]) -> Int
             -> (CVec, RNG, [(Int,Bool)])
    leftStep (_, [], _)     _ = error "sampleAll: empty RNG"
    leftStep (v, r:rs, acc) p =
      let Site x0 x1 = sites psi ! p
          u0  = x0 #> v
          u1  = x1 #> v
          w0  = realPart (dot u0 u0)   -- dot conjugates left arg => ‖u‖²
          w1  = realPart (dot u1 u1)
          tot = w0 + w1
      in if tot < tolP
           then error "sampleAll: prob ~ 0"
           else let bit = r * tot >= w0
                    v'  = if bit then u1 else u0
                    q   = phys2log psi ! p
                in (v', rs, (q, bit) : acc)

-- Helper functions for MPS
bondDimensions :: MPS -> V.Vector Int
bondDimensions mps = bondDim <$> (sites mps)
  where
    bondDim (Site a0 _) = rows a0


maxBondDimension :: MPS -> Int
maxBondDimension mps = maximum . V.toList . bondDimensions $ mps
