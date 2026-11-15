{-# LANGUAGE ScopedTypeVariables #-}
-- | The comparison of the massiv state-vector backend (Old) with its
--   replacement (New) that was run before the old one was retired. The
--   results are in attic/README.md. Not built.
--
--   To run it again, put StatevectorMassiv.hs back into the library as
--   @QPP.Semantics.StatevectorMassiv@ (attic/README.md says how) and add
--
--   > executable StatevectorBench
--   >   default-language: GHC2021
--   >   hs-source-dirs: attic, test
--   >   main-is: StatevectorBench.hs
--   >   other-modules: Generators
--   >   build-depends: base, hmatrix, qpp, time, QuickCheck, tasty-hunit, tasty-quickcheck
--   >   ghc-options: -O1 -threaded -rtsopts
--
--   to qpp.cabal. Then, for example,
--
--   > cabal run StatevectorBench -- sweep 60 +RTS -N1
--   > cabal run StatevectorBench -- programs 60 +RTS -N1
--   > cabal run StatevectorBench -- small +RTS -N1
--   > cabal run StatevectorBench -- one old qft 20 +RTS -s
module Main where

import QPP
import QPP.Semantics (Convertible(..))
import QPP.Semantics.Matrix (CMat)
import qualified QPP.Semantics.Matrix            as M
import qualified QPP.Semantics.Statevector       as New
import qualified QPP.Semantics.StatevectorMassiv as Old
import Programs.QFT (qft)
import Generators (genQOpAt)

import qualified Numeric.LinearAlgebra as H
import Test.QuickCheck.Gen (Gen, unGen, choose)
import Test.QuickCheck.Random (mkQCGen)
import Data.Time.Clock
import System.Timeout
import System.IO
import Control.Exception (evaluate)
import Control.Monad
import System.Environment
import Text.Printf
import Data.Complex
import Data.IORef

------------------------------------------------------------------------
-- helpers
------------------------------------------------------------------------

-- | Hide a value from the optimiser, so that repeated runs are not shared.
opaque :: a -> IO a
opaque x = newIORef x >>= readIORef
{-# NOINLINE opaque #-}

diffM :: CMat -> CMat -> Double
diffM a b = H.sumElements (H.cmap magnitude (a - b))

-- | A deterministic stream of samples in [0,1).
rngFrom :: Int -> [Double]
rngFrom seed = map (\x -> fromIntegral (x `mod` 1000003) / 1000003) (tail (iterate step seed))
  where step x = (x * 1103515245 + 12345) `mod` 2147483648

sample :: Int -> Gen a -> a
sample seed g = unGen g (mkQCGen seed) 30

zeros :: Int -> [Int]
zeros n = replicate n 0

-- | Wall-clock seconds of one run, or Nothing if it exceeds the limit.
timed :: Int -> IO a -> IO (Maybe (a, Double))
timed limitSec act = do
  t0 <- getCurrentTime
  r  <- timeout (limitSec * 1000000) act
  t1 <- getCurrentTime
  pure (fmap (\x -> (x, realToFrac (diffUTCTime t1 t0))) r)

-- | Best of three runs when the first one is short; a single run otherwise.
best :: Int -> IO a -> IO (Maybe (a, Double))
best limitSec act = do
  r1 <- timed limitSec act
  case r1 of
    Just (x, t) | t < 2 -> do
      more <- replicateM 2 (timed limitSec act)
      pure (Just (x, minimum (t : [ t' | Just (_, t') <- more ])))
    _ -> pure r1

data Backend = BOld | BNew | BMat deriving (Eq, Show)

-- | Apply an operator to |0...0> and force every amplitude. The result is
--   returned as a column for comparison.
run :: Backend -> QOp -> Int -> IO CMat
run b op0 n0 = do
  op <- opaque op0
  n  <- opaque n0
  case b of
    BNew -> do let st = New.apply (New.evalOp op) (New.ket (zeros n))
               _ <- evaluate (H.norm_2 (H.flatten st))
               pure st
    BMat -> do let st = M.apply (M.evalOp op) (M.ket (zeros n))
               _ <- evaluate (H.norm_2 (H.flatten st))
               pure st
    BOld -> do let st = Old.apply (Old.evalOp op) (Old.ket (zeros n))
               _ <- evaluate (norm st)
               -- the conversion is outside the measured work only for small n
               pure (if n <= 14 then to st else H.fromLists [[0]])

fmt :: Maybe Double -> String
fmt Nothing  = "    >limit"
fmt (Just t) = printf "%10.4f" t

------------------------------------------------------------------------
-- circuit families
------------------------------------------------------------------------
at :: Int -> Int -> QOp -> QOp           -- place a 1-qubit op at qubit k of n
at n k g = k <@ g @> (n-k-1)

hTensor, hComposed, rxLayer, zzLayer, ghz, bitrev, toffolis :: Int -> QOp
hTensor n   = foldr1 (⊗) (replicate n H)
hComposed n = foldr1 (∘) [ at n k H | k <- [0..n-1] ]
rxLayer n   = foldr1 (∘) [ at n k (R X (1/3)) | k <- [0..n-1] ]
zzLayer n   = foldr1 (∘) [ k <@ R (Z ⊗ Z) (1/5) @> (n-k-2) | k <- [0..n-2] ]
ghz n       = foldr1 (∘) (reverse (at n 0 H : [ k <@ C X @> (n-k-2) | k <- [0..n-2] ]))
bitrev n    = Permute [n-1, n-2 .. 0]
-- an X with n-1 controls, after one H per qubit so that every amplitude is non-zero
toffolis n  = iterate C X !! (n-1) ∘ hTensor n

rzChain :: Int -> Int -> QOp             -- m rotations in a row on qubit 0 of n
rzChain n m = foldr1 (∘) (replicate m (R Z (1/8) @> (n-1)))

trotter :: Int -> Int -> QOp             -- m Trotter steps of a transverse-field Ising chain
trotter n m = foldr1 (∘) (replicate m (rxLayer n ∘ zzLayer n))

-- | key, title, parameters, and the (qubits, operator) of a parameter
families :: [(String, String, [Int], Int -> (Int, QOp))]
families =
  [ ("qft",       "QFT on n qubits",                             [4,6..22], \n -> (n, qft n))
  , ("ghz",       "GHZ preparation on n qubits",                 [4,6..22], \n -> (n, ghz n))
  , ("trotter",   "Ising Trotter, 3 steps, n qubits",            [4,6..20], \n -> (n, trotter n 3))
  , ("htensor",   "H (x) H (x) ... (x) H on n qubits",           [2,4..22], \n -> (n, hTensor n))
  , ("hcomposed", "one H per qubit, composed, n qubits",         [2,4..22], \n -> (n, hComposed n))
  , ("rxlayer",   "one R X per qubit, composed, n qubits",       [2,4..22], \n -> (n, rxLayer n))
  , ("rzchain",   "m successive R Z on qubit 0 of 10 qubits",    [2,4..30], \m -> (10, rzChain 10 m))
  , ("bitrev",    "bit-reversal Permute on n qubits",            [4,8..24], \n -> (n, bitrev n))
  , ("toffolis",  "X with n-1 controls, after H on every qubit", [4,8..20], \n -> (n, toffolis n))
  ]

-- | Sweep one family; a backend that exceeds the limit is not run on larger parameters.
sweep :: Int -> (String, String, [Int], Int -> (Int, QOp)) -> IO ()
sweep limit (key, title, params, build) = do
  printf "\n== %s: %s  (seconds, limit %d s) ==\n" key title limit
  printf "%6s %6s %10s %10s %10s %10s %10s\n"
         ("param" :: String) ("qubits" :: String) ("old" :: String) ("new" :: String)
         ("old/new" :: String) ("matrix" :: String) ("|old-new|" :: String)
  oldDead <- newIORef False
  newDead <- newIORef False
  matDead <- newIORef False
  forM_ params $ \p -> do
    let (n, op) = build p
    nd <- readIORef newDead
    rn <- if nd then pure Nothing else best limit (run BNew op n)
    when (maybe True (const False) rn) (writeIORef newDead True)
    od <- readIORef oldDead
    ro <- if od then pure Nothing else best limit (run BOld op n)
    when (maybe True (const False) ro) (writeIORef oldDead True)
    md <- readIORef matDead
    rm <- if md || n > 10 then pure Nothing else best limit (run BMat op n)
    when (n <= 10 && maybe True (const False) rm) (writeIORef matDead True)
    let ratio = case (ro, rn) of
                  (Just (_, to'), Just (_, tn)) -> printf "%10.1f" (to' / tn)
                  _                             -> "          " :: String
        d = case (ro, rn) of
              (Just (so, _), Just (sn, _)) | n <= 14 -> printf "%10.1e" (diffM so sn)
              _                                      -> "" :: String
    printf "%6d %6d %s %s %s %s %s\n" p n (fmt (snd <$> ro)) (fmt (snd <$> rn)) ratio
           (if n <= 10 then fmt (snd <$> rm) else "          ") d
    hFlush stdout

------------------------------------------------------------------------
-- programs with measurement
------------------------------------------------------------------------
progs :: Int -> IO ()
progs limit = do
  printf "\n== programs: GHZ, then Measure every qubit, then Initialize every qubit to 0  (seconds, limit %d s) ==\n" limit
  printf "%6s %10s %10s %10s %s\n" ("qubits" :: String) ("old" :: String) ("new" :: String) ("old/new" :: String) ("  same outcomes" :: String)
  oldDead <- newIORef False
  forM_ [4,8,12,16,20] $ \n0 -> do
    let prog0 = [Unitary (ghz n0), Measure [0 .. n0-1], Initialize [0 .. n0-1] (replicate n0 False)]
        runNew = do n <- opaque n0; prog <- opaque prog0
                    let (st, outs, _) = New.evalProg prog (New.ket (zeros n)) (rngFrom 7)
                    _ <- evaluate (H.norm_2 (H.flatten st)); pure outs
        runOld = do n <- opaque n0; prog <- opaque prog0
                    let (st, outs, _) = Old.evalProg prog (Old.ket (zeros n)) (rngFrom 7)
                    _ <- evaluate (norm st); pure outs
    rn <- best limit runNew
    od <- readIORef oldDead
    ro <- if od then pure Nothing else best limit runOld
    when (maybe True (const False) ro) (writeIORef oldDead True)
    let ratio = case (ro, rn) of
                  (Just (_, to'), Just (_, tn)) -> printf "%10.1f" (to' / tn)
                  _                             -> "          " :: String
        same  = case (ro, rn) of
                  (Just (oo, _), Just (on, _)) -> if oo == on then "  yes" else "  NO" :: String
                  _                            -> ""
    printf "%6d %s %s %s %s\n" n0 (fmt (snd <$> ro)) (fmt (snd <$> rn)) ratio same
    hFlush stdout

------------------------------------------------------------------------
-- many small operators, as in the tests and the exercises
------------------------------------------------------------------------
small :: IO ()
small = do
  printf "\n== small random operators (depth <= 4), each applied to |0...0>; total seconds for all of them ==\n"
  printf "%6s %6s %10s %10s %10s %10s %10s\n" ("qubits" :: String) ("count" :: String) ("old" :: String)
         ("new" :: String) ("old/new" :: String) ("matrix" :: String) ("max diff" :: String)
  forM_ [(1,2000),(2,2000),(3,2000),(4,1000),(5,1000),(6,500),(8,200)] $ \(n, count :: Int) -> do
    let ops = [ sample (1000 * n + i) (do d <- choose (1, 4); genQOpAt n d) | i <- [1 .. count] ]
    _ <- evaluate (sum (map (length . showOp) ops))
    let total b = do
          r <- best 300 (forM ops (\op -> run b op n))
          pure (fmap (\(sts, t) -> (sts, t)) r)
    rn <- total BNew
    ro <- total BOld
    rm <- total BMat
    let worst = case (ro, rn) of
                  (Just (so, _), Just (sn, _)) -> printf "%10.1e" (maximum (zipWith diffM so sn))
                  _                            -> "" :: String
        ratio = case (ro, rn) of
                  (Just (_, to'), Just (_, tn)) -> printf "%10.1f" (to' / tn)
                  _                             -> "          " :: String
    printf "%6d %6d %s %s %s %s %s\n" n count (fmt (snd <$> ro)) (fmt (snd <$> rn)) ratio (fmt (snd <$> rm)) worst
    hFlush stdout
  printf "\n== the same operators converted to their matrices with `to (evalOp op)`, as the test suites do ==\n"
  printf "%6s %6s %10s %10s %10s %10s\n" ("qubits" :: String) ("count" :: String) ("old" :: String)
         ("new" :: String) ("old/new" :: String) ("max diff" :: String)
  forM_ [(2,1000),(4,500),(6,100)] $ \(n, count :: Int) -> do
    let ops = [ sample (2000 * n + i) (do d <- choose (1, 4); genQOpAt n d) | i <- [1 .. count] ]
    _ <- evaluate (sum (map (length . showOp) ops))
    let conv :: (QOp -> CMat) -> IO [CMat]
        conv f = forM ops $ \op0 -> do op <- opaque op0
                                       let m = f op
                                       _ <- evaluate (H.norm_2 (H.flatten m)); pure m
    rn <- best 300 (conv (\op -> to (New.evalOp op)))
    ro <- best 300 (conv (\op -> to (Old.evalOp op)))
    let worst = case (ro, rn) of
                  (Just (so, _), Just (sn, _)) -> printf "%10.1e" (maximum (zipWith diffM so sn))
                  _                            -> "" :: String
        ratio = case (ro, rn) of
                  (Just (_, to'), Just (_, tn)) -> printf "%10.1f" (to' / tn)
                  _                             -> "          " :: String
    printf "%6d %6d %s %s %s %s\n" n count (fmt (snd <$> ro)) (fmt (snd <$> rn)) ratio worst
    hFlush stdout

------------------------------------------------------------------------
main :: IO ()
main = do
  hSetBuffering stdout LineBuffering
  args <- getArgs
  case args of
    ["sweep", limit]      -> mapM_ (sweep (read limit)) families
    ["sweep", limit, key] -> mapM_ (sweep (read limit)) [ f | f@(k, _, _, _) <- families, k == key ]
    ["programs", limit]   -> progs (read limit)
    ["small"]             -> small
    -- one run, for `+RTS -s`
    ["one", backend, key, p] -> do
      let b = case backend of { "old" -> BOld; "new" -> BNew; _ -> BMat }
          (n, op) = head [ build (read p) | (k, _, _, build) <- families, k == key ]
      r <- timed 3600 (run b op n)
      printf "%s %s %s: %s s\n" backend key p (fmt (snd <$> r))
    _ -> putStrLn "usage: StatevectorBench (sweep LIMIT [FAMILY] | programs LIMIT | small | one (old|new|matrix) FAMILY PARAM)"
