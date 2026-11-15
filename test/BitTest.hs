-- | Tests of `Bit`, the classical bit of "QPP.Bit": its arithmetic modulo 2,
--   how it prints, its conversions, and the literals that are not bits.
module Main where

import QPP
import qualified QPP.Semantics.Matrix as MS

import Control.Exception (ErrorCall(..), evaluate, try)
import Numeric.LinearAlgebra (sumElements)

import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck

import Generators ()     -- the Arbitrary instance of Bit

bits :: [Bit]
bits = [0, 1]

-- | The message of the error that evaluating x raises, if it raises one.
errorOf :: a -> IO (Maybe String)
errorOf x = either (\(ErrorCall m) -> Just m) (const Nothing) <$> try (evaluate x)

main :: IO ()
main = defaultMain $ testGroup "Bit"
  [ testCase "+ is exclusive or, * is and, and negation changes nothing" $ do
      [ a + b | a <- bits, b <- bits ] @?= [0, 1, 1, 0]
      [ a * b | a <- bits, b <- bits ] @?= [0, 0, 0, 1]
      [ a - b | a <- bits, b <- bits ] @?= [0, 1, 1, 0]
      map negate bits @?= bits
  , testCase "bits print as 0 and 1" $ do
      show bits @?= "[0,1]"
      show (Initialize [0, 2] [0, 1]) @?= "Initialize [0,2] [0,1]"
  , testCase "conversions from Bool and Int, and to Int" $ do
      map Bit [False, True] @?= bits
      map fromIntegral [0, 1 :: Int] @?= bits
      map fromEnum bits @?= [0, 1]
      [minBound .. maxBound] @?= bits
  , testCase "a number other than 0 or 1 is not a bit" $ do
      errorOf (2 :: Bit) >>= (@?= Just "2 is not a bit")
      errorOf (fromIntegral (7 :: Int) :: Bit) >>= (@?= Just "7 is not a bit")
      errorOf (sumElements (MS.ket [0, 2])) >>= (@?= Just "2 is not a bit")
  , testProperty "the sum of bits is the parity of the number of ones" $ \bs ->
      sum bs === Bit (odd (length (filter (== 1) bs)))
  ]
