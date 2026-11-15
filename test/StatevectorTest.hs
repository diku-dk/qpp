-- | Properties of the state-vector semantics, "QPP.Semantics.Statevector".
--
--   The semantics has one rule per constructor of `QOp`: @applyOp e@ acts on
--   the first qubits of a state and carries any further qubits along. The
--   suite checks it in six ways.
--
--   * __Matrix semantics, 0 to 8 qubits.__ The matrix of @evalOp e@ is the
--     matrix of "QPP.Semantics.Matrix", and on a state with @t@ further qubits
--     @applyOp e ψ = (M(e) ⊗ I) ψ@. Random expressions of every arity, and
--     every constructor in fixed positions at every size.
--   * __The rules.__ Each rule as it is written, on the blocks and halves of a
--     state. No matrix of a whole operator is formed here.
--   * __Algebra.__ Linearity, norm, adjoints, shapes, several states at once.
--   * __The operator interface.__ Evaluated operators combined with `<>`, `⊗`,
--     `⊕`, `adj` and the vector-space operations have the matrices that the
--     matrix backend gets from the same combinations, and the instances agree
--     with the rules of `applyOp`.
--   * __Programs and measurement.__ The same outcomes and final state as the
--     matrix backend on a shared supply of samples; a measurement against the
--     amplitudes themselves.
--   * __Worked examples__ of the lecture notes' section on this semantics.
--
--   Rotation axes are any Hermitian expression, not only the Pauli strings of
--   `genQOpAt`. Budgets are per group, for the default of 100 QuickCheck tests,
--   and scale with the command line:
--   @cabal run StatevectorTest -- --quickcheck-tests 1000@ runs ten times as
--   many cases.
module Main where

import QPP
import QPP.Semantics (Convertible(..))
import QPP.Semantics.Matrix (CMat)
import qualified QPP.Semantics.Matrix      as MS
import qualified QPP.Semantics.Statevector as SV
import Programs.QFT (qft)

import Numeric.LinearAlgebra
  ( rows, flatten, reshape, toList, toRows, toColumns, fromColumns
  , asColumn, vjoin, kronecker, ident, scale, takeRows, dropRows, (><), norm_2 )
import qualified Numeric.LinearAlgebra as H

import Control.Applicative (liftA3)
import Control.Exception (ErrorCall, evaluate, try)
import Data.Bits (testBit)
import Data.Complex
import Data.List (isInfixOf)

import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck hiding (scale, (><))     -- hmatrix's here

import Generators

------------------------------------------------------------------------
-- Amplitudes
------------------------------------------------------------------------

-- | The state-vector semantics, on any column of amplitudes.
sv :: QOp -> CMat -> CMat
sv = SV.applyOp

-- | r × c entries from the unit square. They are not normalised: the rules
--   are linear, and recursive calls act on parts of a state anyway.
genAmps :: Int -> Int -> Gen CMat
genAmps r c = (r >< c) <$> vectorOf (r * c) (liftA2 (:+) unit unit)
  where unit = choose (-1, 1)

-- | A column of 2^q amplitudes.
genColumn :: Int -> Gen CMat
genColumn q = genAmps (2 ^ q) 1

-- | A column for an expression of a qubits: a, a+1 or a+2 qubits, at most 8.
genFor :: Int -> Gen CMat
genFor a = do
  t <- choose (0, min 2 (8 - a))
  genColumn (a + t)

-- | A state of norm 1 on q qubits.
genState :: Int -> Gen CMat
genState q = do
  ψ <- genColumn q `suchThat` ((> 0.1) . norm_2 . flatten)
  pure (scale (1 / (norm_2 (flatten ψ) :+ 0)) ψ)

column :: [ComplexT] -> CMat
column zs = (length zs >< 1) zs

-- | The 2^n blocks of a column. Block x holds, in order, the amplitudes whose
--   first n qubits read x.
blocks :: Int -> CMat -> [CMat]
blocks n ψ = map asColumn (toRows (reshape (rows ψ `div` 2 ^ n) (flatten ψ)))

-- | Blocks one below the other, as one column.
stack :: [CMat] -> CMat
stack = asColumn . vjoin . map flatten

halves :: CMat -> (CMat, CMat)
halves ψ = (takeRows h ψ, dropRows h ψ)
  where h = rows ψ `div` 2

-- | ⟨a, b⟩, conjugate-linear in the first argument.
innerC :: CMat -> CMat -> ComplexT
innerC a b = sum (zipWith (\x y -> conjugate x * y) (toList (flatten a)) (toList (flatten b)))

-- | The same shape and, up to rounding, the same entries.
(=~) :: CMat -> CMat -> Property
a =~ b
  | H.size a /= H.size b =
      counterexample ("shapes " ++ show (H.size a) ++ " and " ++ show (H.size b)) False
  | otherwise = a ~~ b

infix 4 =~

-- | k tests at the default budget of 100; scales with @--quickcheck-tests@.
budget :: Int -> TestTree -> TestTree
budget k = adjustOption (\(QuickCheckTests n) -> QuickCheckTests (max 1 (n * k `div` 100)))

------------------------------------------------------------------------
-- Expressions
------------------------------------------------------------------------

-- | A well-formed expression on exactly n qubits, of depth at most d. Unlike
--   `genQOpAt`, it also rotates about Hermitian axes that are not Pauli strings.
genExpr :: Int -> Int -> Gen QOp
genExpr n d
  | d <= 0    = genQOpBase n
  | otherwise = frequency $
      [ (2, genQOpBase n)
      , (3, do k <- choose (0, n)
               liftA2 Tensor (genExpr k (d-1)) (genExpr (n-k) (d-1)))
      , (3, liftA2 Compose (genExpr n (d-1)) (genExpr n (d-1)))
      , (1, Adjoint <$> genExpr n (d-1))
      , (2, liftA2 R (genHermitian n (d-1)) genRat)
      ] ++
      [ (2, C <$> genExpr (n-1) (d-1))                                    | n >= 1 ] ++
      [ (2, liftA2 DirectSum (genExpr (n-1) (d-1)) (genExpr (n-1) (d-1))) | n >= 1 ]

-- | An expression of exactly n qubits and depth at most 4.
genExprOf :: Int -> Gen QOp
genExprOf n = choose (0, 4) >>= genExpr n

-- | An expression on exactly n qubits whose matrix is Hermitian, and unitary
--   like every well-formed expression: a rotation axis.
genHermitian :: Int -> Int -> Gen QOp
genHermitian 0 _ = elements [Id 0, Phase 0, Phase 1]       -- the scalars 1 and -1
genHermitian n d
  | d <= 0    = base
  | otherwise = frequency
      [ (3, base)
      , (2, do k <- choose (0, n)
               liftA2 Tensor (genHermitian k (d-1)) (genHermitian (n-k) (d-1)))
      , (2, C <$> genHermitian (n-1) (d-1))
      , (2, liftA2 DirectSum (genHermitian (n-1) (d-1)) (genHermitian (n-1) (d-1)))
        -- a conjugate u† h u of a Hermitian h
      , (1, do u <- genExpr n (d-1)
               h <- genHermitian n (d-1)
               pure (Compose (Adjoint u) (Compose h u)))
      ]
  where
    base = oneof $ [ genPauli n, pure (Id n), Permute <$> genInvolution n ]
                ++ [ elements [X, Y, Z, H] | n == 1 ]

-- | A permutation of [0..n-1] that is its own inverse: some disjoint swaps.
genInvolution :: Int -> Gen [Int]
genInvolution n = do
  xs <- shuffle [0 .. n-1]
  m  <- choose (0, n `div` 2)
  let swaps        = take m (pairs xs)
      partner i    = head ([ b | (a, b) <- swaps, a == i ] ++ [ a | (a, b) <- swaps, b == i ] ++ [i])
  pure (map partner [0 .. n-1])
  where pairs (a:b:rest) = (a, b) : pairs rest
        pairs _          = []

-- | Whether every rotation in the expression has a Hermitian axis or the angle 0.
axesHermitian :: QOp -> Bool
axesHermitian e = case e of
  R ax θ        -> (θ == 0 || hermitian ax) && axesHermitian ax
  C u           -> axesHermitian u
  Adjoint u     -> axesHermitian u
  Tensor u v    -> axesHermitian u && axesHermitian v
  DirectSum u v -> axesHermitian u && axesHermitian v
  Compose u v   -> axesHermitian u && axesHermitian v
  _             -> True
  where hermitian ax = let m = MS.evalOp ax in mDiff m (H.tr m) < tol

-- | `shrinkQOp` may cut a rotation axis down to a part that is not Hermitian;
--   such candidates are not well-formed expressions and are dropped.
shrinkExpr :: Int -> QOp -> [QOp]
shrinkExpr n = filter axesHermitian . shrinkQOp n

-- | A property of every random expression of exactly n qubits.
forExpr :: Testable prop => Int -> (QOp -> prop) -> Property
forExpr n = forAllShrinkShow (genExprOf n) (shrinkExpr n) showOp

-- | An expression of a qubits (0 to 5) together with a column of a + t qubits
--   (t from 0 to 3).
data Applied = Applied Int Int QOp CMat

instance Show Applied where
  show (Applied a t e ψ) =
    showOp e ++ "  (" ++ show a ++ " qubits) on " ++ show (a + t) ++ " qubits: "
             ++ show (toList (flatten ψ))

instance Arbitrary Applied where
  arbitrary = do
    a <- choose (0, 5)
    t <- choose (0, 3)
    liftA2 (Applied a t) (genExprOf a) (genColumn (a + t))
  shrink (Applied a t e ψ) =
       [ Applied a (t-1) e (takeRows (rows ψ `div` 2) ψ) | t > 0 ]
    ++ [ Applied a t e' ψ | e' <- shrinkExpr a e ]

------------------------------------------------------------------------
-- Against the matrix semantics, 0 to 8 qubits
------------------------------------------------------------------------

-- | The operator matrix @to (evalOp e)@ is M(e).
prop_matrix :: Int -> Property
prop_matrix n = forExpr n $ \e ->
  (to (SV.evalOp e) :: CMat) =~ MS.evalOp e

-- | On a column with t further qubits, S(e) ψ = (M(e) ⊗ I) ψ.
prop_carries :: Int -> Int -> Property
prop_carries a t = forExpr a $ \e ->
  forAll (genColumn (a + t)) $ \ψ ->
    sv e ψ =~ (kronecker (MS.evalOp e) (ident (2 ^ t)) H.<> ψ)

-- | Every constructor in fixed positions on exactly n qubits.
fixedCases :: Int -> [(String, QOp)]
fixedCases n =
     [ ("Id", Id n), ("Phase 1/3 ⊗ Id", Phase (1/3) ⊗ Id n), ("Id ⊗ Phase -1/4", Id n ⊗ Phase (-1/4)) ]
  ++ [ (name ++ " on qubit " ++ show k, k <@ g @> (n-k-1))
     | (name, g) <- [("X", X), ("Y", Y), ("Z", Z), ("H", H), ("SX", SX)], k <- [0 .. n-1] ]
  ++ [ ("R " ++ name ++ " 2/7 on qubit " ++ show k, k <@ R g (2/7) @> (n-k-1))
     | (name, g) <- [("X", X), ("Y", Y), ("Z", Z), ("H", H)], k <- [0 .. n-1] ]
  ++ [ ("C X on qubits " ++ show k ++ " and " ++ show (k+1), k <@ C X @> (n-k-2)) | k <- [0 .. n-2] ]
  ++ [ ("rotate the qubits by " ++ show s, Permute ([s .. n-1] ++ [0 .. s-1])) | s <- [1 .. n-1] ]
  ++ [ ("reverse the qubits", Permute [n-1, n-2 .. 0]) ]
  ++ [ ("qft ⊗ Adjoint qft, split after qubit " ++ show (k-1), qft k ⊗ Adjoint (qft (n-k))) | k <- [1 .. n-1] ]
  ++ concat
     [ [ ("X under " ++ show (n-1) ++ " controls", iterate C X !! (n-1))
       , ("H ⊗ ... ⊗ H", foldr1 (⊗) (replicate n H))
       , ("R (Z ⊗ ... ⊗ Z) 1/5", R (foldr1 (⊗) (replicate n Z)) (1/5))
       , ("R (reverse the qubits) 1/3", R (Permute [n-1, n-2 .. 0]) (1/3))
       , ("R (X under controls) -3/8", R (iterate C X !! (n-1)) (-3/8))
       , ("qft", qft n)
       , ("Adjoint qft", Adjoint (qft n))
       , ("qft ∘ Adjoint qft", qft n ∘ Adjoint (qft n))
       , ("C qft", C (qft (n-1)))
       , ("qft ⊕ Adjoint qft", DirectSum (qft (n-1)) (Adjoint (qft (n-1))))
       ]
     | n >= 2 ]

matrixTests :: TestTree
matrixTests = testGroup "Matrix semantics, 0 to 8 qubits"
  [ budget 100 $ testGroup "to (evalOp e) = M(e), random expressions"
      [ testProperty (show n ++ " qubits") (prop_matrix n) | n <- [0 .. 8] ]
  , budget 20 $ testGroup "S(e) ψ = (M(e) ⊗ I) ψ with t further qubits"
      [ testProperty (show a ++ " + " ++ show t ++ " qubits") (prop_carries a t)
      | a <- [0 .. 8], t <- [0 .. 8 - a] ]
  , testGroup "Every constructor in fixed positions"
      [ testCase (show n ++ " qubits: " ++ show (length (fixedCases n)) ++ " expressions") $
          mapM_ (\(name, e) -> assertClose (show n ++ " qubits, " ++ name)
                                           (MS.evalOp e) (to (SV.evalOp e)))
                (fixedCases n)
      | n <- [0 .. 8] ]
  ]

------------------------------------------------------------------------
-- The rules, as they are written
------------------------------------------------------------------------

prop_rule_id :: Property
prop_rule_id =
  forAll (choose (0, 6)) $ \n -> forAll (genFor n) $ \ψ ->
    sv (Id n) ψ =~ ψ

prop_rule_phase :: Property
prop_rule_phase =
  forAll genRat $ \θ -> forAll (genFor 0) $ \ψ ->
    sv (Phase θ) ψ =~ scale (cis (pi * fromRational θ)) ψ

-- | A gate combines the two halves with the entries of its 2×2 matrix.
prop_rule_gate :: Property
prop_rule_gate =
  forAll (elements [X, Y, Z, H, SX]) $ \g -> forAll (genFor 1) $ \ψ ->
    let m        = MS.evalOp g
        g_ i j   = scale (m `H.atIndex` (i, j))
        (ψ0, ψ1) = halves ψ
    in sv g ψ =~ stack [ g_ 0 0 ψ0 + g_ 0 1 ψ1, g_ 1 0 ψ0 + g_ 1 1 ψ1 ]

-- | Each block at b₀…bₙ₋₁ moves to b_{k₀}…b_{kₙ₋₁}.
prop_rule_permute :: Property
prop_rule_permute =
  forAll (choose (0, 6)) $ \n -> forAll (genPermutation n) $ \ks -> forAll (genFor n) $ \ψ ->
    let before   = blocks n ψ
        moved    = blocks n (sv (Permute ks) ψ)
        bits x   = [ testBit x (n - 1 - j) | j <- [0 .. n-1] ]
        index bs = sum [ 2 ^ (n - 1 - j) | (j, True) <- zip [0 :: Int ..] bs ]
        target x = index [ bits x !! k | k <- ks ]
    in conjoin [ (moved !! target x) =~ (before !! x) | x <- [0 .. 2 ^ n - 1 :: Int] ]

prop_rule_rotation :: Property
prop_rule_rotation =
  forAll gen $ \(ax, θ, ψ) ->
    let α = pi * fromRational θ / 2
    in sv (R ax θ) ψ =~ scale (cos α :+ 0) ψ - scale (0 :+ sin α) (sv ax ψ)
  where gen = do n  <- choose (0, 5)
                 ax <- genHermitian n =<< choose (0, 3)
                 θ  <- genRat
                 ψ  <- genFor n
                 pure (ax, θ, ψ)

prop_rule_control :: Property
prop_rule_control =
  forAll gen $ \(u, ψ) ->
    let (ψ0, ψ1) = halves ψ
    in sv (C u) ψ =~ stack [ψ0, sv u ψ1]
  where gen = do n <- choose (0, 5)
                 liftA2 (,) (genExprOf n) (genFor (n + 1))

-- | S(v) on each block for the qubits of u, then S(u) on the result.
prop_rule_tensor :: Property
prop_rule_tensor =
  forAll gen $ \(u, v, ψ) ->
    let χ = stack [ sv v b | b <- blocks (op_qubits u) ψ ]
    in sv (Tensor u v) ψ =~ sv u χ
  where gen = do n <- choose (0, 3)
                 m <- choose (0, 3)
                 liftA3 (,,) (genExprOf n) (genExprOf m) (genFor (n + m))

prop_rule_directsum :: Property
prop_rule_directsum =
  forAll gen $ \(u, v, ψ) ->
    let (ψ0, ψ1) = halves ψ
    in sv (DirectSum u v) ψ =~ stack [sv u ψ0, sv v ψ1]
  where gen = do n <- choose (0, 5)
                 liftA3 (,,) (genExprOf n) (genExprOf n) (genFor (n + 1))

prop_rule_compose :: Property
prop_rule_compose =
  forAll gen $ \(u, v, ψ) ->
    sv (Compose u v) ψ =~ sv u (sv v ψ)
  where gen = do n <- choose (0, 5)
                 liftA3 (,,) (genExprOf n) (genExprOf n) (genFor n)

prop_rule_adjoint :: Applied -> Property
prop_rule_adjoint (Applied _ _ u ψ) =
  sv (Adjoint u) ψ =~ sv (dagger u) ψ

ruleTests :: TestTree
ruleTests = budget 100 $ testGroup "The rules"
  [ testProperty "Id n: ψ"                                    prop_rule_id
  , testProperty "Phase θ: e^{iπθ} ψ"                          prop_rule_phase
  , testProperty "X, Y, Z, H, SX: G on the two halves"         prop_rule_gate
  , testProperty "Permute k: block b₀…bₙ₋₁ moves to b_{k₀}…"   prop_rule_permute
  , testProperty "R u θ: cos α ψ − i sin α S(u) ψ"             prop_rule_rotation
  , testProperty "C u: S(u) on the lower half"                 prop_rule_control
  , testProperty "u ⊗ v: S(v) on each block, then S(u)"        prop_rule_tensor
  , testProperty "u ⊕ v: S(u) and S(v) on the halves"          prop_rule_directsum
  , testProperty "u ∘ v: S(u) (S(v) ψ)"                        prop_rule_compose
  , testProperty "Adjoint u: S(dagger u) ψ"                    prop_rule_adjoint
  ]

------------------------------------------------------------------------
-- Algebra
------------------------------------------------------------------------

prop_linear :: Applied -> Property
prop_linear (Applied a t e ψ) =
  forAll (genColumn (a + t)) $ \φ -> forAll (liftA2 (:+) unit unit) $ \z ->
    sv e (scale z ψ + φ) =~ scale z (sv e ψ) + sv e φ
  where unit = choose (-1, 1)

prop_norm :: Applied -> Property
prop_norm (Applied _ _ e ψ) =
  let normIn  = norm_2 (flatten ψ)
      normOut = norm_2 (flatten (sv e ψ))
  in counterexample (show normIn ++ " became " ++ show normOut) (abs (normOut - normIn) < tol)

prop_inverse :: Applied -> Property
prop_inverse (Applied _ _ e ψ) =
  sv (Adjoint e) (sv e ψ) =~ ψ .&&. sv e (sv (Adjoint e) ψ) =~ ψ

-- | ⟨φ, S(e) ψ⟩ = ⟨S(Adjoint e) φ, ψ⟩: the adjoint, with no matrix involved.
prop_adjoint_inner :: Applied -> Property
prop_adjoint_inner (Applied a t e ψ) =
  forAll (genColumn (a + t)) $ \φ ->
    let lhs = innerC φ (sv e ψ)
        rhs = innerC (sv (Adjoint e) φ) ψ
    in counterexample (show lhs ++ " and " ++ show rhs) (magnitude (lhs - rhs) < tol)

-- | `applyOp` returns the shape it was given, and the amplitudes do not depend
--   on that shape: any number of columns that divides the further amplitudes.
prop_shape :: Applied -> Property
prop_shape (Applied _ t e ψ) =
  forAll (elements [ 2 ^ s | s <- [0 .. t] ]) $ \(c :: Int) ->
    let φ      = reshape c (flatten ψ)
        result = sv e φ
    in H.size result === H.size φ .&&. asColumn (flatten result) =~ sv e ψ

-- | Columns side by side are transformed one by one, also when their number
--   is not a power of two.
prop_columns :: Property
prop_columns =
  forAll gen $ \(e, states) ->
    let u = SV.evalOp e
    in SV.apply u states
         =~ fromColumns [ flatten (SV.apply u (asColumn s)) | s <- toColumns states ]
  where gen = do n <- choose (0, 5)
                 p <- choose (1, 5)
                 liftA2 (,) (genExprOf n) (genAmps (2 ^ n) p)

prop_to_from :: Property
prop_to_from =
  forAll (choose (0, 5 :: Int)) $ \n -> forAll (genAmps (2 ^ n) (2 ^ n)) $ \m ->
    (to (from m :: SV.OpT) :: CMat) =~ m

prop_arity :: Applied -> Property
prop_arity (Applied a _ e _) = n_qubits (SV.evalOp e) === a

algebraTests :: TestTree
algebraTests = budget 100 $ testGroup "Algebra"
  [ testProperty "S(e) is linear"                          prop_linear
  , testProperty "S(e) keeps the norm"                     prop_norm
  , testProperty "S(Adjoint e) undoes S(e), and back"      prop_inverse
  , testProperty "⟨φ, S(e) ψ⟩ = ⟨S(Adjoint e) φ, ψ⟩"       prop_adjoint_inner
  , testProperty "the result has the argument's shape"     prop_shape
  , testProperty "several states side by side"             prop_columns
  , testProperty "to (from m) = m"                         prop_to_from
  , testProperty "n_qubits (evalOp e) is the arity of e"   prop_arity
  , testCase "apply rejects a state of another size" $ do
      r <- try (evaluate (SV.apply (SV.evalOp (C X)) (SV.ket [0, 0, 0])))
      case r of
        Left (err :: ErrorCall) ->
          assertBool (show err) ("2-qubit operator cannot act on 8 amplitudes" `isInfixOf` show err)
        Right _ -> assertFailure "a 2-qubit operator was applied to 3 qubits"
  ]

------------------------------------------------------------------------
-- The operator interface of evaluated operators
------------------------------------------------------------------------

-- | The interface does on evaluated operators what it does on their matrices.
--   Every operator is `from` of its matrix, so random matrices cover them
--   all: a and b on n qubits, c on m. The last line adds one further qubit.
prop_interface :: Property
prop_interface =
  forAll gen $ \(a, b, c, z, ψ) ->
    let u   = from a :: SV.OpT
        v   = from b
        w   = from c
        mat = to :: SV.OpT -> CMat
    in conjoin
         [ mat (u ⊗ w) =~ a ⊗ c
         , mat (u ⊕ v) =~ a ⊕ b
         , mat (u ∘ v) =~ a ∘ b
         , mat (adj u) =~ adj a
         , mat (z .* u .+ v .- u) =~ z .* a .+ b .- a
         , counterexample "inner" (magnitude (inner u v - inner a b) < tol)
         , counterexample "norm"  (abs (norm u - norm a) < tol)
         , SV.runOp (u ⊗ w) ψ =~ (kronecker (a ⊗ c) (ident 2) H.<> ψ)
         ]
  where gen = do n <- choose (0, 3)
                 m <- choose (0, 2)
                 a <- genAmps (2 ^ n) (2 ^ n)
                 b <- genAmps (2 ^ n) (2 ^ n)
                 c <- genAmps (2 ^ m) (2 ^ m)
                 z <- liftA2 (:+) (choose (-1, 1)) (choose (-1, 1))
                 ψ <- genColumn (n + m + 1)
                 pure (a, b, c, z, ψ)

-- | The instances follow the rules of `applyOp`: combining two evaluated
--   expressions gives the evaluated combination.
prop_interface_rules :: Property
prop_interface_rules =
  forAll gen $ \(u, v, w) ->
    conjoin [ mat (SV.evalOp u ⊗ SV.evalOp w) =~ mat (SV.evalOp (u ⊗ w))
            , mat (SV.evalOp u ⊕ SV.evalOp v) =~ mat (SV.evalOp (u ⊕ v))
            , mat (SV.evalOp u ∘ SV.evalOp v) =~ mat (SV.evalOp (u ∘ v))
            , mat (adj (SV.evalOp u))         =~ mat (SV.evalOp (adj u)) ]
  where gen = do n <- choose (0, 3)
                 m <- choose (0, 3)
                 liftA3 (,,) (genExprOf n) (genExprOf n) (genExprOf m)
        mat = to :: SV.OpT -> CMat

interfaceTests :: TestTree
interfaceTests = budget 100 $ testGroup "The operator interface"
  [ testProperty "⊗, ⊕, ∘, adj, .*, .+, .-, inner, norm as on the matrices"  prop_interface
  , testProperty "evalOp u ⊗ evalOp v = evalOp (u ⊗ v); so for ⊕, ∘ and adj" prop_interface_rules
  , testCase "the line of chapter 2: evalOp (C X) ∘ (evalOp H ⊗ evalOp I)" $
      assertClose "matrix" (to (SV.evalOp (C X) ∘ (SV.evalOp H ⊗ SV.evalOp I)))
                           (MS.evalOp (C X ∘ (H ⊗ I)))
  , testCase "operators on different numbers of qubits are not composed" $ do
      r <- try (evaluate (n_qubits (SV.evalOp X <> SV.evalOp (C X))))
      case r of
        Left (err :: ErrorCall) ->
          assertBool (show err) ("cannot be composed" `isInfixOf` show err)
        Right n -> assertFailure ("composed, on " ++ show n ++ " qubits")
  ]

------------------------------------------------------------------------
-- Programs and measurement
------------------------------------------------------------------------

-- | A deterministic supply of samples in (0,1). The offset 1/π keeps them away
--   from 0 and from the probabilities that occur exactly, such as 1/2 and 1/4:
--   a sample equal to a probability would make the outcome depend on rounding.
samples :: Int -> RNG
samples seed = map (\x -> (fromIntegral (x `mod` 999983) + 1 / pi) / 999983) (drop 1 (iterate next seed))
  where next x = (x * 1103515245 + 12345) `mod` 2147483648

-- | Up to six steps on q ≥ 1 qubits: unitaries, measurements and initialisations.
genProgram :: Int -> Gen Program
genProgram q = do
  len <- choose (1, 6)
  vectorOf len $ frequency
    [ (3, Unitary <$> (genExpr q =<< choose (0, 3)))
    , (2, Measure <$> genQubits)
    , (1, do ks <- genQubits
             Initialize ks <$> vectorOf (length ks) arbitrary)
    ]
  where genQubits = do c <- choose (1, q)
                       take c <$> shuffle [0 .. q-1]

-- | Same outcomes, final state and remaining samples as the matrix backend.
prop_program :: Int -> Property
prop_program q =
  forAll (vectorOf q arbitrary) $ \bits -> forAll (genProgram q) $ \prog ->
  forAll (choose (1, 1000000)) $ \seed ->
    let (ψm, om, rm) = MS.evalProg prog (MS.ket bits) (samples seed)
        (ψs, os, rs) = SV.evalProg prog (SV.ket bits) (samples seed)
    in os === om .&&. ψs =~ ψm .&&. take 3 rs === take 3 rm

-- | The amplitudes in which qubit k of q reads c, the others set to 0.
keep :: Int -> Int -> Bit -> CMat -> CMat
keep q k c ψ = column [ if Bit (testBit x (q - 1 - k)) == c then z else 0
                      | (x, z) <- zip [0 :: Int ..] (toList (flatten ψ)) ]

-- | A state on 1 to 6 qubits and one of its qubits.
genMeasured :: Gen (Int, Int, CMat)
genMeasured = do
  q <- choose (1, 6)
  liftA2 ((,,) q) (choose (0, q-1)) (genState q)

-- | measure1 against the amplitudes: the outcome is 0 exactly when the sample
--   lies below the probability of 0, and the state is the normalised branch.
prop_measure :: Property
prop_measure =
  forAll genMeasured $ \(q, k, ψ) -> forAll (choose (0, 0.999999)) $ \r ->
    let p0            = (\n -> n * n) (norm_2 (flatten (keep q k 0 ψ)))
        c             = Bit (r >= p0)
        branch        = keep q k c ψ
        expected      = scale (1 / (norm_2 (flatten branch) :+ 0)) branch
        (ψ', os, rest) = SV.measure1 (ψ, [], r : samples 1) k
    in abs (r - p0) > 1e-9 ==>
         os === [c] .&&. ψ' =~ expected .&&. take 2 rest === take 2 (samples 1)

-- | A second measurement of the same qubit repeats the outcome and changes nothing.
prop_measure_twice :: Property
prop_measure_twice =
  forAll genMeasured $ \(_, k, ψ) -> forAll (choose (1, 1000000)) $ \seed ->
    let (ψ1, o1, rest1) = SV.measure1 (ψ, [], samples seed) k
        (ψ2, o2, _)     = SV.measure1 (ψ1, o1, rest1) k
    in o2 === o1 ++ o1 .&&. ψ2 =~ ψ1

-- | The two projections of a qubit are diagonal, complementary and idempotent.
prop_projections :: Property
prop_projections =
  forAll genMeasured $ \(q, k, ψ) ->
    let p c = SV.apply (SV.measureProjection q k c)
    in conjoin [ p 0 ψ =~ keep q k 0 ψ
               , p 1 ψ =~ keep q k 1 ψ
               , p 0 ψ + p 1 ψ =~ ψ
               , p 0 (p 0 ψ) =~ p 0 ψ
               , p 1 (p 0 ψ) =~ scale 0 ψ ]

-- | After Initialize ks vs, measuring ks gives vs, and Initialize itself
--   records no outcome and uses one sample per qubit.
prop_initialize :: Property
prop_initialize =
  forAll genMeasured $ \(q, _, ψ) -> forAll (genInit q) $ \(ks, vs) ->
  forAll (choose (1, 1000000)) $ \seed ->
    let (ψ1, o1, rest1) = SV.evalProg [Initialize ks vs] ψ (samples seed)
        (_,  o2, _)     = SV.evalProg [Measure ks] ψ1 rest1
    in conjoin [ o1 === []
               , o2 === vs
               , take 2 rest1 === take 2 (drop (length ks) (samples seed))
               , counterexample "norm" (abs (norm_2 (flatten ψ1) - 1) < tol) ]
  where genInit q = do c  <- choose (1, q)
                       ks <- take c <$> shuffle [0 .. q-1]
                       liftA2 (,) (pure ks) (vectorOf c arbitrary)

programTests :: TestTree
programTests = testGroup "Programs and measurement"
  [ budget 25 $ testGroup "evalProg agrees with the matrix backend"
      [ testProperty (show q ++ " qubits") (prop_program q) | q <- [1 .. 8] ]
  , budget 100 $ testGroup "Measurement"
      [ testProperty "measure1 against the amplitudes"           prop_measure
      , testProperty "measuring twice repeats the outcome"       prop_measure_twice
      , testProperty "the two projections of a qubit"            prop_projections
      , testProperty "Initialize sets the qubits and records nothing" prop_initialize
      ]
  ]

------------------------------------------------------------------------
-- Worked examples of the lecture notes
------------------------------------------------------------------------

exampleTests :: TestTree
exampleTests = testGroup "Worked examples"
  [ testCase "H ⊗ C X on the amplitudes 1 to 8" $
      assertClose "amplitudes"
        (sv (H ⊗ C X) oneToEight)
        (scale (1 / sqrt 2) (column [6, 8, 12, 10, -4, -4, -4, -4]))
  , testCase "H ⊗ C X on |010> is (|011> + |111>)/√2" $
      assertClose "state"
        (SV.apply (SV.evalOp (H ⊗ C X)) (SV.ket [0, 1, 0]))
        (scale (1 / sqrt 2) (SV.ket [0, 1, 1] + SV.ket [1, 1, 1]))
  , testCase "C X ⊗ H, the arities exchanged, on the amplitudes 1 to 8" $
      assertClose "amplitudes"
        (sv (C X ⊗ H) oneToEight)
        (scale (1 / sqrt 2) (column [3, -1, 7, -1, 15, -1, 11, -1]))
  , testCase "Measure [2,0] on (|000> + |101>)/√2 with samples 3/4, 1/4" $ do
      let (ψ, os, _) = SV.evalProg [Measure [2, 0]] start (3/4 : 1/4 : samples 1)
      assertClose "state" ψ (SV.ket [1, 0, 1])
      os @?= [1, 1]
  , testCase "Initialize [2,0] [0,0] on the same state and samples" $ do
      let (ψ, os, _) = SV.evalProg [Initialize [2, 0] [0, 0]] start (3/4 : 1/4 : samples 1)
      assertClose "state" ψ (SV.ket [0, 0, 0])
      os @?= []
    -- Chapter 2's check of the Deutsch program, by its four oracles. The
    -- amplitudes of the other outcome must cancel to exactly 0, or the sample
    -- 0.0 selects it.
  , testCase "Deutsch: the measured bit is f(0) ⊕ f(1), for the samples 0.0 and 0.999" $
      sequence_
        [ assertEqual (showOp oracle ++ ", sample " ++ show r) [balanced] os
        | (oracle, balanced) <- [(Id 2, 0), (Id 1 ⊗ X, 0), (C X, 1), ((Id 1 ⊗ X) ∘ C X, 1)]
        , r <- [0.0, 0.999]
        , let (_, os, _) = SV.evalProg [Unitary ((H ⊗ H) >: oracle >: (H ⊗ Id 1)), Measure [0]]
                                       (SV.ket [0, 1]) (repeat r) ]
  ]
  where start      = scale (1 / sqrt 2) (SV.ket [0, 0, 0] + SV.ket [1, 0, 1])
        oneToEight = column (map fromIntegral [1 .. 8 :: Int])

------------------------------------------------------------------------

main :: IO ()
main = defaultMain $ testGroup "StatevectorTest"
  [ matrixTests, ruleTests, algebraTests, interfaceTests, programTests, exampleTests ]
