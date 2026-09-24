module Main where
import QubitOperators
import MatrixSemantics
import PrettyOp
import PrettyMatrix

Rz = R Z

{-| Single-qubit Mach-Zender Interferometer without which-path detector -}
mzi_circuit :: RealT -> Program
mzi_circuit θ = [ H <> Rz θ <> H, Measure [0] ]

{-|
 Mach–Zehnder Interferometer with non-destructive which-path detector.

 Two-qubit model:
      - path qubit p
      - which-path ancilla a

    Circuit:
      |0>_p -- H -- Rz(phi) --●-- H -- measure p
                              |
      |0>_a ------------------X--------
                 (CNOT = perfect which-path marker)
-}
mziwp_circuit :: RealT -> Program
mziwp_circuit θ = [
    H  ⊗   I <>  -- First beam splitter
    Rz θ ⊗ I <>  -- Relative phase between the two arms (adjustable with mirror distance)
    C      X <>  -- Non-destructive which-path detector (e.g. polarization tagging)
    H  ⊗   I,    -- Second beam splitter
    Measure [0]]

measure1 :: Int ->  StateT -> Double -> (Int,StateT) -- Measure ks psi = (outcomes,collapsed_state)
measure1 _ _ [] = "Random numbers exhausted."
measure1 k ψ r  = let
        n  = ilog2 (rows ψ)
        p0 = measureProjection n k 0
        p1 = measureProjection n k 1

        prob0   = norm (p0<>ψ)
        outcome = if r < p0 then 0 else 1
        collapsed_state = if r < p0 then p0<>ψ else p1<>ψ
    in
        (outcome, collapsed_state)

measure :: StateT -> [Int] -> RNG -> (([Int], StateT), RNG)
measure ψ ks rng = let 
    (rs, rng') = splitAt (length ks) rng
    
    apply_measure1 :: ([Int], StateT) -> (Int, Double) -> ([Int], StateT)
    apply_measure1 (outs, state) (k, r) = let
            (outcome, collapsed_state) = measure1 k state r
        in
            (outs ++ [outcome], collapsed_state)
    in
        (foldl apply_measure1 ([], ψ) (zip ks rs), rng')
    
-- TODO: Random monad?
outcomes :: Program -> Int -> RNG -> [Int]
outcomes experiment repetitions rng = let
    n       = prog_width experiment    -- Number of qubits
    p0      = measureProjection n 0 0 -- Projection operator for measuring Q0 = 0
    ψ0      = ket $ replicate 0 n
    run_exp = evalProg experiment ψ0
    evalShot rng = let 
            (ψ_final,rng') = run_exp rng 
        in case rng' of 
            [] -> error "No more random numbers. This never happens."
             r:rng'' -> (fst $ measure1 ψ_final 0 r, rng'')




main :: IO()
main = do

