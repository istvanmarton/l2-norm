-- optimized version of L2.hs
--
-- Optimizations:
--  - vector norms are cached (~ 2x speedup)
--  - not all submatrix norms are calculated (~ 2x speedup)

import qualified Data.Vector.Unboxed as V
import System.Environment

type Vec = V.Vector Int

instance Num Vec where
  a + b = V.zipWith (+) a b

type Norm1 = Int

norm1 :: Vec -> Norm1
norm1 v = V.sum (V.map abs v)

type VecWithNorm = (Vec, Norm1)

withNorm1 :: Vec -> VecWithNorm
withNorm1 v = (v, norm1 v)

type L2Norm = Int

type MatWithNorms = [(L2Norm, Vec)]

f :: VecWithNorm -> VecWithNorm -> MatWithNorms -> L2Norm -> L2Norm
f (va, na) (vb, nb) m i = case m of
  [] -> max (na + nb) i
  (n, v): m'
    | i >= n + na + nb -> i
    | otherwise -> f (withNorm1(v + va)) (vb, nb) m' (f (va, na) (withNorm1(v + vb)) m' i)

f' v m i = f (withNorm1(V.replicate (V.length v) 0)) (withNorm1 v) m i

g [] = (0, [])
g (v: m) = (j, (j, v): mn)  where
  (i, mn) = g m
  j = f' v mn i

conc m1 m2@((n, _): _) = map (\(a, b) -> (a + n, b)) m1 ++ m2

l2 :: [Vec] -> L2Norm
l2 (v: m) = f' v (mn1 `conc` mn2) (max n1 n2)
  where
    (m1, m2) = splitAt (length m `div` 2) m

    (n1, mn1) = g m1
    (n2, mn2) = g m2

main :: IO ()
main = do
    [f] <- getArgs
    s <- readFile f
    print . l2 . map (V.fromList . map read . words) . filter (not . null) $ lines s
