---- L2 norm of a matrix
--
-- compilation:
--
--   ghc -O2 L2
--
-- usage example:
--
--   ./L2 test.mat
--
import qualified Data.Vector.Unboxed as V
import System.Environment


-- vector of 64 bit integers
type Vec = V.Vector(Int)

-- addition of two vectors
instance Num Vec where
  a + b = V.zipWith (+) a b

-- Manhattan norm (non-negative integer)
type Norm1 = Int

-- calculate Manhattan norm
norm1 :: Vec -> Norm1
norm1(v) = V.sum(V.map abs v)

---- L2 norm of a matrix
--
-- l2(m)  ==  max_{ma, mb: ma \/ mb == m} (norm1(sum(ma)) + norm1(sum(mb)))
--
-- where (ma \/ mb) is the union of sets of vectors ma and mb.
--
type L2Norm = Int

---- matrix annotated with the l2 norms of its submatrices
--
-- Each matrix row v_j is annotated with the l2 norm of [v_j, v_(j+1), ..., v_n]
-- where [v_1, v_2, ..., v_n] is the original matrix.
--
-- For example, if the original matrix is 
--   [v1, v2, v3] 
-- then the annotated matrix is
--   [(n1, v1), (n2, v2), (n3, v3)]
-- where
--   n1 = l2([v1, v2, v3])
--   n2 = l2([v2, v3])
--   n3 = l2([v3])
--
type MatWithNorms = [(L2Norm, Vec)]


---- the main algorithm
--
-- Semantics of f is given by the following equation:
--
-- f(va, vb, m, i)
--    ==
-- max(i, max_{ma, mb: ma \/ mb == m} (norm1(a + sum(ma)) + norm1(b + sum(mb))))
--
f :: (Vec, Vec, MatWithNorms, L2Norm) -> L2Norm
--
-- base case: the matrix has zero rows
--
f(va, vb, [], i) = max i (norm1(va) + norm1(vb))
--
-- recursive case: 
--   n is the l2 norm of the matrix
--   v is the first row of the matrix
--   m is the rest of the matrix
--
f(va, vb, (n, v): m, i)
    | i >= n + norm1(va) + norm1(vb) = i
    | otherwise                      = f(v + va, vb, m, f(va, v + vb, m, i))


---- annotate a matrix with the l2 norm of its submatrices
--
-- The result's first element is the l2 norm of the whole matrix.
--
g :: [Vec] -> (L2Norm, MatWithNorms)
g([]) = (0, [])
g(v: m) = (j, (j, v): mn)
  where
    (i, mn) = g(m)
    j = f(V.replicate (V.length v) 0, v, mn, i)


-- calculation of l2 norm only
l2 :: [Vec] -> L2Norm
l2 = fst . g

-- command line interface
main :: IO ()
main = do
    args <- getArgs
    case args of
      [name] -> do
        s <- readFile(name)
        (print . l2 . map (V.fromList . map read . words) . filter (not . null) . lines)(s)
      _ -> error "usage: L2 <filename>"
