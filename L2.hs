-- optimized version of L2.hs
--
-- Optimizations:
--  - vector norms are cached
--  - not all submatrix norms are calculated
--  - parallel execution

import System.Environment
import qualified Data.Vector.Unboxed as V
import Control.Monad
import Control.Monad.Par.IO
import Control.Monad.IO.Class
import Control.Monad.Par.Class
import Data.IORef
import Control.Concurrent

type Vec = V.Vector Int

(.+.) :: Vec -> Vec -> Vec
a .+. b = V.zipWith (+) a b

type Norm1 = Int

norm1 :: Vec -> Norm1
norm1 v = V.sum (V.map abs v)

type VecWithNorm = (Vec, Norm1)

withNorm1 :: Vec -> VecWithNorm
withNorm1 v = (v, norm1 v)

type L2Norm = Int

type Mat = [Vec]
type MatWithNorms = [(L2Norm, Vec)]

f :: VecWithNorm -> VecWithNorm -> MatWithNorms -> L2Norm -> L2Norm
f (va, na) (vb, nb) m i = case m of
  [] -> max (na + nb) i
  (n, v): m'
    | i >= n + na + nb -> i
    | otherwise -> f (withNorm1(v .+. va)) (vb, nb) m' (f (va, na) (withNorm1(v .+. vb)) m' i)

fPar :: IORef L2Norm -> VecWithNorm -> VecWithNorm -> MatWithNorms -> ParIO L2Norm
fPar best (va, na) (vb, nb) m = case m of
  (n, v): m' | n == maxBound -> do
    c1 <- spawn_ (fPar best (withNorm1 $ v .+. va) (vb, nb) m')
    c2 <- spawn_ (fPar best (va, na) (withNorm1 $ v .+. vb) m')
    liftM2 max (get c1) (get c2)
  _ -> do
    i <- liftIO $ readIORef best
    let i' = f (va, na) (vb, nb) m i
    i' `seq` liftIO (atomicModifyIORef best (\i -> (max i' i, ())))
    pure i'

{-# INLINE f' #-}
f' :: L2Norm -> Vec -> MatWithNorms -> IO L2Norm
f' guess v vs = do
  best <- newIORef guess
  runParIO $ fPar best (withNorm1 $ V.replicate (V.length v) 0) (withNorm1 v) vs

g :: Int -> Mat -> IO MatWithNorms
g _ [] = pure []
g l (v: vs) = do
  m <- g (l-1) vs
  let m' = zipWith (\i (a, b) -> (if i < (length vs + 1) `div` 4 then maxBound else a, b)) [0..] m
  i <- if l > 0 then pure maxBound else f' 0 v m'
  pure ((i, v): m)

l2 :: Int -> Mat -> IO L2Norm
l2 guess (v:vs) = do
  m <- g ((length vs + 1) `div` 4) vs
  f' guess v m

main :: IO ()
main = do
  args <- getArgs
  case args of
    [guess, f] -> do
      s <- readFile f
      let m = (map (V.fromList . map read . words) . takeWhile (not . null) . lines) s
      i <- l2 (read guess) m
      print i
