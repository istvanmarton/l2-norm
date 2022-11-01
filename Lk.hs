{-# language BangPatterns, ViewPatterns #-}

import Control.Concurrent
import Control.Monad
import Control.Monad.Par.IO
import Control.Monad.IO.Class
import Control.Monad.Par.Class
import Data.IORef
import Data.List
import Data.Maybe
import qualified Data.Vector.Unboxed as V
import Options.Applicative


type Vec = V.Vector Int

{-# INLINE zero #-}
zero :: Vec -> Vec
zero v = V.replicate (V.length v) 0

(.+.) :: Vec -> Vec -> Vec
a .+. b = V.zipWith (+) a b

type Norm1 = Int

norm1 :: Vec -> Norm1
norm1 v = V.sum (V.map abs v)

type VecWithNorm = (Vec, Norm1)
type VecsWithNorm = [VecWithNorm]

{-# INLINE withNorm1 #-}
withNorm1 :: Vec -> VecWithNorm
withNorm1 v = (v, norm1 v)

type L2Norm = Int
type Witness = [VecsWithNorm]
type L2Witness = (L2Norm, Witness)

{-# INLINE maxWitness #-}
maxWitness :: L2Witness -> L2Witness -> L2Witness
maxWitness w1 w2 = if fst w1 < fst w2 then w2 else w1

maxWitness' :: L2Witness -> L2Witness -> L2Witness
maxWitness' w1 w2 = case compare (fst w1) (fst w2) of
  LT -> w2
  GT -> w1
  EQ -> if null (snd w1) then w2 else w1

type Mat = [Vec]
type MatWithNorms = [(L2Norm, Vec)]
type Order = Int

summ :: VecsWithNorm -> L2Norm
summ [] = 0
summ ((_, n): _) = n

{-# INLINE withN #-}
withN v [] = (v, norm1 v): []
withN v vs@((_, n): _) = (v, n + norm1 v): vs

{-# INLINE cons #-}
cons (va, n) vs rs = (va, n - summ vs + summ rs): rs

fD :: Witness -> VecsWithNorm -> VecsWithNorm -> MatWithNorms -> L2Witness -> L2Witness
fD w rs vs m i = case m of
  [] -> maxWitness (summ rs + summ vs, w) i
  (n, v): m'
    | fst i >= n + summ rs + summ vs -> i
    | otherwise -> fL w rs vs v m' (fR w rs vs v m' i)

fR w rs (vv@(va, n): vs) v m i
   | n > 0 = fD (rs: w) rs (withN(v .+. va) vs) m (fR w (cons vv vs rs) vs v m i)
   | otherwise = fD (rs: w) rs (withNorm1(v): vs) m i
fR w _ _ v m i = i

fL w (vv@(va, _): rs) vs v m i = fL w rs (cons vv rs vs) v m (fD (rs: w) rs (withN(v .+. va) vs) m i)
fL w _ _ v m i = i

fDP :: Witness -> VecsWithNorm -> VecsWithNorm -> MatWithNorms -> IORef L2Witness -> ParIO L2Witness
fDP w rs vs m i = case m of
  [] -> pure (summ rs + summ vs, w)
  ((== maxBound) -> True, v): m' -> fLP w rs vs v m' i `merge` fRP w rs vs v m' i
  _ -> do
    i_ <- liftIO $ readIORef i
    let i' = fD w rs vs m i_
    i' `seq` liftIO (atomicModifyIORef i (\i_ -> (maxWitness i' i_, ())))
    pure i'

fRP w rs (vv@(va, n): vs) v m i
  | n > 0 = fDP (rs: w) rs (withN(v .+. va) vs) m i `merge` fRP w (cons vv vs rs) vs v m i
  | otherwise = fDP (rs: w) rs (withNorm1(v): vs) m i
fRP w _ _ v m i = pure (0, [])

fLP w (vv@(va, _): rs) vs v m i = fLP w rs (cons vv rs vs) v m i `merge` fDP (rs: w) rs (withN(v .+. va) vs) m i
fLP w _ _ v m i = pure (0, [])

merge a b = do
  x <- spawn a
  y <- spawn b
  liftM2 maxWitness' (get x) (get y)

{-# INLINE f' #-}
f' :: Order -> L2Norm -> Vec -> MatWithNorms -> IO L2Witness
f' order guess v vs = do
  best <- newIORef (guess, [])
  runParIO $ fDP [[]] [] (withNorm1 v: replicate (order-1) (withNorm1 (V.replicate (V.length v) 0))) vs best

g :: Order -> Int -> Mat -> IO MatWithNorms
g _ _ [] = pure []
g order l (v: vs) = do
  m <- g order (l-1) vs
  let m' = zipWith (\i (a, b) -> (if i < (length vs + 1) `div` 4 then maxBound else a, b)) [0..] m
  i <- if l > 0 then pure maxBound else fst <$> f' order 0 v m'
  pure ((i, v): m)

l2 :: Int -> Order -> L2Norm -> Mat -> IO L2Witness
l2 depth order guess (v:vs) = do
  m <- g order depth vs
  f' order guess v m

compute :: Maybe Int -> Order -> L2Norm -> String -> IO ()
compute depth order guess f = do
  s <- readFile f
  let m = (map (V.fromList :: [Int] -> Vec) . filter (not . null) . map (map read . words) . lines) s
  (i, w_) <- l2 (fromMaybe (length m `div` 4) depth) order guess m

  let
    w = map length (reverse w_)
    -- `norm` is redifined here for better compiler optimizations
    norm l = V.sum (V.map abs (foldl (V.zipWith (+)) (V.replicate (V.length (head m)) 0) l))
    rows = zip w m
    i' = sum [norm [r | (j', r) <- rows, j' == j] | j <- [0..order-1]]

  if i == i' then putStr $ unlines $
      [ "Row numbers in the partitions:" ]
      ++ [ "  " ++ unwords [show i | (j', i) <- zip w [1..], j' == j] | j <- [0..order-1]] ++
      [ "L" ++ show order ++ " norm:"
      , "  " ++ show i'
      ]
    else error $ unlines
      [ "!!! ERROR !!!"
      , "Computed L" ++ show order ++ " norm:"
      , "  " ++ show i
      , "Recalculated norm:"
      , "  " ++ show i'
      , "Maybe the guessed value was too high, try with a lower guessed value."
      ]


main :: IO ()
main = join (execParser opts)
 where
    opts = info (helper <*> options)
      ( fullDesc
     <> progDesc "Partition the rows of a matrix in n disjoint sets such that the sum of their summed Manhattan norm is maximal."
      )

    options :: Parser (IO ())
    options = compute
      <$> optional (option auto $ short 'd' <> long "depth" <> metavar "NAT" <> help "parallel depth - default is height/4")
      <*> (fromMaybe 2 <$> optional (option auto $ short 'o' <> long "order" <> metavar "NAT" <> help "order - default is 2" <> completeWith ["2"]))
      <*> (fromMaybe 0 <$> optional (option auto $ short 'g' <> long "guessed" <> metavar "NAT" <> help "guessed result - default is 0" <> completeWith ["0"]))
      <*> (argument str (metavar "FILE" <> action "filename"))
