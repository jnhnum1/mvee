module Main(main) where

import Control.Monad
import Control.Monad.State.Strict

import Data.Ellipsoid
import Data.Ellipsoid.Examples.Car

import Graphics.Gnuplot.Simple

import Numeric.Container
import System.Random

-- blargh
list2tup :: (Int, Int) -> [a] -> (a, a)
list2tup (i, j) coords = (coords !! i, coords !! j)

main :: IO ()
main = do
  let 
    initEllipsoid = 
      (toVector $ mkCar (mkLoc 25 7.5) (mkLoc 15 7.5) 0,
      scale 0.5 $ ident 5)
  let randPtGen = randPtIn initEllipsoid
  stdGen <- getStdGen
  forM_ [8000] $ \numPts -> do
      let pts = evalState (replicateM numPts randPtGen) stdGen
          afterPts = map (toVector . stepCar 5 . fromVector) pts
      print afterPts
      let
          afterEllipsoid = mvee 1e-3 afterPts
          tups = (list2tup (0, 1) . toList) `map` afterPts
      print afterEllipsoid
      plotDots [] tups
