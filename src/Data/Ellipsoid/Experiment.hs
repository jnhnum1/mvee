module Main(main) where
-- module Data.Ellipsoid.Experiment where

import Control.Monad
import Control.Monad.State.Strict

import Data.Accessor
import Data.Ellipsoid
import Data.Ellipsoid.Examples.Car

import Graphics.Rendering.Chart hiding (Vector, Point)
import Graphics.Rendering.Chart.Gtk

import Numeric.Container
import System.Random

-- blargh
list2tup :: (Int, Int) -> [a] -> (a, a)
list2tup (i, j) coords = (coords !! i, coords !! j)

unitSample :: Int -- ^ the number of samples per dimensions
           -> Int -- ^ the number of dimensions
           -> [Point]
unitSample k d = 
  let choices = [-1, -1 + 2 / fromIntegral (k - 1)..1] 
    in map fromList $ filter withinUnit $ replicateM d choices
  where
    withinUnit xs = sum (map (\x -> x^2) xs) <= 1.0001

latticeSample :: Int -- ^ the number of samples per dimension (in the unit
                     -- ellipsoid)
              -> Ellipsoid -- the ellipsoid which you are sampling from
              -> [Point]
latticeSample k ell = 
  let dim = ellipsoidDim ell
      unitSamples = unitSample k dim
      transformation = transformFromUnit ell -- inverse of l^t
  in map transformation unitSamples

makeConvergencePlot :: Ellipsoid 
                    -> (Point -> Point)
                    -> [Int] -- ^ numbers of samples per dimension
                    -> PlotPoints Int Double
makeConvergencePlot ell trans sampleNums =
  let points = [(i, latticeSample i ell) | i <- sampleNums]
      volumes = [(i, volumeFactor $ mvee 1e-3 pts) | (i, pts) <- points]
  in
    plot_points_values ^= volumes $ defaultPlotPoints

makeRandConvergencePlot :: (RandomGen g, MonadState g m)
                    => Ellipsoid 
                    -- ^ the initial ellipsoid from which to sample
                    -> (Point -> Point) 
                    -- ^ the transformation to apply to the sample points
                    -> [Int]
                    -> Int
                    -> m (PlotLines Int Double)

makeRandConvergencePlot ell trans numPtses reps = do
  ptses <- replicateM reps genPoints
  return $! plot_lines_values ^= ptses
         $ defaultPlotLines
  where 
    genPoints = do
      ellipsoids <- mapM (sampleImage ell trans) numPtses
      let volumePts :: [(Int, Double)]
          volumePts = zip numPtses (map volumeFactor ellipsoids)
      return $! volumePts

main :: IO ()
main = do
  let 
    initEllipsoid = 
      (toVector $ mkCar (mkLoc 25 7.5) (mkLoc 15 7.5) 0,
      scale 0.5 $ ident 5)
  stdGen <- getStdGen
  let 
      transform = (toVector . stepCar 5 . fromVector)
      plot = plot_points_title ^= "MVEE volume" $ 
          makeConvergencePlot initEllipsoid transform [1..9]
      layout = layout1_title ^= "Sampling Convergence" 
             $ layout1_plots ^= [Left (toPlot plot)]
             $ defaultLayout1
  renderableToWindow (toRenderable layout) 640 480
