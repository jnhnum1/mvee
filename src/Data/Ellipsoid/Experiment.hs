module Main(main) where
-- module Data.Ellipsoid.Experiment where

import Control.Lens hiding ((^=))

import Control.Monad
import Control.Monad.State.Strict

import Data.Accessor ((^=))
import Data.Ellipsoid
import Data.Ellipsoid.Examples.Car

import Data.Colour
import Data.Colour.Names

import Debug.Trace

import Graphics.Rendering.Chart hiding (Vector, Point)
import Graphics.Rendering.Chart.Gtk 
import Numeric.Container
import Numeric.LinearAlgebra.Algorithms (inv)
import System.Random

list2tup :: (Int, Int) -> [a] -> (a, a)
list2tup (i, j) coords = (coords !! i, coords !! j)


pointPlotConfig :: PlotPoints x y
pointPlotConfig = plot_points_style ^= exes 5 1 (opaque blue) $ defaultPlotPoints


unitSample :: Int -- ^ the number of samples per dimensions
           -> Int -- ^ the number of dimensions
           -> [Point]
unitSample k d = 
  let choices = [-1, -1 + 2 / fromIntegral (k - 1)..1] 
    in map fromList $ filter withinUnit $ replicateM d choices
  where
    withinUnit xs = sum (map (\x -> x^2) xs) <= 1.0001

skewUnitSample :: Int -- ^ the number of samples per dimensions
           -> Int -- ^ the number of dimensions
           -> [Point]
skewUnitSample k d = 
  let choices = [-1, -1 + 2 / fromIntegral (k - 1)..1] 
    in map (fromList . skew) $ replicateM d choices
  where
    skew xs = 
      let alpha = minimum [abs (1 / x) | x <- xs]
          alphanorm = sqrt $ sum [alpha^2 * x^2 | x <- xs] 
      in [x / alphanorm | x <- xs]

latticeSample :: Int -- ^ the number of samples per dimension (in the unit
                     -- ellipsoid)
              -> Ellipsoid -- the ellipsoid which you are sampling from
              -> [Point]
latticeSample k ell = 
  let dim = ellipsoidDim ell
      unitSamples = unitSample k dim
      transformation = transformFromUnit ell -- inverse of l^t
  in map transformation unitSamples

skewSample :: Int -- ^ the number of samples per dimension (in the unit
                     -- ellipsoid)
              -> Ellipsoid -- the ellipsoid which you are sampling from
              -> [Point]
skewSample k ell = 
  let dim = ellipsoidDim ell
      unitSamples = skewUnitSample k dim
      transformation = transformFromUnit ell -- inverse of l^t
  in map transformation unitSamples

makeSamples :: Int -- ^ number samples per dimension
           -> Ellipsoid
           -> (Point -> Point)
           -> ImageSamples
makeSamples k ell trans = 
  map trans (latticeSample k ell)

scatterPlot :: ImageSamples -> (Int, Int) -> PlotPoints Double Double
scatterPlot pts (i1, i2) =
  let pts2d = map (list2tup (i1, i2) . toList) pts
  in plot_points_values ^= pts2d $ pointPlotConfig

pointPlot :: [(Int, Double)] -> PlotPoints Int Double
pointPlot pts = plot_points_values ^= pts $ pointPlotConfig
  
makeConvergencePlot :: [(Int, ImageSamples)]
                    -> PlotPoints Int Double
makeConvergencePlot samples = 
  let points = mapped._2 %~ volumeFactor.mvee $ samples
  in
    pointPlot points

buildMapM :: Monad m => (a -> m b) -> [a] -> m [(a,b)]
buildMapM f xs = liftM (zip xs) (mapM f xs)

makeRandConvergencePlot :: (RandomGen g, MonadState g m)
                    => Ellipsoid 
                    -- ^ the initial ellipsoid from which to sample
                    -> (Point -> Point) 
                    -- ^ the transformation to apply to the sample points
                    -> [Int]
                    -> m (PlotLines Int Double)
makeRandConvergencePlot ell trans numPtses = do
  pts <- genPoints
  return $! plot_lines_values ^= [pts]
         $ defaultPlotLines
  where 
    genPoints = do
      samples <- buildMapM (randomSample ell) numPtses 
      let 
          imageSamples = mapped . _2 . mapped %~ trans $ samples
          volumePts = mapped . _2 %~ volumeFactor . mvee $ imageSamples
      return $! volumePts

main :: IO ()
main = do
  let 
    initEllipsoid = 
      (toVector $ mkCar (mkLoc 25 7.5) (mkLoc 15 7.5) 0,
        scale 0.5 $ ident 5, inv $ scale 0.5 $ ident 5)
    transform = (toVector . stepCar 5 . fromVector)
    samples = flip map [3..9] $ (\k -> 
        (k, map transform $ skewSample k initEllipsoid))
        --(k, map transform $ latticeSample k initEllipsoid))
    layout = layout1_title ^= "Skew Image Points (car.backposition)"
           $ layout1_plots ^= [Left (toPlot $ makeConvergencePlot samples)]
           $ defaultLayout1
  renderableToWindow (toRenderable layout)
  forM_ samples $ \(k, sample) -> do
    let
      plot = scatterPlot sample (0, 1)
    renderableToWindow (toRenderable layout) 640 480
