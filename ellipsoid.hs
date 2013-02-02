import Control.Monad
import Control.Monad.State.Class
import Control.Monad.State.Strict

import Data.Functor
import Data.Functor.Identity
import Data.Packed.Matrix

import Graphics.Gnuplot.Simple

import Foreign.Storable

import Numeric.Container
import Numeric.LinearAlgebra.Algorithms
import Numeric.LinearAlgebra.Util

import System.Random

import Prelude hiding (max)

type Point = Vector Double

-- the ellipsoid (c, A) is the set of points x such that
-- (x - c)'A(x - c) <= 1
type Ellipsoid = (Point, Matrix Double)

-- TODO make tail-recursive
findMax :: Ord a => [a] -> (Int, a)
findMax [x] = (0, x)
findMax (x:xs) = 
  let (i, m) = findMax xs in
    if x > m
      then (0, x)
      else (1 + i, m)

-- > iterateCollect f x
-- > [[x], [f x, x], [f f x, f x, x], ... ]
iterateCollect :: (a -> a) -> a -> [[a]]
iterateCollect f x =
  let y = [x] : (map addF y)
  in y
    where addF xs@(x:_) = (f x) : xs

iterateUntil :: ([a] -> Bool) -> (a -> a) -> a -> a
iterateUntil pred f x =
  head $ head $ dropWhile (not . pred) (iterateCollect f x)

-- In the absence of the ST monad, can use this to update an index of a vector
-- according to a given transformation function.
modifyVec :: (Storable a) => Int -> (a->a) -> Vector a -> Vector a
modifyVec i f v =
  mapVectorWithIndex g v
    where 
      g j vj =
        if i == j 
          then f vj
          else vj

-- mvee eps pts approximately computes the minimum volume ellipsoid containing
-- all the points, where eps is a small number used for convergence testing.
mvee :: Double -> [Point] -> Ellipsoid
mvee eps pts =
  let d :: Num a => a -- stupid monomorphism restriction
      d = fromIntegral $ dim $ head pts

      numPts :: Num a => a -- stupid monomorphism restriction
      numPts = fromIntegral $ length pts

      p :: Matrix Double
      p = fromColumns pts

      q :: Matrix Double
      q = fromBlocks [[p], [ones 1 numPts]]

      q' :: Matrix Double
      q' = trans q

      u0 :: Vector Double
      u0 = scale (1 / numPts) (numPts |> (repeat 1))

      withinThreshold :: [Vector Double] -> Bool
      withinThreshold [] = False
      withinThreshold [x] = False
      withinThreshold (x : y : _) = norm (x - y) < eps

      updateU :: Vector Double -> Vector Double
      updateU u =
        let x = q <> diag u <> q'
            m = takeDiag $ q' <> inv x <> q
            (i, max) = findMax $ toList m
            step_size = (max - d - 1)  / ((d + 1) * (max - 1))
        in modifyVec i (+step_size) (scale (1 - step_size) u)
  in 
    let u = iterateUntil withinThreshold updateU u0
        pu = p <> u
        puMat = asColumn pu
    in (pu, (1 / d) * inv (p <> (diag u) <> (trans p) - puMat <> (trans puMat))) -- ugly :(

newtype VolumeOrd = VolumeOrd {getEllipsoid :: Ellipsoid} deriving (Eq)
instance Ord VolumeOrd where
  compare v1 v2 =
    let (_, mat1) = getEllipsoid v1
        (_, mat2) = getEllipsoid v2
    in compare (det mat2) (det mat1)

-- what is the volume of the ellipsoid relative to the unit n-sphere
volumeFactor :: Ellipsoid -> Double
volumeFactor (_, mat) = 1 / sqrt (det mat)

-- A version of randomR which is polymorphic over a monad stack with a MonadState
-- component storing state of type RandomGen.
randomRM :: (Random a, RandomGen g, MonadState g m) => 
  (a, a) -> m a
randomRM bounds = do
  gen <- get
  let (result, newgen) = randomR bounds gen
  put newgen
  return result

-- takes a dimension and returns a monadic computation of a random point in the unit sphere of tha dimension.
-- NOTE this is not currently uniform sampling from the ball.
randPtInUnit :: (RandomGen g, MonadState g m) => Int -> m Point
randPtInUnit n = 
  (fromList . snd) `liftM` helper n
  where
    helper 0 = return (0, [])
    helper n = do
      (normSoFarSq, coordsSoFar) <- helper (n - 1)
      let coordBound = sqrt (1 - normSoFarSq)
      nextCoord <- randomRM (-coordBound, coordBound)
      return (nextCoord^2 + normSoFarSq, nextCoord : coordsSoFar)

-- Note that a unit ball is characterized by (x - c)'(x - c) <= 1
-- An ellipsoid which is a transformation by A of the unit ball is 
-- thus represented by the following:
-- (A^-1 x)'(A^-1 x) <= 1
-- x'(A^-t A^-1) x <= 1
--
-- We represent ellipsoids as the middle part, so by getting the cholesky 
-- factorization, we get A^-1, whose inverse is then A.
-- So we can generate random points from our ellipsoid by generating them 
-- from the unit ball, and then applying A.
randPtIn :: (RandomGen g, MonadState g m) => Ellipsoid -> m Point
randPtIn (center, mat) = let
  upperInv = inv (chol mat) in
    do
      unitPt <- randPtInUnit (dim center)
      return $ (upperInv <> unitPt) + center

repeatM :: (Monad m) => Int -> m a -> m [a]
repeatM n m =
  sequence $ replicate n m

-- blargh
list2tup :: [a] -> (a, a)
list2tup [x, y] = (x, y)

main = do
  let points = fromList <$> [[1,0],[0,1],[-1, 0]]
      ellipsoid = mvee 1e-4 points
  print ellipsoid
  let randPtGen = randPtIn ellipsoid
  stdGen <- getStdGen
  let pts = evalState (repeatM 10000 randPtGen) stdGen
      tups = (list2tup . toList) `map` pts 
  plotDots [] tups 
