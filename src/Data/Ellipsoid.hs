{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies, FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances, UndecidableInstances #-}

module Data.Ellipsoid(
  Point,
  Ellipsoid,
  randPtIn,
  mvee,
  volumeFactor,
  VolumeOrd) where 

import Control.Monad
import Control.Monad.State.Class

import Data.List(tails, find)
import Data.Maybe(fromJust)
import Data.Packed.Matrix

import Debug.Trace

import Foreign.Storable (Storable)

import Numeric.Container hiding (find, Mul(..))
import qualified Numeric.Container as C
import Numeric.LinearAlgebra.Algorithms (chol,inv,det,pnorm,NormType(..))
import Numeric.LinearAlgebra.Util

import System.Random

import Prelude hiding (max)

type Point = Vector Double

-- the ellipsoid (c, A) is the set of points x such that
-- (x - c)'A(x - c) <= 1
type Ellipsoid = (Point, Matrix Double)

-- sigh, hacking around the poorly designed hmatrix api
class Mul a b c | a b -> c where
  infixl 7 <>
  (<>) :: a -> b -> c

instance (Product e) => Mul (Matrix e) (Matrix e) (Matrix e) where
  (<>) = (C.<>)

instance (Product e) => Mul (Matrix e) (Vector e) (Vector e) where
  (<>) = (C.<>)

instance (Product e) => Mul (Vector e) (Matrix e) (Vector e) where
  (<>) = (C.<>)

newtype Diag a = Diag [a]

-- instance (Num a) => Mul (Diag a) (Diag a) (Diag a) where
  -- (Diag d) <> (Diag f) = Diag $ zipWith (*) d f

-- instance (Container Vector a) => Mul (Diag a) (Matrix a) (Matrix a) where
  -- (Diag d) <> m = fromRows $ zipWith scale d (toRows m)

instance (Container Vector a) => Mul (Matrix a) (Diag a) (Matrix a) where
  m <> (Diag d) = fromColumns $ zipWith scale d (toColumns m)

-- TODO make tail-recursive
findMax :: Ord a => [a] -> (Int, a)
findMax [x] = (0, x)
findMax (x:xs) = 
  let (i, m) = findMax xs in
    if x > m
      then (0, x)
      else (1 + i, m)

iterateUntil :: ([a] -> Bool) -> (a -> a) -> a -> a
iterateUntil predicate f x =
  head . fromJust . find predicate . tails $ iterate f x 

-- In the absence of the ST monad, can use this to update an index of a vector
-- according to a given transformation function.
modifyVec :: (Storable a) => Int -> (a->a) -> Vector a -> Vector a
modifyVec i f =
  mapVectorWithIndex g
    where 
      g j vj =
        if i == j 
          then f vj
          else vj

-- diagMult a b = takeDiag (a <> b)
diagMult :: Product e => Matrix e -> Matrix e -> Vector e
diagMult a b = 
  fromList $ zipWith dot (toRows a) (toColumns b)

-- mvee eps pts approximately computes the minimum volume ellipsoid containing
-- all the points, where eps is a small number used for convergence testing.
-- TODO there is still a big memory leak in here somewhere. search and destroy.
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
      u0 = scale (1 / numPts) (numPts |> repeat 1)

      withinThreshold :: [Vector Double] -> Bool
      withinThreshold [] = False
      withinThreshold [_] = False
      withinThreshold (x : y : _) = 
        let
          diff = norm (x - y)
        in trace (show diff) (diff < eps)

      updateU :: Vector Double -> Vector Double
      updateU u =
        let x = q <> Diag (toList u) <> q'
            m = diagMult q' (inv x <> q)
            (i, max) = findMax $ toList m
            step_size = (max - d - 1)  / ((d + 1) * (max - 1))
        in modifyVec i (+step_size) (scale (1 - step_size) u)
  in 
    let u = iterateUntil withinThreshold updateU u0
        pu = p <> u
        puMat = asColumn pu
    in (pu, scale (1 / d) $ inv $ p <> Diag (toList u) <> trans p -
                                  puMat <> trans puMat)

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
randPtInUnit n {- dimension -} = 
  do
    p <- fromList `liftM` replicateM n (randomRM (-1, 1))
    if pnorm Frobenius p <= 1
      then return p
      else randPtInUnit n

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
  symmetrized = scale 0.5 (mat + trans mat)
  upperInv = inv (chol symmetrized) in
    do
      unitPt <- randPtInUnit (dim center)
      return $ (upperInv <> unitPt) + center

