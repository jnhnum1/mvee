{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies, FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances, UndecidableInstances #-}

module Data.Ellipsoid(
  Point,
  Ellipsoid,
  mvee,
  project,
  volumeFactor,
  randomSample,
  ellipsoidDim,
  transformFromUnit,
  ImageSamples,
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
import Numeric.LinearAlgebra.Algorithms hiding (i, eps)
import Numeric.LinearAlgebra.Util

import System.Random

import Prelude hiding (max)

type Point = Vector Double
type ImageSamples = [Point]
type Index = Int

-- the ellipsoid (c, L^t, L^-t) is the set of points x such that
-- (x - c)'LL^T(x - c) <= 1
-- equivalently
-- (c + L^-t w), for ||w|| <= 1
-- equivalently
-- {x | ||L^t(x - c)|| <= 1}
type Ellipsoid = (Point, Matrix Double, Matrix Double)

ellipsoidDim :: Ellipsoid -> Int
ellipsoidDim (_, mat, _) = rows mat

transformFromUnit :: Ellipsoid -> Vector Double -> Vector Double
transformFromUnit (c, _, invmat) v = 
  c + invmat <> v

lq :: Field t => Matrix t -> (Matrix t, Matrix t)
lq mat = 
  let (q, r) = qr (trans mat)
  in (trans r, trans q)
-- see https://tcg.mae.cornell.edu/pubs/Pope_FDA_08.pdf
-- Produces an ellipsoid projected onto the given axes
project :: Ellipsoid -> [Vector Double] -> Ellipsoid
project (c, lt, ltinv) basis = 
  let t = fromColumns basis
      t' = trans t
      c' = t' <> c
      (u, s) = leftSV $ t' <> ltinv
      (l', _) = lq $ u <> inv (diag s)
  in
    (c', l', inv l')

-- sigh, hacking around the poorly designed hmatrix api
-- this class replaces the hmatrix's Mul class.
class Mul a b c | a b -> c where
  infixl 7 <>
  (<>) :: a -> b -> c

instance (Product e) => Mul (Matrix e) (Matrix e) (Matrix e) where
  (<>) = (C.<>)

instance (Product e) => Mul (Matrix e) (Vector e) (Vector e) where
  (<>) = (C.<>)

instance (Product e) => Mul (Vector e) (Matrix e) (Vector e) where
  (<>) = (C.<>)

-- This type represents a diagonal matrix, but more space- and
-- time-efficiently.
newtype Diag a = Diag [a]

-- Technically I don't need the first two of these, but they're included for
-- completeness
instance (Num a) => Mul (Diag a) (Diag a) (Diag a) where
  (Diag d) <> (Diag f) = Diag $ zipWith (*) d f

instance (Container Vector a) => Mul (Diag a) (Matrix a) (Matrix a) where
  (Diag d) <> m = fromRows $ zipWith scale d (toRows m)

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

-- This takes an almost-symmetric matrix, and makes it symmetric.  This is
-- necessary if, for example, we want to take the cholesky factorization.
symmetrized :: Matrix Double -> Matrix Double
symmetrized mat = scale 0.5 $ mat + trans mat

-- diagMult a b = takeDiag (a <> b)
-- this is a more efficient way of computing the diagonal of a product of two
-- matrices.  It is particularly more efficient when, as below, we are
-- multiplying a (large x few) matrix and a (few x large) matrix.
diagMult :: Product e => Matrix e -> Matrix e -> Vector e
diagMult a b = 
  fromList $ zipWith dot (toRows a) (toColumns b)

-- mvee eps pts approximately computes the minimum volume ellipsoid containing
-- all the points, where eps is a small number used for convergence testing.
-- This algorithm comes from
-- http://stackoverflow.com/questions/1768197/bounding-ellipse/1768440#1768440
mvee' :: Double -> [Point] -> Ellipsoid
mvee' eps pts =
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
        norm (x - y) < eps

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
        lt = chol $ symmetrized $ scale (1 / d) $ inv $ 
                  p <> Diag (toList u) <> trans p - puMat <> trans puMat
    in trace "mvee done" $ (pu, lt, inv lt)

mvee :: [Point] -> Ellipsoid
mvee = mvee' 1e-3

newtype VolumeOrd = VolumeOrd {getEllipsoid :: Ellipsoid} deriving (Eq)
instance Ord VolumeOrd where
  compare (VolumeOrd ell1) (VolumeOrd ell2) = 
    compare (volumeFactor ell1) (volumeFactor ell2)

-- what is the volume of the ellipsoid relative to the unit n-sphere
volumeFactor :: Ellipsoid -> Double
volumeFactor (_, mat, _) = 1 / sqrt (det $ trans mat <> mat)

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
randPtIn (center, mat, upperInv) = 
  do
    unitPt <- randPtInUnit (dim center)
    return $ (upperInv <> unitPt) + center

randomSample :: (RandomGen g, MonadState g m) => 
             Ellipsoid -- ^ sampled ellipsoid
             -> Int -- ^ the number of sample points
             -> m ImageSamples
randomSample init numPts = replicateM numPts (randPtIn init)
