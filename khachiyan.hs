import Data.Functor
import Data.Packed.Matrix

import Foreign.Storable

import Numeric.Container
import Numeric.LinearAlgebra.Algorithms
import Numeric.LinearAlgebra.Util

import Prelude hiding (max)

type Point = Vector Double
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


main = do
  print $ mvee 1e-5 $ fromList <$> [[1,0],[0,1],[-1,0],[0,-1]]
