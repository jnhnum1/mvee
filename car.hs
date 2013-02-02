{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies #-}

import Data.Packed.Vector

class Vectorizable a s | a -> s where
  toVector :: a -> Vector s
  fromVector :: Vector s -> a

data Loc = Loc {x :: Double, y :: Double} deriving (Show)
data Car = Car {front :: Loc, back :: Loc, angle :: Double} deriving (Show)

instance Vectorizable Car Double where
  toVector (Car f b t) = fromList [(x f), (y f), (x b), (y b), t]
  fromVector v = 
    let
      [xf, yf, xb, yb, t] = toList v
    in Car (Loc xf yf) (Loc xb yb) t

main = print $ (Car (Loc 1 2) (Loc 3 4) 5)
