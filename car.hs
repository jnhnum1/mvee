{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies #-}

import Control.Lens

import Data.Packed.Vector

class Vectorizable a s | a -> s where
  toVector :: a -> Vector s
  fromVector :: Vector s -> a

data Loc = Loc {_x :: Double, _y :: Double} deriving (Show)
makeLenses ''Loc

data Car = Car {_front :: Loc, _back :: Loc, _angle :: Double} deriving (Show)
makeLenses ''Car

instance Vectorizable Car Double where
  toVector (Car f b t) = fromList [(x f), (y f), (x b), (y b), t]
  fromVector v = 
    let
      [xf, yf, xb, yb, t] = toList v
    in Car (Loc xf yf) (Loc xb yb) t

stepCar :: Double -> Double -> Car -> Car
stepCar v dt (Car f b theta) = 
  let
    yl = (y f) - (y b)
    xl = (x f) - (x b)
    h = sqrt $ yl * yl + xl * xl
    coa = xl / h
    sia = yl / h
    dy = v * ((sin theta) * coa + (cos theta) * sia)
    dx = v * ((cos theta) * coa - (sin theta) * sia)

main = print $ (Car (Loc 1 2) (Loc 3 4) 5)
