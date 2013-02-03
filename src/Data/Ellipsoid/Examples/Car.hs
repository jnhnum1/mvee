{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies, TemplateHaskell  #-}
module Data.Ellipsoid.Examples.Car (
  Vectorizable(..),
  Loc,
  mkLoc,
  Car,
  mkCar,
  stepCar) where

import Control.Lens
import Control.Lens.TH

import Data.Packed.Vector

-- sigh, most recent version of lens package isn't on hackage
infixl 1 &
(&) :: a -> (a -> b) -> b
a & f = f a

class Vectorizable a s | a -> s where
  toVector :: a -> Vector s
  fromVector :: Vector s -> a

data Loc = Loc {_x :: Double, _y :: Double} deriving (Show)
makeLenses ''Loc

mkLoc x y = Loc x y

data Car = Car {_front :: Loc, _back :: Loc, _angle :: Double} deriving (Show)
makeLenses ''Car

mkCar f b t = Car f b t

instance Vectorizable Car Double where
  toVector (Car f b t) = fromList [f^.x, f^.y, b^.x, b^.y, t]
  fromVector v = 
    let
      [xf, yf, xb, yb, t] = toList v
    in Car (Loc xf yf) (Loc xb yb) t

stepCar :: Double -> Car -> Car
stepCar v car = 
  let
    yl = car^.front.y - car^.back.y
    xl = car^.front.x - car^.back.x
    theta = car^.angle
    h = sqrt $ yl * yl + xl * xl
    coa = xl / h
    sia = yl / h
    dy = v * ((sin theta) * coa + (cos theta) * sia)
    dx = v * ((cos theta) * coa - (sin theta) * sia)
    tt = (dx + xl) * coa + (dy + yl) * sia
    q = dx * coa + xl * coa + dy * sia + yl * sia -
      0.5 * sqrt (4*tt^2 - 4*(dx^2 + 2*dx*xl + dy^2 + 2*dy*yl))
  in car & front.x +~ dx & front.y +~ dy & back.x +~ q * coa 
         & back.y +~ q * sia

-- main = print $ toVector `map` iterate (stepCar 5) (Car (Loc 25 7.5) (Loc 15 7.5) 0.5)
