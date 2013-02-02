{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies, TemplateHaskell  #-}

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

data Car = Car {_front :: Loc, _back :: Loc, _angle :: Double} deriving (Show)
makeLenses ''Car

instance Vectorizable Car Double where
  toVector (Car f b t) = fromList [(_x f), (_y f), (_x b), (_y b), t]
  fromVector v = 
    let
      [xf, yf, xb, yb, t] = toList v
    in Car (Loc xf yf) (Loc xb yb) t

stepCar :: Double -> Double -> Car -> Car
stepCar v dt car = 
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

main = print $ (Car (Loc 1 2) (Loc 3 4) 5)
