{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE PatternGuards #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}

{-# OPTIONS_GHC -Wall #-}
module Numeric.LevMar (
  Params
, Samples
, Model
, Jacobian
, Info        (..)
, StopReason  (..)
, Options     (..)
, Constraints (..)
, defaultOpts
, levmar_p1m1_approx
, levmar_approx
) where

import Numeric.LinearAlgebra
import Data.Vector.Storable (Vector)
import qualified Data.Vector.Storable as VS

import Control.Applicative
import Control.Arrow
import Control.Monad
import Data.Data
import Data.Maybe
import Data.Monoid (Monoid(..))

import Debug.Trace

type Params r = Data.Vector.Storable.Vector r
type Samples r = Data.Vector.Storable.Vector r

{-| A functional relation describing measurements represented as a function
from a vector of parameters to a vector of expected samples.

 * Ensure that the length @m@ of the parameter vector equals the length of the
   initial parameter vector in 'levmar'.

 * Ensure that the length @n@ of the output sample vector equals the length of
   the sample vector in 'levmar'.

 * Ensure that the length @n@ of the output sample vector vector is bigger than or
   equal to the length @m@ of the parameter vector.
-}
type Model r = Params r -> Samples r

{-| The <http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant jacobian>
of the 'Model' function. Expressed as a function from a vector of
parameters to a matrix which for each expected sample describes the
partial derivatives of the parameters.

 * Ensure that the length @m@ of the parameter vector equals the length of the
   initial parameter vector in 'levmar'.

 * Ensure that the output matrix has the dimension @n><m@ where @n@ is the
   number of samples and @m@ is the number of parameters.
-}
type Jacobian r = Params r -> Matrix r

-- | Minimization options
data Options r =
    Opts { optScaleInitMu      :: !r -- ^ Scale factor for initial /mu/.
         , optStopNormInfJacTe :: !r -- ^ Stopping thresholds for @||J^T e||_inf@.
         , optStopNorm2Dp      :: !r -- ^ Stopping thresholds for @||Dp||_2@.
         , optStopNorm2E       :: !r -- ^ Stopping thresholds for @||e||_2@.
         , optDelta            :: !r -- ^ Step used in the difference
                                    -- approximation to the Jacobian. If
                                    -- @optDelta<0@, the Jacobian is approximated
                                    -- with central differences which are more
                                    -- accurate (but slower!)  compared to the
                                    -- forward differences employed by default.
         } deriving (Eq, Ord, Read, Show, Data, Typeable)

-- | Default minimization options
defaultOpts :: Fractional r => Options r
defaultOpts = Opts { optScaleInitMu      = 1e-3
                   , optStopNormInfJacTe = 1e-17
                   , optStopNorm2Dp      = 1e-17
                   , optStopNorm2E       = 1e-17
                   , optDelta            = 1e-6
                   }

-- | Ensure that these vectors have the same length as the number of parameters.
data Constraints r = Constraints
    { lowerBounds       :: !(Maybe (Params r))            -- ^ Optional lower bounds
    , upperBounds       :: !(Maybe (Params r))            -- ^ Optional upper bounds
    , weights           :: !(Maybe (Params r))            -- ^ Optional weights
    , linearConstraints :: !(Maybe (LinearConstraints r)) -- ^ Optional linear constraints
    } deriving (Read, Show, Typeable)

-- calculate the appropriate projection function into the boxed constraints
projectBox :: (VS.Storable r, Ord r) => Constraints r -> Params r -> Params r
projectBox Constraints{..} = ubProj . lbProj
  where
    lbProj = maybe id (VS.zipWith max) lowerBounds
    ubProj = maybe id (VS.zipWith min) upperBounds

-- | Linear constraints consisting of a constraints matrix, @k><m@ and
--   a right hand constraints vector, of length @k@ where @m@ is the number of
--   parameters and @k@ is the number of constraints.
type LinearConstraints r = (Matrix r, Vector r)

-- | * 'mempty' is defined as a 'Constraints' where all fields are 'Nothing'.
--
--   * 'mappend' merges two 'Constraints' by taking the first non-'Nothing' value
--     for each field.
instance Monoid (Constraints r) where
    mempty = Constraints Nothing Nothing Nothing Nothing
    mappend (Constraints lb1 ub1 w1 l1)
            (Constraints lb2 ub2 w2 l2) = Constraints (lb1 `mplus` lb2)
                                                      (ub1 `mplus` ub2)
                                                      (w1  `mplus` w2)
                                                      (l1  `mplus` l2)


-- | Information regarding the minimization.
data Info r = Info
  { infNorm2initE      :: !r          -- ^ @||e||_2@             at initial parameters.
  , infNorm2E          :: !r          -- ^ @||e||_2@             at estimated parameters.
  , infNormInfJacTe    :: !r          -- ^ @||J^T e||_inf@       at estimated parameters.
  , infNorm2Dp         :: !r          -- ^ @||Dp||_2@            at estimated parameters.
  , infMuDivMax        :: !r          -- ^ @\mu/max[J^T J]_ii ]@ at estimated parameters.
  , infNumIter         :: !Int        -- ^ Number of iterations.
  , infStopReason      :: !StopReason -- ^ Reason for terminating.
  , infNumFuncEvals    :: !Int        -- ^ Number of function evaluations.
  , infNumJacobEvals   :: !Int        -- ^ Number of jacobian evaluations.
  , infNumLinSysSolved :: !Int        -- ^ Number of linear systems solved,
                                     --   i.e. attempts for reducing error.
  } deriving (Eq, Ord, Read, Show, Data, Typeable)

-- | Reason for terminating.
data StopReason
  = SmallGradient  -- ^ Stopped because of small gradient @J^T e@.
  | SmallDp        -- ^ Stopped because of small Dp.
  | MaxIterations  -- ^ Stopped because maximum iterations was reached.
  | SingularMatrix -- ^ Stopped because of singular matrix. Restart from current
                   --   estimated parameters with increased 'optScaleInitMu'.
  | SmallestError  -- ^ Stopped because no further error reduction is
                   --   possible. Restart with increased 'optScaleInitMu'.
  | SmallNorm2E    -- ^ Stopped because of small @||e||_2@.
  | InvalidValues  -- ^ Stopped because model function returned invalid values
                   --   (i.e. NaN or Inf). This is a user error.
    deriving (Eq, Ord, Read, Show, Data, Typeable, Enum)


data LevMarError
    = LevMarError                    -- ^ Generic error (not one of the others)
    | LapackError                    -- ^ A call to a lapack subroutine failed
                                     --   in the underlying C levmar library.
    | FailedBoxCheck                 -- ^ At least one lower bound exceeds the
                                     --   upper one.
    | MemoryAllocationFailure        -- ^ A call to @malloc@ failed in the
                                     --   underlying C levmar library.
    | ConstraintMatrixRowsGtCols     -- ^ The matrix of constraints cannot have
                                     --   more rows than columns.
    | ConstraintMatrixNotFullRowRank -- ^ Constraints matrix is not of full row
                                     --   rank.
    | TooFewMeasurements             -- ^ Cannot solve a problem with fewer
                                     --   measurements than unknowns.  In case
                                     --   linear constraints are provided, this
                                     --   error is also returned when the number
                                     --   of measurements is smaller than the
                                     --   number of unknowns minus the number of
                                     --   equality constraints.
      deriving (Eq, Ord, Read, Show, Data, Typeable)

incrStep :: Info k -> Info k
incrStep info = info { infNumIter = infNumIter info +1 }

incrFunc :: Info k -> Info k
incrFunc info = info {infNumFuncEvals = infNumFuncEvals info + 1 }

incrLin :: Info k -> Info k
incrLin info = info {infNumLinSysSolved = infNumLinSysSolved info + 1 }

incrJacEval :: Info k -> Info k
incrJacEval info = info {infNumJacobEvals = infNumJacobEvals info + 1 }

setMuDivMax :: k -> Info k -> Info k
setMuDivMax x info = info{infMuDivMax = x }

levmar_p1m1_approx :: (Fractional n, Ord n) => (n->n) -> n -> n -> Int -> Options n -> (n, Info n)
levmar_p1m1_approx model p0 sample iters Opts{..} = outerloop info0 stop0 v0
  where
    j (pM,yM) (pN,yN) = (yN-yM) / (pN-pM)
    p1 = p0*(1+1e-4)
    y0 = model p0
    y1 = model p1  
    j0 = j (p0,y0) (p1,y1)
    a0 = j0*j0
    ep0 = sample-y0
    g0 = j0*ep0
    u0 = optScaleInitMu*a0
    v0 = (p0,y0,a0,ep0,g0,u0,2)
    stop0 = SmallGradient <$ guard (abs g0 <= optStopNormInfJacTe)
    info0 = Info
        { infNorm2initE = g0*g0
        , infNorm2E = g0
        , infNormInfJacTe = 0
        , infNorm2Dp      = abs (optDelta)
        , infMuDivMax     = u0/a0
        , infNumIter      = 0
        , infStopReason   = InvalidValues
        , infNumFuncEvals = 2  -- for y0,y1, which we need for g0
        , infNumJacobEvals = 0
        , infNumLinSysSolved = 0
        }

    outerloop info@Info{infNumIter} doStop v@(!pLast,!_yLast,!_a',!_ep,!_g,!_mu,!_uStep)
      | Just infStopReason <- doStop = (pLast, info{infStopReason})
      | infNumIter >= iters = (pLast, info{infStopReason=MaxIterations})
      | otherwise  = case innerstep info v of
                        (vNext,info', stop') -> outerloop info' stop' vNext

    innerstep info v@(!pLast, !yLast, !a', !ep, !g, !mu, !uStep) = case (abs dp <= dpThreshold,rho > 0) of
        (True, _)     -> (v, info{infNorm2Dp = dp}, Just SmallDp)
        -- reject this dp
        (False,False) -> ((pLast,yLast,a',ep,g,mu*uStep,2*uStep),incrLin $ incrFunc info,Nothing)
        (False, True) ->
            let !jNew  = j (pLast,yLast) (pNew,yNew)
                !aNew  = jNew*jNew       -- hessian
                !gNew  = jNew*epNew
                rhoTmp = 2*rho-1
                !uNew  = mu*max (1/3) (1- rhoTmp*rhoTmp*rhoTmp)
                !doStop = msum
                          [ SmallNorm2E   <$ guard (norm2Ep <= optStopNorm2E)
                          , SmallGradient <$ guard (abs gNew <= optStopNormInfJacTe)
                          ]
                !vNext = (pNew,yNew,aNew,epNew,gNew,uNew,2)
                !info' = incrLin . incrStep . incrFunc $ info
                          { infNorm2E = norm2Ep
                          , infNormInfJacTe = abs gNew
                          , infNorm2Dp = norm2dp
                          , infMuDivMax = mu/aNew
                          }
            in (vNext, info', doStop)
      where
          dpThreshold = optStopNorm2Dp * pLast*pLast
          dp    = g / (a'+mu)
          pNew  = pLast+dp
          yNew  = model pNew
          epNew = sample-yNew
          norm2Ep = epNew*epNew
          norm2dp = dp*dp
          rho = 2 * (infNorm2E info - norm2Ep) / (dp*(mu*dp+g)) - 1
{-# INLINEABLE levmar_p1m1_approx #-}

-- TODO: replace with super-efficient implementation
-- | numerically estimate the transposed Jacobian using forward differentiation
jacobianTForward :: (Element r, Fractional r) => (Model r) -> Params r -> Samples r -> r -> Matrix r
jacobianTForward f x !y0 delta = fromRows $ map dY [0..VS.length x-1]
  where
    dY i = let dx = x VS.// [(i, xi+dxi)]
               xi  = x VS.! i
               dxi = xi*delta
               dxiInv = 1/dxi
           in VS.zipWith (\dy y -> (dy-y)*dxiInv) (f dx) y0

broydens_approx
    :: (Element r, Fractional r, Product r, Container Vector r, Num (Vector r))
    => Matrix r -> Samples r -> Samples r -> Params r -> Matrix r
broydens_approx jLast yLast yNew dp =
    jLast + ((asColumn $ scale (1/(dot dp dp)) (yNew - yLast - (jLast `mXv` dp)))
              `mXm` asRow dp)

levmar_approx
    :: (Show n, Container Vector n, RealOf n ~ n, Product n, Ord n, Num (Vector n), Field n)
    => (Samples n-> Params n)
    -> Maybe (Jacobian n)
    -> Samples n
    -> Params n
    -> Int
    -> Options n
    -> Constraints n
    -> (Samples n, Info n)
levmar_approx model jac'm p0 sample iters opts constr@Constraints{..} =
    case linearConstraints of
      Nothing -> maybe (levmar_approx_boxed model) (levmar_boxed model) jac'm
                   p0 sample iters opts (projectBox constr)
      Just (constraintA,constraintB) ->
        first reconstrain $ maybe (levmar_approx_boxed model')
                                  (levmar_boxed model' . reconstrainJ) jac'm
                                  p0' sample iters opts project
        -- for linear equality constraints, transform the system into a
        -- related, unconstrained system via variable elimination, then call the solver on that.
        --
        -- any solution to Ax = b (linear constraints) can be expressed as
        -- x = Q1 R^-T P^T b + Q2 x2
        --
        -- so we search for an unconstrained x2.
        --
        -- For LEC with box constraints, combining these two is equivalent to
        -- solving a system of linear equations at each iteration.  levmar (c)
        -- instead implements penalty terms for values outside the box, but we're
        -- going to do the linear solve.
        where
            project x = case msum [lowerBounds, upperBounds] of
                          Nothing -> x
                          Just _  -> let constrX = reconstrain x
                                         constrX' = projectBox constr constrX
                                     in if constrX == constrX'
                                        then x
                                        else case toColumns $ linearSolveLS q2 (asColumn $ constrX' - constFactor) of
                                            [x'2] -> x'2
                                            _     -> error "levmar_approx: internal error in projection"
            m = rows constraintA
            aT = trans constraintA
            (_,_,perm,_) = lu constraintA
            (q,r') = qr (aT `mXm` perm)
            q1 = takeColumns m q
            q2 = dropColumns m q
            r = takeRows m r'
            invTrans = inv . trans
            constFactor = q1 `mXv` (invTrans r `mXv` (trans perm `mXv` constraintB))
            -- get the constrained x for a given x2
            reconstrain x = constFactor + (q2 `mXv` x)
            reconstrainJ j = (`mXm` q2) . j . reconstrain
            model' = model . reconstrain
            p0' = project p0Unproj

            -- for the initial x2, find the solutions to the
            -- linear system Q2 x2 = x - constFactor
            [p0Unproj] = toColumns $ linearSolveLS q2 (asColumn $ p0 - constFactor)
{-# INLINEABLE levmar_approx #-}

-- main levmar function with approximate jacobian.  Supports box-constrained
-- systems only.
levmar_approx_boxed
    :: forall n. (Show n, Container Vector n, RealOf n ~ n, Product n, Ord n, Num (Vector n), Field n)
    => (Samples n-> Params n)
    -> Samples n
    -> Params n
    -> Int
    -> Options n
    -> (Samples n -> Samples n) -- projection function for box constraints
    -> (Samples n, Info n)
levmar_approx_boxed model p0 sample iters Opts{..} proj = outerloop info0 stop0 v0
  where
    y0 = model p0
    j0 = jacobianTForward model p0 y0 1e-4
    a0 = calcHess j0
    ep0 = sample - y0
    g0 = calcG j0 ep0
    innersize = rows a0
    calcHess j = j `mXm` trans j
    calcG j ep = j `mXv` ep

    u0 = optScaleInitMu * VS.maximum (takeDiag a0)
    v0 = (p0,y0,j0,a0,ep0,g0,u0,2,False)
    stop0 :: Maybe StopReason
    stop0 = SmallGradient <$ guard (normInf g0 <= optStopNormInfJacTe)
    info0 :: Info n
    norm2' v = let v' = norm2 v in v'*v'
    info0 = Info
        { infNorm2initE   = norm2' ep0
        , infNorm2E       = norm2' ep0
        , infNormInfJacTe = normInf g0
        , infNorm2Dp      = 0
        , infMuDivMax     = u0 / normInf (takeDiag a0)
        , infNumIter      = 0
        , infStopReason   = InvalidValues
        , infNumFuncEvals = 1+VS.length p0  -- for initial jacobian
        , infNumJacobEvals = 1
        , infNumLinSysSolved = 0
        }
    outerloop info@Info{infNumIter} doStop v@(!pLast,!_yLast,!_jLast,!_a',!_ep,!_g,!_mu,!_uStep,_lastAccepted)
      | Just infStopReason <- doStop = (pLast, info{infStopReason})
      | infNumIter >= iters = (pLast, info{infStopReason=MaxIterations})
      | otherwise  = case innerstep info v of
                        (vNext,info', stop') -> outerloop info' stop' vNext

    innerstep info v@(!pLast, !yLast, !jLast, !a', !ep, !g, !mu, !uStep, !lastAccepted) =
      case (norm2dp <= dpThreshold,rho > 0) of
        (True, _)     -> (v, info{infNorm2Dp = normInf dp}, Just SmallDp)

        -- reject this dp
        (False,False) ->
            let (!jNext,!aNext,!gNext,!nextAccepted) =
                  case (fullJRecalc, lastAccepted || uphill) of
                      -- after accepting a dP, update the jacobian estimate
                      -- and related terms
                      -- also update the estimate if we appear to be going
                      -- uphill
                      (True, _) -> (j, calcHess j, calcG j ep, False)
                          where
                              j = jacobianTForward model pLast yLast
                                          (1e-4 * VS.minimum (VS.map abs dp))
                      (False, True) -> (j, calcHess j, calcG j ep, lastAccepted)
                          where
                              j = trans $ broydens_approx (trans jLast) yLast yNew dp
                      (False, False) -> (jLast, a', g, False)
            in ( (pLast,yLast,jNext,aNext,ep,gNext, mu*uStep,2*uStep,nextAccepted)
               , incrJac . setMuDivMax muDivMax . incrLin . incrStep $ incrFunc info
               , Nothing)

        -- accept this dp
        (False, True) ->
            let !jNew  = if fullJRecalc
                            then jacobianTForward model pNew yNew (1e-4 * VS.minimum (VS.map abs dp))
                            else trans $ broydens_approx (trans jLast) yLast yNew dp
                !aNew  = calcHess jNew
                !gNew  = calcG jNew epNew
                normInfGNew = normInf gNew
                rhoTmp = 2*rho-1
                !uNew  = mu*max (1/3) (1 - rhoTmp*rhoTmp*rhoTmp)
                !doStop = msum
                          [ SmallNorm2E   <$ guard (norm2Ep <= optStopNorm2E)
                          , SmallGradient <$ guard (normInfGNew <= optStopNormInfJacTe)
                          ]
                !vNext = (pNew,yNew,jNew,aNew,epNew,gNew,uNew,2,not fullJRecalc)
                !info' = incrJac . setMuDivMax muDivMax . incrLin . incrStep . incrFunc $ info
                          { infNorm2E = norm2Ep
                          , infNormInfJacTe = normInfGNew
                          , infNorm2Dp = norm2dp
                          }
            in (vNext, info', doStop)
      where
          dpThreshold = optStopNorm2Dp * norm2' pLast
          -- we use (A+mu*I) for the left factor
          -- the Marquardt step apparently is supposed to scale by the max
          -- value of A or something, but no example code seems to implement
          -- that...
          dp    = head . toColumns $ linearSolve
                      (a'+diagRect 0 (VS.replicate innersize mu) innersize innersize)
                      (asColumn g)
          pNew  = proj $ pLast+dp
          yNew  = model pNew
          epNew = sample-yNew
          norm2Ep = norm2' epNew :: RealOf n
          norm2dp = norm2' dp
          rho = 2 * (infNorm2E info -  norm2Ep) / (dot dp (scale mu dp + g)) - 1
          muDivMax = mu*uStep / VS.maximum (VS.map abs $ takeDiag a')

          fullJRecalc = (1+infNumIter info) `rem` (max 10 (cols j0)) == 0
          uphill = infNorm2E info - norm2Ep > 0
          -- increment the info if we need to calculate a full Jacobian
          -- estimation
          incrJac = if fullJRecalc then (!! VS.length p0) . iterate incrFunc . incrJacEval else id
{-# INLINEABLE levmar_approx_boxed #-}

-- main levmar function with definite jacobian.  Supports box-constrained
-- systems only.
levmar_boxed
    :: forall n. (Show n, Container Vector n, RealOf n ~ n, Product n, Ord n, Num (Vector n), Field n)
    => (Samples n -> Params n)
    -> Jacobian n
    -> Samples n
    -> Params n
    -> Int
    -> Options n
    -> (Samples n -> Samples n) -- projection function for box constraints
    -> (Samples n, Info n)
levmar_boxed model jac p0 sample iters Opts{..} proj = outerloop info0 stop0 v0
  where
    y0 = model p0
    j0 = jac p0
    jT0 = trans j0
    a0 = calcHess jT0 j0
    ep0 = sample - y0
    g0 = calcG jT0 ep0
    innersize = rows a0
    calcHess jT j = jT `mXm` j
    calcG jT ep = jT `mXv` ep

    u0 = optScaleInitMu * VS.maximum (takeDiag a0)
    v0 = (p0,a0,g0,u0,2)
    stop0 :: Maybe StopReason
    stop0 = SmallGradient <$ guard (normInf g0 <= optStopNormInfJacTe)
    info0 :: Info n
    norm2' v = let v' = norm2 v in v'*v'
    info0 = Info
        { infNorm2initE   = norm2' ep0
        , infNorm2E       = norm2' ep0
        , infNormInfJacTe = normInf g0
        , infNorm2Dp      = 0
        , infMuDivMax     = u0 / normInf (takeDiag a0)
        , infNumIter      = 0
        , infStopReason   = InvalidValues
        , infNumFuncEvals = 1+VS.length p0  -- for initial jacobian
        , infNumJacobEvals = 1
        , infNumLinSysSolved = 0
        }
    outerloop info@Info{infNumIter} doStop v@(!pLast,!_a',!_g,!_mu,!_uStep)
      | Just infStopReason <- doStop = (pLast, info{infStopReason})
      | infNumIter >= iters = (pLast, info{infStopReason=MaxIterations})
      | otherwise  = case innerstep info v of
                        (vNext,info', stop') -> outerloop info' stop' vNext

    innerstep info v@(!pLast, !a', !g, !mu, !uStep) =
      case (norm2dp <= dpThreshold,rho > 0) of
        (True, _)     -> (v, info{infNorm2Dp = normInf dp}, Just SmallDp)

        -- reject this dp
        (False,False) ->
            ( (pLast,a',g, mu*uStep,2*uStep)
            , setMuDivMax muDivMax . incrLin . incrStep $ incrFunc info
            , Nothing)

        -- accept this dp
        (False, True) ->
            let !jNew  = jac pNew
                !jT    = trans jNew
                !aNew  = calcHess jT jNew
                !gNew  = calcG jT epNew
                normInfGNew = normInf gNew
                rhoTmp = 2*rho-1
                !uNew  = mu*max (1/3) (1 - rhoTmp*rhoTmp*rhoTmp)
                !doStop = msum
                          [ SmallNorm2E   <$ guard (norm2Ep <= optStopNorm2E)
                          , SmallGradient <$ guard (normInfGNew <= optStopNormInfJacTe)
                          ]
                !vNext = (pNew,aNew,gNew,uNew,2)
                !info' = incrJacEval . setMuDivMax muDivMax . incrLin . incrStep . incrFunc $ info
                          { infNorm2E = norm2Ep
                          , infNormInfJacTe = normInfGNew
                          , infNorm2Dp = norm2dp
                          }
            in (vNext, info', doStop)
      where
          dpThreshold = optStopNorm2Dp * norm2' pLast
          -- we use (A+mu*I) for the left factor
          -- the Marquardt step apparently is supposed to scale by the max
          -- value of A or something, but no example code seems to implement
          -- that...
          dp    = head . toColumns $ linearSolve
                      (a'+diagRect 0 (VS.replicate innersize mu) innersize innersize)
                      (asColumn g)
          pNew  = proj $ pLast+dp
          yNew  = model pNew
          epNew = sample-yNew
          norm2Ep = norm2' epNew :: RealOf n
          norm2dp = norm2' dp
          rho = 2 * (infNorm2E info -  norm2Ep) / (dot dp (scale mu dp + g)) - 1
          muDivMax = mu*uStep / VS.maximum (VS.map abs $ takeDiag a')
{-# INLINEABLE levmar_boxed #-}

y0x [j,k] = j*j+2*k
y1x [j,_k] = 1/j
y v = let xs = VS.toList v in VS.fromList [y0x xs, y1x xs]
jy' [j,k] = fromLists [[2*j, 2],[-1/(j*j),0]]
jy = jy' . VS.toList
constrB = fromList [9::Double]
constrA = (><)  1 2 [4, -1 :: Double]
constr = mempty { linearConstraints = Just (constrA,constrB) } :: Constraints Double
