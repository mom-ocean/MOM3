c====================== include file "scalar.h" ========================
c
c     various scalar quantities:
c
c     dtts   = time step for density & tracers (in seconds)
c     dtuv   = time step for baroclinic velocity (in seconds)
c     dtsf   = time step for barotropic velocity (in seconds)
c     c2dtts = 2*dtts
c     c2dtuv = 2*dtuv
c     c2dtsf = 2*dtsf
c     smooth = smoothing coefficient for robert time filter
c     acor   = (>0, 0) = (implicit, explicit) treatment of coriolis
c               term for internal and external modes.
c     rho0   = mean density for Boussinesq approximation
c     rho0r  = 1/rho0 
c     omega  = earth`s rotation rate (radians/sec)
c     radius = earth`s radius (cm)
c     grav   = earth`s gravitational acceleration (cm/sec**2)
c     cdbot  = bottom drag coefficient
c     ncon   = number of  passes through convective code in tracer
c     gcor   = time centering for coriolis term
c
c     taux0  = constant zonal windstress (dynes/cm**2) for idealized
c              equatorial studies
c     tauy0  = constant meridional windstress (dynes/cm**2) for
c              idealized equatorial studies
c
      common /scalar/ dtts, dtuv, dtsf, c2dtts, c2dtuv, c2dtsf, acor
      common /scalar/ smooth
      common /scalar/ rho0, rho0r, omega, radius, grav, cdbot, gcor
      common /scalar/ taux0, tauy0
      common /scalri/ ncon
c
c     bring in non dimensional constants
c
#include "ndcon.h"
c
