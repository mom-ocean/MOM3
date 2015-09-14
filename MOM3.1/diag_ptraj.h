c====================== include file "diag_ptraj.h" ====================
c
c     inputs:
c
c     nptraj  = number of particle trajectories. particles are
c               assumed to be of zero mass, neutrally buoyant and move
c               with the local three dimensional flow during each time
c               step.
c
c               initially, particles are randomly distributed within
c               the volume given by:
c
c     ptslon  = starting longitude for initial particle distribution
c     ptelon  = ending longitude for initial particle distribution
c     ptslat  = starting latitude for initial particle distribution
c     ptelat  = ending latitude for initial particle distribution
c     ptsdpt  = starting depth for initial particle distribution
c     ptedpt  = ending depth for initial particle distribution
c
c     outputs:
c
c     pxyz    = particle coordinates. index (1,2,3) is for particle 
c               (longitude, latitude, depth).
c
c     pijk    = the particle is bounded by the volume with vertices
c               given by the eight nearest surrounding model grid points
c               on the "xu","yu", and "z" grids. index (1,2,3) locates
c               the (longitude, latitude, depth) index of the deepest 
c               northeast corner of this bounding volume.
c
c     em      = matrix of deformation rates for calculation of
c               lyapunov exponents
c
c     ptdone  = boolean for limiting multiple passes on trajectories
c               within one time step
c
c     initpt  = boolean for initializing particle positions
c               (t,f) = (initialize, do not initialize)
c     
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c             
      parameter (nptraj = 200)
      integer pijk
      logical ptdone, initpt
c
      common /cptrji/ pijk(3,nptraj)
      common /cptraj/ pxyz(3,nptraj)
      common /cptraj/ ptslon, ptelon, ptslat, ptelat, ptsdpt, ptedpt
#ifdef lyapunov
      common /cptraj/ em(2,2,nptraj)
#endif
      common /cptra2/ ptdone(nptraj), initpt
