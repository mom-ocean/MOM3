c====================== include file "emode.h" =========================
c
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
c
c     surface pressure/free surface time centering parameters:
c
c     alph,gam = parameters to center the coriolis and surface pressure 
c                gradients in time (leapfrog)
c
c     theta    =parameters to center the coriolis and surface pressure 
c                gradients in time (mixing ts )
c
c     apgr     = is = to alpha/theta leapfrog/mixing ts
c
# ifdef implicit_free_surface
c
c     to set free surface time centering parameters...
c     note: for proper time centering of coriolis
c     and pressure gradients alph, gam should
c     satisfy gam = 1 -2*alph.  for stability
c     alph should be > 1/4.  recommended values
c     are alph = gam = 1/3.
# endif
      common /emode/  alph, gam, theta, apgr
#endif
c
c     mxscan  = max number of allowable scans for poisson solvers
c     mscan   = actual number of scans taken by poisson solvers
c     tolrsf  = tolerance for stream function calculation. 
c               the solution is halted when it is within "tolrsf"
c               of the "true" solution assuming geometric convergence.
c     tolrsp  = tolerance for surface pressure calculation
c               the solution is halted when it is within "tolrsp"
c               of the "true" solution assuming geometric convergence.
c     tolrfs  = tolerance for implicit free surface calculation
c               the solution is halted when it is within "tolrfs"
c               of the "true" solution assuming geometric convergence.
c     esterr  = estimated maximum error in elliptic solver assuming
c               geometric convergence
c    
#ifndef explicit_free_surface      
      common /emode/  mxscan, mscan, tolrsf, tolrsp, tolrfs
      common /emode/  esterr
#endif
c
c     nisle = number of land masses
c     nippts= number of land mass perimeter points
c     iperm = "i" coordinate for the land mass perimeter point
c     jperm = "j" coordinate for the land mass perimeter point
c     iofs  = offset for indexing into the land mass perimeter points
c     imask = controls whether calculations get done on perimeters
c     set mask for land mass perimeters on which to perform calculations
c     imask(-n) = .false.  [no equations ever on dry land mass n]
c     imask(0)  = .true.   [equations at all mid ocean points]
c     imask(n)  = .true./.false [controls whether there will be
c                                equations on the ocean perimeter of
c                                land mass n]
c     note: land mass 1 is the northwest-most land mass
c     for the numbering of the other landmasses, see generated map(i,j)
      logical imask
      common /emodei/ nippts(mnisle), iofs(mnisle), iperm(maxipp)
      common /emodei/ jperm(maxipp), nisle, imain
      common /emodel/ imask (-mnisle:mnisle)
#ifdef zero_island_flow
c
c     make net transport between "land_mass_a" and "land_mass_b" = zero.
c     land_mass_a = any land mass number from array map(i,jrow)
c     land_mass_b = any other land mass number from array map(i,jrow)
c     ipermt      = temporary array needed for re-organizing perimeters
c
      common /emodei/ land_mass_a, land_mass_b
      common /emodei/ ipermt(maxipp,2)
#endif
c
      character*16 variable
      logical converged
      common /emodel/ converged
      common /emodec/ variable
c
