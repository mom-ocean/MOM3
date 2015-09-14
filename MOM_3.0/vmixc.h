c====================== include file "vmixc.h" =========================
c
c         vertical mixing coefficients and related variables
c
c     kappa_h = constant vertical diffusion coefficient (cm**2/sec)
c     kappa_m = constant vertical viscosity coefficient (cm**2/sec)
c
c     visc_cbu  = viscosity coeff at bottom of U cell (cm**2/sec) 
c     diff_cbt  = diffusion coeff at bottom of T cell (cm**2/sec)
c     visc_cbu_limit = largest allowable "visc_cbu" 
c     diff_cbt_limit = largest allowable "diff_cbt" 
c     aidif = coefficient for implicit time differencing for
c             vertical diffusion. aidif=1 gives the fully implicit
c             case. aidif=0 gives the fully explicit case
c             note: not used unless "implicitvmix" or "isoneutralmix"
c                   is enabled
c
c     author:   r.c.pacanowski     e-mail  rcp@gfdl.gov
c=======================================================================
c
c
c     set parameter nsmix to allow for double diffusion
c
#ifdef kppvmix
      parameter ( nsmix=2 ) 
#else
      parameter ( nsmix=1 )
#endif
      real kappa_h,  kappa_m
      common /vmixr0/ visc_cbu_limit, diff_cbt_limit, aidif
      common /vmixr0/ kappa_h, kappa_m
      common /vmixr1/ visc_cbu(imt,km,2:jmw-1)
      common /vmixr1/ diff_cbt(imt,km,2:jmw-1,nsmix)
c
#if defined bryan_lewis_vertical
      common /vmixr0/ Ahv(km)
#endif
c
#ifdef ppvmix
c
c     variables for pacanowski-philander vertical diffusion
c
c     fricmx = max vertical mixing coefficient
c     wndmix = min vertical mixing in level 1 to simulate wind mixing
c     diff_cbt_back = background "diff_cbt"
c     visc_cbu_back = background "visc_cbu"
c
      common /vmixr0/ wndmix, fricmx, diff_cbt_back, visc_cbu_back
#endif
c
#ifdef kppvmix
c
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c     Define various parameters and common blocks for kmix vertical-
c     mixing scheme; used in "kmixs.F" subroutines
c-----------------------------------------------------------------------
c
c
c     jwtype    = Jerlov water type        (1 to 5)
c     ghats     = nonlocal transport           (s/cm^2)
c     hbl       = boundary layer depth (cm)
c     bfsfc     = surface buoyancy forcing     (cm^2/s^3)
c     ws        = scalar velocity scale
c     wm        = momentum velocity scale
c     caseA     = = 1 in case A; =0 in case B
c     stable    = = 1 in stable forcing; =0 in unstable
c     dkm1      = boundary layer difs at kbl-1 level
c     blmc      = boundary layer mixing coefficients
c     sigma     = normalized depth (d / hbl)
c     rhosfc    = potential density of sfc layer (g/cm^3)
c     talpha    = d(rho)/ d(pot.temperature)   (g/cm^3/C)
c     sbeta     = d(rho)/ d(salinity)        (g/cm^3/PSU)
c     alphaDT   = alpha * DT  across interfaces  (g/cm^3)
c     betaDS    = beta  * DS  across interfaces  (g/cm^3)
c     ustar     = surface friction velocity        (cm/s)
c     Bo        = surface turb buoy. forcing   (cm^2/s^3)
c     Bosol     = radiative buoy forcing       (cm^2/s^3)
c     dbloc     = local delta buoy at interfaces (cm/s^2)
c     dbsfc     = delta buoy w/ respect to sfc   (cm/s^2)
c     kbl       = index of first grid level below hbl
c
c-----------------------------------------------------------------------
c
      common /vmixr3/ ws(imt,2:jmw-1), wm(imt,2:jmw-1)
      common /vmixr3/ Bo(imt,2:jmw-1), Bosol(imt,2:jmw-1)
      common /vmixr3/ jwtype(imt,2:jmw-1)
      common /vmixr3/ dbloc(imt,km,2:jmw-1), dbsfc(imt,km,2:jmw-1)
      common /vmixr3/ bfsfc(imt,2:jmw-1), ustar(imt,2:jmw-1)
      common /vmixr3/ caseA(imt,2:jmw-1), stable(imt,2:jmw-1)
      common /vmixr3/ sigma(imt,2:jmw-1), rhosfc(imt,2:jmw-1)
      common /vmixr3/ talpha(imt,km,2:jmw-1), sbeta(imt,km,2:jmw-1)
      common /vmixr3/ alphaDT(imt,km,2:jmw-1),betaDS(imt,km,2:jmw-1)
      common /vmixr3/ dkm1(imt,2:jmw-1,3), blmc(imt,km,2:jmw-1,3)
      common /vmixr3/ dVsq(imt,km,jmw), Rib(imt,2,2:jmw-1)
      common /vmixr3/ kbl(imt,2:jmw-1), hbl(imt,2:jmw-1)
      common /vmixr3/ ghats(imt,km,2:jmw-1)
      common /vmixr3/ gat1(imt,2:jmw-1,3), dat1(imt,2:jmw-1,3)
c
c-----------------------------------------------------------------------
c     parameters for several subroutines
c
c     epsilon = nondimensional extent of the surface layer = 0.1
c     vonk    = von Karman`s constant                      = 0.4
c     conc1,conam,concm,conc2,zetam,conas,concs,conc3,zetas
c             = scalar coefficients
c     lri     = logical switch for shear instability mixing
c     ldd     = logical switch for double-diffusive mixing
c-----------------------------------------------------------------------
c
      logical lri,ldd
      common /kmixcom/ epsilon, vonk, conc1, conam, concm, conc2
     &,                zetam, conas, concs, conc3, zetas, lri, ldd
c
c-----------------------------------------------------------------------
c     parameters for subroutine "bldepth"
c
c
c     to compute depth of boundary layer:
c
c     Ricr    = critical bulk Richardson Number            = 0.3
c     cekman  = coefficient for ekman depth                = 0.7 
c     cmonob  = coefficient for Monin-Obukhov depth        = 1.0
c     concv   = ratio of interior buoyancy frequency to 
c               buoyancy frequency at entrainment depth    = 1.8
c     hbf     = fraction of bounadry layer depth to 
c               which absorbed solar radiation 
c               contributes to surface buoyancy forcing    = 1.0
c     Vtc     = non-dimensional coefficient for velocity
c               scale of turbulant velocity shear
c               (=function of concv,concs,epsilon,vonk,Ricr)
c-----------------------------------------------------------------------
c
      common /kmixcbd/ Ricr, cekman, cmonob, concv, hbf, Vtc
c
c-----------------------------------------------------------------------
c     parameters and common arrays for subroutines "kmixinit" 
c     and "wscale"
c
c
c     to compute turbulent velocity scales:
c
c     nni     = number of values for zehat in the look up table
c     nnj     = number of values for ustar in the look up table
c
c     wmt     = lookup table for wm, the turbulent velocity scale 
c               for momentum
c     wst     = lookup table for ws, the turbulent velocity scale 
c               for scalars
c     deltaz  = delta zehat in table
c     deltau  = delta ustar in table
c     zmin    = minimum limit for zehat in table (cm3/s3)
c     zmax    = maximum limit for zehat in table
c     umin    = minimum limit for ustar in table (cm/s)
c     umax    = maximum limit for ustar in table
c-----------------------------------------------------------------------
c
      parameter (nni = 890, nnj = 480)
c
      common /kmixcws/ wmt(0:nni+1,0:nnj+1), wst(0:nni+1,0:nnj+1),
     &                 deltaz, deltau, zmin, zmax, umin, umax
c
c-----------------------------------------------------------------------
c     parameters for subroutine "ri_iwmix"
c
c
c     to compute vertical mixing coefficients below boundary layer:
c
c     Riinfty = local Richardson Number limit 
c               for shear instability                      = 0.7
c     visc_cbu_limit = viscosity max due to shear instability     
c     diff_cbt_limit = diffusivity ..  
c     visc_cbu_iw    = viscosity background due to internal waves  
c     diff_cbt_iw    = diffusivity ..     
c     visc_con_limit = viscosity due to convective instability
c     diff_con_limit = diffusivity ..
c     visc_cbt       = viscosity coeff at bottom of "T" cell
c-----------------------------------------------------------------------
c
      common /kmixcri/ Riinfty
     &,                 visc_cbu_iw, diff_cbt_iw
     &,                 visc_con_limit, diff_con_limit
c
c-----------------------------------------------------------------------
c     parameters for subroutine "ddmix"
c
c
c     to compute additional diffusivity due to double diffusion:
c
c     Rrho0   = limit for double diffusive density ratio
c     dsfmax  = maximum diffusivity in case of salt fingering (m2/s)
c-----------------------------------------------------------------------
c
      common /kmixcdd/ Rrho0, dsfmax
c
c-----------------------------------------------------------------------
c     parameters for subroutine "blmix"
c
c
c     to compute mixing within boundary layer:
c
c     cstar   = proportionality coefficient for nonlocal transport
c     cg      = non-dimensional coefficient for counter-gradient term
c-----------------------------------------------------------------------
c
      common /kmixcbm/ cstar, cg
#endif
