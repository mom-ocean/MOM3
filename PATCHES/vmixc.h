!====================== include file "vmixc.h" =========================
!
!         vertical mixing coefficients and related variables
!
!     kappa_h = constant vertical diffusion coefficient (cm**2/sec)
!     kappa_m = constant vertical viscosity coefficient (cm**2/sec)
!
!     visc_cbu  = viscosity coeff at bottom of U cell (cm**2/sec) 
!     diff_cbt  = diffusion coeff at bottom of T cell (cm**2/sec)
!     visc_cbu_limit = largest allowable "visc_cbu" 
!     diff_cbt_limit = largest allowable "diff_cbt" 
!     aidif = coefficient for implicit time differencing for
!             vertical diffusion. aidif=1 gives the fully implicit
!             case. aidif=0 gives the fully explicit case
!             note: not used unless "implicitvmix" or "isoneutralmix"
!                   is enabled
!
!     author:   r.c.pacanowski     e-mail  rcp@gfdl.gov
!=======================================================================
!
!
!     set parameter nsmix to allow for double diffusion
!
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
!
#if defined bryan_lewis_vertical
      common /vmixr0/ Ahv(km)
#endif
!
#ifdef ppvmix
!
!     variables for pacanowski-philander vertical diffusion
!
!     fricmx = max vertical mixing coefficient
!     wndmix = min vertical mixing in level 1 to simulate wind mixing
!     diff_cbt_back = background "diff_cbt"
!     visc_cbu_back = background "visc_cbu"
!
      common /vmixr0/ wndmix, fricmx, diff_cbt_back, visc_cbu_back
#endif
!
#ifdef kppvmix
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!     Define various parameters and common blocks for kmix vertical-
!     mixing scheme; used in "kmixs.F" subroutines
!-----------------------------------------------------------------------
!
!
!     jwtype    = Jerlov water type        (1 to 5)
!     ghats     = nonlocal transport           (s/cm^2)
!     hbl       = boundary layer depth (cm)
!     bfsfc     = surface buoyancy forcing     (cm^2/s^3)
!     ws        = scalar velocity scale
!     wm        = momentum velocity scale
!     caseA     = = 1 in case A; =0 in case B
!     stable    = = 1 in stable forcing; =0 in unstable
!     dkm1      = boundary layer difs at kbl-1 level
!     blmc      = boundary layer mixing coefficients
!     sigma     = normalized depth (d / hbl)
!     rhosfc    = potential density of sfc layer (g/cm^3)
!     talpha    = d(rho)/ d(pot.temperature)   (g/cm^3/C)
!     sbeta     = d(rho)/ d(salinity)        (g/cm^3/PSU)
!     alphaDT   = alpha * DT  across interfaces  (g/cm^3)
!     betaDS    = beta  * DS  across interfaces  (g/cm^3)
!     ustar     = surface friction velocity        (cm/s)
!     Bo        = surface turb buoy. forcing   (cm^2/s^3)
!     Bosol     = radiative buoy forcing       (cm^2/s^3)
!     dbloc     = local delta buoy at interfaces (cm/s^2)
!     dbsfc     = delta buoy w/ respect to sfc   (cm/s^2)
!     kbl       = index of first grid level below hbl
!
!-----------------------------------------------------------------------
!
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
!
!-----------------------------------------------------------------------
!     parameters for several subroutines
!
!     epsilon = nondimensional extent of the surface layer = 0.1
!     vonk    = von Karman`s constant                      = 0.4
!     conc1,conam,concm,conc2,zetam,conas,concs,conc3,zetas
!             = scalar coefficients
!     lri     = logical switch for shear instability mixing
!     ldd     = logical switch for double-diffusive mixing
!-----------------------------------------------------------------------
!
      logical lri,ldd
      common /kmixcom/ epsilon, vonk, conc1, conam, concm, conc2
     &,                zetam, conas, concs, conc3, zetas, lri, ldd
!
!-----------------------------------------------------------------------
!     parameters for subroutine "bldepth"
!
!
!     to compute depth of boundary layer:
!
!     Ricr    = critical bulk Richardson Number            = 0.3
!     cekman  = coefficient for ekman depth                = 0.7 
!     cmonob  = coefficient for Monin-Obukhov depth        = 1.0
!     concv   = ratio of interior buoyancy frequency to 
!               buoyancy frequency at entrainment depth    = 1.8
!     hbf     = fraction of bounadry layer depth to 
!               which absorbed solar radiation 
!               contributes to surface buoyancy forcing    = 1.0
!     Vtc     = non-dimensional coefficient for velocity
!               scale of turbulant velocity shear
!               (=function of concv,concs,epsilon,vonk,Ricr)
!-----------------------------------------------------------------------
!
      common /kmixcbd/ Ricr, cekman, cmonob, concv, hbf, Vtc
!
!-----------------------------------------------------------------------
!     parameters and common arrays for subroutines "kmixinit" 
!     and "wscale"
!
!
!     to compute turbulent velocity scales:
!
!     nni     = number of values for zehat in the look up table
!     nnj     = number of values for ustar in the look up table
!
!     wmt     = lookup table for wm, the turbulent velocity scale 
!               for momentum
!     wst     = lookup table for ws, the turbulent velocity scale 
!               for scalars
!     deltaz  = delta zehat in table
!     deltau  = delta ustar in table
!     zmin    = minimum limit for zehat in table (cm3/s3)
!     zmax    = maximum limit for zehat in table
!     umin    = minimum limit for ustar in table (cm/s)
!     umax    = maximum limit for ustar in table
!-----------------------------------------------------------------------
!
      parameter (nni = 890, nnj = 480)
!
      common /kmixcws/ wmt(0:nni+1,0:nnj+1), wst(0:nni+1,0:nnj+1),
     &                 deltaz, deltau, zmin, zmax, umin, umax
!
!-----------------------------------------------------------------------
!     parameters for subroutine "ri_iwmix"
!
!
!     to compute vertical mixing coefficients below boundary layer:
!
!     Riinfty = local Richardson Number limit 
!               for shear instability                      = 0.7
!     visc_cbu_limit = viscosity max due to shear instability     
!     diff_cbt_limit = diffusivity ..  
!     visc_cbu_iw    = viscosity background due to internal waves  
!     diff_cbt_iw    = diffusivity ..     
!     visc_con_limit = viscosity due to convective instability
!     diff_con_limit = diffusivity ..
!     visc_cbt       = viscosity coeff at bottom of "T" cell
!-----------------------------------------------------------------------
!
      common /kmixcri/ Riinfty
     &,                 visc_cbu_iw, diff_cbt_iw
     &,                 visc_con_limit, diff_con_limit
!
!-----------------------------------------------------------------------
!     parameters for subroutine "ddmix"
!
!
!     to compute additional diffusivity due to double diffusion:
!
!     Rrho0   = limit for double diffusive density ratio
!     dsfmax  = maximum diffusivity in case of salt fingering (m2/s)
!-----------------------------------------------------------------------
!
      common /kmixcdd/ Rrho0, dsfmax
!
!-----------------------------------------------------------------------
!     parameters for subroutine "blmix"
!
!
!     to compute mixing within boundary layer:
!
!     cstar   = proportionality coefficient for nonlocal transport
!     cg      = non-dimensional coefficient for counter-gradient term
!-----------------------------------------------------------------------
!
      common /kmixcbm/ cstar, cg
#endif
