c====================== include file "diag.h" ==========================
c
c
c     variables used for computing diagnostics:
c
c     tcella  = T cell surface area cm**2 (entire ocean)
c     ucella  = U cell surface area cm**2 (entire ocean)
c     tcellv  = T cell volume cm**3 (entire ocean)
c     ucellv  = U cell volume cm**3 (entire ocean)
c     tcellv_j= T cell volume cm**3 (latitude row j)
c     ucellv_j= U cell volume cm**3 (latitude row j)
c
      common /cdiag/ tcellv, ucellv
      common /cdiag/ tcella(km), ucella(km)
      common /cdiag/ tcellv_j(jmt), ucellv_j(jmt)

#ifdef time_step_monitor
c
c     ektot    = "total" kinetic energy per unit volume at "tau". units
c                ergs/cm**3 = dyn/cm**2 = g/cm/sec**2 = 10**-7 J/cm**3.
c                ektot is the "total" ke in the sense that it considers
c                both the internal and external modes summed over the
c                entire ocean volume. The contributions of
c                vertical motions are neglected on the basis of scaling
c                arguments (i.e., w**2 << u**2 + v**2).
c     dtabs    = absolute value of rate of change of tracer per unit
c                volume centered at "tau"
c     tbar     = first moment of tracer at time "tau"
c     travar   = variance = second moment of tracer about mean at "tau"
c     tratot   = total mass of tracer at time "tau"
c
      common /cdiag/ ektot(0:km,jmt)
      common /cdiag/ dtabs(0:km,nt,jmt)
      common /cdiag/ travar(0:km,nt,jmt)
      common /cdiag/ tbar(0:km,nt,jmt)
      common /cdiag/ tratot(nt)
# if defined explicit_free_surface 
c     etaana(2) = height excess of the total model area
c     etaana(3) = kinetic energy in the barotropic mode
c     etaana(4) = potential energy in the barotropic mode
c     (analysis according Killworth, Stainforth, Webb and Paterson)
c
c     etaana(1) = total volume rate input from fresh water (Sv)
c     volume_tot= total accumulated volume from sff (cm3)
c     eta_mix   = viscosity(cm2/sec) for eta_laplacian
c     neta      = index of etat which is filtered with eta_laplacian
c 
      common /emoder/ etaana(4), volume_tot
c     
# endif
#endif
#ifdef energy_analysis
c
c     engint   = volume averaged internal mode energy integral
c                components
c     engext   = volume averaged external mode energy integral
c                components
c     buoy     = volume averaged buoyancy
c     plicin   = internal mode imbalance
c     plicex   = external mode imbalance
c     buoerr   = buoyancy error
c     enleak   = non-linear error
c
c     tcerr    = maximum "t" cell continuity error
c     ucerr    = maximum "u" cell continuity error
c     itcerr   = "i" index corresponding to "tcerr"
c     jtcerr   = "jrow" index corresponding to "tcerr"
c     ktcerr   = "k" index corresponding to "tcerr"
c     iucerr   = "i" index corresponding to "ucerr"
c     jucerr   = "jrow" index corresponding to "ucerr"
c     kucerr   = "k" index corresponding to "ucerr"
c
c     wtbot    = maximum "adv_vbt" error at ocean bottom
c     iwtbot   = "i" index corresponding to "wtbot"
c     jwtbot   = "jrow" index corresponding to "wtbot"
c     kwtbot   = "k" index corresponding to "wtbot"
c     wubot    = maximum "adv_vbu" at ocean bottom
c     iwubot   = "i" index corresponding to "wubot"
c     jwubot   = "jrow" index corresponding to "wubot"
c     kwubot   = "k" index corresponding to "wubot"
c
c     wtlev    = zonally integrated adv_vbt for each level
c     wulev    = zonally integrated adv_vbu for each level
c
      common /cdiag/ buoy(0:km,jmt), engint(0:km,8,jmt), engext(8,jmt)
      common /cdiag/ tcerr(jmt), ucerr(jmt)
      common /cdiag/ plicin, plicex, buoerr, enleak  
      common /cdiag/ wtbot(jmt), wubot(jmt) 
      common /cdiag/ wtlev(km,0:jmt), wulev(km,0:jmt)
      common /cdiagi/ itcerr(jmt), jtcerr(jmt), ktcerr(jmt) 
      common /cdiagi/ iucerr(jmt), jucerr(jmt), kucerr(jmt) 
      common /cdiagi/ iwtbot(jmt), jwtbot(jmt), kwtbot(jmt)
      common /cdiagi/ iwubot(jmt), jwubot(jmt), kwubot(jmt)
#endif
#ifdef term_balances
# include "diag_termbal.h"
#endif
#ifdef gyre_components
c
c     ttn      = northward transport of tracer components
c
c     ttn2     = northward transport of tracers for ocean basins
c                  (.,.,.,0)       Global
c                  (.,.,.,1:nhreg) Ocean basins
c                also,
c                  (6,.,.,.) total transport due to advection
c                            (including gm if -Dgm_advect) 
c                  (7,.,.,.) total transport due to diffusion 
c                            (including skew flux if -Dgm_skew)
c                  (8,.,.,.) total transport
      common /gyres/ ttn(8,jmt,ntmin2)
      common /gyres/ ttn2(6:8,jmt,nt,0:nhreg)
#endif
#ifdef meridional_overturning
c
c     vmsf     = overturning streamfunction with depth as vert coord
c     sigAmsf  = overturning streamfunction with sigAvert as vert coord
c     sigBmsf  = overturning streamfunction with sigBvert as vert coord
c
      common /cdiag/ vmsf(jmt,0:km)
# ifdef meridional_overturning_density
      common /cdiag/ sigAmsf(jmt,km),sigAvert(km)
      common /cdiag/ sigBmsf(jmt,km),sigBvert(km)
# endif
# ifdef gent_mcwilliams
      common /cdiag/ vmsf_gm(jmt,0:km)
# endif
#endif

#ifdef meridional_overturning_basin
c
c     vmsf     = vertical_meridional stream function by basin
c              basin #0  = global domain (all basins together)
c              other basins may be setup by the researcher
c              (i.e. basin 1 could be the atlantic, 2->the pacific, etc)
c              (default is nbasin = 1 which treats the entire ocean
c                    as basin #1. (redundant with basin #0))         
c
c     nbasin = number of basins for calculation of meridional transport
c              basin #0  = global domain (all basins together)
c              other basins may be setup by the researcher
c              (i.e. basin 1 could be the atlantic, 2->the pacific, etc)
c              (default is nbasin = 1 which treats the entire ocean
c                    as basin #1. (redundant with basin #0))         
c
#if !defined test_case_A
      parameter (nbasin=1)
#else
      parameter (nbasin=2)
#endif
      common /cbasin/ basinmask(imt,jmt,0:nbasin)
      common /cdiag/ vmsf(jmt,0:km,0:nbasin)
# ifdef gent_mcwilliams
      common /cdiag/ vmsf_gm(jmt,0:km,0:nbasin)
# endif
#endif

#ifdef show_zonal_mean_of_sbc
c
c     zmsmf    = zonal mean surface momentum flux
c     zmstf    = zonal mean surface tracer flux
c     zmsm     = zonal mean surface momentum
c     zmst     = zonal mean surface tracers
c     zmsff    = zonal mean surface fresh water flux (Sv)
c     zmau     = surface area weighting for "u" latitudes
c     zmat     = surface area weighting for "t" latitudes
c
      common /cdiag/ zmsmf(jmt,2), zmstf(jmt,nt)
      common /cdiag/ zmau(jmt), zmat(jmt)
      common /cdiag/ zmsm(jmt,2), zmst(jmt,nt)
      common /cdiag/ zmsff(jmt)
#endif
#ifdef tracer_yz
c
c     tyz(,,,1)  = zonal mean of tracer T
c     tyz(,,,2)  = zonal mean of d(T)/dt
c     tyz(,,,3)  = zonal mean of advection of T
c     tyz(,,,4)  = zonal mean of diffusion of T
c     tyz(,,,5)  = zonal mean of source of T
c
      common /cdiag/ tyz(jmt,km,nt,5)
#endif
#ifdef save_convection
c
c     convect0  = temperature before explicit convection
c     convect1  = time rate of change of temperature due to 
c                 explicit convection. set to epsln over land
c                 for use in identifying land cells
c
      common /exconv/ excnv0(imt,km,2:jmw-1)
      common /exconv/ excnv1(imt,km,2:jmw-1)
#endif
#ifdef save_mixing_coeff
c
c     ce = ceoff on east face of cell (1 is for momentum, 2 for tracers)
c     cn = ceoff on north face of cell (1 is for momentum,2 for tracers)
c     cb = ceoff on bottom face of cell(1 is for momentum,2 for tracers)
c
      common /mixcoef/ ce(imt,km,2:jmw-1,2)
      common /mixcoef/ cn(imt,km,2:jmw-1,2)
      common /mixcoef/ cb(imt,km,2:jmw-1,2)
#endif
#ifdef tracer_averages
c
c     common for tracer averages within horizontal regions
c
c
c     sumbt = volume weighted tracer sum for a given region
c     avgbt = volume weighted tracer average for a given region
c     sumbk = volume weighted tracer sum for a region at a k level
c     avgbk = volume weighted tracer average for a region at a k level
c     sumgt = total global volume weighted tracer sum
c     avggt = total global volume weighted tracer average
c     sumgk = total global volume weighted tracer sum at a k level
c     avggk = total global volume weighted tracer average at a k level
c     sumbf = area weighted tracer flux sum for a given region
c     avgbf = average area weighted tracer flux for a region
c     sumgf = global area weighted tracer flux sum
c     avggf = global average area weighted tracer flux
c
c
      common /tavgr/ sumbk(nhreg,km,nt), sumbt(nhreg,nt), sumgk(km,nt)
      common /tavgr/ sumgt(nt), sumbf(nhreg,nt), sumgf(nt)
      common /tavgr/ avgbk(nhreg,km,nt), avgbt(nhreg,nt), avggk(km,nt)
      common /tavgr/ avggt(nt), avgbf(nhreg,nt), avggf(nt)
#endif
#ifdef save_fct
      common /fctdiag/ fct_minus_ctr(imt,km,2:jmw-1)
      common /fctdiag/ flux_diffX(imt,km,2:jmw-1)
      common /fctdiag/ flux_diffY(imt,km,2:jmw-1)
      common /fctdiag/ flux_diffZ(imt,km,2:jmw-1)
#endif
#if defined meridional_tracer_budget
      
c
c                      meridional tracer balance
c
c     to calculate the meridional tracer balance, each term in the
c     tracer equation is integrated in depth,longitude and time for each
c     latitude. the time integration is over "tmbint"(see switch.h) days
c     for the calculation to make sense, the longitudinal
c     integration must go around the world (in the case of an open
c     latitude {eg: 60 deg S}) or between land masses (in the case of
c     a closed latitude {eg: 30 deg S from south america to africa})
c     the domain may be divided into arbitrary basins (atlantic, indian
c     pacific) but for certain latitudes, results from basins may have
c     to be combined to satisfy the above conditions. 
c
c
c     ntmbb    = number of basins over which the tracer meridional
c                balance will be calculated. the test case assumes 1
c                basin as defined in "setocn.F"
c     tstor    = average latitudinal storage of tracer
c     tdiv     = average latitudinal divergence of tracer
c     tflux    = average latitudinal surface flux of tracer
c     tdif     = average latitudinal diffusion of tracer
c     tsorc    = average latitudinal additional source of tracer
c     smdvol   = ocean volume of basin at latitude
c     numtmb   = number of time steps over which the terms have been
c                accumulated
      parameter (ntmbb=1)
      common /tmbr/ tstor(jmt,nt,0:ntmbb), tdiv(jmt,nt,0:ntmbb)
      common /tmbr/ tflux(jmt,nt,0:ntmbb), tdif(jmt,nt,0:ntmbb)
      common /tmbr/ tsorc(jmt,nt,0:ntmbb), smdvol(jmt,0:ntmbb)
      common /tmbi/ numtmb
#endif
c-----------------------------------------------------------------------
c       for opening output unints using mpp_open. 
c
c       MPP_ASCII: formatted read/write
c       MPP_NETCDF saves 32 bit netCDF file 
c       MPP_IEEE32 saves 32 bit IEEE file 
c       MPP_NATIVE saves data in same form as computed
c-----------------------------------------------------------------------
#if defined use_netCDF
        integer, parameter :: OUTPUT_FORM = MPP_NETCDF
#else
        integer, parameter :: OUTPUT_FORM = MPP_IEEE32
#endif




