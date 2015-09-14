!====================== include file "diag.h" ==========================
!
!
!     variables used for computing diagnostics:
!
!     tcella  = T cell surface area cm**2 (entire ocean)
!     ucella  = U cell surface area cm**2 (entire ocean)
!     tcellv  = T cell volume cm**3 (entire ocean)
!     ucellv  = U cell volume cm**3 (entire ocean)
!     tcellv_j= T cell volume cm**3 (latitude row j)
!     ucellv_j= U cell volume cm**3 (latitude row j)
!
      common /cdiag/ tcellv, ucellv
      common /cdiag/ tcella(km), ucella(km)
      common /cdiag/ tcellv_j(jmt), ucellv_j(jmt)

#ifdef time_step_monitor
!
!     ektot    = "total" kinetic energy per unit volume at "tau". units
!                ergs/cm**3 = dyn/cm**2 = g/cm/sec**2 = 10**-7 J/cm**3.
!                ektot is the "total" ke in the sense that it considers
!                both the internal and external modes summed over the
!                entire ocean volume. The contributions of
!                vertical motions are neglected on the basis of scaling
!                arguments (i.e., w**2 << u**2 + v**2).
!     dtabs    = absolute value of rate of change of tracer per unit
!                volume centered at "tau"
!     tbar     = first moment of tracer at time "tau"
!     travar   = variance = second moment of tracer about mean at "tau"
!     tratot   = total mass of tracer at time "tau"
!
      common /cdiag/ ektot(0:km,jmt)
      common /cdiag/ dtabs(0:km,nt,jmt)
      common /cdiag/ travar(0:km,nt,jmt)
      common /cdiag/ tbar(0:km,nt,jmt)
      common /cdiag/ tratot(nt)
# if defined explicit_free_surface 
!     etaana(2) = height excess of the total model area
!     etaana(3) = kinetic energy in the barotropic mode
!     etaana(4) = potential energy in the barotropic mode
!     (analysis according Killworth, Stainforth, Webb and Paterson)
!
!     etaana(1) = total volume rate input from fresh water (Sv)
!     volume_tot= total accumulated volume from sff (cm3)
!     eta_mix   = viscosity(cm2/sec) for eta_laplacian
!     neta      = index of etat which is filtered with eta_laplacian
! 
      common /emoder/ etaana(4), volume_tot
      common /emoder/ volume_tot_pe, volume_tot_leap_pe
      common /emoder/ flux_sffs, flux_taum1_sffs
!     
# endif
#endif
#ifdef energy_analysis
!
!     engint   = volume averaged internal mode energy integral
!                components
!     engext   = volume averaged external mode energy integral
!                components
!     buoy     = volume averaged buoyancy
!     plicin   = internal mode imbalance
!     plicex   = external mode imbalance
!     buoerr   = buoyancy error
!     enleak   = non-linear error
!
!     tcerr    = maximum "t" cell continuity error
!     ucerr    = maximum "u" cell continuity error
!     itcerr   = "i" index corresponding to "tcerr"
!     jtcerr   = "jrow" index corresponding to "tcerr"
!     ktcerr   = "k" index corresponding to "tcerr"
!     iucerr   = "i" index corresponding to "ucerr"
!     jucerr   = "jrow" index corresponding to "ucerr"
!     kucerr   = "k" index corresponding to "ucerr"
!
!     wtbot    = maximum "adv_vbt" error at ocean bottom
!     iwtbot   = "i" index corresponding to "wtbot"
!     jwtbot   = "jrow" index corresponding to "wtbot"
!     kwtbot   = "k" index corresponding to "wtbot"
!     wubot    = maximum "adv_vbu" at ocean bottom
!     iwubot   = "i" index corresponding to "wubot"
!     jwubot   = "jrow" index corresponding to "wubot"
!     kwubot   = "k" index corresponding to "wubot"
!
!     wtlev    = zonally integrated adv_vbt for each level
!     wulev    = zonally integrated adv_vbu for each level
!
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
!
!     ttn      = northward transport of tracer components
!
!     ttn2     = northward transport of tracers for ocean basins
!                  (.,.,.,0)       Global
!                  (.,.,.,1:nhreg) Ocean basins
!                also,
!                  (6,.,.,.) total transport due to advection
!                            (including gm if -Dgm_advect) 
!                  (7,.,.,.) total transport due to diffusion 
!                            (including skew flux if -Dgm_skew)
!                  (8,.,.,.) total transport
      common /gyres/ ttn(8,jmt,ntmin2)
      common /gyres/ ttn2(6:8,jmt,nt,0:nhreg)
#endif
#ifdef meridional_overturning
!
!     vmsf     = overturning streamfunction with depth as vert coord
!     sigAmsf  = overturning streamfunction with sigAvert as vert coord
!     sigBmsf  = overturning streamfunction with sigBvert as vert coord
!
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
!
!     vmsf     = vertical_meridional stream function by basin
!              basin #0  = global domain (all basins together)
!              other basins may be setup by the researcher
!              (i.e. basin 1 could be the atlantic, 2->the pacific, etc)
!              (default is nbasin = 1 which treats the entire ocean
!                    as basin #1. (redundant with basin #0))         
!
!     nbasin = number of basins for calculation of meridional transport
!              basin #0  = global domain (all basins together)
!              other basins may be setup by the researcher
!              (i.e. basin 1 could be the atlantic, 2->the pacific, etc)
!              (default is nbasin = 1 which treats the entire ocean
!                    as basin #1. (redundant with basin #0))         
!
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
!
!     zmsmf    = zonal mean surface momentum flux
!     zmstf    = zonal mean surface tracer flux
!     zmsm     = zonal mean surface momentum
!     zmst     = zonal mean surface tracers
!     zmsff    = zonal mean surface fresh water flux (Sv)
!     zmau     = surface area weighting for "u" latitudes
!     zmat     = surface area weighting for "t" latitudes
!
      common /cdiag/ zmsmf(jmt,2), zmstf(jmt,nt)
      common /cdiag/ zmau(jmt), zmat(jmt)
      common /cdiag/ zmsm(jmt,2), zmst(jmt,nt)
      common /cdiag/ zmsff(jmt)
#endif
#ifdef tracer_yz
!
!     tyz(,,,1)  = zonal mean of tracer T
!     tyz(,,,2)  = zonal mean of d(T)/dt
!     tyz(,,,3)  = zonal mean of advection of T
!     tyz(,,,4)  = zonal mean of diffusion of T
!     tyz(,,,5)  = zonal mean of source of T
!
      common /cdiag/ tyz(jmt,km,nt,5)
#endif
#ifdef save_convection
!
!     convect0  = temperature before explicit convection
!     convect1  = time rate of change of temperature due to 
!                 explicit convection. set to epsln over land
!                 for use in identifying land cells
!
      common /exconv/ excnv0(imt,km,2:jmw-1)
      common /exconv/ excnv1(imt,km,2:jmw-1)
#endif
#ifdef save_mixing_coeff
!
!     ce = ceoff on east face of cell (1 is for momentum, 2 for tracers)
!     cn = ceoff on north face of cell (1 is for momentum,2 for tracers)
!     cb = ceoff on bottom face of cell(1 is for momentum,2 for tracers)
!
      common /mixcoef/ ce(imt,km,2:jmw-1,2)
      common /mixcoef/ cn(imt,km,2:jmw-1,2)
      common /mixcoef/ cb(imt,km,2:jmw-1,2)
#endif
#ifdef tracer_averages
!
!     common for tracer averages within horizontal regions
!
!
!     sumbt = volume weighted tracer sum for a given region
!     avgbt = volume weighted tracer average for a given region
!     sumbk = volume weighted tracer sum for a region at a k level
!     avgbk = volume weighted tracer average for a region at a k level
!     sumgt = total global volume weighted tracer sum
!     avggt = total global volume weighted tracer average
!     sumgk = total global volume weighted tracer sum at a k level
!     avggk = total global volume weighted tracer average at a k level
!     sumbf = area weighted tracer flux sum for a given region
!     avgbf = average area weighted tracer flux for a region
!     sumgf = global area weighted tracer flux sum
!     avggf = global average area weighted tracer flux
!
!
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
      
!
!                      meridional tracer balance
!
!     to calculate the meridional tracer balance, each term in the
!     tracer equation is integrated in depth,longitude and time for each
!     latitude. the time integration is over "tmbint"(see switch.h) days
!     for the calculation to make sense, the longitudinal
!     integration must go around the world (in the case of an open
!     latitude {eg: 60 deg S}) or between land masses (in the case of
!     a closed latitude {eg: 30 deg S from south america to africa})
!     the domain may be divided into arbitrary basins (atlantic, indian
!     pacific) but for certain latitudes, results from basins may have
!     to be combined to satisfy the above conditions. 
!
!
!     ntmbb    = number of basins over which the tracer meridional
!                balance will be calculated. the test case assumes 1
!                basin as defined in "setocn.F"
!     tstor    = average latitudinal storage of tracer
!     tdiv     = average latitudinal divergence of tracer
!     tflux    = average latitudinal surface flux of tracer
!     tdif     = average latitudinal diffusion of tracer
!     tsorc    = average latitudinal additional source of tracer
!     smdvol   = ocean volume of basin at latitude
!     numtmb   = number of time steps over which the terms have been
!                accumulated
      parameter (ntmbb=1)
      common /tmbr/ tstor(jmt,nt,0:ntmbb), tdiv(jmt,nt,0:ntmbb)
      common /tmbr/ tflux(jmt,nt,0:ntmbb), tdif(jmt,nt,0:ntmbb)
      common /tmbr/ tsorc(jmt,nt,0:ntmbb), smdvol(jmt,0:ntmbb)
      common /tmbi/ numtmb
#endif
!-----------------------------------------------------------------------
!       for opening output unints using mpp_open. 
!
!       MPP_ASCII: formatted read/write
!       MPP_NETCDF saves 32 bit netCDF file 
!       MPP_IEEE32 saves 32 bit IEEE file 
!       MPP_NATIVE saves data in same form as computed
!-----------------------------------------------------------------------
#if defined use_netCDF
        integer, parameter :: OUTPUT_FORM = MPP_NETCDF
#else
        integer, parameter :: OUTPUT_FORM = MPP_IEEE32
#endif




