      module datamod
!
!=======================================================================
!     ramdisk3d = ramdrive area for 3-D data
!     ramdisk2d = ramdrive area for 2-D data
!
!     author: R.C.Pacanowski   e-mail -> rcp@gfdl.gov
!=======================================================================
!
      real, dimension(:,:,:), allocatable :: ramdisk3d
      real, dimension(:,:), allocatable :: ramdisk2d
      integer :: pe, npes
      end module




      module memory_window
!
!=======================================================================
!     storage needed to solve the 3-D baroclinic and tracer equations
!
!     author: R.C.Pacanowski   e-mail -> rcp@gfdl.gov
!=======================================================================
!
#include "derived_options.h"
#if defined new_memory_window
!     taum1 = tau-1 time level for variables in MW (memory window)
!     tau   = tau   time level for variables in MW
!     taup1 = tau+1 time level for variables in MW
!
      integer taum1, tau, taup1
!
!     prognostic variables u = velocity, t = tracer
!
!     u(i,k,j,n,tau) = total velocity in cm/sec where:
!      i   = index for longitude
!      k   = index for depth
!      j   = index for latitude row within MW
!      n   = component (1 = zonal, 2 = meridional)
!      tau = (-1,0,1) = time level (tau-1, tau, tau+1)
!      (only internal modes are kept on disk and at tau+1 in the MW)
!
!     t(i,k,j,n,tau) = tracer where:
!      i   = index for longitude
!      k   = index for depth
!      j   = index for latitude row within MW
!      n   = component (1 = potential temperature deg C, 2 = salinity.
!            units = deviation from 0.035 grams of salt/cm**3 of water)
!            if nt > 2 then other tracers are allowed.
!      tau = (-1,0,1) = time level (tau-1, tau, tau+1)
!
      real, dimension(:,:,:,:,:), allocatable :: u, t
!
!     Naming convention: adv = advective, diff = diffusive, f = flux
!                        v = velocity, n = northern face, t = T-cell
!                        b = bottom face, e = eastern face, u = U-cell
!
      real, dimension(:,:,:), allocatable :: adv_vet, adv_vnt, adv_vbt
      real, dimension(:,:,:), allocatable :: adv_veu, adv_vnu, adv_vbu
      real, dimension(:,:,:), allocatable :: adv_fe, adv_fn, adv_fb
      real, dimension(:,:,:), allocatable :: diff_fe, diff_fn, diff_fb
      real, dimension(:,:,:,:), allocatable :: grad_p
!
!     rho = density
!     tmask = (1,0) on (ocean,land) T-cell
!     umask = (1,0) on (ocean,land) U-cell
!
      real, dimension(:,:,:), allocatable :: rho, tmask, umask
!
!     smf = surface momentum flux, bmf = bottom momentum flux
!     stf = surface tracer flux,   btf = bottom tracer flux
!     units = dyn/cm**2, cal/cm**2/s, and (g salt)/cm**2/s
!     smf_turb = turbulent part of surface momentum flux
!     stf_turb = turbulent part of the surface tracer flux
!     
      real, dimension(:,:,:), allocatable :: smf,  bmf, stf, btf
      real, dimension(:,:,:), allocatable :: smf_turb,  stf_turb
!
!     naming convention: cst = cosine on grid point within T-cell
!                        csu = cosine on grid point within U-cell
!                        dxt = delta x for a T-cell on equator
!                        dxu = delta x for a U-cell on equator
!                        r   = reciprocal
!                        ah  = horizontal diffusion coeff
!                        am  = horizontal viscosity coeff
!                        met = metric
!
      real, dimension(:,:), allocatable :: cstdxtr,  cstdxur
      real, dimension(:,:), allocatable :: cstdxt2r, ah_cstdxur
      real, dimension(:,:), allocatable :: csudxtr,  csudxur
      real, dimension(:,:), allocatable :: csudxu2r, am_csudxtr
      real, dimension(:,:), allocatable :: advmet
#ifdef pressure_gradient_average
      real, dimension(:,:,:), allocatable :: rhotaum1, rhotaup1
      real, dimension(:,:,:), allocatable :: rhotilde
#endif
!
#ifdef partial_cell
!     dh = delta thickness, dhwt = vertical distance between T points 
      real, dimension(:,:,:), allocatable :: dht, dhu, dhwt
      real, dimension(:,:,:), allocatable :: pc_sink
#endif
!
#if defined fourth_order_tracer_advection || defined quicker
!     adv_f4n = 4th order advective flux
      real, dimension(:,:,:,:), allocatable :: adv_f4n
#endif
!
#ifdef velocity_horz_biharmonic
!     delsq = del**2 of horizontal velocity 
      real, dimension(:,:,:,:), allocatable :: del2_vel
#endif
#ifdef tracer_horz_biharmonic
!     delsq = del**2 of tracers
      real, dimension(:,:,:,:), allocatable :: del2_tracer
#endif
!
#if defined isoneutralmix
      real, dimension(:,:,:), allocatable :: diff_fbiso
#endif
!
#ifdef source_term
      real, dimension(:,:,:), allocatable :: source
#endif
!
#if defined implicitvmix || defined isoneutralmix
      real, dimension(:,:,:), allocatable :: zzi
#endif
!
#if defined linearized_advection || defined equatorial_thermocline
      real, dimension(:,:), allocatable :: tbarz
#endif
!
#ifdef fct
!
!     anti_fe  = antidiffusive flux across the eastern face of a T cell
!     anti_fn  = antidiffusive flux across the northern face of a T cell
!     anti_fb  = antidiffusive flux across the bottom face of a T cell
!     R_plusY  = ratio of maximal feasible to maximal possible change
!                of tracer T in subroutine tracer.F, N-S dim delimiter
!     R_minusY = ratio of minimal feasible to minimal possible change
!     t_lo     = low order solution for tracer
!     R_plus3 = ratio of maximal feasible to maximal possible change
!               of tracer T in subroutine tracer.F, 3-d delimiter
!     R_minus3= ratio of minimal feasible to minimal possible change
!               of tracer T in subroutine tracer.F, 3-d delimiter
!
      real, dimension(:,:,:,:), allocatable :: anti_fe, anti_fn, anti_fb
      real, dimension(:,:,:,:), allocatable :: R_plusY, R_minusY
# ifdef fct_3d
      real, dimension(:,:,:,:), allocatable :: R_plus3, R_minus3
# endif
# ifdef fct_dlm2
      real, dimension(:,:,:,:), allocatable :: t_lo
# endif
#endif
!
#else
      real dummy_real
#endif
      end module






      module arrays_2d
!
!=======================================================================
!     storage needed to solve 2-D barotropic equation plus some
!     additional 2-D arrays
!
!     author: R.C.Pacanowski   e-mail -> rcp@gfdl.gov
!=======================================================================
!
!     kmt = number of vertical boxes over T points
!     kmu = number of vertical boxes over U points
      integer, dimension(:,:), allocatable :: kmt, kmu
!
!     psi   = stream function (,,1) is for tau; (,,2) is for tau-1
!     zu    = vertically averaged forcing from momentum equations
!             (,,1) is zonal and (,,2) is meridional component
!     ztd   = curl of "zu" for the stream function equation
!     ptd   = time change of stream function
!     h     = depth over U points
!     hr    = reciprocal depth over U points
!     res   = residual from elliptic solver
!
!     map   = land mass map distinguishing, ocean, land, and perimiters
!
      integer, dimension(:,:), allocatable :: map
      real, dimension(:,:), allocatable :: res, hr, h
      real, dimension(:,:,:), allocatable :: zu
#if defined partial_cell
!
!     htp   = depth (cm) from surface to the bottom of partial T-cells
!     hup   = depth (cm) from surface to the bottom of partial U-cells
!
      real, dimension(:,:), allocatable :: htp, hup
#endif
#ifndef explicit_free_surface
      real, dimension(:,:), allocatable :: ptd
#endif
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
!
!     ps    = surface pressure (,,1) is for tau; (,,2) is for tau-1
!     divf  = barotropic divergence of uncorrected ubar & vbar (rhs of 
!             surface pressure eqn)
!     ubar  = barotropic velocity defined on U point for "tau" 
!             (,,1) is zonal and (,,2) is meridional velocity
!     ubarm1= barotropic velocity defined on U point for "tau-1" 
!             (,,1) is zonal and (,,2) is meridional velocity
      real, dimension(:,:,:), allocatable :: ps, uhat, ubar, ubarm1
      real, dimension(:,:), allocatable :: pguess, divf
#endif
#if defined explicit_free_surface
!     ps    = surface pressure (,,1) is for tau, (,,2) is for tau-1 
!     etat  = surface height on T-cell. 
!             (,,1) is for leap-frog tau and taup1,
!             (,,2) is for time mean tau and tau-1, 
!             (,,3) is for time mean tau and taup1
!     etau  = surface height on U-cells at tau.  diagnosed from etat
!     depthu = time varying depth of a U column; updated on baroclinic times.
!     depthur= 1/depthu
!     thicktau  =top level thickness at tau; from leap-frog etat 
!     emsk  = surface land/sea mask on T-cells
!     umsk  = surface land/sea mask on U-cells
!     ubar  = barotropic velocity defined on U point for "tau" 
!             (,,1) is zonal and (,,2) is meridional velocity
!     ubarm1= barotropic velocity defined on U point for "tau-1" 
!             (,,1) is zonal and (,,2) is meridional velocity
!     convU = convergence of vertically integrated flow (cm/sec)
!     rhosurf = density in surface level.
!     surftracer = tracer in surface level at tau and taup1.
!     psiu, psiv are approximate streamfunctions for depthu*ubar 
      real, dimension(:,:,:), allocatable :: ps
      real, dimension(:,:,:), allocatable :: etat
      real, dimension(:,:),   allocatable :: etau, depthu, depthur
      real, dimension(:,:),   allocatable :: thicktau
      real, dimension(:,:),   allocatable :: emsk, umsk
      real, dimension(:,:,:), allocatable :: ubar, ubarm1
      real, dimension(:,:),   allocatable :: convU
      real, dimension(:,:),   allocatable :: rhosurf
# if defined time_step_monitor && defined explicit_free_surface_nonlinear
      real, dimension(:,:,:,:), allocatable :: surftracer
# endif
# if defined snapshots || defined time_averages 
      real, dimension(:,:),   allocatable :: psiu, psiv
# endif
# ifdef explicit_fresh_water_flux
!     sfft = surface freshwater flux on t-cells
!     sffu = surface freshwater flux on u-cells
      real, dimension(:,:), allocatable :: sfft
      real, dimension(:,:), allocatable :: sffu
#  ifdef need_sfft_expl
      real, dimension(:,:), allocatable :: sfft_expl
#  endif
#  if !defined simple_sbc
!     strf = turbulent and non-neutral part of surface tracer flux on t-cells
!     (In the case of fresh water flux a neutral tracer flux 
!     does not change surface tracer concentration, T_1 = T_w)
      real, dimension(:,:,:), allocatable :: strf
#  endif
# endif
# if defined explhmix
!     bhd, chd, dhd, ghd, hhd
!            - coefficients for the horizontal friction
!              gives a minor correction
!              only available with explicit_eb or explicit_efb
      real, dimension(:), allocatable :: bhd, chd, dhd, ghd, hhd 
# endif
# if defined restnosp
      real, dimension(:,:), allocatable :: ptd
# endif
#endif
!
#ifdef stream_function
      real, dimension(:,:,:), allocatable :: psi
      real, dimension(:,:), allocatable :: ztd
#endif
#if defined rigid_lid_surface_pressure || defined implicit_free_surface || defined stream_function
      real, dimension(:,:,:,:), allocatable :: cf

#endif
#ifdef neptune
!     unep  = neptune velocity
!     pnep  = neptune streamfunction
      real, dimension(:,:,:), allocatable :: unep
      real, dimension(:,:), allocatable :: pnep
#endif
#if defined tracer_averages || defined term_balances || defined gyre_components
!     mskhr = mask field defining regions in the horizontal
!              (eg: mskhr(i,j) = n indicates point (i,j) is in the
!               "nth" horizontal region   where n=1..nhreg)
!               The "mskhr" masks are used in "diagt1" when
!               computing volume weighted tracer averages and in
!               "baroclinic.F", "tracer.F" and "diag2.F" when computing
!               term balances for tracers and momentum.
      integer, dimension(:,:), allocatable :: mskhr
#endif
#ifdef meridional_tracer_budget
!     msktmb   = mask for basin numbers (1..ntmbb. 0 is reserved for
!                sum of all basins)
      integer, dimension(:,:), allocatable :: msktmb
#endif
#ifdef kppvmix
!     hblt = boundary layer depth, on "t-grid"        (cm)
!     hmxl = mixed layer depth on "t-grid", where the
!            buoyancy difference from the surface
!            reaches a threshold, dbsfchmxl (set in
      real, dimension(:,:), allocatable :: hblt, hmxl
#endif
!     cori   = coriolis parameter for velocity component "n"
#if defined rot_grid
      real, dimension(:,:,:), allocatable :: cori
#else
      real, dimension(:,:), allocatable :: cori
#endif
#ifdef explicit_eta_dpdc
!     varables needed for mass conservation of the delplus-delcross filter.     
!     The fields are allocated and defined in subroutine init_surface.
!     ncoast   = number of coastal land points
!     icoast   = vector of the zonal index of coastal land points
!     jcoast   = vector of the meridional index of coastal land points
!     splr     = weight of the surrounding ocean points in 
!                "plus" directions 
!     The temporary fields are used to save memory.
!     neta      = index of etat which is filtered with eta_dpdc if
!     eta_mix   = dimensionless strength of the eta_dpdc filter
!                 explicit_eta_dpdc_once is on
! 
      integer                         :: ncoast, neta
      real                            :: eta_mix
      real, dimension(:), allocatable :: splrtmp
      real, dimension(:), allocatable :: splr
      real, dimension(:), allocatable :: icoasttmp
      real, dimension(:), allocatable :: icoast
      real, dimension(:), allocatable :: jcoasttmp
      real, dimension(:), allocatable :: jcoast
#endif
#ifdef explicit_eta_laplacian
!     neta      = index of etat which is filtered with eta_laplacian
!     eta_mix   = viscosity(cm2/s) for eta_laplacian
! 
      integer                         :: neta
      real                            :: eta_mix
#endif

      end module
