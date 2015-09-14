!======================= include file "mw.h" ===========================
!
!     Refer to the manual for a description of the memory window (MW)
!
!     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
!=======================================================================
!
#if !defined new_memory_window
!     taum1 = tau-1 time level for variables in MW
!     tau   = tau   time level for variables in MW
!     taup1 = tau+1 time level for variables in MW
!
      integer taum1, tau, taup1
      common /mwi/ taum1, tau, taup1
!
!-----------------------------------------------------------------------
!     MW arrays for prognostic equations:
!-----------------------------------------------------------------------
!
!     u(i,k,j,n,tau) = total velocity where:
!      i   = index for longitude
!      k   = index for depth
!      j   = index for latitude row within MW
!      n   = component (1 = zonal, 2 = meridional)
!      tau = time level (tau-1, tau, tau+1)
!             
!      (only internal modes are on disk and at tau+1 in the MW)
!
!     t(i,k,j,n,tau) = tracer where:
!      i   = index for longitude
!      k   = index for depth
!      j   = index for latitude row within MW
!      n   = component (1 = temperature, 2 = salinity)
!            if nt > 2 then other tracers are allowed.
!      tau = time level (tau-1, tau, tau+1)
!
!     note: temperature is potential temperature in degrees Celsius and 
!           salinity is in "model units", the deviation from 0.035 grams
!           of salt/cm**3 of water, or, assuming a water density of
!           1 gram/cm**3, the deviation from 0.035 g of salt/g of water.
!           one can convert model units to the common units of parts
!           per thousand (ppt) by adding 0.035 grams/cm**3 to the model
!           units and then multiplying by 1000. 
!     
      common /mw/ u(imt,km,jmw,2,-1:1), t(imt,km,jmw,nt,-1:1)
!
!-----------------------------------------------------------------------
!     MW arrays for diagnostic equations and workspace:
!-----------------------------------------------------------------------
!
!     diagnostic advective velocities are in units of cm/sec
!
!     adv_vet  = advective velocity on the eastern face of a T cell
!     adv_vnt  = advective velocity on the northern face of a T cell
!     adv_veu  = advective velocity on the eastern face of a U cell
!     adv_vnu  = advective velocity on the northern face of a U cell
!     adv_vbt  = advective velocity on the bottom face of a T cell
!     adv_vbu  = advective velocity on the bottom face of a U cell
!
!     rho      = density at center of a T cell in units of gm/cm**3
!                note: there is an arbitrary constant which is only a
!                function of depth in "rho". It is related to 
!                subtracting a reference level density for purposes of
!                accuracy.
!
!     pressure = hydrostatic pressure at the grid point within a T cell
!     grad_p   = hydrostatic pressure gradient for U cell. There are
!                two components: (1,2) is for (dp/dx, dp/dy)
!
      common /mw/ adv_vet(imt,km,2:jmw), adv_vnt(imt,km,1:jmw)
      common /mw/ adv_veu(imt,km,2:jmw-1)
      common /mw/ adv_vnu(imt,km,1:jmw-1)
      common /mw/ adv_vbt(imt,0:km,2:jmw)
#ifdef trajectories
      common /mw/ adv_vbu(imt,0:km,1:jmw-1)
#else
      common /mw/ adv_vbu(imt,0:km,2:jmw-1)
#endif
      common /mw/ rho(imt,km,1:jmw)
#ifdef pressure_gradient_average
      common /mw/ rhotaum1(imt,km,1:jmw)
      common /mw/ rhotaup1(imt,km,1:jmw)
      common /mw/ rhotilde(imt,km,1:jmw)
#endif
      common /mw/ pressure(imt,km,1:jmw)
      common /mw/ grad_p(imt,km,1:jmw-1,2)
!
!     tmask = tracer cell land/sea mask   = (0.0, 1.0) on (land, sea)
!     umask = velocity cell land/sea mask = (0.0, 1.0) on (land, sea)
!
      common /mw/ tmask(imt,km,1:jmw), umask(imt,km,1:jmw)
!
#ifdef partial_cell
      common /mw/ dht(imt,km,jmw+1), dhu(imt,km,jmw)
      common /mw/ dhwt(imt,0:km,jmw)
      common /mw/ pc_sink(imt,km,2:jmw-1)
      common /mw/ ztp(imt,km,jmw)
#endif
!
#if defined fourth_order_tracer_advection || defined quicker
!
!     adv_f4n = 4th order advective flux
!               
      common /mw/ adv_f4n(imt,km,1:jmw-1,nt)
#endif
!
#if defined velocity_horz_biharmonic
!
!     delsq for horz velocity components
!
      common /mw/ del2_vel(imt,km,2:jmw-1,2)
#endif
#if defined tracer_horz_biharmonic
!
!     delsq for tracers 
!
      common /mw/ del2_tracer(imt,km,2:jmw-1,nt)
#endif
!
!
!     these workspace arrays are recalculated for each component of the
!     equations so do not have to be moved as the MW moves northward. 
!
!     adv_fe   = advective flux across the eastern face of a cell
!     adv_fn   = advective flux across the northern face of a cell
!                (removed in most cases and put directly into the 
!                 statement functions for speed optimization.)
!     adv_fb   = advective flux across the bottom face of a cell
!
!     diff_fe  = diffusive flux across the eastern face of a cell
!     diff_fn  = diffusive flux across the northern face of a cell
!     diff_fb  = diffusive flux across the bottom face of a cell
!     source   = source term
!
      common /mw/ adv_fe(imt,km,2:jmw-1)
      common /mw/ adv_fn(imt,km,1:jmw-1)
      common /mw/ adv_fb(imt,0:km,2:jmw-1)
!
      common /mw/ diff_fe(imt,km,2:jmw-1)
      common /mw/ diff_fn(imt,km,1:jmw-1)
      common /mw/ diff_fb(imt,0:km,2:jmw-1)
#if defined isoneutralmix
      common /mw/ diff_fbiso(imt,0:km,2:jmw-1)
#endif
!
#ifdef source_term
      common /mw/ source(imt,km,2:jmw-1)
#endif
!
#if defined implicitvmix || defined isoneutralmix
      common /mw/ zzi(imt,km,2:jmw-1)
#endif
!
!     these grid factors are for optimizations
!     
      common /mw/ cstdxtr(imt,jmw), cstdxur(imt,jmw)
      common /mw/ cstdxt2r(imt,jmw),ah_cstdxur(imt,jmw)
      common /mw/ csudxtr(imt,jmw), csudxur(imt,jmw)
      common /mw/ csudxu2r(imt,jmw),am_csudxtr(imt,jmw)
!
!     advmet = coeff for metric advection.
!
      common /advec/ advmet(jmt,2)
!
!
#if defined linearized_advection || defined equatorial_thermocline
!
!     initial temperature profile independent of latitude and longitude
!
      common /profile/ tbarz(km,nt)
#endif
!
!     smf = surface momentum flux
!          1 => zonal wind stress (dynes/cm**2)
!          2 => meridional wind stress (dynes/cm**2)
!     bmf = bottom momentum flux
!          1 => zonal bottom drag (dynes/cm**2)
!          2 => meridional bottom drag (dynes/cm**2)
!     stf = surface tracer flux
!          1 => surface heat flux (cal/cm**2/sec = cm*degC/sec = ly/sec)
!          2 => surface salt flux (grams of salt/cm**2/sec)
!     stf_turb 
!         = the turbulent part of surface tracer flux
!     btf = bottom tracer flux (for consistency but normally zero!)
!          1 => bottom heat flux (cal/cm**2/sec = cm*degC/sec = ly/sec)
!          2 => bottom salt flux (grams of salt/cm**2/sec)
!
# ifdef air_pressure
      common /mw/ smf(imt,1:jmw,3),  bmf(imt,1:jmw,2)
# else
      common /mw/ smf(imt,1:jmw,2),  bmf(imt,1:jmw,2)
# endif
      common /mw/ stf(imt,1:jmw,nt), btf(imt,1:jmw,nt)
      common /mw/ stf_turb(imt,1:jmw,nt)
!
#ifdef tracer_advection_fct
!               "antidiffusive" flux as in Zalesak, 1989( see FCTstm.h)
!               same for R+, R-
!     anti_fe  = antidiffusive flux across the eastern face of a T cell
!     anti_fn  = antidiffusive flux across the northern face of a T cell
!     anti_fb  = antidiffusive flux across the bottom face of a T cell
!     R_plusY  = ratio of maximal feasible to maximal possible change
!                of tracer T in subroutine tracer.F, N-S dim delimiter
!     R_minusY = ratio of minimal feasible to minimal possible change
!                of tracer T in subroutine tracer.F, N-S dim delimiter
# ifdef fct_dlm2
!     t_lo     = low order solution for tracer
# endif
# ifdef fct_3d
!     R_plus3 = ratio of maximal feasible to maximal possible change
!               of tracer T in subroutine tracer.F, 3-d delimiter
!     R_minus3= ratio of minimal feasible to minimal possible change
!               of tracer T in subroutine tracer.F, 3-d delimiter
# endif
      common /mw/ anti_fe(imt,km,2:jmw-1,nt)
      common /mw/ anti_fn(imt,km,1:jmw-1,nt)
      common /mw/ anti_fb(imt,0:km,2:jmw-1,nt)
      common /mw/ R_plusY(imt,km,1:jmw-1,nt)
      common /mw/ R_minusY(imt,km,1:jmw-1,nt)
# ifdef fct_3d
      common /mw/ R_plus3(imt,km,1:jmw-1,nt)
      common /mw/ R_minus3(imt,km,1:jmw-1,nt)
# endif
# ifdef fct_dlm2
      common /mw/ t_lo(imt,km,2:jmw-1,nt)
# endif
#endif
!
#endif
