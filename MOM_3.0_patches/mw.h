c======================= include file "mw.h" ===========================
c
c     Refer to the manual for a description of the memory window (MW)
c
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
#if !defined new_memory_window
c     taum1 = tau-1 time level for variables in MW
c     tau   = tau   time level for variables in MW
c     taup1 = tau+1 time level for variables in MW
c
      integer taum1, tau, taup1
      common /mwi/ taum1, tau, taup1
c
c-----------------------------------------------------------------------
c     MW arrays for prognostic equations:
c-----------------------------------------------------------------------
c
c     u(i,k,j,n,tau) = total velocity where:
c      i   = index for longitude
c      k   = index for depth
c      j   = index for latitude row within MW
c      n   = component (1 = zonal, 2 = meridional)
c      tau = time level (tau-1, tau, tau+1)
c             
c      (only internal modes are on disk and at tau+1 in the MW)
c
c     t(i,k,j,n,tau) = tracer where:
c      i   = index for longitude
c      k   = index for depth
c      j   = index for latitude row within MW
c      n   = component (1 = temperature, 2 = salinity)
c            if nt > 2 then other tracers are allowed.
c      tau = time level (tau-1, tau, tau+1)
c
c     note: temperature is potential temperature in degrees Celsius and 
c           salinity is in "model units", the deviation from 0.035 grams
c           of salt/cm**3 of water, or, assuming a water density of
c           1 gram/cm**3, the deviation from 0.035 g of salt/g of water.
c           one can convert model units to the common units of parts
c           per thousand (ppt) by adding 0.035 grams/cm**3 to the model
c           units and then multiplying by 1000. 
c     
      common /mw/ u(imt,km,jmw,2,-1:1), t(imt,km,jmw,nt,-1:1)
c
c-----------------------------------------------------------------------
c     MW arrays for diagnostic equations and workspace:
c-----------------------------------------------------------------------
c
c     diagnostic advective velocities are in units of cm/sec
c
c     adv_vet  = advective velocity on the eastern face of a T cell
c     adv_vnt  = advective velocity on the northern face of a T cell
c     adv_veu  = advective velocity on the eastern face of a U cell
c     adv_vnu  = advective velocity on the northern face of a U cell
c     adv_vbt  = advective velocity on the bottom face of a T cell
c     adv_vbu  = advective velocity on the bottom face of a U cell
c
c     rho      = density at center of a T cell in units of gm/cm**3
c                note: there is an arbitrary constant which is only a
c                function of depth in "rho". It is related to 
c                subtracting a reference level density for purposes of
c                accuracy.
c
c     pressure = hydrostatic pressure at the grid point within a T cell
c     grad_p   = hydrostatic pressure gradient for U cell. There are
c                two components: (1,2) is for (dp/dx, dp/dy)
c
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
c
c     tmask = tracer cell land/sea mask   = (0.0, 1.0) on (land, sea)
c     umask = velocity cell land/sea mask = (0.0, 1.0) on (land, sea)
c
      common /mw/ tmask(imt,km,1:jmw), umask(imt,km,1:jmw)
c
#ifdef partial_cell
      common /mw/ dht(imt,km,jmw+1), dhu(imt,km,jmw)
      common /mw/ dhwt(imt,0:km,jmw)
      common /mw/ pc_sink(imt,km,2:jmw-1)
      common /mw/ ztp(imt,km,jmw)
#endif
c
#if defined fourth_order_tracer_advection || defined quicker
c
c     adv_f4n = 4th order advective flux
c               
      common /mw/ adv_f4n(imt,km,1:jmw-1,nt)
#endif
c
#if defined velocity_horz_biharmonic
c
c     delsq for horz velocity components
c
      common /mw/ del2_vel(imt,km,2:jmw-1,2)
#endif
#if defined tracer_horz_biharmonic
c
c     delsq for tracers 
c
      common /mw/ del2_tracer(imt,km,2:jmw-1,nt)
#endif
c
c
c     these workspace arrays are recalculated for each component of the
c     equations so do not have to be moved as the MW moves northward. 
c
c     adv_fe   = advective flux across the eastern face of a cell
c     adv_fn   = advective flux across the northern face of a cell
c                (removed in most cases and put directly into the 
c                 statement functions for speed optimization.)
c     adv_fb   = advective flux across the bottom face of a cell
c
c     diff_fe  = diffusive flux across the eastern face of a cell
c     diff_fn  = diffusive flux across the northern face of a cell
c     diff_fb  = diffusive flux across the bottom face of a cell
c     source   = source term
c
      common /mw/ adv_fe(imt,km,2:jmw-1)
      common /mw/ adv_fn(imt,km,1:jmw-1)
      common /mw/ adv_fb(imt,0:km,2:jmw-1)
c
      common /mw/ diff_fe(imt,km,2:jmw-1)
      common /mw/ diff_fn(imt,km,1:jmw-1)
      common /mw/ diff_fb(imt,0:km,2:jmw-1)
#if defined isoneutralmix
      common /mw/ diff_fbiso(imt,0:km,2:jmw-1)
#endif
c
#ifdef source_term
      common /mw/ source(imt,km,2:jmw-1)
#endif
c
#if defined implicitvmix || defined isoneutralmix
      common /mw/ zzi(imt,km,2:jmw-1)
#endif
c
c     these grid factors are for optimizations
c     
      common /mw/ cstdxtr(imt,jmw), cstdxur(imt,jmw)
      common /mw/ cstdxt2r(imt,jmw),ah_cstdxur(imt,jmw)
      common /mw/ csudxtr(imt,jmw), csudxur(imt,jmw)
      common /mw/ csudxu2r(imt,jmw),am_csudxtr(imt,jmw)
c
c     advmet = coeff for metric advection.
c
      common /advec/ advmet(jmt,2)
c
c
#if defined linearized_advection || defined equatorial_thermocline
c
c     initial temperature profile independent of latitude and longitude
c
      common /profile/ tbarz(km,nt)
#endif
c
c     smf = surface momentum flux
c          1 => zonal wind stress (dynes/cm**2)
c          2 => meridional wind stress (dynes/cm**2)
c     bmf = bottom momentum flux
c          1 => zonal bottom drag (dynes/cm**2)
c          2 => meridional bottom drag (dynes/cm**2)
c     stf = surface tracer flux
c          1 => surface heat flux (cal/cm**2/sec = cm*degC/sec = ly/sec)
c          2 => surface salt flux (grams of salt/cm**2/sec)
c     stf_turb 
!         = the turbulent part of surface tracer flux
c     btf = bottom tracer flux (for consistency but normally zero!)
c          1 => bottom heat flux (cal/cm**2/sec = cm*degC/sec = ly/sec)
c          2 => bottom salt flux (grams of salt/cm**2/sec)
c
# ifdef air_pressure
      common /mw/ smf(imt,1:jmw,3),  bmf(imt,1:jmw,2)
# else
      common /mw/ smf(imt,1:jmw,2),  bmf(imt,1:jmw,2)
# endif
      common /mw/ stf(imt,1:jmw,nt), btf(imt,1:jmw,nt)
      common /mw/ stf_turb(imt,1:jmw,nt)
c
#ifdef tracer_advection_fct
c               "antidiffusive" flux as in Zalesak, 1989( see FCTstm.h)
c               same for R+, R-
c     anti_fe  = antidiffusive flux across the eastern face of a T cell
c     anti_fn  = antidiffusive flux across the northern face of a T cell
c     anti_fb  = antidiffusive flux across the bottom face of a T cell
c     R_plusY  = ratio of maximal feasible to maximal possible change
c                of tracer T in subroutine tracer.F, N-S dim delimiter
c     R_minusY = ratio of minimal feasible to minimal possible change
c                of tracer T in subroutine tracer.F, N-S dim delimiter
# ifdef fct_dlm2
c     t_lo     = low order solution for tracer
# endif
# ifdef fct_3d
c     R_plus3 = ratio of maximal feasible to maximal possible change
c               of tracer T in subroutine tracer.F, 3-d delimiter
c     R_minus3= ratio of minimal feasible to minimal possible change
c               of tracer T in subroutine tracer.F, 3-d delimiter
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
c
#endif
