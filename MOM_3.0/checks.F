      subroutine checks (errorc, vmixset, htr_mixset, hvel_mixset)
      use mom3_mpp_mod
      use arrays_2d
      use memory_window
      use sbc_info
      use sbc_ocean
      logical errorc, vmixset, htr_mixset, hvel_mixset
#include "param.h"
#include "accel.h"
#include "coord.h"
#include "derived_options.h"
#include "grdvar.h"
#include "hmixc.h"
#include "index.h"
#include "iounit.h"
#ifdef isoneutralmix
# include "isopyc.h"
#endif
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "taskrows.h"
#include "vmixc.h"
#ifdef save_xbts
# include "diag_xbt.h"
#endif
      integer(kind=8) :: icg, jcg
      common / checks_common / icg, jcg, cgmax, gridmin, sml
c
c-----------------------------------------------------------------------
c     do consistency checks before allowing MOM to continue
c         
c     author:   R. C. Pacanowski     e-mail=> rcp@gfdl.gov
c    (mpp items)   v. balaji         e-mail=> vb@gfdl.gov
c-----------------------------------------------------------------------
c
      write (stdout,'(/,20x,a,/,8x,a,a,/)')
     & 'G E N E R A L    C O N S I S T E N C Y    C H E C K S'
     &,'(A "rudimentary intelligence" analysis of all specifications'
     &,' and their interactions)'
c
#ifdef monthly_averages
        write (stdout,'(/,(1x,a))')
     & '==> Warning: option "monthly_averages" is enabled. This option'
     &,'             will over-ride the specified averaging period for'
     &,'             any diagnostic which averages data in time.      '
#endif
c
c     explicit free surface option checks
c
# if defined explicit_free_surface_linear && !defined explicit_free_surface
      write(stdout,9000)  
     &'==> Error: explicit_free_surface_linear must also have        '
     &,'          explicit_free_surface enabled.                     '
        errorc = .true.
# endif
#if defined explicit_free_surface_nonlinear 
      write(stdout,9000)
     & '==> Warning: Option explicit_free_surface_nonlinear is enabled.'
     &,'             dzt(k=1)+eta must be positive in order to       '
     &,'             preclude vanishingly small top cell thicknesses.'
     &,'             If dzt(k=1)+eta < 0, the model will not conserve'
     &,'             first moments and may blow-up.                  '
      if(dtuv .ne. dtts) then
      write(stdout,9000)
     & '==> Warning: Option explicit_free_surface_nonlinear is enabled'
     &,'             yet dtts is not equal to dtuv.  In order to      '
     &,'             have a conservative model, set dtts=dtuv.        '
      endif
# if !defined explicit_free_surface
      write(stdout,9000)  
     &'==> Error: explicit_free_surface_nonlinear must also have     '
     &,'          explicit_free_surface enabled.                     '
        errorc = .true.
# endif
# if !defined partial_cell
      write(stdout,9000)  
     &'==> Error: explicit_free_surface_nonlinear must also have     '
     &,'          partial_cell enabled.  If wish to run free surface '
     &,'          without partial cells, enable the option           '
     &,'          -Dexplicit_free_surface_linear, which will provide '
     &,'          a nonconservative form of the free surface method. '
        errorc = .true.
# endif
#endif
#if defined explicit_free_surface_linear
        write (stdout,9000)
     & '==> Warning: a variable thickness for level k=1 implied by    '
     &,'             option explicit_free_surface has not been        '
     &,'             accounted for in the depth dependent equations.  '  
     &,'             Consequently, the model is not conservative.     '
     &,'             Disable -Dexplicit_free_surface_linear and enable'
     &,'             -Dexplicit_free_surface_nonlinear to get a       '
     &,'             variable thickness and to be conservative.       '
#endif
#if defined explicit_free_surface_linear && defined explicit_free_surface_nonlinear 
      write(stdout,9000)  
     &'==> Error: Options explicit_free_surface_nonlinear and         '
     &,'          explicit_free_surface_linear cannot both be enabled.'
        errorc = .true.
#endif
#if defined explicit_free_surface && defined time_averaged_surface
# if !defined explicit_free_surface_linear && !defined explicit_free_surface_nonlinear 
      write(stdout,9000)  
     &'==> Error: Option explicit_free_surface must also have either '
     &,'          explicit_free_surface_linear or                    '
     &,'          explicit_free_surface_nonlinear enabled.           '
        errorc = .true.
# endif       
#endif
#if defined time_averaged_surface && defined explicit_dpdc
      write(stdout,9000)  
     &'==> Error: Option "explicit_dpdc" is not available with option' 
     &,'          "time_averaged_surface"                            '
        errorc = .true.
#endif

#ifdef explicit_free_surface
# if !defined time_averaged_surface
        if (dtuv .ne. dtts) then
          write (stdout,9000)
     & '==> Warning: when dtts > dtuv, the explicit free surface      '
     &,'             may blow up unless options "explicit_eb" and     '
     &,'             "explicit_efb" are disabled.                     '
        endif
#  if defined fourfil || defined firfil
	write (stdout,9000)
     & '==> Error: polar filtering with option fourfil or firfil is '
     &,'           not advisable without using the time averaged    '
     &,'           version of the explicit free surface.    '
	  errorc = .true.
#  endif
# endif
# if !defined explicit_efb && !defined explicit_eb
	  write (stdout,9000)
     & '==> Warning: because options explicit_eb or explicit_efb have '
     &,'             not been enabled, the explicit free surface '
     &,'             solution will be time averaged over two baroclinic'
     &,'             timesteps to keep the barotropic solution defined'
     &,'             at the same time as the baroclinic solution. This'             
     &,'             will allow stretching of tracer timesteps relative'
     &,'             to baroclinic timesteps.'
# endif
# if defined explhmix && !defined velocity_horz_mix_const
	  write (stdout,9000)
     & '==> Error: option explhmix works only with                   '
     &,'    option "velocity_horz_mix_const".  Disable "explhmix"    '
	  errorc = .true.
# endif
# if defined explhmix && !defined velocity_horz_laplacian
	  write (stdout,9000)
     & '==> Error: option "explhmix" only works with                 '
     &,'     only "velocity_horz_laplacian".                         '
	  errorc = .true.
# endif
#endif
c
c     fresh water flux option checks
c
#if defined explicit_fresh_water_flux
      ifresh = index_of_sbc ('fresh wtr')
      if (ifresh .eq. 0) then
        write(stdout,*)  
     &  '==> Error:  Option explicit_fresh_water_flux is enabled but '
     &, '            the fresh water S.B.C. has not been setup. '
        errorc = .true.
      endif
      i1 = index_of_sbc ('fresh wtr')
      i2 = index_of_sbc ('saltflux')
      if (i1 .eq. i2 .and. i1 .ne. 0) then
        write(stdout,*)  
     &  '==> Error:  cannot have both a saltflux and fresh water S.B.C'
     &, '            for tracer n=2. Configure for one or the other. '
        errorc = .true.
      endif
# if !defined explicit_free_surface
      write(stdout,9000)  
     &'==> Error:  Option explicit_fresh_water_flux requires option '
     &,'           explicit_free_surface be enabled.                '
      errorc = .true.
# endif
# if defined explicit_free_surface_linear
      write(stdout,9000)  
     &'==> Warning: Adding fresh water with option                  '
     &,'            "explicit_free_surface_linear" will not conserve'
     &,'            tracers.  Switch to option                      '
     &,'            "explicit_free_surface_nonlinear" to do so.     '
# endif
# if !defined salinity_psu
      write(stdout,9000)  
     &'==> Error:  Option explicit_fresh_water_flux also requires   '
     &,'           option salinity_psu to be enabled.               '
      errorc = .true.
# endif
#endif
#if defined river_inflow && !defined explicit_fresh_water_flux
      write(stdout,*)  
     & '==> Error:  Option river_inflow has not been coded to work    '
     &,'            without option explicit_fresh_water_flux.         '
      errorc = .true.
#endif
#if defined salinity_psu && defined coupled
      write(stdout,9000)  
     &'==> Error:  Option "coupled" was not changed to work with    '
     &,'           option salinity_psu              .               '
      errorc = .true.
#endif
c
c
#if defined symmetry && defined extend_isotropic_grid
        write (stdout,'(/,(1x,a))')
     & '==> Error: option extend_isotropic_grid  is not allowed with   '
     &,'           option symmetry.                               '
        errorc = .true.
#endif
c
c       options which have not been parallelized  
c
#if defined parallel_1d
# if defined crayio
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "crayio" yet.                           '
        errorc = .true.
# endif
# if defined fio
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "fio" yet.                              '
        errorc = .true.
# endif
# if defined coupled
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "coupled" yet.                          '
        errorc = .true.
# endif
# if defined gyre_components
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "gyre_components" yet.                  '
        errorc = .true.
# endif
# if defined meridional_overturning
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "meridional_overturning" yet.           '
        errorc = .true.
# endif
# if defined save_density_terms 
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "save_density_terms" yet.           '
        errorc = .true.
# endif
# if defined term_balances
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "term_balances" yet.                    '
        errorc = .true.
# endif
# if defined tracer_averages
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "tracer_averages" yet.                  '
        errorc = .true.
# endif
# if defined show_zonal_mean_of_sbc
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "show_zonal_mean_of_sbc" yet.           '
        errorc = .true.
# endif
# if defined diagnostic_surf_height
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "diagnostic_surf_height".               '
        errorc = .true.
# endif
# if defined meridional_tracer_budget
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "meridional_tracer_budget" yet.         '
        errorc = .true.
# endif
# if defined trajectories
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "trajectories" yet.                     '
        errorc = .true.
# endif
# if defined stability_tests
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "stability_tests" yet.                  '
        errorc = .true.
# endif
# if defined matrix_sections
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "matrix_sections" yet.                  '
        errorc = .true.
# endif
# if defined tracer_yz
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "tracer_yz" yet.                        '
        errorc = .true.
# endif
# if defined cross_flow
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "cross_flow" yet.                       '
        errorc = .true.
# endif
# if defined density_terms 
        write (stdout,'(/,(1x,a))')
     & '==> Error: option parallel_1d has has not been extended to work'
     &,'           with option "density_terms" yet.                    '
        errorc = .true.
# endif
#endif
#if defined air_pressure
# if !defined explicit_free_surface
        write (stdout,'(/,(1x,a))')
     & '==> Error:  option air_pressure has only been coded for option'
     &,'            explicit_free_surface.                            '
        errorc = .true.
# endif
#endif
#if !defined parallel_1d
      if (num_processors .ne. 1) then
        write (stdout,'(1x,a)')
     & '==> Error: option parallel_1d must be enabled if               '
     &,'           num_processors > 1                                  '
        write (stdout,*) ' num_processors = ',num_processors
        errorc = .true.
      endif
#endif
#if defined stream_function
      if (num_processors .ne. 1) then
        write (stdout,'(1x,a)')
     & '==> Error: option stream_function only works on one processor'
        write (stdout,*) ' num_processors = ',num_processors
        errorc = .true.
      endif
# if defined parallel_1d
      write (stdout,'(1x,a)')
     & '==> Error: option parallel_1d has not been extended to work    '
     &,'           with option stream_function                         '
      errorc = .true.
# endif
#endif
#if defined rigid_lid_surface_pressure && defined parallel_1d
      write (stdout,'(1x,a)')
     & '==> Error: option parallel_1d has not been extended to work    '
     &,'           with option rigid_lid_surface_pressure              '
      errorc = .true.
#endif
      if (imt .lt. 3) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  parameter "imt" less than 3 is not allowed'        
        errorc = .true.
      endif
c
      if (jmt .lt. 4) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  parameter "jmt" less than 4 is not allowed'        
        errorc = .true.
      endif
      if (num_processors .gt. 1 .and. num_processors .le. jmt-2) then
        write (stdout,'(/,1x,a,i4,a)')
     & '==> Warning:  "num_processors" =',num_processors
     &,' implies multi-tasking (parallel processing)'
      elseif (num_processors .eq. 1) then
        write (stdout,'(/,(1x,a))')
     & '==> Warning:  "num_processors" = 1 implies unitasking          '
      elseif (num_processors .gt. jmt-2 .or. num_processors .lt. 1) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  "num_processors" cannot be > jmt-2 or < 1          '
        write (stdout,*)' num_processors = ',num_processors
        write (stdout,*) 
     &  'num_processors is set via namelist "procs" in run script '
	errorc = .true.
      endif
c
#ifdef xlandmix
# ifndef source_term
      write (stdout,9000)
     & '==> Error: ifdef "xlandmix" is enabled but "source_term"     '
     &,'             is not. "source_term" needs to be enabled  '
      errorc = .true.
# endif
#endif
c
c     more options not checked out in parallel
c
#if defined parallel_1d
# if defined xlandmix
      write (stdout,'(/,(1x,a))')
     & '==> Error:  option "xlandmix" has not been checked out for use '
     &,'            with option "parallel_1d".                         '
      errorc = .true.
# endif
# if defined pressure_gradient_average
      write (stdout,'(/,(1x,a))')
     & '==> Error: option "pressure_gradient_average" needs to have row'
     &,'           jrow=jecomp updated via communication calls. This'
     &,'           has not been done yet.'
	errorc = .true.
# endif
# ifdef equivalence_mw
      write (stdout,'(/,(1x,a))')
     & '==> Error:  "equivalence_mw" is not allowed with option       '
     &,'            "parallel_1d"                                     '
      errorc = .true.
# endif       
# if defined timing
      if (num_processors .gt. 1) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  option "timing" is not allowed with option        '
     &,'            "parallel_1d"                                     '
        errorc = .true.
      endif
# endif       
#endif
c
c     some misc option checks
c
#if defined save_density_terms 
        write (stdout,'(/,(1x,a))')
     & '==> Warning: option save_density_terms has not been run    '
     &,'             much at GFDL.  It may be buggy.               '
#endif
#ifdef pressure_gradient_average
      if (.not. eb) then
        write (stdout,'(/,(1x,a))')
     & '==> Warning: option "pressure_gradient_average" may need      '
     &,'             euler backward mixing timesteps to suppress the  '
     &,'             computational modes.                             '
      endif
# if !defined max_window
      write (stdout,'(/,(1x,a))')
     & '==> Error: option "pressure_gradient_average" requires option '
     &,'            max_window.                                       '
      errorc = .true.
# endif
# if defined pressure_depth
      write (stdout,'(/,(1x,a))')
     & '==> Error: option "pressure_gradient_average" has not been    '
     &,'           checked out with option "pressure_depth"           '
      errorc = .true.
# endif
# ifdef robert_time_filter
      write (stdout,'(/,(1x,a))')
     & '==> Error: "robert_time_filter" has not been implemented for  '
     &,'            option "pressure_gradient_average" yet.           '
      errorc = .true.
# endif
      write (stdout,'(/,(1x,a))')
     & '==> Error: "pressure_gradient_average" has not been checked   '
     &,'            out yet. Remove this stop condition to continue.  '
      errorc = .true.
#endif
c
#ifdef linearized_advection
# if defined tracer_advection_fct
      write (stdout,'(/,(1x,a))')
     & '==> Error: "tracer_advection_fct" is incompatible with '
     &,'             "linearized_advection", Cannot have both. '
        errorc = .true.
# endif
# if defined fourth_order_tracer_advection
      write (stdout,'(/,(1x,a))')
     & '==> Error:  "fourth_order_tracer_advection" is incompatible   '
     &,'            with "linearized_advection".                      '
     &,'            Can`t have both.                                  '
        errorc = .true.
# endif
# if defined quicker
      write (stdout,'(/,(1x,a))')
     & '==> Error:  "quicker" is incompatible                         '
     &,'            with "linearized_advection".                      '
     &,'            Can`t have both.                                  '
        errorc = .true.
# endif
# if defined levitus_ic || defined idealized_ic
      write (stdout,'(/,(1x,a))')
     & '==> Error:  "levitus_ic" or "idealized_ic" cannot be          '
     &,'            be enabled when using "linearized_advection".     '
        errorc = .true.
# endif
# if !defined linearized_density
      write (stdout,'(/,(1x,a))')
     & '==> Error:  "linearized_density" must be enabled when using   '
     &,'            "linearized_advection".                           '
        errorc = .true.
# endif
#else
# if !defined obctest && !defined obctest2 && !defined obctest3
#  if !defined levitus_ic && !defined idealized_ic
#   if !defined equatorial_thermocline && !defined read_ic
      write (stdout,'(/,(1x,a))')
     & '==> Error:  either "levitus_ic", "idealized_ic", "read_ic" or '
     &,'            "equatorial_thermocline" must be enabled to       '
     &,'             specify the type of initial conditions           '
        errorc = .true.
#   endif
#  endif
# endif
c
c       Check for obsolete options
c
#if defined xbts 
        write (stdout,9000)
     & '==> Error: option "xbts" has been changed to "save_xbts".     '
        errorc = .true.
#endif
#if defined netcdf 
        write (stdout,9000)
     & '==> Error: option "netcdf" has been changed to "use_netCDF".  '
        errorc = .true.
#endif
#if defined isopycmix
        write (stdout,9000)
     & '==> Error: isopycmix is not an option. Use -Dredi_diffusion   '
     &,'           refer to the manual for other isoneutral options.  '
        errorc = .true.
#endif
#if defined merid_by_basin
        write (stdout,9000)
     & '==> Error: Option "merid_by_basin" has been changed to        '
     &,'           "meridional_overturning_basin"                     '
        errorc = .true.
#endif
#if defined consthmix
        write (stdout,9000)
     & '==> Error: consthmix is not an option. Use instead            '
     &,'    tracer_horz_mix_const and/or velocity_horz_mix_const    '
        errorc = .true.
#endif
#if defined smagnlmix
        write (stdout,9000)
     & '==> Error: smagnlmix is not an option. Use instead            '
     &,'    tracer_horz_mix_smag and/or velocity_horz_mix_smag      '
        errorc = .true.
#endif
#if defined varhmix
        write (stdout,9000)
     & '==> Error: varhmix is not an option. Use instead              '
     &,'    tracer_horz_mix_var and/or velocity_horz_mix_var        '
        errorc = .true.
#endif
c
c Check initial condition options
c
      n = 0
#if defined obctest || defined obctest2 || defined obctest3
      n = n + 1
#endif
#if defined levitus_ic
      n = n + 1
#endif
#if defined read_ic
      n = n + 1
#endif
#if defined idealized_ic
      n = n + 1
#endif
#if defined equatorial_thermocline
      n = n + 1
#endif
      if (n .ne. 1) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  Either "levitus_ic" or "idealized_ic" or "read_ic"'
     &,'            or "equatorial_thermocline" must be enabled.      '
     &,'             choose only one for initial conditions.          '
        errorc = .true.
      endif
#endif
c
c     Check boundary condition options
c
#if defined minimize_sbc_memory && defined time_mean_sbc_data
      write (stdout,'(/,(1x,a))')
     & '==> Error:  "minimize_sbc_memory" and "time_mean_sbc_data"    '
     &,'            are incompatible. "minimize_sbc_memory" is        '
     &,'            intended for option "time_varying_sbc_data".      '
        errorc = .true.
#endif
c
      n = 0
#if defined time_mean_sbc_data
      n = n + 1
#endif
#if defined time_varying_sbc_data
      n = n + 1
#endif
#if defined simple_sbc
      n = n + 1
#endif
#if defined coupled
      n = n + 1
#endif
      if (n .ne. 1) then
        write (stdout,'(/,(1x,a))')
     & '==> Error:  Either "simple_sbc" or "time_mean_sbc_data" or    '
     &,'            "time_varying_sbc_data" or "coupled" must be      '
     &,'             enabled. choose one surface boundary condition   '
        errorc = .true.
      endif
#if defined constant_taux && !defined simple_sbc
      write (stdout,'(/,(1x,a))')
     & '==> Error: "constant_taux"  only works with "simple_sbc".     '
     &,'           it has not been extended to work with other SBCs   '
      errorc = .true.
#endif
#if defined constant_tauy && !defined simple_sbc
      write (stdout,'(/,(1x,a))')
     & '==> Error: "constant_tauy"  only works with "simple_sbc".     '
     &,'           it has not been extended to work with other SBCs   '
      errorc = .true.
#endif
c
c     Check lateral mixing options
c 
#if defined tracer_horz_mix_const 
      if (htr_mixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "tracer_horz_mix_const"  cannot be enabled        '
     &,'    because another horizontal mixing scheme has been enabled.'
        errorc     = .true.
      else
        htr_mixset = .true.
      endif
#endif
#if defined tracer_horz_mix_var
      if (htr_mixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "tracer_horz_mix_var"  cannot be enabled        '
     &,'    because another horizontal mixing scheme has been enabled.'
        errorc     = .true.
      else
        htr_mixset = .true.
      endif
#endif
#if defined tracer_horz_mix_var && defined isoneutralmix
        write (stdout,'(/,(1x,a))')
     & '==> Error: "tracer_horz_mix_var" cannot be enabled         '
     &,'    along with "isoneutralmix".'
        errorc     = .true.
#endif
#if defined tracer_horz_mix_smag && !defined velocity_horz_mix_smag
        write (stdout,'(/,(1x,a))')
     & '==> Error: "tracer_horz_mix_smag"  cannot be enabled        '
     &,'    unless "velocity_horz_mix_smag" is also enabled.        '
        errorc     = .true.
#endif
#if defined tracer_horz_mix_var && !defined velocity_horz_mix_smag && !defined velocity_horz_mix_var
        write (stdout,'(/,(1x,a))')
     & '==> Error: "tracer_horz_mix_var"  cannot be enabled           '
     &,'    unless "velocity_horz_mix_smag" or "velocity_horz_mix_var"'
     &,'    are also enabled.                                         '
        errorc     = .true.
#endif
#if defined bryan_lewis_horizontal
      if (htr_mixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "bryan_lewis_horizontal"  cannot be enabled        '
     &,'    because another horizontal mixing scheme has been enabled.'
        errorc     = .true.
      else
        htr_mixset = .true.
      endif
#endif
      if (.not.htr_mixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: No horizontal tracer diffusivity scheme enabled.  '
        errorc = .true.
      endif
#if defined velocity_horz_mix_const 
      if (hvel_mixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "velocity_horz_mix_const"  cannot be enabled      '
     &,'    because another horizontal mixing scheme has been enabled.'
        errorc     = .true.
      else
        hvel_mixset = .true.
      endif
#endif
#if defined velocity_horz_mix_var
      if (hvel_mixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "velocity_horz_mix_var"  cannot be enabled        '
     &,'    because another horizontal mixing scheme has been enabled.'
        errorc     = .true.
      else
        hvel_mixset = .true.
      endif
#endif
#if defined am_cosine && !defined velocity_horz_mix_var
      write (stdout,'(/,(1x,a))')
     & '==> Error: "velocity_horz_mix_var" must be on with "am_cosine"'
      errorc = .true.
#endif
#if defined am_cosinep1
      write (stdout,'(/,(1x,a))')
     & '==> Error:"am_cosinep1" is not an option anymore. '
      errorc = .true.
#endif
#if defined am_taper_highlats
      write (stdout,'(/,(1x,a))')
     & '==> Error:"am_taper_highlats" is not an option anymore. '
      errorc = .true.
#endif
      if (.not.hvel_mixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: No horizontal velocity viscosity scheme enabled.  '
        errorc = .true.
      endif
#ifdef neptune
# if defined partial_cell
        write (stdout,9000)
     & '==> Error:  "partial_cell", has not been extended to work with'
     &,'             "neptune".                 '
        errorc = .true.
# endif
#endif
c
c     Check vertical mixing options
c
#ifdef constvmix
      if (vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "constvmix"  cannot be enabled because another     '
     &,'            vertical mixing scheme has been enabled           '
        errorc = .true.
      else
        vmixset = .true.
      endif
#endif
#if defined kppmix && !defined fill_isolated_cells
        write (stdout,'(/,(1x,a))')
     & '==> Error: "kppmix" cannot work in diffusion-only grid points'
     &,'          use "fill_isolated_cells" to eliminate such points.'
        errorc = .true.
#endif
#ifdef bryan_lewis_vertical
      if (.not. vmixset) then
        write (stdout,'(/,(1x,a/1x,a/1x,a/1x,a/1x,a/1x,a))')
     & '==> Error: "bryan_lewis_vertical" tracer diffusion requires   '
     &,'           "constvmix", "ppvmix", "tcvmix", or "kppvmix".     '
     &,'           If "constvmix" is enabled, "bryan_lewis_vertical"  '
     &,'           coefficients replace kappa_h. If "ppvmix","tcvmix",'
     &,'           or "kppvmix" is enabled, "bryan_lewis_vertical"    '
     &,'           coefficients are added to the predicted coeffs.    '
        errorc = .true.
      endif
#endif
      if (.not.vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: No vertical mixing scheme has been enabled  '
        errorc = .true.
      endif
#ifdef tcvmix
      if (vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "tcvmix" cannot be enabled because another         '
     &,'            vertical mixing scheme has been enabled           '
        errorc = .true.
      else
        vmixset = .true.
      endif
#endif
      if (.not.vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: there is no vertical mixing scheme enabled         '
        errorc = .true.
      endif

#if defined partial_cell
	write (stdout,9000)
     & '==> Warning: option partial_cell uses the UNESCO equation of '
     &,'             state for hydrostatic pressure gradients.       '
#endif
#if defined bbl_ag
	write (stdout,9000)
     & '==> Warning: option bbl_ag uses the UNESCO equation of       '
     &,'             state for hydrostatic pressure gradients.       '
#endif
c
c       Options for memory management of diagnostics 
c
#if defined diagnostics_lomem && defined diagnostics_himem
          write(stdout,9000)
     & '==> Error: both "diagnostics_lomem" and                      '
     &,'           "diagnostics_himem" are enabled. Can only use one.'
        errorc = .true.        
#endif
#if defined diagnostics_lomem
          write(stdout,9000)
     & '==> Note: "diagnostics_lomem" is enabled.                    '
     &,'          Diagnostics of 3d fields will be buffered within   '
     &,'          the memory window.  This approach saves memory, but'
     &,'          is considerably slower than diagnostics_himem.     '
#endif
#if defined diagnostics_himem
          write(stdout,9000)
     & '==> Note: "diagnostics_himem" is enabled.                    '
     &,'          Diagnostics of 3d fields will be buffered over the '
     &,'          full model.   This approach requires much memory,  ' 
     &,'          but is considerably faster than diagnostics_lomem. '
#endif
c
c       Check for incompatible options
c        
#if defined fourfil && defined firfil
	write (stdout,9000)
     & '==> Error: either use option fourfil or firfil but not both.  '
	  errorc = .true.
#endif
#if defined rigid_lid_surface_pressure
        write (stdout,9000)
     & '==> Warning: when using the "rigid_lid_surface_pressure" method'
     &,'             a checkerboard may develop in "ps". Using option '
     &,'             "remove_ps_checkerboard" may help to remove it.'
#endif
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
# if defined sf_5_point
        write (stdout,9000)
     & '==> Warning:  "sf_5_point" is ignored when using the'
     &,'              "implicit_free_surface" or '
     &,'              "rigid_lid_surface_pressure" option.'
     &,'              conjugate gradient with 9 point numerics is used.'
# endif
#endif
c
#if defined implicit_free_surface
# if defined fourfil || defined firfil
        write (stdout,9000)
     & '==> Error: the implicit free surface equations are not         '
     &,'           compatible with options firfil or fourfil           '
     &,'           remove the stop condition to continue               '
        errorc = .true.  
# endif
#endif
#if !defined implicit_free_surface && !defined explicit_free_surface
# if !defined rigid_lid_surface_pressure && !defined stream_function 
        write (stdout,9000)
     & '==> Error: enable "stream_function", "implicit_free_surface"   '
     &,'           "explicit_free_surface" or the                      '
     &,'           "rigid_lid_surface_pressure" option                 '
        errorc = .true.  
# endif
#endif
c
#ifdef stream_function
# if !defined sf_5_point && !defined sf_9_point
        write (stdout,9000)
     & '==> Error: choose "sf_5_point" or "sf_9_point" when using'
     &,'           the stream_function method'
        errorc = .true.  
# endif
# if defined sf_5_point && defined sf_9_point
        write (stdout,9000)
     & '==> Error: choose "sf_5_point" or "sf_9_point" when using'
     &,'           the stream_function method (not both)'
        errorc = .true.  
# endif
#endif
#if !defined stream_function && defined zero_island_flow
        write (stdout,9000)
     & '==> Error:  "zero_island_flow" requires "stream_function"      '
     &,'            to be enabled                                      '
        errorc = .true.  
#endif
#if defined parallel_1d
        write (stdout,9000)
     & '==> Warning:  "parallel_1d" is enabled.                        '
        write (stdout,*) '              Note:', num_processors
     &, ' processors were requested.'
c
c     look for potential load imbalance 
c
      ncr = jecomp3(1)-jscomp3(1)+1
      i   = 0
      do n=2,num_processors
        if (ncr .ne. (jecomp3(n)-jscomp3(n)+1)) i = 1
      enddo
      if (i .eq. 1) then
        write (stdout,*)
     & '==> Warning: There is a potential load imbalance if the number '
     &,'             of computed rows is not the same on all processors'
        do n=1,num_processors
          ncr = jecomp3(n)-jscomp3(n)+1
          write (stdout,*) 'Processor #',n,' computes ',ncr,' rows'
        enddo
      endif
#endif
#if defined sixth_order_window
      if (jmw .ne. 7 .and. jmw .ne. window_size) then
        write (stdout,9000)
     & '==> Error: jmw is set wrong for a sixth_order_window'
     &,'           only the following sizes are allowed:'
        write (stdout,*)'          jmw=7 or jmw=',window_size
        write (stdout,*)'          you have set jmw=',jmw
	write (stdout,*)'          reset jmw to ',window_size
	write (stdout,*)'          in size.h'
        errorc = .true.  
      else
        write (stdout,9000)
     & '==> Note: a sixth_order_window is being used.'
      endif
#endif
#if defined fourth_order_window
      if (jmw .ne. 5 .and. jmw .ne. window_size) then
        write (stdout,9000)
     & '==> Error: jmw is set wrong for a fourth_order_window'
     &,'           only the following sizes are allowed:'
        write (stdout,*)'          jmw=5 or jmw=',window_size
        write (stdout,*)'          you have set jmw=',jmw
	write (stdout,*)'          reset jmw to ',window_size
	write (stdout,*)'          in size.h'
        errorc = .true.  
      else
        write (stdout,9000)
     & '==> Note: a fourth_order_window is being used.'
      endif
#endif
#if defined second_order_window
      if (jmw .ne. 3 .and. jmw .ne. window_size) then
        write (stdout,9000)
     & '==> Error: jmw is set wrong for a second_order_window'
     &,'           only the following sizes are allowed:'
        write (stdout,*)'          jmw=3 or jmw=',window_size
        write (stdout,*)'          you have set jmw=',jmw
	write (stdout,*)'          reset jmw to ',window_size
	write (stdout,*)'          in size.h'
        errorc = .true.  
      else
        write (stdout,9000)
     & '==> Note: a second_order_window is being used.'
      endif
#endif
#if defined max_window
      write (stdout,9000)  '==> Warning: The MW is open all the way ' 
#endif
#if defined pressure_depth && !defined max_window
      write (stdout,9000)
     & '==> Error: option pressure_depth requires option max_window'
        errorc = .true.  
#endif
#if defined pressure_depth
# if !defined partial_cell && !defined bbl_ag
      write (stdout,9000)
     & '==> Error: option pressure_depth requires option partial_cell'
     &,'           or option bbl_ag'
        errorc = .true.  
# endif
#endif
c
      if (nkflds .lt. 2) then
        write (stdout,9000)
     & '==> Error: "nkflds" must be at least 2                        '
        write (stdout,*)'          nkflds is set = ',nkflds
        errorc = .true.  
      endif
c
#ifdef equivalence_mw
        write (stdout,9000)
     & '==> Warning: The MW is equivalenced to the poisson solver     '
     &,'           workspace to reduce memory requirements. Verify    ' 
     &,'           that results are correct by comparing to a case    ' 
     &,'           with the "equivalence_mw" option removed.          ' 
      if (jmw .eq. jmt) then
        write (stdout,9000)
     & '==> Error: The MW cannot be equivalenced to the poisson solver'
     &,'           workspace but when "jmw" = "jmt"'
     &,'           remove the "equivalence_mw" option!' 
        errorc = .true.
      endif
#endif
#ifdef read_my_kmt
        write (stdout,9000)
     & '==> Warning: The topography is being imported. Are you sure   '
     &,'             that the "kmt" field matches the vertical grid?  ' 
     &,'             i.e. Are the grid points within levels given by  '
     &,'             k = kmt(i,jrow) really defined at depth "zt(k)"? ' 
#endif
c
#ifdef restorst
      if (damp_factor(1) .eq. c0 .or. damp_factor(2) .eq. c0) then
        write (stdout,9000)
     & '==> Error: the damping time scale "dampts" must be set > 0    '
     &,'           when enabling option "restorst"                    '
        errorc = .true.
      endif
      if (damp_factor(1) .ne. damp_factor(2)) then
        write (stdout,9000)
     & '==> Warning: are you aware that the damping time scale "dampts"'
     &,'           is set differently for T and S in option "restorst"?'
     &,'           (they can be set the same in the namelist)          '
      endif
#else
        write (stdout,9000)
     & '==> Warning: the damping time scale "dampts" is > zero but    '
     &,'           the "restorst" option is not enabled               '
#endif
#if !defined ramdrive && !defined crayio && !defined fio
# if !defined ssread_sswrite
        write (stdout,9000)
     & '==> Error: either "ramdrive", crayio", "fio", or              '
     &,'           "ssread_sswrite" option must be turned on.         '
      errorc = .true.
# endif
#endif
#if defined ramdrive && defined crayio
        write (stdout,9000)
     & '==> Error:  "ramdrive" and "crayio" are incompatible          '
     &,'             choose one or the other option.                  '
        errorc = .true.
#endif
#if defined ramdrive && defined fio
        write (stdout,9000)
     & '==> Error:  "ramdrive" and "fio" are incompatible             '
     & '             choose one or the other option.                  '
        errorc = .true.
#endif
#if defined ramdrive && defined ssread_sswrite
        write (stdout,9000)
     & '==> Error:  "ramdrive" and "ssread_sswrite" are incompatible  '
     & '             choose one or the other option.                  '
        errorc = .true.
#endif
#if defined crayio && defined fio
        write (stdout,9000)
     & '==> Error:  "crayio" and "fio" are incompatible               '
     &,'             choose one or the other option.                  '
        errorc = .true.
#endif
#if defined ssread_sswrite && defined fio
        write (stdout,9000)
     & '==> Error:  "ssread_sswrite" and "fio" are incompatible       '
     &,'             choose one or the other option.                  '
        errorc = .true.
#endif
#if defined crayio && defined ssread_sswrite
        write (stdout,9000)
     & '==> Error:  "crayio" and "ssread_sswrite" are incompatible    '
     &,'             choose one or the other option.                  '
        errorc = .true.
#endif
#if defined crayio && !defined _CRAY
        write (stdout,9000)
     & '==> Error:  You must be on a Cray (where macro _CRAY is       '
     &,'             enabled) for "crayio" to work                    '
        errorc = .true.
#endif
#if defined ssread_sswrite && !defined CRAYPVP
        write (stdout,9000)
     & '==> Error:  You must be on a Cray PVP machine                 '
     &,'                     for "ssread_sswrite" to work             '
        errorc = .true.
#endif
#if defined ssread_sswrite
        write (stdout,9000)
     & '==> Warning:  "ssread_sswrite" implies CRAY SSD is available  '
     &,'             (SSD is Solid State Disk)                        '
#endif
c
#if defined fourth_order_tracer_advection
        write (stdout,9000)
     & '==> Warning: "fourth_order_tracer_advection" has a more       '
     &,'             stringent stability requirement than the standard'
     &,'             second order advection. If problems arise, the   '
     &,'             time step may have to be reduced. Refer to the   '
     &,'             manual                                           '
        write (stdout,9000)
     & '==> Warning: "fourth_order_tracer_advection" assumes constant '
     &,'             grid spacing.                                    '
# if defined partial_cell
        write (stdout,9000)
     & '==> Error:  "partial_cell", has not been extended to work with'
     &,'             "fourth_order_tracer_advection".                 '
        errorc = .true.
# endif
#endif
c
#if defined cross_flow && !defined isoneutralmix
        write (stdout,9000)
     & '==> Error: isoneutralmix is needed for  cross_flow           '
        errorc = .true.
#endif
#if defined save_ncdiff && !defined nonconst_diffusivity
        write (stdout,9000)
     & '==> Error: "nonconst_diffusivity" is needed for "save_ncdiff" '
        errorc = .true.
#endif
#if save_density_terms && !defined redi_diffusion
        write (stdout,9000)
     &'  ==> Error: redi_diffusion is needed for                      '
     &,' save_density_terms.  Set ahiso = 0 if not wanting            '
     &,' any nontrivial diffusive fluxes, but wanting to diagnose the '
     &,' potential density terms.'
        errorc = .true.
#endif
c
#ifdef redi_diffusion
# if defined tracer_horz_mix_const || defined tracer_horz_mix_var
        write (stdout,9000)
     & '==> Note: since redi_diffusion was specified for tracers,     '
     &,'    kappa_h and Ah will be used as background mixing coeffs   '
# endif
# if defined tracer_horz_biharmonic
        write (stdout,9000)
     & '==> Error:  tracer_horz_biharmonic is incompatible with    '
     &,'    redi_diffusion.                                        '
        errorc = .true.
# endif
#endif
#if !defined redi_diffusion && defined gm_advect
        write (stdout,9000)
     & '==> Error:"gm_advect" without "redi_diffusion" is not supported'
        errorc = .true.
#endif
c
c
#ifdef obc
# if defined partial_cell
        write (stdout,9000)
     & '==> Error:  "partial_cell", has not been extended to work with'
     &,'             "open boundary conditions".                 '
        errorc = .true.
# endif
# ifdef symmetry
        write (stdout,9000)
     &  '==> Error:  "obc" and "symmetry" together is not compatible  '
     &, '             choose one or the other option.                 '
        errorc = .true.
# endif
# if defined tracer_horz_biharmonic || defined velocity_horz_biharmonic
        write (stdout,9000)
     &  '==> Warning:  At the open boundaries, laplacian mixing is used'
# endif
# if defined explicit_free_surface && !defined orlanski
        write (stdout,9000)
     &  '==> Note:  The use of "obc" together with free surface is   '
     & ,'           tested for "orlanski" only. Without "orlanski" it'
     & ,'           is not fully tested out.                         '
#  if !defined time_averaged_surface
        write (stdout,9000)
     &  '==> Error:  The use of "obc" together with free surface is  '
     & ,'            prepared only for use with                      '
     & ,'            "time_averaged_surface".                        '
        errorc = .true.
#  endif
# endif
# ifdef cyclic 
#  if defined obc_west || defined obc_east
        write (stdout,9000)
     &  '==> Error:  "obc_west" or "obc_east" and "cyclic" together is'
     &, '             not compatible. Choose one or the other option. '
        errorc = .true.
#  endif
# endif
#endif
c
#if defined orlanski && !defined obc
        write (stdout,9000)
     &  '==> Error:  For Orlanski radiation condition you have to     '
     &, '            choose an open boundary.                         '
        errorc = .true.
#endif
c
#if defined obcparameter 
# ifndef obc
        write (stdout,9000)
     &  '==> Error:  For writing obcparameter into snapshots you have '
     &, '            to choose open boundaries.                       '
        errorc = .true.
# endif
# ifndef snapshots
        write (stdout,9000)
     &  '==> Error:  For writing obcparameter into snapshots you have '
     &, '            to choose "snapshots".                           '
        errorc = .true.
# endif
#endif
c
c
#if defined implicit_free_surface && defined stream_function
        write (stdout,9000)
     & '==> Error:  "implicit_free_surface" and "stream_function" are '
     &,'             incompatible. choose one or the other option.    '
        errorc = .true.
#endif
#if defined rigid_lid_surface_pressure && defined stream_function
        write (stdout,9000)
     & '==> Error: "rigid_lid_surface_pressure" and "stream_function" '
     &,'           are incompatible. choose one or the other option.  '
        errorc = .true.
#endif
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
# ifdef diagnostic_surf_height
        write (stdout,9000)
     & '==> Warning:  "diagnostic_surf_height" does not apply under'
     &,'              the "rigid_lid_surface_pressure" or'
     &,'              "implicit_free_surface" options.'
# endif
#endif
#ifdef stream_function
# ifdef sf_9_point
        write (stdout,9000)
     & '==> Warning: although "sf_9_point" conserves total energy     '
     &,'             it has a 2 grid point null space.                '
# endif
# ifdef sf_5_point
        write (stdout,9000)
     & '==> Warning: although "sf_5_point" has no null space, it does '
     &,'             not conserve total energy.                       '
# endif
#endif
c
#if defined free_slip
        write (stdout,9000)
     & '==> Error: option "free_slip" has been removed because it can'
     &,'           not be made to work for arbitrary geometry        '
        errorc = .true.
#endif
#if defined redi_diffusion && defined tracer_horz_mix_const
        if ((ah+ahisop) .gt. 1.e11) then
          write (stdout,9000)
     & '==> Error: "ahisop"+"ah" is too large for the                 '
     &,'             "redi_diffusion" mixing option                   '
          errorc = .true.
        endif
#endif
c
        if (dtsf .le. c0) then
          write (stdout,9000)
     & '==> Error: need to set the external mode time step "dtsf"     '
          errorc = .true.
          dtsf = epsln
	endif
c
        if (dtuv .le. c0) then
          write (stdout,9000)
     & '==> Error: need to set the internal mode time step "dtuv"     '
          errorc = .true.
          dtuv = epsln
	endif
c
        if (dtts .le. c0) then
          write (stdout,9000)
     & '==> Error: need to set the density time step "dtts"           '
          errorc = .true.
          dtts = epsln
	endif
c
        critv = 1.e-6
        if (mod(rundays,dtts*secday) .gt. critv) then
          t1 = nint(rundays/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: there must be an integral number of density time '
     &,'    steps within "rundays" (the integration time).            '
          write (stdout,*) '               (changed "rundays" from '
     &,   rundays,' days to ', t1,' days to insure this condition)'
          rundays = t1
        endif
c
c-----------------------------------------------------------------------
c       time step checks
c-----------------------------------------------------------------------
c
#if defined damp_inertial_oscillation
      if (acor .lt. p5) then
        write (stdout,9000)
     & '==> Error: "acor" must be between 0.5 & 1.0 when              '
     &,'           option damp_inertial_oscillation is enabled.       '
          errorc = .true.
      endif
# ifdef tracer_advection_fct
      write (stdout,9000)
     & '==> Error: non-linear instability will result if             '
     &,'           "tracer_advection_fct" is used with               '
     &,'           "damp_inertial_oscillation"                       '
          errorc = .true.
# endif
#else
      if (acor .ne. 0) then
        write (stdout,9000)
     & '==> Error: "acor" must=0 when option damp_inertial_oscillation'
     &,'           is not enabled.  '
          errorc = .true.
      else
c
c       check for marginally resolved inertial oscillation
c
# ifdef rot_grid
        fmax = epsln
        do jrow=jscomp,jecomp
          do i=2,imtm1
	    if (kmu(i,jrow) .ne. 0) then
              fmax = max(fmax,abs(cori(i,jrow,1)))
	    endif
          enddo
        enddo
# else
        fmax = epsln
        do jrow=jscomp,jecomp
	  do i=2,imtm1
	    if (kmu(i,jrow) .ne. 0) then
              fmax = max(fmax,abs(c2*omega*sine(jrow)))
	    endif
	  enddo
        enddo
# endif
        dtuvmax = (1.0/6.0)*(c2*pi)/fmax
        write (stdout,'(/a,/,11x,a,f8.0,a)')
     &' ==> Note: 6 timesteps/(min inertial period) implies a maximum'
     &,'dtuv of ',dtuvmax,' sec.'
        if (dtuv .gt. dtuvmax) then
          write (stdout,9000)
     & '==> Error: the inertial oscillation is not resolved. reduce    '
     &,'           "dtuv" or use option "damp_inertial_oscillation"    '
     &,'            or remove stop condition to continue.'
          errorc = .true.
        endif
      endif
#endif
c
c       CFL time step limit for internal gravity waves
c
        cgint = 300.0
	dtcg  = 1.e10
	icg   = 1
	jcg   = jmt/2
        dxumin = 1.0e20
	do jrow=jscomp,jecomp
	  do i=2,imt-1
	    if (kmu(i,jrow) .ne. 0) then
	      dxumin = min(dxumin,dyu(jrow),csu(jrow)*dxu(i))
	      dti = 0.5*dxumin/cgint
	      if (dti .lt. dtcg) then
	        dtcg = dti
	        icg  = i
	        jcg  = jrow
	      endif
	    endif
	  enddo
	enddo
	dtcg = nint(dtcg)
	write (stdout,*)
     & '==> Note: based on a maximum internal gravity wavespeed of'
     &, nint(cgint/100),' m/sec,'
        write (stdout,*)
     & '          the CFL restriction implies a maximum dtuv =',dtcg
     &,' sec'
        write (stdout,*)
     &'          at ','(i,jrow) = (',icg,',',jcg,') on the U grid where'
        write (stdout,*)
     &'          the minimum grid spacing = ~',nint(dxumin/1.e5),' km.'
        write (stdout,*)
     &'          The specified dtuv = ',nint(dtuv),'sec.'
#if defined fourfil || defined firfil
        flmax = max(filter_reflat_n,-filter_reflat_s)
        dxfil = cos(flmax*pi/180.0)*dxt(2)
	dtf   = nint(0.5*dxfil/cgint)
        dtcg  = max(dtcg,dtf)
	write (stdout,*)
     & '   (since polar filtering with option fourfil or firfil is'
	write (stdout,*)
     & '    enabled, the maximum dtuv can be extended to ~',dtcg,' sec)'
        write (stdout,*)
        dxmax = c0
	dxmin = 1.e10
        do i=2,imt-1
	  dxmax = max(dxu(i),dxmax)
	  dxmin = min(dxu(i),dxmin)
	enddo
	if (dxmax .ne. dxmin) then
	  write (stdout,9000)
     & '==> Error: Cannot have variable resolution in zonal direction'
     &,'           with polar filering (fourfil or firfil)           '
          write (stdout,*)'           dxmax=',dxmax,' dxmin=',dxmin
          errorc = .true.
	endif
#endif
c
c       CFL horizontal advective restriction for tracers
c
        itmin  = imt/2
        jtmin  = jmt/2
        dxtmin = min(dyt(jmt/2),cst(jmt/2)*dxt(imt/2))
	do jrow=jscomp,jecomp
	  do i=2,imt-1
	    if (kmt(i,jrow) .ne. 0) then
	      dxtmin2 = min(dyt(jrow),cst(jrow)*dxt(i))
              if(dxtmin2 .lt. dxtmin) then
                dxtmin = dxtmin2
                itmin  = i
                jtmin  = jrow
              endif
            endif
	  enddo
	enddo
c
        dt_adv = dtts
	c_adv = nint(0.5*dxtmin/dt_adv)
	write (stdout,*)
     & '==> Note: based on a minimum T-cell grid length of approx'
     &, nint(dxtmin/1.e5),' km'
        write (stdout,*)
     & '          at the (i,j) grid point (',itmin,',',jtmin
     &,  ') and a tracer timestep dtts =',nint(dtts),' sec,'
        write (stdout,*)
     & '          the CFL restriction implies a maximum horizontal '
        write (stdout,*)
     & '          advection velocity of ',nint(c_adv),' cm/sec.'
        write (stdout,*)
     & '          advection velocity greater than this is unstable'
#if defined fourfil || defined firfil
        flmax  = max(filter_reflat_n,-filter_reflat_s)
        dxfil  = cos(flmax*pi/180.0)*dxt(2)
	dcf    = nint(0.5*dxfil/dt_adv)
        dcf = max(c_adv,dcf)
	write (stdout,*)
     & '   (since polar filtering with option fourfil or firfil is'
	write (stdout,*)
     & '    enabled the maximum advective velocity is ~',dcf,' cm/sec.)'
        write (stdout,*)
#endif
c
c       CFL vertical advective restriction for tracers
c
        dt_adv = dtts
	dztmin = dzt(1)
	do k=2,km
	  if (dztmin .gt. dzt(k)) then dztmin = dzt(k)
	enddo
	c_adv = 0.5*dztmin/dt_adv
	write (stdout,*)
     & '==> Note: based on a minimum T-cell grid thickness of approx'
     &, nint(dztmin/1.e2),' m'
        write (stdout,*)
     & '          and a tracer timestep dtts =',nint(dtts),' sec,'
        write (stdout,*)
     & '          the CFL restriction implies a maximum vertical '
        write (stdout,'(a,es10.3,a)')
     & '          advection velocity of ',c_adv,' cm/sec.'
        write (stdout,*)
     & '          advection velocity greater than this is unstable'
#ifdef stream_function
        write (stdout,9000)
     & '==> Warning: Large topographic gradients will reduce the max   '
     &,'             allowable timestep when using the stream function.'
     &,'             Free surface methods do not have this restriction.'
#endif
	if (dtuv .gt. dtcg) then
	  write (stdout,9000)
     &'==> Error: reduce dtuv or remove this stop condition.'
          errorc = .true.
        endif
#ifdef explicit_free_surface
	if (abs(nint(dtuv/dtsf)*dtsf-dtuv) .gt. critv) then
	  write (stdout,9000)
     & '==> Warning: for explicit free surface mod(dtuv,dtsf)=0 must  '
     &,'             be satisfied. It is not. Original setting: '
	  write (stdout,*) 'dtuv =  ',dtuv,' sec'
	  write (stdout,*) 'dtsf =  ',dtsf,' sec'
	  write (stdout,9000) '             New setting:               '
          dtsf = dtuv/(int(dtuv/dtsf)+1)
	  write (stdout,*) 'dtsf =  ',dtsf,' sec'
	endif
        if( abs(dtuv/dtsf/2 - nint(dtuv/dtsf/2)) .gt. critv ) then        
          write(*,9000)
     & '==> Error: for explicit free surface dtuv/dtsf must be         '
     & ,'          an even integer to maintain proper centering        '
     & ,'          of barotropic time on the baroclinic time.          '
     & ,'          Reset dtsf and/or dtuv to satisfy this constraint.  '
	  errorc = .true.
        endif 
c
        icg   = 1
	jcg   = 1
	cgmax = 1.0
	dtcgmin = 1.e10
	gridmin = dxu(2)
        do jrow=jscomp,jecomp
	  do i=2,imt-1
	    if (kmu(i,jrow) .gt. 0) then
	      depth = zw(kmu(i,jrow))
	      cgrav = sqrt(grav*depth)
	      gridsp = min(dyu(jrow), csu(jrow)*dxu(i))
	      dtcg = 0.5*(gridsp)/cgrav
	      if (dtcg .lt. dtcgmin) then
	        gridmin = gridsp
                dtcgmin = dtcg
	        cgmax = cgrav
	        icg   = i
	        jcg   = jrow
	      endif
	    endif
	  enddo
	enddo
	dtcgmin = nint(dtcgmin)
	cgmax   = nint(cgmax)
# ifdef parallel_1d
! tp: this code does not work correcly on T90 but does on T3E
! tp: Balaji should look at this
        sml = dtcgmin
        call mpp_max(sml)
        if( pe.EQ.0 )then
            if( sml.NE.dtcgmin )
     &           call mpp_transmit( icg, 4, NULL_PE, icg, 4, ANY_PE  )
        else
            if( sml.EQ.dtcgmin )
     &           call mpp_transmit( icg, 4,       0, icg, 4, NULL_PE )
        endif
        call mpp_sync()
# endif
	write (stdout,9000) '==> Note: for the explicit free surface.'
        write (stdout,*)
     & '          The fastest gravity wave speed =  ',cgmax,'cm/sec'
     &,' at i=',icg,' and jrow=',jcg,' on the U grid.'
        write (stdout,*)
     & '          where the grid length = ~',nint(gridmin/1.e5),' km.'
	write (stdout,*)
     & '          For stability, dtsf must be less than ',dtcgmin,' sec'
	write (stdout,*)
     & '          The specified dtuv =  ',dtuv,' sec and dtsf = '
     &, dtsf,' sec. '
	write (stdout,*)
# if defined time_averaged_surface
     & '          This implies ',2*nint(dtuv/dtsf),' (due to averaging)'
# else
     & '          This implies ',nint(dtuv/dtsf)
# endif
     &,' barotropic time steps'
        write (stdout,*) 
     & '          per baroclinic time step.'
# if defined explicit_polar_filter 
	dtf_fs   = (dxfil/gridmin)*dtcgmin
	write (stdout,*)
     & '   (since polar filtering with option fourfil or firfil'
     &,'    is enabled, '
	write (stdout,*)
     & '    the maximum dtsf can be extended to ~',nint(dtf_fs),' sec.)'
        write (stdout,*)
	if (dtsf .gt. dtf_fs) then
	  write (stdout,9000)'==> Error: for explicit free surface.  '
     &,'          dtsf > time step allowed by fastest wave speed'
	  errorc = .true.
	endif
# else
	if (dtsf .gt. dtcgmin) then
	  write (stdout,9000)'==> Error: for explicit free surface.  '
     &,'          dtsf > time step allowed by fastest wave speed'
	  errorc = .true.
	endif
#  if defined fourfil || defined firfil
	write (stdout,9000)
     & '==> Note: No polar filtering applied to barotropic fields.   '
     &,'          enable option "explicit_polar_filter" to do so. But'
     &,'          note that for many configurations, polar filtering '
     &,'          of these fields can be more expensive than         '
     &,'          resolving the dynamics.                            '
#  endif
# endif
#endif
c
        if (itmb) then
          write (stdout,9000)
     & '==> Warning: "itmb" is set to "true". set it to "false" in    '
     &,'            subsequent runs to prevent the time independent   '
     &,'            basin mask from being written more than once. This'
     &,'            reduces the size of the diagnostic file.          '
	endif
c
        if (itrmb) then
          write (stdout,9000)
     & '==> Warning: "itrmb" is set to "true". set it to "false" in   '
     &,'            subsequent runs to prevent the time independent   '
     &,'            region masks from being written more than once.   '
     &,'            This reduces the size of the diagnostic file.     '
	endif
c
        if (itavg) then
          write (stdout,9000)
     & '==> Warning: "itavg" is set to "true". set it to "false" in   '
     &,'            subsequent runs to prevent the time independent   '
     &,'            region masks from being written more than once.   '
     &,'            This reduces the size of the diagnostic file.     '
	endif
#ifdef trajectories
        if (trajint .eq. c0) then
          write (stdout,9000)
     & '==> Warning: "trajectories" is set but "trajint"=0. this means'
     &,'            particle positions will be written every time step'
        endif
        if (trajint .lt. c0) then
          write (stdout,9000)
     & '==> Error: "trajectories" is set but "trajint" < 0. this means'
     &,'            particle positions will not be written.           '
          errorc = .true.
        endif
# ifndef lyapunov
          write (stdout,9000)
     & '==> Warning: "trajectories" is enabled but "lyapunov" is not. '
# endif
#else
# ifdef lyapunov
          write (stdout,9000)
     & '==> Warning: "trajectories" is not enabled but "lyapunov" is. '
# endif
#endif
#ifdef meridional_tracer_budget
        if (mod(tmbint,dtts*secday) .gt. critv .and. tmbint .gt. c0)
     &    then
          t1 = nint(tmbint/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: "tmbint" does not contain an integral number of  '
     &,'    density time steps "dtts".                                '
          write (stdout,*) '               (changed "tmbint" from '
     &,   tmbint,' days to ', t1,' days to insure this condition)'
          tmbint = t1
        endif
        if (tmbint .eq. c0) then
          write (stdout,9000)
     & '==> Error: averaging interval "tmbint"  must be specified when'
     &,'   the "meridional_tracer_budget" option is enabled           '
          errorc = .true.
        endif
        if (tmbint .gt. tmbper) then
          write (stdout,9000)
     & '==> Warning: the interval "tmbint" exceeds the averaging'
     &,'          period "tmbper" for option "meridional_tracer_budget"'
        endif
        if (tmbint .lt. tmbper) then
          write (stdout,9000)
     & '==> Warning: averaging period "tmbper" exceeds interval'
     &,'            "tmbint". Setting tmbper = tmbint          '
          tmbper = tmbint
	endif
        if (tmbper .eq. c0) then
          write (stdout,9000)
     & '==> Warning: the averaging period "tmbper" is zero. The '
     &,'            average will be over only one time step!    '
        endif
#else
      if (tmbint .gt. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging interval "tmbint" is > zero but the'
     &,'             the "meridional_tracer_budget" option is not on. '
      endif
#endif
#ifdef time_averages
# if defined time_averages_disk
          write (stdout,9000)
     & '==> Error: "time_averages_disk" is no longer available.       '
          errorc = .true.
# endif 
        if (mod(timavgint,dtts*secday).gt.critv .and. timavgint .gt. c0)
     &    then
          t1 = nint(timavgint/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: "timavgint" does not contain an integral number of'
     &,'    density time steps "dtts".                                '
          write (stdout,*) '               (changed "timavgint" from '
     &,   timavgint,' days to ', t1,' days to insure this condition)'
          timavgint = t1
        endif
        if (timavgint .eq. c0) then
          write (stdout,9000)
     & '==> Warning: averaging interval "timavgint" = 0.0 implies '
     &,'   no averaging when option "time_averages" is enabled    '
        endif
        if (timavgint .gt. timavgper) then
          write (stdout,9000)
     & '==> Warning: the interval "timavgint" exceeds the averaging'
     &,'             period "timavgper" for option "time_averages" '
        endif
        if (timavgint .lt. timavgper) then
          write (stdout,9000)
     & '==> Warning: averaging period "timavgper" exceeds interval'
     &,'            "timavgint". Setting timavgper = timavgint    '
          timavgper = timavgint
        endif
        if (timavgper .eq. c0) then
          write (stdout,9000)
     & '==> Warning: the averaging period "timavgper" is zero. The '
     &,'            average will be over only one time step!       '
        endif
#else
      if (timavgint .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging period "timavgint"  is > zero but  '
     &,'           the "time_averages" option is not enabled          '
      endif
#endif
#ifdef save_xbts
        if (mod(xbtint,dtts*secday) .gt. critv .and. xbtint .gt. c0)
     &    then
          t1 = nint(xbtint/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: "xbtint" does not contain an integral number of  '
     &,'    density time steps "dtts".                                '
          write (stdout,*) '               (changed "xbtint" from '
     &,   xbtint,' days to ', t1,' days to insure this condition)'
          xbtint = t1
        endif
        if (xbtint .eq. c0) then
          write (stdout,9000)
     & '==> Warning: averaging interval "xbtint" = 0.0 implies no '
     &,'   averaging when option "save_xbts" is enabled           '
        endif
        if (numxbt .eq. 0) then
          write (stdout,9000)
     &   '==> Error: the number of XBT stations detected is zero      '
     &,  '           at least one station must be specified           '
          errorc = .true.
        endif
        if (xbtint .gt. xbtper) then
          write (stdout,9000)
     & '==> Warning: the interval "xbtint" exceeds the  averaging'
     &,'             period "xbtper" for option "save_xbts"      '
        endif
        if (xbtint .lt. xbtper) then
          write (stdout,9000)
     & '==> Warning: averaging period "xbtper" exceeds interval'
     &,'            "xbtint". Setting xbtper = xbtint          '
          xbtper = xbtint
        endif
        if (xbtper .eq. c0) then
          write (stdout,9000)
     & '==> Warning: the averaging period "xbtper" is zero. The '
     &,'            average will be over only one time step!    '
        endif
#else
      if (xbtint .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging interval "xbtint"  is > zero but   '
     &,'             the "save_xbts" option is not enabled            '
      endif
#endif
#ifdef diagnostic_surf_height
        if (mod(dspint,dtts*secday) .gt. critv .and. dspint .gt. c0)
     &    then
          t1 = nint(dspint/(dtts*secday))*dtts*secday
          write (stdout,9000)
     & '==> Warning: "dspint" does not contain an integral number of  '
     &,'    density time steps "dtts".                                '
          write (stdout,*) '               (changed "dspint" from '
     &,   dspint,' days to ', t1,' days to insure this condition)'
          dspint = t1
        endif
        if (dspint .eq. c0) then
          write (stdout,9000)
     & '==> Warning: averaging interval "dspint" = 0.0 implies no '
     &,'   averaging when option "diagnostic_surf_height" is enabled'
        endif
        if (dspint .gt. dspper) then
          write (stdout,9000)
     & '==> Warning: the interval "dspint" exceeds the averaging '
     &,'            period "dspper" for option "diagnostic_surf_height"'
        endif
        if (dspint .lt. dspper) then
          write (stdout,9000)
     & '==> Warning: averaging period "dspper" exceeds interval'
     &,'            "dspint". Setting dspper = dspint          '
          dspper = dspint
	endif
        if (dspper .eq. c0) then
          write (stdout,9000)
     & '==> Warning: the averaging period "dspper" is zero. The '
     &,'            average will be over only one time step!    '
        endif
#else
      if (dspint .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging interval "dspint"  is > zero but'
     &,'             option "diagnostic_surf_height" is not enabled'
      endif
#endif
c
      if ((dtuv .ne. dtsf) .or. (dtuv .ne. dtts)) then
        write (stdout,9000)
     & '==> Warning: use of unequal time steps implies the transient  '
     &,'             response may be inaccurate.                      '
      endif
c
c     check for mixing coefficients larger than stability permits
c
      dymin  = dyt(2)
      dxmin  = dxt(2)
      dzmin  = dzt(1)
      xlmax  = dtxcel(1)
#ifdef isoneutralmix
      fzimax = fzisop(1)
#endif
      do jrow=2,jmtm1
        dymin  = min(dymin,dyt(jrow))
      enddo
      do i=2,imtm1
        dxmin  = min(dxmin,dxt(i))
      enddo
      do k=2,km
        xlmax  = max(xlmax,dtxcel(k))
        dzmin  = min(dzmin,dzt(k))
#ifdef isoneutralmix
        fzimax = max(fzimax,fzisop(k))
#endif
      enddo
c
      if (xlmax .gt. c1) then
        write (stdout,9000)
     & '==> Warning: use of accelerated time steps implies the        '
     &,'             transient response is unimportant and multiple   '
     &,'             equilibria do not exist. stability tests will    '
     &,'             use "dtts" multiplied by the maximum "dtxcel"    '
      endif
c
c Check that size of viscosity is not too large
c
#if defined velocity_horz_laplacian && !defined velocity_horz_mix_smag
      dtxl = dtts*xlmax
      num = 0
      do j=2,jmtm1
        dxdymn = c1/(c1/(dxmin*cst(j))**2 + c1/dyt(j)**2)
        ah_crit = p25*dxdymn/dtxl
        if (ah .gt. ah_crit) then
          num = num + 1
          if (num .eq. 1) write (stdout,9000)
     & '==> Warning: lateral diffusive criteria exceeded for "ah".    '
     &,'              use a smaller "dtts", "dtxcel", or "ah"         '
          write (stdout,'(a48,f6.2,a5,i3,a,e14.7)') ' at latitude '
     &,   yt(j),',  j=',j,', the critical value of "ah" =',ah_crit
        endif
      enddo
      num = 0
      do j=2,jmtm2
        dxdymn = c1/(c1/(dxmin*csu(j))**2 + c1/dyu(j)**2)
        am_crit = p25*dxdymn/(dtuv + epsln)
        amtest  = am
# if defined am_cosine
        amtest = am*csu(j)
	write (stdout,'(a,i4,a,f6.2,a,1pe14.7)')
     & 'Lateral viscosity coefficient for jrow = ',j,' lat=',yu(j)
     &,  ' is am*cos(jrow) = ',amtest
# endif
# if defined am_cosinep1
        amtest = am*(csu(j)+c1)
	write (stdout,'(a,i4,a,f6.2,a,1pe14.7)')
     & 'Lateral viscosity coefficient for jrow = ',j,' lat=',yu(j)
     &,  ' is am*(cos(jrow)+1) = ',amtest
# endif
# if defined am_taper_highlats
        amtaper_min=am/10.
        amtaper_start_lat=60.
        amtaper_end_lat=87.
        amtaper = am
        if(abs(yu(j)) .gt. amtaper_start_lat) then
           amtaper = am -
     $          (am-amtaper_min)*(amtaper_start_lat-abs(yu(j)))**2
     $          /(amtaper_start_lat-amtaper_end_lat)**2
        endif
        if (amtaper .lt. amtaper_min) amtaper = amtaper_min
	write (stdout,'(a,i4,a,f6.2,a,1pe14.7)')
     & 'Lateral viscosity coefficient for j = ',j,' lat=',yu(j)
     &,  ' is = ',amtaper
        amtest = amtaper 
# endif
        if (amtest .gt. am_crit) then
          num = num + 1
          if (num .eq. 1) write (stdout,9000)
     & '==> Warning: lateral diffusive criteria exceeded for "am".    '
     &,'             use a smaller "dtuv" or "am"                     '
          write (stdout,'(a48,f6.2,a5,i3,a,e14.7)') ' at latitude '
     &,   yu(j),',  j=',j,', the critical value of "am" =',am_crit
        endif
      enddo
      num = 0
      do j=2,jmtm2
        cbeta = 2.28e-13
        am_munk = cbeta*csu(j)*(dxmin*csu(j)*sqrt(3.0)/pi)**3
        if (am .le. am_munk) then
          num = num + 1
          if (num .eq. 1) write (stdout,9000)
     & '==> Warning: Munk boundary layer unresolved. use a            '
     &,'             larger "am" (showing first 50 violations)        '
          if (num .le. 50) then
            write (stdout,'(a48,f6.2,a5,i3,a,e14.7)') ' at latitude '
     &,     yu(j),',  j=',j,', the Munk value of "am" =',am_munk
          endif
        endif
      enddo
#endif
#ifdef sponges
# ifndef source_term
      write (stdout,9000)
     & '==> Error: ifdef "sponges" is enabled but "source_term"     '
     &,'             is not. "source_term" needs to be enabled  '
      errorc = .true.
# endif
#endif
#ifdef shortwave
# ifndef source_term
      write (stdout,9000)
     & '==> Error: ifdef "shortwave" is enabled but "source_term"     '
     &,'             is not. "source_term" needs to be enabled  '
      errorc = .true.
# endif
# ifndef coupled
      write (stdout,9000)
     & '==> Warning: ifdef "shortwave" is enabled but "coupled"       '
     &,'             is not. Is the shortwave component of radiation  '
     &,'             defined ?                                        '
# endif
#else
      if (dzt(1) .lt. 20.0e2) then
        write (stdout,9000)
     & '==> Warning: if shallow mixed layers develop, then enabling   '
     &,'             ifdef "shortwave" may help to deepen them. note  '
     &,'             that either you or the atmosphere must provide   '
     &,'             the solar short wave as a boundary condition.    '
      endif
#endif
#ifdef constvmix
      do k=1,km
        if ((dtts*dtxcel(k)*kappa_h)/dzt(k)**2 .ge. p25) then
          write (stdout,9000)
     & '==> Warning: vertical diffusive criteria exceeded on "kappa_h"'
     &,'             use a smaller "dtts", "dtxcel", and/or "kappa_h" '
         write (stdout,'(a48,i3)') ' at level =',k
        endif
      enddo
      if ((dtuv*kappa_m)/dzmin**2 .ge. p25) then
        write (stdout,9000)
     & '==> Warning: vertical diffusive criteria exceeded on "kappa_m"'
     &,'           use a smaller "dtuv" and/or "kappa_m"              '
      endif
#endif
#ifdef implicitvmix
      write (stdout,9000)
     & '==> Warning: "implicitvmix" turns convective adjustment off.  '
     &,'              It will be handled by the vertical mixing       '
#else
# if defined oldconvect
      write (stdout,9000)
     & '==> Warning: the old convective scheme is enabled. Verify     '
     &,'              that "ncon" is set to what is wanted:           '
      write (stdout,'(32x,a,i3)') '"ncon" is currently set to ',ncon
# else
      write (stdout,9000)
     & '==> Warning: the full convective scheme is enabled.           '
     &,'             it will ignore "ncon" and remove all instability '
#  if defined fullconvect
      write (stdout,9000)
     & '==> Warning: the full convective scheme is the default.       '
     &,'            option -Dfullconvect is no longer needed. Remove it'
#  endif
# endif
#endif
#if defined tcvmix && !defined implicitvmix
      do k=1,km
        if ((dtts*dtxcel(k)*diff_cbt_limit)/dzt(k)**2 .ge. p25) then
          write (stdout,9000)
     & '==> Warning: vertical diffusive criteria exceeded for         '
     &,'            "diff_cbt_limit".  use a smaller "dtts", "dtxcel" '
     &,'            ,and/or "diff_cbt_limit"                          '
          write (stdout,'(a48,i3)') ' at level =',k
        endif
      enddo
      if ((dtuv*visc_cbu_limit)/dzmin**2 .ge. p25) then
        write (stdout,9000)
     & '==> Warning: vertical diffusive criteria exceeded for         '
     &,'             "visc_cbu_limit". use a smaller "dtuv" or        '
     &,'             "visc_cbu_limit"                                 '
      endif
#endif
c
c     check range of implicit factors
c
#if !defined implicitvmix && !defined redi_diffusion
      if (aidif .ne. c0) then
        write (stdout,9000)
     & '==> Warning: "aidif" was not = 0.0  but "implicitvmix" or     '
     &,'             "redi_diffusion" is not defined, so vertical     '
     &,'              mixing will be fully explicit.                  '
      endif
#endif
#if (defined velocity_horz_mix_const || defined am_cosine) && defined stream_function
c
c-----------------------------------------------------------------------
c     search for topographic instabilities (based  on the  work of
c     Peter Killworth  ...  eqn 11 from ocean modeling nov 1987)
c-----------------------------------------------------------------------
c
      num = 50
      do j=jscomp,jecomp
        dysq = dyt(j)**2
        do i=2,imtm1
          if (kmu(i+1,j-1) .ne. 0 .and. kmu(i+1,j) .ne. 0) then
            dxsq = (dxt(i)*cst(j))**2
# if defined am_cosine
            clix = am*csu(j)*dtuv/dxsq
# else
            clix = am*dtuv/dxsq
# endif
            h1   = zw(kmu(i+1,j-1))
            h2   = zw(kmu(i+1,j))
            hx   = (8.0*h1*h2/(h1+h2)**2 + dxsq/dysq)/(4.0 + dxsq/dysq)
            if (clix .ge. hx .and. num .ge. 0) then
              num = num - 1
              write (stdout,'(1x,a56,a31,i4,a1,i4,a1)')
     &        '==>   Warning: Killworth topographic roughness condition'
     &,       ' exceeded at location (i,jrow) = (',i+1,',',j,')'
              if (num .eq. 0) then
                write (stdout,9000)
     &         '==> Warning: msgs terminated after 50 cases were found '
              endif
            endif
          endif
        enddo
      enddo
#endif
c
c     verify that the domain boundary is valid
c
#if defined cyclic && defined solid_walls
        write (stdout,9000)
     & '==> Error: Both "cyclic" and "solid_walls cannot be enabled' 
        errorc = .true.
#endif
#if !defined cyclic && !defined solid_walls && !defined obc
        write (stdout,9000)
     & '==> Error: Either "cyclic", "solid_walls 
     &  or "obc" must be enabled'
        errorc = .true.
#endif
#ifdef solid_walls 
# if !defined obc_west && !defined obc_east
      jw = 0
      je = 0
      do j=jscomp,jecomp
        if (kmt(1,j) .ne. 0) jw = j
	if (kmt(imt,j) .ne. 0) je = j
      enddo
      if (je .ne. 0 .or. jw .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along   '
     &,'           the eastern or western boundary. This is not       '
     &,'           correct when enabling option "solid_walls"         '
        if (jw .ne. 0) write (stdout,*) '           at i=1 and j=',jw
        if (je .ne. 0) write (stdout,*) '           at i=imt and j=',je
      endif 
# endif
#endif
#if !defined obc_south && !defined obc_north
      in = 0
      is = 0
      do i=1,imt
        if (kmt(i,jstask) .ne. 0 .and. jstask .eq. 1) is = i
	if (kmt(i,jetask) .ne. 0 .and. jetask .eq. jmt) in = i
      enddo
      if (is .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along   '
     &,'           the southern boundary.                             '
        write (stdout,*) '           at j=1 and i=',is
      endif 
# ifndef symmetry
      if (in .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along   '
     &,'           the northern boundary.                             '
        write (stdout,*) '           at j=jmt and i=',in
      endif 
# endif
#endif
c
c     verify that each ocean point is at least 2 levels deep
c
      numk = 0
      do jrow=jscomp,jecomp
        do i=1,imt
	  if (kmt(i,jrow) .eq. 1) then
            numk = numk + 1
	    errorc = .true.
	    write (stdout,*)
     &	            ' Error: kmt(',i,',',jrow,') = 1 is not allowed'
	  endif
        enddo
      enddo
      if (numk .ne. 0) then	
        write (stdout,9000)
     & '==> Error: "kmt" must be at least 2 levels deep at all ocean  '
     &,'           points.                                            '
      endif
c
#if !defined quicker && defined ncar_upwind3
      write (stdout,9000)
     & '==> Error: "ncar_upwind3" only works when option              '
     &,'           option "quicker"  is enabled                       '
      errorc = .true.
#endif

#if defined save_fct && !defined tracer_advection_fct
      write (stdout,9000)
     & '==> Error:"save_fct" is enabled without "tracer_advection_fct"'
      errorc = .true.
#endif
#ifdef tracer_advection_fct
# if defined partial_cell
        write (stdout,9000)
     & '==> Warning: "partial_cell", has not been extended to work with'
     &,'             "tracer_advection_fct".                           '
c        errorc = .true.
# endif
# if  defined fourth_order_tracer_advection
        write (stdout,9000)
     & '==> Error: "tracer_advection_fct" is not compatible with      '
     &,'           "fourth_order_tracer_advection"                    '
     &,'           Can`t specify more than one tracer advection scheme'
        errorc = .true.
# endif
# if  defined quicker
        write (stdout,9000)
     & '==> Error: "tracer_advection_fct" is not compatible with      '
     &,'           "quicker"                                          '
     &,'           Can`t specify more than one tracer advection scheme'
        errorc = .true.
# endif
# if defined fct_dlm1 && defined fct_dlm2
      write (stdout,9000)
     & '==> Error: more than one fct delimiter type is specified      '
     &,'           specify either fct_dlm1 or fct_dlm2                '
     &,'           Default: using fct_dlm1                            '
        errorc = .true.
# endif
# if !defined fct_dlm1 && !defined fct_dlm2
      write (stdout,9000)
     & '==> Warning: fct delimiter type is not specified              '
     &,'             specify either fct_dlm1 or fct_dlm2              '
     &,'    Default: using fct_dlm1                                   '
# endif
#else
# if defined fct_dlm1 || defined fct_dlm2 || defined fct_3d
      write (stdout,9000)
     & '==> Warning: main fct option is not specified                 '
     &,'             minor fct options will be ignored                '
# endif
#endif
c
#if defined no_velocity_mixing_type
      write (stdout,9000)
     & '==> Error: Either option "velocity_horz_biharmonic" or option'
     &,'          "velocity_horz_laplacian" must be enabled.'
        errorc = .true.
#endif
#if defined no_tracer_mixing_type
      write (stdout,9000)
     & '==> Error: Either option "tracer_horz_biharmonic" or option'
     &,'          "tracer_horz_laplacian" must be enabled.'
        errorc = .true.
#endif
c
#if defined velocity_horz_biharmonic && defined second_order_window
      write (stdout,9000)
     & '==> Error: "fourth_order_window" should be enabled with option'
     &,'          "velocity_horz_biharmonic" in "derived_options.h" but'
     &,'          it is not. '
        errorc = .true.
#endif
#if defined tracer_horz_biharmonic && defined second_order_window
      write (stdout,9000)
     & '==> Error: "fourth_order_window" should be enabled with option'
     &,'           "tracer_horz_biharmonic" in "derived_options.h" but'
     &,'           it is not. '
        errorc = .true.
#endif
#if defined tracer_advection_fct && defined second_order_window
      write (stdout,9000)
     & '==> Error: "fourth_order_window" should be enabled with option'
     &,'           "tracer_advection_fct" in                '
     &,'           "derived_options.h" but it is not.       '
        errorc = .true.
#endif
#if defined fct
      write (stdout,9000)
     & '==> Error: "fct" has been renamed "tracer_advection_fct" '
     &, '   Change name in script to get flux corrected transport'
        errorc = .true.
#endif

#if defined quicker && defined second_order_window
      write (stdout,9000)
     & '==> Error: "fourth_order_window" should be enabled with option'
     &,'           "quicker". in "derived_options.h" but it is not.   '
        errorc = .true.
#endif
#if defined fourth_order_tracer_advection
# if defined second_order_window
      write (stdout,9000)
     & '==> Error: "fourth_order_window" should be enabled with option'
     &,'           "fourth_order_tracer_advection". in                '
     &,'            "derived_options.h" but it is not.                '
        errorc = .true.
# endif
#endif
#if defined partial_cell
# if defined term_balances && !defined max_window
        write (stdout,9000)
     & '==> Error:  "partial_cell", requires option "max_window" to   '
     &,'             work with option "term_balances"                 '
        errorc = .true.
# endif
# if defined save_xbts && !defined max_window
        write (stdout,9000)
     & '==> Error:  "partial_cell", requires option "max_window" to   '
     &,'             work with option "save_xbts"                     '
        errorc = .true.
# endif
#endif
#if defined meridional_overturning && defined meridional_overturning_basin
         write (stdout,9000)
     & '==> Error:  meridional_overturning is redundant with'
     &,' meridional_overturning_basin'
        errorc = .true.
# endif
#ifdef parallel_1d
        if (ALL(kmt(1:imt,jscomp:jecomp) .eq. 0)) then
            write( stdout,'(2a,3i6)' )
     &    '==> Warning: There are no ocean KMT points in'
     &,   '            this parallel domain. Problem domain pn,js,je='
     &,           pn, jscomp, jecomp
        endif
        if (ALL(kmu(1:imt,jscomp:jecomp) .eq. 0)) then
            write( stdout,'(2a,3i6)' )
     &    '==> Warning: There are no ocean KMU points in'
     &,   '            this parallel domain. Problem domain pn,js,je='
     &,           pn, jscomp, jecomp
        endif
# ifndef robert_time_filter
        write( stdout,'(2a)' )
     &       '==> Error: Option "parallel_1d" requires option',
     &       '"robert_time_filter".'
        errorc = .TRUE.
# endif
# if !defined use_libSMA && !defined use_libMPI
        write( stdout,'(2a)' )
     &       '==> Error: Option "parallel_1d" requires option',
     &       '"use_libSMA or use_libMPI".'
        errorc = .TRUE.
# endif
#endif
        
      write (stdout,'(/,20x,a,/)')
     &         ' E N D    C O N S I S T E N C Y    C H E C K S'
#ifdef parallel_1d
      call mpp_sync()
#endif
      if (errorc) then
        write (stdout,'(/a/a)')
     &  'An error flag has been detected in "checks"'
     &,' Use "grep -i error " on this printout file to find ALL errors.'
        call abort()
      endif
c
      return
9000  format (/,(1x,a))
      end
