      subroutine set_time_switches
c
c=======================================================================
c
c               S E T   E V E N T   S W I T C H E S
c
c     Subroutine set_time_switches sets switches located in "switch.h"
c     that trigger periodically recurring events in the model
c     such as diagnostics and end-of-run. All switches are set either
c     .true. or .false. every time step.
c
c     input:
c     1) the time manager module (time_manager_mod)
c     2) switch.h file containing all switch intervals in units of days.
c     3) the following time structures: dt_time, model_time,
c        time_since_base
c
c     output:
c
c     a logical switch corresponding to each switch interval plus
c     additional logical switches for end-of-month, end-of-run, etc.
c
c     author:       R.C. Pacanowski     e-mail=>  rcp@gfdl.gov
c=======================================================================
c
      use time_manager_mod
      type(time_type) ::  target_time, time
c
#include "derived_options.h"
#include "stdunits.h"
#include "switch.h"
#include "tmngr.h"
      integer days, months, years, seconds, minutes, hours
c
c-----------------------------------------------------------------------
c     save previous stamp and generate a new stamp
c-----------------------------------------------------------------------
c
      pstamp = stamp
c
      call get_date (model_time, years, months, days
     &,                          hours, minutes, seconds)
      call set_stamp (stamp, years, months, days
     &,                          hours, minutes, seconds)
# if defined separate_diag_files
      xstamp = stamp
# else
      xstamp = stop_stamp
# endif
c
c-----------------------------------------------------------------------
c     compute relative time since Initial Condition
c-----------------------------------------------------------------------
c
      call get_time (model_time - ic_time, isec, iday)
      relday = iday + isec/86400.
      prelyr = relyr
      relyr  = relday/yrlen
c
c-----------------------------------------------------------------------
c     is it the first time step of a run ?
c-----------------------------------------------------------------------
c
      if (model_time == start_time) then
        first = .true.
      else
        first = .false.
      endif
c
c-----------------------------------------------------------------------
c     is it within 1/2 time step of the end of a day ?
c-----------------------------------------------------------------------
c
      eoday = set_switch (1.0, model_time, dt_time)
c
c-----------------------------------------------------------------------
c     is it within 1/2 time step of the middle of a month ?
c-----------------------------------------------------------------------
c
      length = days_in_month(model_time)
      time = set_time(0, length)/2
      target_time = set_date(years, months, 1) + time
      call get_time (target_time, isec, iday)
      dayint = iday + isec/86400.0
      midmon = set_switch (dayint, model_time, dt_time)
c
c-----------------------------------------------------------------------
c     is it within 1/2 time step of the end of the month ?
c-----------------------------------------------------------------------
c
      if (days .gt. 15) then
        length = days_in_month(model_time)
        target_time = set_date(years, months, 1) + set_time(0, length)
      else
        target_time = set_date(years, months, 1)
      endif
      call get_time (target_time, isec, iday)
c
c     if target_time = 1/1/1900 then iday=0. disallow this case 
c
      if (iday .eq. 0 .and. isec .eq. 0) iday = 30
      dayint = iday + isec/86400.
      eomon = set_switch (dayint, model_time, dt_time)
      if (midmon .and. eomon) midmon = .false.
c
c-----------------------------------------------------------------------
c     is it within 1/2 time step of the end of the year ?
c-----------------------------------------------------------------------
c
      time = set_date(years, 12, 1)
      if (months .gt. 6) then
        length = days_in_month(time)
        target_time = time + set_time(0, length)
      else
        target_time = time
      endif
      call get_time (target_time, isec, iday)
      dayint = iday + isec/86400.
      eoyear = set_switch (dayint, model_time, dt_time)
c
c-----------------------------------------------------------------------
c     is it the last time step of the run ?
c-----------------------------------------------------------------------
c
      if (model_time .ge. stop_time) then
        eorun = .true.
      else
        eorun = .false.
      endif
c
c-----------------------------------------------------------------------
c     is time mixing to be done now ?
c-----------------------------------------------------------------------
c
#if defined robert_time_filter
      leapfrog = .true.
#else
      if (nmix .eq. 0 .or. nmix .eq. 1) then
        leapfrog = .true.
      else
        if (mod(itt,nmix) .eq. 1) then
          leapfrog = .false.
        else
          leapfrog = .true.
        endif
      endif
#endif
c
#ifdef tracer_averages
c
c-----------------------------------------------------------------------
c     are regional tracer averages to be done at this time ?
c-----------------------------------------------------------------------
c
      tavgts = set_switch (tavgint, time_since_base, dt_time)
#endif
c
c
#if defined cross_flow
c
c-----------------------------------------------------------------------
c     is saving of diapycnal and isopycnal components of flow to be done
c     at this time ?
c-----------------------------------------------------------------------
c
      crossts = set_switch (crossint, time_since_base, dt_time)
#endif
#if defined save_fct
c
c-----------------------------------------------------------------------
c     is saving of difference between FCT and leapfrog advection to be
c     done at this time ?
c-----------------------------------------------------------------------
c
      fctts = set_switch (fctint, time_since_base, dt_time)
#endif
#if defined save_pressure 
c
c-----------------------------------------------------------------------
c     is saving of pressure to be done at this time ?
c-----------------------------------------------------------------------
c
      pressts = set_switch (pressint, time_since_base, dt_time)
#endif
#if defined save_mixing_coeff
c
c-----------------------------------------------------------------------
c     is saving of estimated mixing coefficients to be done at
c     this time ?
c-----------------------------------------------------------------------
c
      cmixts = set_switch (cmixint, time_since_base, dt_time)
#endif
#ifdef save_convection
c
c-----------------------------------------------------------------------
c     is saving of temperature change due to explicit convection
c     to be done at this time ?
c-----------------------------------------------------------------------
c
      exconvts = set_switch (exconvint, time_since_base, dt_time)
#endif
#ifdef energy_analysis
c
c-----------------------------------------------------------------------
c     are global energetics to be done at this time ?
c-----------------------------------------------------------------------
c
      glents = set_switch (glenint, time_since_base, dt_time)
#endif
c
#ifdef term_balances
c
c-----------------------------------------------------------------------
c     are tracer & momentum term balances to be done at this time ?
c-----------------------------------------------------------------------
c
      trmbts = set_switch (trmbint, time_since_base, dt_time)
#endif
c
#ifdef matrix_sections
c
c-----------------------------------------------------------------------
c     are x-z printouts to be done at this time ?
c-----------------------------------------------------------------------
c
      prxzts =  set_switch (prxzint, time_since_base, dt_time)
#endif
c
#ifdef show_external_mode
c
c-----------------------------------------------------------------------
c     is the external mode printout to be done at this time ?
c-----------------------------------------------------------------------
c
      extts = set_switch (extint, time_since_base, dt_time)
#endif
c
#ifdef tracer_yz
c
c-----------------------------------------------------------------------
c     is time to write the zonally averaged tracer components?
c-----------------------------------------------------------------------
c
      tyzts = set_switch (tyzint, time_since_base, dt_time)
#endif
c
#if defined meridional_overturning || defined meridional_overturning_basin
c
c-----------------------------------------------------------------------
c     is time to write the vertical_meridional stream function?
c-----------------------------------------------------------------------
c
      vmsfts = set_switch (vmsfint, time_since_base, dt_time)
#endif
c
#ifdef gyre_components
c
c-----------------------------------------------------------------------
c     are the gyre components to be done at this time ?
c-----------------------------------------------------------------------
c
      gyrets = set_switch (gyreint, time_since_base, dt_time)
#endif
c
#ifdef time_step_monitor
c
c-----------------------------------------------------------------------
c     is printing of time step integrals to be done at this time ?
c-----------------------------------------------------------------------
c
      tsits = set_switch (tsiint, time_since_base, dt_time)
#endif
c
#ifdef snapshots
c
c-----------------------------------------------------------------------
c     is time to write a snapshot of the MOM data ?
c-----------------------------------------------------------------------
c
      snapts = set_switch (snapint, time_since_base, dt_time)
#endif
c
#ifdef save_density_terms
c
c-----------------------------------------------------------------------
c     is time to write out local potential density terms ?
c-----------------------------------------------------------------------
c
      rhots = set_switch (rhoint, time_since_base, dt_time)
#endif
c
#if defined hl_diffusivity || defined vmhs_diffusivity
c
c-----------------------------------------------------------------------
c     Time to compute a new nonconstant diffusivity. 
c     Determined as a multiple of itt
c-----------------------------------------------------------------------
c
       Ndiffint = 1 + nint(diffint)
       diffts   = mod(itt,Ndiffint) .eq. 0
#endif
c
#if defined save_ncdiff
c
c-----------------------------------------------------------------------
c     Time to map arrays for nonconstant diffusivities
c-----------------------------------------------------------------------
c
      diffsnapts = set_switch (diffsnapint, time_since_base, dt_time)
#endif
c
#ifdef trajectories
c
c-----------------------------------------------------------------------
c     is it time to write particle trajectories?
c     they must be written at the end of the run for restarting.
c-----------------------------------------------------------------------
c
      trajts = set_switch (trajint, time_since_base, dt_time)
      trajts = trajts .or. eorun
#endif
c
#ifdef stability_tests
c
c-----------------------------------------------------------------------
c     is it time to do a stability diagnosis ?
c-----------------------------------------------------------------------
c
      stabts = set_switch (stabint, time_since_base, dt_time)
#endif
c
#ifdef show_zonal_mean_of_sbc
c
c-----------------------------------------------------------------------
c      is it time to do a boundary condition diagnosis ?
c-----------------------------------------------------------------------
c
      zmbcts = set_switch (zmbcint, time_since_base, dt_time)
#endif
c
#ifdef diagnostic_surf_height
c
c-----------------------------------------------------------------------
c     is it time to accumulate the forcing or "average and write" the
c      diagnostic surface height?
c-----------------------------------------------------------------------
c
# if defined monthly_averages
      dspts    = eomon .or. eorun
      dspperts = .true.
# else
      dspts = set_switch (dspint, time_since_base, dt_time)
      idays = int(dspper)
      iseconds = (dspper - idays)*86400
      time = set_time (iseconds, idays)
      dspperts = set_switch (dspint, time_since_base + time/2, time)
# endif
#endif
c
#ifdef time_averages
c
c-----------------------------------------------------------------------
c     is it time to accumulate or "average and write" time mean data?
c-----------------------------------------------------------------------
c
# if defined monthly_averages
      timavgts    = eomon .or. eorun
      timavgperts = .true.
# else
      timavgts = set_switch (timavgint, time_since_base, dt_time)
      idays = int(timavgper)
      iseconds = (timavgper - idays)*86400
      time = set_time (iseconds, idays)
      timavgperts = set_switch (timavgint, time_since_base + time/2
     &,                         time)
# endif
#endif
c
#ifdef save_xbts
c
c-----------------------------------------------------------------------
c     is it time to accumulate or "average and write" XBT data?
c-----------------------------------------------------------------------
c
# if defined monthly_averages
      xbtts    = eomon .or. eorun
      xbtperts = .true.
# else
      xbtts = set_switch (xbtint, time_since_base, dt_time)
      idays = int(xbtper)
      iseconds = (xbtper - idays)*86400
      time = set_time (iseconds, idays)
      xbtperts = set_switch (xbtint, time_since_base + time/2, time)
# endif
#endif
c
#ifdef meridional_tracer_budget
c
c-----------------------------------------------------------------------
c     is it time to accumulate or "average and write" the meridional 
c     tracer budget?
c-----------------------------------------------------------------------
c
# if defined monthly_averages
      tmbts    = eomon .or. eorun
      tmbperts = .true.
# else
      tmbts = set_switch (tmbint, time_since_base, dt_time)
      idays = int(tmbper)
      iseconds = (tmbper - idays)*86400
      time = set_time (iseconds, idays)
      tmbperts = set_switch (tmbint, time_since_base + time/2, time)
# endif
#endif
c
      return
      end

