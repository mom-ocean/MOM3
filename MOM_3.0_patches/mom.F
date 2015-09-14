      subroutine mom
c
c=======================================================================
c
c                        GFDL Modular Ocean Model
c
c     A primitive equation ocean model developed by researchers at  
c         the Geophysical Fluid Dynamics Laboratory /NOAA.
c                                
c     For details, refer to the manual.
c
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c    (mpp items)  v. balaji     e-mail=> vb@gfdl.gov
c=======================================================================
c
      use arrays_2d
      use datamod
      use memory_window

      use mom3_mpp_mod
#ifdef debug_parallel
      use chksum_mod
#endif
c
      external blkdta
c
#include "param.h"
#include "emode.h"
#include "iounit.h"
#ifdef obc
# include "cobc.h"
#endif
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "taskrows.h"
#include "tmngr.h"
#if defined max_window && defined parallel_1d
      real,dimension(imt*km,jmw,nt) :: tmp_var_t
      real,dimension(imt*km,jmw,2) :: tmp_var_u
      pointer (ptr_t,tmp_var_t)
      pointer (ptr_u,tmp_var_u)
#endif
      save cps, itt1st
c
#ifdef new_timing
      common /timers_vb / iitk(99)
#endif
      logical goto_2nd_pass
#ifdef timing
      call tic ('driver', 'ocean')
      call tic ('ocean', 'misc')
#endif
c     
c-----------------------------------------------------------------------
c     integrate equations for one time step
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c     update timestep counter, set corresponding model time, and set
c     time dependent logical switches which determine program flow.
c-----------------------------------------------------------------------
c
      itt = itt + 1
      call tmngr (dtts)
      call set_time_switches
c
c-----------------------------------------------------------------------
c     read system clock to calculate compute_time/grid_point/timestep 
c     (neglecting overhead of "setocn")
c-----------------------------------------------------------------------
c
      if (first) then
	itt1st = itt
        cps    = cpusec()
      endif
c
c-----------------------------------------------------------------------
c     initialize diagnostic variables
c-----------------------------------------------------------------------
c
#ifdef new_timing
      call system_clock(iit0)
#endif
      call diag_initial (jscomp, jecomp)
#ifdef new_timing
      call system_clock(iit1)
      iitk(5) = iitk(5) + iit1 - iit0
      iit0 = iit1
#endif
c
c-----------------------------------------------------------------------
c     adjust various quantities for leapfrog/mixing timesteps
c
c     leapfrog----------> h(tau+1) = h(tau-1) + 2dt*F(tau)
c                          
c     forward-----------> tau-1 <= tau
c                           h(tau+1) = h(tau-1) + dt*F(tau)
c
c     euler backward:     tau-1 <= tau
c       euler1----------> h(tau` ) = h(tau-1) + dt*F(tau)
c       euler2----------> h(tau+1) = h(tau-1) + dt*F(tau`)
c-----------------------------------------------------------------------
c
      if (leapfrog) then
c
c       normal leapfrog time step
c
        euler1  = .false.
	euler2  = .false.
	forward = .false.
	eots    = .true.
c
        c2dtts  = c2*dtts
        c2dtuv  = c2*dtuv
        c2dtsf  = c2*dtsf
      else
c
c       mixing time step (forward step or euler backward step)
c
        if (eb) then
          euler1  = .true.
	  euler2  = .false.
	  forward = .false.
	  eots    = .false.
	else
          euler1  = .false.
	  euler2  = .false.
	  forward = .true.
	  eots    = .true.
	endif
c
        c2dtts = dtts
        c2dtuv = dtuv
        c2dtsf = dtsf
c
#ifdef stream_function
        do jrow=1,jmt
          do i=1,imt
            psi(i,jrow,2) = psi(i,jrow,1)
          enddo
	enddo
#endif
      endif
c
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
        if (leapfrog) then
	  apgr = alph
        else
	  apgr = theta
        do jrow=jstask,jetask
            do i=1,imt
              ubarm1(i,jrow,1) = ubar(i,jrow,1)
	      ubarm1(i,jrow,2) = ubar(i,jrow,2)
            enddo
	  enddo
        endif
#endif
#ifdef explicit_free_surface
        if (.not.leapfrog) then
          do jrow=jstask,jetask
            do i=1,imt
              ubarm1(i,jrow,1) = ubar(i,jrow,1)
              ubarm1(i,jrow,2) = ubar(i,jrow,2)
            enddo
          enddo
        endif
#endif
c
#ifdef implicit_free_surface
      if (leapfrog) then
        if (acor .ne. c0) acor = alph
      else
        if (acor .ne. c0) acor = theta
      endif
#endif
c
c-----------------------------------------------------------------------
c     set time centering "gcor" for coriolis term
c-----------------------------------------------------------------------
c
      if (acor .eq. c0) then
        gcor = c1
      elseif (acor .ne. c0) then
	gcor = c0
      endif
c
#if defined implicit_free_surface
c
c     set "gcor" for leapfrog or mixing time steps
c
      if (leapfrog) then
	if (acor .eq. c0) then
	  gcor = c1
	else
	  gcor = gam
	endif
      else
        gcor = c0
      endif
#endif
c
c-----------------------------------------------------------------------
c     update pointers to tau-1, tau, & tau+1 data on disk.
c     for latitude rows they point to latdisk
c     for 2D fields they point to records on kflds
c-----------------------------------------------------------------------
c
      taum1disk = mod(itt+1,2) + 1
      taudisk   = mod(itt  ,2) + 1
      taup1disk = taum1disk
c
c-----------------------------------------------------------------------
c     update pointers (indices) to tau-1, tau, & tau+1 data in the MW
c-----------------------------------------------------------------------
c
#if defined max_window
c
c       rotate time levels instead of moving data
c
        taum1 = mod(itt+0,3) - 1
        tau   = mod(itt+1,3) - 1
        taup1 = mod(itt+2,3) - 1
#else
c
c       time levels are being held constant in time.
c
#endif
c
c=======================================================================
c
c               SOLVE THE BAROCLINIC AND TRACER EQUATIONS
c
c     Latitude rows are solved using a MW (memory window) approach.
c     Refer to the manual for  details. 
c     
c=======================================================================
c
      goto_2nd_pass = .true.
      do while (goto_2nd_pass)
c
#ifdef trace_indices
      if (leapfrog) write (stdout,'(/1x,a)')
     &   'Starting a --->Leapfrog timestep'
      if (forward) write (stdout,'(/1x,a)')
     &   'Starting a --->Forward timestep'
      if (euler1) write (stdout,'(/1x,a)')
     &   'Starting the first step of an ---> Euler backward timestep'
      if (euler2) write (stdout,'(/1x,a)')
     &   'Starting the second step of an --->Euler backward tstep'
      write (stdout,'(//1x,7(a,i4))')
     & "=> In mom: itt=",itt," disk indices: taum1disk =",taum1disk
     &," taudisk=",taudisk," taup1disk=",taup1disk
     &,", memory window indices: taum1=",taum1," tau=",tau
     &," taup1 ",taup1
#endif
#if !defined ramdrive
c
c-----------------------------------------------------------------------
c     establish a queue for "look ahead" reading "tau-1" and "tau" data
c-----------------------------------------------------------------------
c
      do jrow=jstask,jetask
	if (leapfrog .or. euler2) then
          call ofind (latdisk(taum1disk), nslab, jrow)
	endif
        call ofind (latdisk(taudisk), nslab, jrow)
      enddo
#endif
c
c-----------------------------------------------------------------------
c     Load the processor`s memory window in preparation for solving
c     tracer and baroclinic equations. (See the manual)
c-----------------------------------------------------------------------
c
#ifdef new_timing
      call system_clock(iit0)
#endif
      do mw=1,num_loads
c
        ncrows = window_size-2*jbuf
	joff = jscomp - jbuf - 1 + (mw-1)*ncrows
        js = (jstask - joff)*(1/mw) + (1 + num_rows_to_copy)*(1-1/mw)
	je = min(window_size,jetask-joff)
c
	is = 2
	ie = imt - 1
c
c       set rows  within the MW for solving prognostic equations
c       (rows are clipped at jrow=2 and jrow=jmt-1)
c
        js_calc = 1+jbuf
        je_calc = ncrows+jbuf
c
c-----------------------------------------------------------------------
c       load prognostic and related variables into the MW
c-----------------------------------------------------------------------
c
	call loadmw (joff, js, je, is, ie, latdisk(taum1disk)
     &,                latdisk(taudisk))
c
c-----------------------------------------------------------------------
c       set vertical boundary conditions for momentum and tracers
c-----------------------------------------------------------------------
c
	call setvbc (joff, js, je, is, ie)
c
#ifdef air_pressure
c-----------------------------------------------------------------------
c       add air pressure
c-----------------------------------------------------------------------
c
        call add_air_pressure (js, je, is, ie)
c
#endif
c-----------------------------------------------------------------------
c       compute pressure gradients
c-----------------------------------------------------------------------
c
        call pressure_gradients (joff, js, je)
c
c-----------------------------------------------------------------------
c       calculate advection velocities for T cells and U cells     
c-----------------------------------------------------------------------
c
        call adv_vel (joff, js, je, is, ie)
c
c-----------------------------------------------------------------------
c       set horizontal mixing coefficients for momentum and tracers
c-----------------------------------------------------------------------
c
        call hmixc (joff, js, je, is, ie)
c
#ifdef isoneutralmix
c
c-----------------------------------------------------------------------
c       calculate isopycnal diffusion tensor components (and
c       gent_mcwilliams advective velocities) for use with tracers
c-----------------------------------------------------------------------
c
        call isopyc (joff, js, je, is, ie)
#endif
c
c-----------------------------------------------------------------------
c       set vertical mixing coefficients for momentum and tracers
c-----------------------------------------------------------------------
c
        call vmixc (joff, js, je, is, ie)
c
#if defined tracer_horz_biharmonic
c
c-----------------------------------------------------------------------
c       calculate del**2 quantities for biharmonic tracer mixing
c-----------------------------------------------------------------------
c
        call delsq_tracer (joff, js, je, is, ie)
#endif
c
#if defined velocity_horz_biharmonic
c
c-----------------------------------------------------------------------
c       calculate del**2 quantities for biharmonic velocity mixing
c-----------------------------------------------------------------------
c
        call delsq_velocity (joff, js, je, is, ie)
#endif
c
c
#if defined obc_south || defined obc_north
# include "mom_obc.inc"
#else 
c
c-----------------------------------------------------------------------
c       compute tracers and internal mode velocities
c-----------------------------------------------------------------------
c
        call tracer     (joff, js_calc, je_calc, is, ie)
        call baroclinic (joff, js_calc, je_calc, is, ie)
#endif
c
c-----------------------------------------------------------------------
c       calculate diagnostics
c-----------------------------------------------------------------------
c
#ifdef new_timing
        call system_clock(iit2)
#endif
        call diag_mw_inside (joff, js_calc, je_calc, is, ie)
#ifdef new_timing
        call system_clock(iit3)
        iitk(5) = iitk(5) + iit3 - iit2
#endif
c
#if defined robert_time_filter
c
c-----------------------------------------------------------------------
c       robert time filter for velocity and tracers
c-----------------------------------------------------------------------
c
        call r_time_filter (joff, js_calc, je_calc, is, ie
     &,                     latdisk(taudisk))
#endif
c
c-----------------------------------------------------------------------
c       write prognostic variables from the MW to disk "tau+1"
c-----------------------------------------------------------------------
c
#if defined max_window
c
c         do nothing since variables are already in "tau+1" MW slots
c
#else
          call putmw (joff, js_calc, je_calc, latdisk(taup1disk))
#endif
c
      enddo
#if defined parallel_1d
c
c-----------------------------------------------------------------------
c     update boundary latitudes when using more than one processor
c-----------------------------------------------------------------------
c
!     update tau and tau+1 disk
# ifdef max_window
      ptr_t = loc(t(1,1,1,1,tau))
      call mpp_update_domains( tmp_var_t, mom_prognostic_domain )
      ptr_t = loc(t(1,1,1,1,taup1))
      call mpp_update_domains( tmp_var_t, mom_prognostic_domain )
      ptr_u = loc(u(1,1,1,1,tau))
      call mpp_update_domains( tmp_var_u, mom_prognostic_domain )
      ptr_u = loc(u(1,1,1,1,taup1))
      call mpp_update_domains( tmp_var_u, mom_prognostic_domain )
# endif
# ifdef ramdrive
      call mpp_update_domains( ramdisk3d, mom_ramdisk_domain )
# endif
#endif
#ifdef new_timing
      call system_clock(iit1)
      iitk(3) = iitk(3) + iit1 - iit0
      iit0 = iit1
#endif
#if defined debug_parallel
#ifdef max_window
      do jrow=jstask,jetask
	joff = jscomp - jbuf - 1
	j = jrow - joff
        write( stdout,'(a,2i6,4z18)' )
     &   'MOM chksum(u,t) on itt, jrow=', itt, jrow,
     &    checksum (u(1,1,j,1,taup1), imt, km),
     &    checksum (u(1,1,j,2,taup1), imt, km),
     &    checksum (t(1,1,j,1,taup1), imt, km),
     &    checksum (t(1,1,j,2,taup1), imt, km)
      enddo
#else
      do jrow=jstask,jetask
        write( stdout,'(a,2i6,2z18)' )
     &  'MOM chksum(ramdisk3d) on itt, jrow=', itt, jrow
     &,  chksum( ramdisk3d(1,jrow,1),size(ramdisk3d,1) )
     &,  chksum( ramdisk3d(1,jrow,2),size(ramdisk3d,1) )
      enddo
#endif
#endif
c
c=======================================================================
c
c               SOLVE THE BAROTROPIC EQUATION
c
c=======================================================================
c
#ifdef new_timing
      call system_clock(iit0)
#endif
      call barotropic (cori, itt, jstask, jetask)
#ifdef new_timing
      call system_clock(iit1)
      iitk(4) = iitk(4) + iit1 - iit0
      iit0 = iit1
#endif
c
c-----------------------------------------------------------------------
c     if this is the 1st pass of an euler backward timestep, set the
c     disk pointers so the proper time levels are read on the 2nd pass 
c     and go back to do the 2nd pass.
c-----------------------------------------------------------------------
c
      if (euler1) then
        eots      = .true.
	euler1    = .false.
	euler2    = .true.
        ntaux     = taum1disk
        taum1disk = taudisk
        taudisk   = taup1disk
        taup1disk = ntaux
      else
        goto_2nd_pass = .false.
      endif
c
      enddo
c
#if defined max_window
      if (euler2) then
c
c       shuffle "tau" and "tau+1" after euler backward to
c       insure data is in the right place for the next timestep
c
	call euler_shuffle
c
c       re-establish correct pointers for this timestep
c
        taum1 = mod(itt+0,3) - 1
        tau   = mod(itt+1,3) - 1
        taup1 = mod(itt+2,3) - 1
      endif
#endif
c
c-----------------------------------------------------------------------
c     output all remaining diagnostics
c-----------------------------------------------------------------------
c
#ifdef new_timing
      call system_clock(iit0)
#endif
      call diag_mw_outside(jscomp, jecomp)

#ifdef new_timing
      call system_clock(iit1)
      iitk(5) = iitk(5) + iit1 - iit0
      iit0 = iit1
#endif
c
c-----------------------------------------------------------------------
c     if it`s the last timestep then clean things up otherwise return
c-----------------------------------------------------------------------
c
      if (eorun) then
c
c       calculate cp_time/grid_point/timestep (neglecting setocn)
c
        ncalls = itt - itt1st + 1
        cpe    = cpusec()
        write (stdout,'(/,1x, a, f9.3, a,1x,a,1pg14.7,a)')
     &    'Timing Summary: ', cpe-cps,  ' cp seconds,'
     &,    '(',(cpe-cps)/(ncalls*imt*jmt*km)
     &,    ' cp seconds per grid point per time step)'
c
c       write a restart
c
        if (restrt) call write_restart (jstask, jetask)
c
        write (stdout,'(1x,a)') 'MOMdone'
        call oclose (kflds)
	call oclose (latdisk(1))
	call oclose (latdisk(2))
      endif
c
#ifdef timing
      call toc ('ocean', 'misc')
      call toc ('driver', 'ocean')
#endif
c
      return
      end


      function cpusec ()
c
c     user supplied function to measure cpu time in seconds
c
      real(kind=4) :: second
      cpusec = second()
      return
      end


