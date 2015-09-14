#include "derived_options.h"
c
c=======================================================================
c
c                           Timing Module
c 
c     This module contains as set of general purpose timing routines
c     which are useful in timing sections of code.
c
c     If compiled with the option -Dtiming, the timers are active.
c     If -Dtiming is not selected, all subroutines are null routines
c     and calls to them need not be removed from source code.
c
c     Timers are organized by group name and subgroup name, both
c     character strings.  For each group timer and subtimers, summary
c     statistics are printed showing total time, number of calls,
c     average time per call, percent of group time, and percent of total
c     time that this timer was active.
c
c     A timer is turned on by a call to the subroutine "tic" and turned
c     off by a call to the subroutine "toc".
c
c          call tic ('group name', 'subgroup name')
c          ...code to be timed....
c          call toc ('group name', 'subgroup name')
c
c     A subroutine "ticr" is provided to resume a previous call to "tic"
c     without increasing the number of calls statistic.  This is useful
c     for timing a part of the code that is separated by some untimed
c     code.
c
c          call tic ('group name', 'subgroup name')
c          ...first part of code....
c          call toc ('group name', 'subgroup name')
c          ...untimed code...
c          call ticr ('group name', 'subgroup name')
c          ...second part of code....
c          call toc ('group name', 'subgroup name')
c
c     Times are cpu times from a call to "tic" to its matched call to
c     "toc".  The function "second" provides the raw cpu times.  If a
c     subroutine call occurs in a section of code being timed, the
c     entire execution time of the subroutine will be included (unless
c     it contains a call to "toc" that timer).
c
c     Any number of timers from different groups may be simultaneously
c     active, and they may be "tic"-ed or "toc"-ed independently.
c     Within a group, at most one timer is active at any time.  If a
c     second timer from the group is "tic"-ed, the first is 
c     automatically "toc"-ed and pushed on a stack for that group.
c     The "tic" time for the second timer is precisely the "toc" time
c     for the first. When the second timer is "toc"-ed, the first timer
c     is popped from the stack and resumed at precisely the time the
c     second is "toc"-ed.
c
c          subroutine sub_name
c            call tic ('sub_name', 'miscellaneous')
c            ...miscellaneous subroutine code...
c            call tic ('sub_name', 'subgroup a')
c              call sub_a  (x, y, z)
c            call toc ('sub_name', 'subgroup a')
c            ...more miscellaneous code...
c            call tic ('sub_name', 'subgroup b')
c              call sub_b  (u, v, w)
c            call toc ('sub_name', 'subgroup b')
c            ...still more miscellaneous code...
c            call toc ('sub_name', 'miscellaneous')
c            return
c            end
c
c   The ('sub_name', 'miscellaneous') timer captures all the time of the
c   subroutine sub_name not specifically timed by the timers('sub_name',
c   'subgroup a') and ('sub_name', 'subgroup b').
c
c   WARNINGS are issued for improper nesting of "tic" and "toc" up to
c   a maximum of max_timer_warnings (in timer.h).  Many improper nesting
c   problems are handled gracefully and the timings may still be useful.
c   Individual WARNING messages indicate the corrective action taken in
c   each case.
c
c   Total time for the run is timed from the call to "inittimers",
c   normally one of the first executable statements of the run, until
c   the call to "showtimers", normally one of the last executable
c   statements of the run.
c
c   Accuracy of the timings depends on the precision of the Fortran
c   readable clock accessed by the function "second".  On machines where
c   the clock returns only integer seconds, the preprocessor (cpp)
c   variable "real_seconds" is turned off, and raw, unadjusted timings
c   that include the overhead of the calls to "tic" and "toc" are
c   printed by "showtimers".  The overhead is relatively small for
c   unnested "tic"-"toc" pairs or innermost nested "tic"-"toc" pairs,
c   but will be somewhat higher for timed that code that itself contains
c   calls to "tic" and "toc".  Also bear in mind that the timing for
c   each individual "tic"-"toc" pair may be in error by a full second
c   because of the limited precision of the clock.  On such machines,
c   it may be necessary to artificially repeat the calculation between
c   "tic" and "toc" to get times long enough that the uncertainty of
c   one second is not significant.
c
c   If greater precision is available to the function "second", the
c   preprocessor variable "real_seconds" is defined, and times shown by
c   "showtimers" are adjusted times, with elaborate compensations made
c   for the cpu time spent in various portions of "tic" and "toc".
c   Calibration loops are executed and timed in both "inittimers" and
c   in "showtimers" to compensate both for combinations of timers within
c   timers and for the increased search time as the list of timer group
c   names and subgroup names grows.  Although in heavily instrumented
c   code, the time spent in "tic" and "toc" may be a third of the total
c   cpu time, experience on the Cray ymp has shown these adjusted times
c   to be quite reliable and accurate.  For your information, the time
c   per "tic"-"toc" pair is also printed by "showtimers", as are
c   warnings about "tic"-ed timers still not "toc"-ed by the end of the
c   run.
c
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
c


#ifdef test_timer
# include "util.F"
      program driver
c
c=======================================================================
c
c     This driver simulates a time step loop for solving a 
c     tracer equation in various  ways.
c     
c
c     To run the driver, execute the "run_timer" script
c
c     author:         r. c. pacanowski    e-mail=> rcp@gfdl.gov
c=======================================================================
# include "stdunits.h"
c
c     allowable range for jmw is 3 to jmt
c
      parameter (imt=92, jmt=60, km=15, jmw=3, nt=2, kmm1=km-1)
      parameter (is=2, ie=imt-1)
      parameter (ncrows = jmw - 2)
      parameter (nmwlat=(jmt-2)/ncrows)
      parameter (maxmw = nmwlat + (jmt-3)/(ncrows*nmwlat))
      common/a/ t(imt,km,jmw), diff_fe(imt,km,jmw), diff_fn(imt,km,jmw)
      common/a/ t2(imt,km,jmw)
      common/a/ adv_fe(imt,km,jmw), adv_fn(imt,km,jmw)
      common/a/ adv_fb(imt,km,jmw), diff_fb(imt,km,jmw)
      common/a/ adv_vet(imt,km,jmw), adv_vnt(imt,km,jmw)
      common/a/ adv_vbt(imt,km,jmw)
      common/a/ adv_vet1(imt,km,jmw), adv_vnt1(imt,km,jmw)
      common/a/ adv_vet2(imt,km,jmw), adv_vnt2(imt,km,jmw)
      common/a/ ta(imt,km,jmw),tb(imt,km,jmw) 
      common/a/ tmask(imt,km,jmw), diff2_n(jmw)
      common/a/ diff2_e(imt,jmw), tdt(km), dztr(km), diff2_b(km)
      common/a/ cstdxtr(imt,jmw), cstdytr(jmw), dzt2r(km)
      common/a/ cstdxt2r(imt,jmw), cstdyt2r(jmw)
      common/a/ diff4_w(imt,jmw)
      common/a/ diff4_e(imt,jmw)
      common/a/ diff4_n(jmt), diff4_s(jmt)
c
c-----------------------------------------------------------------------
c     various forms for  diffusion
c-----------------------------------------------------------------------
c
c     flux form using temp storage for diffusive fluxes
c
      DIFF_Tx(i,k,j) = (diff_fe(i,  k,j)*tmask(i+1,k,j)
     &                - diff_fe(i-1,k,j)*tmask(i-1,k,j))*cstdxtr(i,j)
      DIFF_Ty(i,k,j) = (diff_fn(i,k,j  )*tmask(i,k,j+1)
     &                - diff_fn(i,k,j-1)*tmask(i,k,j-1))*cstdytr(j)
      DIFF_Tz(i,k,j) = (diff_fb(i,k-1,j) - diff_fb(i,k,j))*dztr(k)
c
c     no temp storage. diffusive fluxes are computed redundantly
c
      T2xx(i,k,j) =
     &       (diff2_e(i,  j)*(t(i+1,k,j) - t(i,  k,j))*tmask(i+1,k,j)
     &      - diff2_e(i-1,j)*(t(i,  k,j) - t(i-1,k,j))*tmask(i-1,k,j)
     &       )*cstdxtr(i,j)
      T2yy(i,k,j) = 
     &         (diff2_n(j  )*(t(i,k,j+1) - t(i,k,j  ))*tmask(i,k,j+1)
     &        - diff2_n(j-1)*(t(i,k,j  ) - t(i,k,j-1))*tmask(i,k,j-1)
     &                                       )*cstdytr(j)
c
c     no temp storage. diffusive fluxes are redundant but optimized
c
      T4xx(i,k,j) =
     &        diff4_e(i,j)*(t(i+1,k,j) - t(i,  k,j))*tmask(i+1,k,j)
     &      - diff4_w(i,j)*(t(i,  k,j) - t(i-1,k,j))*tmask(i-1,k,j)
      T4yy(i,k,j) = 
     &         diff4_n(j)*(t(i,k,j+1) - t(i,k,j  ))*tmask(i,k,j+1)
     &       - diff4_s(j)*(t(i,k,j  ) - t(i,k,j-1))*tmask(i,k,j-1)
c
c-----------------------------------------------------------------------
c     various forms for advection
c-----------------------------------------------------------------------
c
c     flux form of advection using temp storage for advective fluxes
c
      ADV_Tx(i,k,j) = (adv_fe(i,k,j) - adv_fe(i-1,k,j))*cstdxt2r(i,j)
      ADV_Ty(i,k,j) = (adv_fn(i,k,j) - adv_fn(i,k,j-1))*cstdyt2r(j)
      ADV_Tz(i,k,j) = (adv_fb(i,k-1,j) - adv_fb(i,k,j))*dzt2r(k)
c
c     no temp storage. advective fluxes are computed redundantly
c
      UT2x(i,k,j) = (adv_vet(i,k,j)*(t(i+1,k,j) + t(i,k,j))
     &            -  adv_vet(i-1,k,j)*(t(i,k,j) + t(i-1,k,j))
     &              )*cstdxt2r(i,j)
      VT2y(i,k,j) = (adv_vnt(i,k,j)*(t(i,k,j+1) + t(i,k,j))
     &             - adv_vnt(i,k,j-1)*(t(i,k,j) + t(i,k,j-1))
     &              )*cstdyt2r(j)
c
c     no temp storage. advective fluxes are redundant but optimized
c
      UT4x(i,k,j) = adv_vet1(i,k,j)*(t(i+1,k,j) + t(i,k,j))
     &           -  adv_vet2(i,k,j)*(t(i,k,j) + t(i-1,k,j))
      VT4y(i,k,j) = adv_vnt1(i,k,j)*(t(i,k,j+1) + t(i,k,j))
     &            - adv_vnt2(i,k,j)*(t(i,k,j) + t(i,k,j-1))
c
      stdout = 6
c
      write (stdout,*) ' '
      write (stdout,*) 
     &' This driver gives timing estimates for various forms of'
      write (stdout,*) 
     &' the tracer equation for memory window size "jmw"'
      write (stdout,*) 
     &' (Note: "jmw" can be changed in the parameter statement)'
      write (stdout,*) ' '
      write (stdout,*) 
     &' The purpose is to show computational efficiency of'
      write (stdout,*) 
     &' these forms as a function of window size on'
      write (stdout,*) ' a variety of computer architectures'
      write (stdout,*) ' '
      write (stdout,*) 
     &' MOM uses the optimum form for the CRAY  taking'
      write (stdout,*) 
     &' into account various window sizes. This may or may'
      write (stdout,*) 
     &' not be optimum on other computer architectures'
      write (stdout,*) ' '
      write (stdout,*) 
     &' No one form may be optimum on all architectures for all'
      write (stdout,*) 
     &' memory window sizes and it may be difficult to predict'
      write (stdout,*) 
     &' which will be fastest on a given architecture.'
      write (stdout,*) 
     &' This is so because it involves compiler optimization,'
      write (stdout,*) 
     &' cache size, vector/scalar speeds, bus speed, etc.' 
      write (stdout,*) 
     &' By inspecting these timing results, the optimum form'
      write (stdout,*) 
     &' for a specific architecture can be identified. It is'
      write (stdout,*) 
     &' then left to the user to decide if changing the form'
      write (stdout,*) 
     &' in MOM to suit his configuration is worthwhile.'
      write (stdout,*) ' '
      write (stdout,*) 
     &' Note that this code only approximates the solving of'
      write (stdout,*) 
     &' tracer equations and is only intended as a guide'
c
      write (stdout,*) '  "jmw" =', jmw
c
      call inittimers
c
c-----------------------------------------------------------------------
c     initialize the variables with garbage
c-----------------------------------------------------------------------
c
      do k=1,km
        tdt(k) = 1.e-20
	dztr(k) = sin(k*2.0)
	dzt2r(k) = 0.5*sin(k*2.0)
	diff2_b(k) = cos(dztr(k))
      enddo
      do j=1,jmt
        diff4_n(j) = sin(j*2.0)
	diff4_s(j) = cos(j*3.1)
      enddo
      do j=1,jmw
	diff2_n(j) = cos(sin(j*0.1))
	cstdytr(j) = cos(6.123*j) 
	cstdyt2r(j) = 0.5*cos(6.123*j) 
        do i=1,imt
	  diff2_e(i,j) = sin(i*j*0.2)
	  cstdxtr(i,j) = sin(i*j*0.111)
	  cstdxt2r(i,j) = 0.5*sin(i*j*0.111)
          diff4_w(i,j) = sin(j*3.4)
	  diff4_e(i,j) = cos(j*5.6)
	enddo
	do k=1,km
	  do i=1,imt
	    t(i,k,j) = sin(0.123*i*j*k)
	    ta(i,k,j) = 0.0
	    tb(i,k,j) = 0.0
	    tmask(i,k,j) = 1.0
	    if (mod(i*j*k,5) .eq. 0) tmask(i,k,j) = 0.0
	    diff_fe(i,k,j) = 0.0
	    diff_fn(i,k,j) = 0.0
	    adv_fe(i,k,j) = 0.0
	    adv_fn(i,k,j) = 0.0
	    adv_fb(i,k,j) = 0.0
	    diff_fb(i,k,j) = 0.0
	    adv_vet(i,k,j) = sin(0.1*i*j*k)
	    adv_vnt(i,k,j) = cos(0.12*i*j*k)
	    adv_vbt(i,k,j) = cos(0.123*i*j*k)
	    adv_vet1(i,k,j) = sin(0.144*i*j*k)
	    adv_vnt1(i,k,j) = cos(0.1244*i*j*k)
	    adv_vet2(i,k,j) = sin(0.155*i*j*k)
	    adv_vnt2(i,k,j) = cos(0.166*i*j*k)
	  enddo
	enddo
      enddo
c
c=======================================================================
c     simulate a time step loop and the flexible memory window of MOM 
c
c     neq      = number of tracer eqns
c     num_reps = makes enough work between tic/toc pairs for accurate
c                timing when machine dependent "second" function only
c                returns time to nearest second (as on SGI)
c=======================================================================
c
      num_timesteps = 1
      num_reps      = 5
      neq           = 2
      do numts=1,num_timesteps
c
        call tic ('test 1','storing all intermediate fluxes')
	do nrep = 1,num_reps
        do mw = 1,maxmw
	  joff   = (mw-1)*ncrows
          js = 2
          je = min(jmw,jmt-joff) -1
	  do n=1,neq
          do j=js,je
	    do k=1,km
	      do i=is-1,ie
	        diff_fe(i,k,j) = diff2_e(i,j)*(t(i+1,k,j) - t(i,k,j))
	        adv_fe(i,k,j) = adv_vet(i,k,j)*(t(i,k,j) + t(i+1,k,j))
              enddo
	    enddo
	  enddo
          do j=js-1,je
	    do k=1,km
	      do i=is,ie
	        diff_fn(i,k,j) = diff2_n(j)*(t(i,k,j+1) - t(i,k,j))
                adv_fn(i,k,j) = adv_vnt(i,k,j)*(t(i,k,j+1) + t(i,k,j))
              enddo
	    enddo
	  enddo
          do j=js,je
            do k=1,km-1
              do i=is,ie
	        adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,j) + t(i,k+1,j))
                diff_fb(i,k,j)  = diff2_b(k)*(t(i,k,j) - t(i,k+1,j))
              enddo
            enddo
          enddo
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        ta(i,k,j) = tb(i,k,j) + tdt(k)*(
     &                 DIFF_Tx(i,k,j) + DIFF_Ty(i,k,j) + DIFF_Tz(i,k,j)
     &                -ADV_Tx(i,k,j) - ADV_Ty(i,k,j) - ADV_Tz(i,k,j))
              enddo
	    enddo
	  enddo
          do j=1,jmw
	    do k=1,km
	      do i=1,imt
	        tb(i,k,j) = t(i,k,j)
	        t(i,k,j)  = ta(i,k,j)
              enddo
	    enddo
	  enddo
        enddo
        enddo
        enddo
        call toc ('test 1', 'storing all intermediate fluxes')
c
c-----------------------------------------------------------------------
c
        call tic ('test 2', 'recalc fluxes in DIFF_Ty')
	do nrep = 1,num_reps
        do mw = 1,maxmw
	  joff   = (mw-1)*ncrows
          js = 2
          je = min(jmw,jmt-joff) -1
	  do n=1,neq
          do j=js,je
	    do k=1,km
	      do i=is-1,ie
	        diff_fe(i,k,j) = diff2_e(i,j)*(t(i+1,k,j) - t(i,k,j))
	        adv_fe(i,k,j) = adv_vet(i,k,j)*(t(i,k,j) + t(i+1,k,j))
              enddo
	    enddo
	  enddo
          do j=js-1,je
	    do k=1,km
	      do i=is,ie
                adv_fn(i,k,j) = adv_vnt(i,k,j)*(t(i,k,j+1) + t(i,k,j))
              enddo
	    enddo
	  enddo 
          do j=js,je
            do k=1,km-1
              do i=is,ie
	        adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,j) + t(i,k+1,j))
                diff_fb(i,k,j)  = diff2_b(k)*(t(i,k,j) - t(i,k+1,j))
              enddo
            enddo
          enddo
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        ta(i,k,j) = tb(i,k,j) + tdt(k)*(
     &                  DIFF_Tx(i,k,j) + T2yy(i,k,j) + DIFF_Tz(i,k,j)
     &                 -ADV_Tx(i,k,j) - ADV_Ty(i,k,j) - ADV_Tz(i,k,j))
              enddo
	    enddo
	  enddo
          do j=1,jmw
	    do k=1,km
	      do i=1,imt
	        tb(i,k,j) = t(i,k,j)
	        t(i,k,j)  = ta(i,k,j)
              enddo
	    enddo
	  enddo
        enddo
	enddo
	enddo
        call toc ('test 2', 'recalc fluxes in DIFF_Ty')
c
c-----------------------------------------------------------------------
c
        call tic ('test 2a', 'recalc using T4yy')
	do nrep = 1,num_reps
        do mw = 1,maxmw
	  joff   = (mw-1)*ncrows
          js = 2
          je = min(jmw,jmt-joff) -1
	  do n=1,neq
          do j=js,je
	    do k=1,km
	      do i=is-1,ie
	        diff_fe(i,k,j) = diff2_e(i,j)*(t(i+1,k,j) - t(i,k,j))
	        adv_fe(i,k,j) = adv_vet(i,k,j)*(t(i,k,j) + t(i+1,k,j))
              enddo
	    enddo
	  enddo
          do j=js-1,je
	    do k=1,km
	      do i=is,ie
                adv_fn(i,k,j) = adv_vnt(i,k,j)*(t(i,k,j+1) + t(i,k,j))
              enddo
	    enddo
	  enddo
          do j=js,je
            do k=1,km-1
              do i=is,ie
	        adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,j) + t(i,k+1,j))
                diff_fb(i,k,j)  = diff2_b(k)*(t(i,k,j) - t(i,k+1,j))
              enddo
            enddo
          enddo
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        ta(i,k,j) = tb(i,k,j) + tdt(k)*(
     &                   DIFF_Tx(i,k,j) + T4yy(i,k,j) + DIFF_Tz(i,k,j)
     &                  -ADV_Tx(i,k,j) - ADV_Ty(i,k,j) - ADV_Tz(i,k,j))
              enddo
	    enddo
	  enddo
          do j=1,jmw
	    do k=1,km
	      do i=1,imt
	        tb(i,k,j) = t(i,k,j)
	        t(i,k,j)  = ta(i,k,j)
              enddo
	    enddo
	  enddo
        enddo
	enddo
	enddo
        call toc ('test 2a', 'recalc using T4yy')
c
c-----------------------------------------------------------------------
c
        call tic ('test 2b', 'recalc using T4yy & VT2y')
	do nrep = 1,num_reps
        do mw = 1,maxmw
	  joff   = (mw-1)*ncrows
          js = 2
          je = min(jmw,jmt-joff) -1
	  do n=1,neq
          do j=js,je
	    do k=1,km
	      do i=is-1,ie
	        diff_fe(i,k,j) = diff2_e(i,j)*(t(i+1,k,j) - t(i,k,j))
	        adv_fe(i,k,j) = adv_vet(i,k,j)*(t(i,k,j) + t(i+1,k,j))
              enddo
	    enddo
	  enddo
          do j=js,je
            do k=1,km-1
              do i=is,ie
	        adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,j) + t(i,k+1,j))
                diff_fb(i,k,j)  = diff2_b(k)*(t(i,k,j) - t(i,k+1,j))
              enddo
            enddo
          enddo
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        ta(i,k,j) = tb(i,k,j) + tdt(k)*(
     &                    DIFF_Tx(i,k,j) + T4yy(i,k,j) + DIFF_Tz(i,k,j)
     &                   -ADV_Tx(i,k,j) - VT2y(i,k,j) - ADV_Tz(i,k,j))
              enddo
	    enddo
	  enddo
          do j=1,jmw
	    do k=1,km
	      do i=1,imt
	        tb(i,k,j) = t(i,k,j)
	        t(i,k,j)  = ta(i,k,j)
              enddo
	    enddo
	  enddo
        enddo
	enddo
	enddo
        call toc ('test 2b', 'recalc using T4yy & VT2y')
c
c-----------------------------------------------------------------------
c
        call tic ('test 2c', 'recalc using T4yy & VT4y')
	do nrep = 1,num_reps
        do mw = 1,maxmw
	  joff   = (mw-1)*ncrows
          js = 2
          je = min(jmw,jmt-joff) -1
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        adv_vnt1(i,k,j) = adv_vnt(i,k,j)*cstdyt2r(j)
	        adv_vnt2(i,k,j) = adv_vnt(i,k,j-1)*cstdyt2r(j)
              enddo
	    enddo
	  enddo
	  do n=1,neq
          do j=js,je
	    do k=1,km
	      do i=is-1,ie
	        diff_fe(i,k,j) = diff2_e(i,j)*(t(i+1,k,j) - t(i,k,j))
	        adv_fe(i,k,j) = adv_vet(i,k,j)*(t(i,k,j) + t(i+1,k,j))
              enddo
	    enddo
	  enddo
          do j=js,je
            do k=1,km-1
              do i=is,ie
	        adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,j) + t(i,k+1,j))
                diff_fb(i,k,j)  = diff2_b(k)*(t(i,k,j) - t(i,k+1,j))
              enddo
            enddo
          enddo
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        ta(i,k,j) = tb(i,k,j) + tdt(k)*(
     &                    DIFF_Tx(i,k,j) + T4yy(i,k,j) + DIFF_Tz(i,k,j)
     &                   -ADV_Tx(i,k,j) - VT4y(i,k,j) - ADV_Tz(i,k,j))
              enddo
	    enddo
	  enddo
          do j=1,jmw
	    do k=1,km
	      do i=1,imt
	        tb(i,k,j) = t(i,k,j)
	        t(i,k,j)  = ta(i,k,j)
              enddo
	    enddo
	  enddo
        enddo
	enddo
	enddo
        call toc ('test 2c', 'recalc using T4yy & VT4y')
c
c-----------------------------------------------------------------------
c
        call tic ('test 3', 'recalc fluxes in DIFF_Tx')
	do nrep = 1,num_reps
        do mw = 1,maxmw
	  joff   = (mw-1)*ncrows
          js = 2
          je = min(jmw,jmt-joff) -1
	  do n=1,neq
          do j=js,je
	    do k=1,km
	      do i=is-1,ie
	        adv_fe(i,k,j) = adv_vet(i,k,j)*(t(i,k,j) + t(i+1,k,j))
              enddo
	    enddo
	  enddo
          do j=js-1,je
	    do k=1,km
	      do i=is,ie
	        diff_fn(i,k,j) = diff2_n(j)*(t(i,k,j+1) - t(i,k,j))
                adv_fn(i,k,j) = adv_vnt(i,k,j)*(t(i,k,j+1) + t(i,k,j))
              enddo
	    enddo
	  enddo
          do j=js,je
            do k=1,km-1
              do i=is,ie
	        adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,j) + t(i,k+1,j))
                diff_fb(i,k,j)  = diff2_b(k)*(t(i,k,j) - t(i,k+1,j))
              enddo
            enddo
          enddo
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        ta(i,k,j) = tb(i,k,j) + tdt(k)*(
     &                   T2xx(i,k,j) + DIFF_Ty(i,k,j) + DIFF_Tz(i,k,j)
     &                  -ADV_Tx(i,k,j) - ADV_Ty(i,k,j) - ADV_Tz(i,k,j))
              enddo
	    enddo
	  enddo
          do j=1,jmw
	    do k=1,km
	      do i=1,imt
	        tb(i,k,j) = t(i,k,j)
	        t(i,k,j)  = ta(i,k,j)
              enddo
	    enddo
	  enddo
        enddo
	enddo
	enddo
        call toc ('test 3', 'recalc fluxes in DIFF_Tx')
c
c-----------------------------------------------------------------------
c
        call tic ('test 3a', 'recalc using T4xx')
	do nrep = 1,num_reps
        do mw = 1,maxmw
	  joff   = (mw-1)*ncrows
          js = 2
          je = min(jmw,jmt-joff) -1
          do j=js,je
            do i=is,ie
              diff4_e(i,j) = diff2_e(i,j)*cstdxtr(i,j)
	      diff4_w(i,j) = diff2_e(i-1,j)*cstdxtr(i,j)
            enddo
          enddo
	  do n=1,neq
          do j=js,je
	    do k=1,km
	      do i=is-1,ie
	        adv_fe(i,k,j) = adv_vet(i,k,j)*(t(i,k,j) + t(i+1,k,j))
              enddo
	    enddo
	  enddo
          do j=js-1,je
	    do k=1,km
	      do i=is,ie
	        diff_fn(i,k,j) = diff2_n(j)*(t(i,k,j+1) - t(i,k,j))
                adv_fn(i,k,j) = adv_vnt(i,k,j)*(t(i,k,j+1) + t(i,k,j))
              enddo
	    enddo
	  enddo
          do j=js,je
            do k=1,km-1
              do i=is,ie
	        adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,j) + t(i,k+1,j))
                diff_fb(i,k,j)  = diff2_b(k)*(t(i,k,j) - t(i,k+1,j))
              enddo
            enddo
          enddo
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        ta(i,k,j) = tb(i,k,j) + tdt(k)*(
     &                  T4xx(i,k,j) + DIFF_Ty(i,k,j) + DIFF_Tz(i,k,j)
     &                 -ADV_Tx(i,k,j) - ADV_Ty(i,k,j) - ADV_Tz(i,k,j))
              enddo
	    enddo
	  enddo
          do j=1,jmw
	    do k=1,km
	      do i=1,imt
	        tb(i,k,j) = t(i,k,j)
	        t(i,k,j)  = ta(i,k,j)
              enddo
	    enddo
	  enddo
        enddo
	enddo
	enddo
        call toc ('test 3a', 'recalc using T4xx')
c
c-----------------------------------------------------------------------
c
        call tic ('test 3b', 'recalc using T4xx & UT2x')
	do nrep = 1,num_reps
        do mw = 1,maxmw
	  joff   = (mw-1)*ncrows
          js = 2
          je = min(jmw,jmt-joff) -1
          do j=js,je
            do i=is,ie
              diff4_e(i,j) = diff2_e(i,j)*cstdxtr(i,j)
	      diff4_w(i,j) = diff2_e(i-1,j)*cstdxtr(i,j)
            enddo
          enddo
	  do n=1,neq
          do j=js-1,je
	    do k=1,km
	      do i=is,ie
	        diff_fn(i,k,j) = diff2_n(j)*(t(i,k,j+1) - t(i,k,j))
                adv_fn(i,k,j) = adv_vnt(i,k,j)*(t(i,k,j+1) + t(i,k,j))
              enddo
	    enddo
	  enddo
          do j=js,je
            do k=1,km-1
              do i=is,ie
	        adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,j) + t(i,k+1,j))
                diff_fb(i,k,j)  = diff2_b(k)*(t(i,k,j) - t(i,k+1,j))
              enddo
            enddo
          enddo
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        ta(i,k,j) = tb(i,k,j) + tdt(k)*(
     &                   T4xx(i,k,j) + DIFF_Ty(i,k,j) + DIFF_Tz(i,k,j)
     &                  -UT2x(i,k,j) - ADV_Ty(i,k,j) - ADV_Tz(i,k,j))
              enddo
	    enddo
	  enddo
          do j=1,jmw
	    do k=1,km
	      do i=1,imt
	        tb(i,k,j) = t(i,k,j)
	        t(i,k,j)  = ta(i,k,j)
              enddo
	    enddo
	  enddo
        enddo
	enddo
	enddo
        call toc ('test 3b', 'recalc using T4xx & UT2x')
c
c-----------------------------------------------------------------------
c
        call tic ('test 4', 'recalc fluxes in DIFF_Tx & DIFF_Ty')
	do nrep = 1,num_reps
        do mw = 1,maxmw
	  joff   = (mw-1)*ncrows
          js = 2
          je = min(jmw,jmt-joff) -1
	  do n=1,neq
          do j=js,je
	    do k=1,km
	      do i=is-1,ie
	        adv_fe(i,k,j) = adv_vet(i,k,j)*(t(i,k,j) + t(i+1,k,j))
              enddo
	    enddo
	  enddo
          do j=js-1,je
	    do k=1,km
	      do i=is,ie
                adv_fn(i,k,j) = adv_vnt(i,k,j)*(t(i,k,j+1) + t(i,k,j))
              enddo
	    enddo
	  enddo
          do j=js,je
            do k=1,km-1
              do i=is,ie
	        adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,j) + t(i,k+1,j))
                diff_fb(i,k,j)  = diff2_b(k)*(t(i,k,j) - t(i,k+1,j))
              enddo
            enddo
          enddo
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        ta(i,k,j) = tb(i,k,j) + tdt(k)*(
     &                   T2xx(i,k,j) + T2yy(i,k,j) + DIFF_Tz(i,k,j)
     &                  -ADV_Tx(i,k,j) - ADV_Ty(i,k,j) - ADV_Tz(i,k,j))
              enddo
	    enddo
	  enddo
          do j=1,jmw
	    do k=1,km
	      do i=1,imt
	        tb(i,k,j) = t(i,k,j)
	        t(i,k,j)  = ta(i,k,j)
              enddo
	    enddo
	  enddo
        enddo
	enddo
	enddo
        call toc ('test 4', 'recalc fluxes in DIFF_Tx & DIFF_Ty')
c
c-----------------------------------------------------------------------
c
        call tic ('test 4a', 'recalc using T4xx & T4yy')
	do nrep = 1,num_reps
        do mw = 1,maxmw
	  joff   = (mw-1)*ncrows
          js = 2
          je = min(jmw,jmt-joff) -1
          do j=js,je
            do i=is,ie
              diff4_e(i,j) = diff2_e(i,j)*cstdxtr(i,j)
	      diff4_w(i,j) = diff2_e(i-1,j)*cstdxtr(i,j)
            enddo
          enddo
	  do n=1,neq
          do j=js,je
	    do k=1,km
	      do i=is-1,ie
	        adv_fe(i,k,j) = adv_vet(i,k,j)*(t(i,k,j) + t(i+1,k,j))
              enddo
	    enddo
	  enddo
          do j=js-1,je
	    do k=1,km
	      do i=is,ie
                adv_fn(i,k,j) = adv_vnt(i,k,j)*(t(i,k,j+1) + t(i,k,j))
              enddo
	    enddo
	  enddo
          do j=js,je
            do k=1,km-1
              do i=is,ie
	        adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,j) + t(i,k+1,j))
                diff_fb(i,k,j)  = diff2_b(k)*(t(i,k,j) - t(i,k+1,j))
              enddo
            enddo
          enddo
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        ta(i,k,j) = tb(i,k,j) + tdt(k)*(
     &                   T4xx(i,k,j) + T4yy(i,k,j) + DIFF_Tz(i,k,j)
     &                  -ADV_Tx(i,k,j) - ADV_Ty(i,k,j) - ADV_Tz(i,k,j))
              enddo
	    enddo
	  enddo
          do j=1,jmw
	    do k=1,km
	      do i=1,imt
	        tb(i,k,j) = t(i,k,j)
	        t(i,k,j)  = ta(i,k,j)
              enddo
	    enddo
	  enddo
        enddo
	enddo
	enddo
        call toc ('test 4a', 'recalc using T4xx & T4yy')
c
c-----------------------------------------------------------------------
c
        call tic ('test 4b', 'recalc T4xx,T4yy,UT2x,VT2y')
	do nrep = 1,num_reps
        do mw = 1,maxmw
	  joff   = (mw-1)*ncrows
          js = 2
          je = min(jmw,jmt-joff) -1
          do j=js,je
            do i=is,ie
              diff4_e(i,j) = diff2_e(i,j)*cstdxtr(i,j)
	      diff4_w(i,j) = diff2_e(i-1,j)*cstdxtr(i,j)
            enddo
          enddo
	  do n=1,neq
          do j=js,je
            do k=1,km-1
              do i=is,ie
	        adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,j) + t(i,k+1,j))
                diff_fb(i,k,j)  = diff2_b(k)*(t(i,k,j) - t(i,k+1,j))
              enddo
            enddo
          enddo
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        ta(i,k,j) = tb(i,k,j) + tdt(k)*(
     &                   T4xx(i,k,j) + T4yy(i,k,j) + DIFF_Tz(i,k,j)
     &                  -UT2x(i,k,j) - VT2y(i,k,j) - ADV_Tz(i,k,j))
              enddo
	    enddo
	  enddo
          do j=1,jmw
	    do k=1,km
	      do i=1,imt
	        tb(i,k,j) = t(i,k,j)
	        t(i,k,j)  = ta(i,k,j)
              enddo
	    enddo
	  enddo
        enddo
	enddo
	enddo
        call toc ('test 4b', 'recalc T4xx,T4yy,UT2x,VT2y')
c
c-----------------------------------------------------------------------
c
        call tic ('test 4c', 'recalc T4xx,T4yy,UT2x,VT4y')
	do nrep = 1,num_reps
        do mw = 1,maxmw
	  joff   = (mw-1)*ncrows
          js = 2
          je = min(jmw,jmt-joff) -1
          do j=js,je
            do i=is,ie
              diff4_e(i,j) = diff2_e(i,j)*cstdxtr(i,j)
	      diff4_w(i,j) = diff2_e(i-1,j)*cstdxtr(i,j)
            enddo
          enddo
	  do n=1,neq
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        adv_vnt1(i,k,j) = adv_vnt(i,k,j)*cstdyt2r(j)
	        adv_vnt2(i,k,j) = adv_vnt(i,k,j-1)*cstdyt2r(j)
              enddo
	    enddo
	  enddo
          do j=js,je
            do k=1,km-1
              do i=is,ie
	        adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,j) + t(i,k+1,j))
                diff_fb(i,k,j)  = diff2_b(k)*(t(i,k,j) - t(i,k+1,j))
              enddo
            enddo
          enddo
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        ta(i,k,j) = tb(i,k,j) + tdt(k)*(
     &                   T4xx(i,k,j) + T4yy(i,k,j) + DIFF_Tz(i,k,j)
     &                  -UT2x(i,k,j) - VT4y(i,k,j) - ADV_Tz(i,k,j))
              enddo
	    enddo
	  enddo
          do j=1,jmw
	    do k=1,km
	      do i=1,imt
	        tb(i,k,j) = t(i,k,j)
	        t(i,k,j)  = ta(i,k,j)
              enddo
	    enddo
	  enddo
        enddo
        enddo
	enddo
        call toc ('test 4c', 'recalc T4xx,T4yy,UT2x,VT4y')
c
c-----------------------------------------------------------------------
c
        call tic ('test 4d', 'recalc T4xx,T4yy,UT4x,VT4y')
	do nrep = 1,num_reps
        do mw = 1,maxmw
	  joff   = (mw-1)*ncrows
          js = 2
          je = min(jmw,jmt-joff) -1
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        adv_vnt1(i,k,j) = adv_vnt(i,k,j)*cstdyt2r(j)
	        adv_vnt2(i,k,j) = adv_vnt(i,k,j-1)*cstdyt2r(j)
	        adv_vet1(i,k,j) = adv_vet(i,k,j)*cstdxt2r(i,j)
	        adv_vet2(i,k,j) = adv_vet(i-1,k,j)*cstdxt2r(i,j)
              enddo
	    enddo
	  enddo
          do j=js,je
            do i=is,ie
              diff4_e(i,j) = diff2_e(i,j)*cstdxtr(i,j)
	      diff4_w(i,j) = diff2_e(i-1,j)*cstdxtr(i,j)
            enddo
          enddo
	  do n=1,neq
          do j=js,je
            do k=1,km-1
              do i=is,ie
	        adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,j) + t(i,k+1,j))
                diff_fb(i,k,j)  = diff2_b(k)*(t(i,k,j) - t(i,k+1,j))
              enddo
            enddo
          enddo
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        ta(i,k,j) = tb(i,k,j) + tdt(k)*(
     &                  T4xx(i,k,j) + T4yy(i,k,j) + DIFF_Tz(i,k,j)
     &                 -UT4x(i,k,j) - VT4y(i,k,j) - ADV_Tz(i,k,j))
              enddo
	    enddo
	  enddo
          do j=1,jmw
	    do k=1,km
	      do i=1,imt
	        tb(i,k,j) = t(i,k,j)
	        t(i,k,j)  = ta(i,k,j)
              enddo
	    enddo
	  enddo
        enddo
        enddo
	enddo
        call toc ('test 4d', 'recalc T4xx,T4yy,UT4x,VT4y')
c
c-----------------------------------------------------------------------
c
        call tic ('test 4e', 'recalc T4xx,T4yy,UT4x,VT2y')
	do nrep = 1,num_reps
        do mw = 1,maxmw
	  joff   = (mw-1)*ncrows
          js = 2
          je = min(jmw,jmt-joff) -1
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        adv_vet1(i,k,j) = adv_vet(i,k,j)*cstdxt2r(i,j)
	        adv_vet2(i,k,j) = adv_vet(i-1,k,j)*cstdxt2r(i,j)
              enddo
	    enddo
	  enddo
          do j=js,je
            do i=is,ie
              diff4_e(i,j) = diff2_e(i,j)*cstdxtr(i,j)
	      diff4_w(i,j) = diff2_e(i-1,j)*cstdxtr(i,j)
            enddo
          enddo
	  do n=1,neq
          do j=js,je
            do k=1,km-1
              do i=is,ie
	        adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,j) + t(i,k+1,j))
                diff_fb(i,k,j)  = diff2_b(k)*(t(i,k,j) - t(i,k+1,j))
              enddo
            enddo
          enddo
          do j=js,je
	    do k=1,km
	      do i=is,ie
	        ta(i,k,j) = tb(i,k,j) + tdt(k)*(
     &                  T4xx(i,k,j) + T4yy(i,k,j) + DIFF_Tz(i,k,j)
     &                 -UT4x(i,k,j) - VT2y(i,k,j) - ADV_Tz(i,k,j))
              enddo
	    enddo
	  enddo
          do j=1,jmw
	    do k=1,km
	      do i=1,imt
	        tb(i,k,j) = t(i,k,j)
	        t(i,k,j)  = ta(i,k,j)
              enddo
	    enddo
	  enddo
        enddo
        enddo
	enddo
        call toc ('test 4e', 'recalc T4xx,T4yy,UT4x,VT2y')
c
c-----------------------------------------------------------------------
c       limit size of numbers to prevent blow up
c
        do mw = 1,maxmw
	  joff   = (mw-1)*ncrows
          js = 2
          je = min(jmw,jmt-joff) -1
          do j=js,je
	    do k=1,km
	      do i=is-1,ie+1
	        t(i,k,j) = sin(t(i,k,j))
	      enddo
	    enddo
	  enddo
        enddo
c
      enddo
      write (stdout,*) t(2,2,2), diff_fn(2,2,2),diff_fe(2,2,2)
      call showtimers
      call abort()
      end

#endif

c#define debug_timer
c#define show_tictocs

c     #define show_tictocs  shows all calls to tic, toc, and ticr
c            and all implied tic/tocs arising from time push/pops

c real_seconds is defined if the cpu clock shows milliseconds or better
# define real_seconds
# ifdef sgi_mipspro
#  undef real_seconds
# endif
      subroutine inittimers
#ifdef timing
c
c=======================================================================
c     initialize timers
c
c         active(i) = j   the one and only active subtimer in the group
c                   = 0   no subtimer in the group active
c
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
c
# include "stdunits.h"
# include "timer.h"
c
# ifndef SGICRAY
      write (stdout,'(/a/a/a/)')
     &  'ERROR: All timings are based on calls to the function "second"'
     &, '       Please modify the "second" function at the bottom of'
     &, '       timer.F to use your system clock'
      if (1 .eq. 1) call abort()
# endif
c
c     initialize clock to detect timing intervals that include midnight
c
      prev_second = 0.0
c
      s = second()
      do i=0,ntimers
        timer(i) = 'unused'
        active(i) = 0
        istackptr(i) = 0
        group_overflow(i) = .false.
        do j=1,nsubtimers
          subtimer(i,j) = 'unused'
          starttime(i,j) = 0.0
          elapsedtime(i,j) = 0.0
          ncalls(i,j) = 0
          nresumes(i,j) = 0
          subtimer_active(i,j) = 0
          nested_calls(i,j) = 0
          weighted_calls(i,j) = 0
          npops(i,j) = 0
          npushes(i,j) = 0
        end do
        subtimer(i,nsubtimers) = 'don''t use'
      end do
      timer(ntimers) = 'don''t use'
      ngroup_overflow = .false.
      ntictoc_errors = 0

c     initialize 'total time' timer

      timer(0) = 'total time'
      subtimer(0,1) = ' '
      ncalls(0,1) = 1
      active(0) = 1
      starttime(0,1) = s

c     initial timing calibration for tic/toc calls
c     time sets of 1000 null tic/toc pairs

      max_timer_warnings = 100
      print_timer_warnings = .false.
      print_tictocs = .false.

      call tic ('tic/toc', '1000 calls + 1000 overheads')
      do i=1,1000
        call tic ('tic/toc', '1000 empty calls')
        call toc ('tic/toc', '1000 empty calls')
      end do
      call toc ('tic/toc', '1000 calls + 1000 overheads')

      call tic ('tic/toc3', '1000 calls + 1000 overheads')
      do i=1,1000
        call tic ('tic/toc', '1000 empty calls3')
        call toc ('tic/toc', '1000 empty calls3')
      end do
      call toc ('tic/toc3', '1000 calls + 1000 overheads')

      call tic ('tic/toc4', '1000 calls + 1000 overheads')
      do i=1,1000
        call tic ('tic/toc', '1000 empty calls4')
        call toc ('tic/toc', '1000 empty calls4')
      end do
      call toc ('tic/toc4', '1000 calls + 1000 overheads')

      call tic ('tic/toc5', '1000 calls + 1000 overheads')
      do i=1,1000
        call tic ('tic/toc', '1000 empty calls5')
        call toc ('tic/toc', '1000 empty calls5')
      end do
      call toc ('tic/toc5', '1000 calls + 1000 overheads')

      call tic ('tic/toc', '1000 calls + 1000 overhead5')
      do i=1,1000
        call tic ('tic/toc', '1000 empty calls5')
        call toc ('tic/toc', '1000 empty calls5')
      end do
      call toc ('tic/toc', '1000 calls + 1000 overhead5')

      call tic ('tic/toc', '1000 calls + 1000 overhead8')
      do i=1,1000
        call tic ('tic/toc', '1000 empty calls8')
        call toc ('tic/toc', '1000 empty calls8')
      end do
      call toc ('tic/toc', '1000 calls + 1000 overhead8')
c
      print_timer_warnings = .true.
      print_tictocs = .true.

#endif
      return
      end


      subroutine tic (name, subname)
#ifdef timing
c
c=======================================================================
c     start a timer
c
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
c
# include "stdunits.h"
# include "timer.h"
      character*(*) name, subname
      character*40 name40, subname40
      logical error

      name40 = name
      subname40 = subname
      call findtimer (name40, i, error)
      if (.not. error) then
        call findsubtimer (subname40, i, j, error)
        if (.not. error) then
#ifdef show_tictocs
          call show_tictoc('tic', i, j, 0.0)
c         show_tictoc will print a time that is slightly too early, but
c         this printout must precede warnings based on this tic.
#endif
          if (active(i) .eq. j) then
            if (print_timer_warnings) then
              write (stdout,'(a,a17,tr1,a32,/,tr12,a)')
     &            'WARNING:    two successive tics for ', timer(i)
     &,            subtimer(i,j)
     &,            'tic ignored'
            end if
            call inc_timer_warn
          else
            call tictoc_nesting ('tic',i,j,s)
c
c           tictoc nesting will return a starttime s for this tic
c
            starttime(i,j) = s
            ncalls(i,j) = ncalls(i,j) + 1
            ncalls(i,1) = ncalls(i,1) + 1
          end if
        else
          group_overflow(i) = .true.
        end if
      else
        ngroup_overflow = .true.
      end if
#endif
      return
      end

c-----------------------------------------------------------------------
c
      subroutine toc (name, subname)
#ifdef timing
c
c=======================================================================
c     stop a timer
c
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
c
# include "stdunits.h"
# include "timer.h"
      character*(*) name, subname
      character*40 name40, subname40
      logical error
      logical keep_going

      s = second()
      name40 = name
      subname40 = subname
      call findtimer (name40, i, error)
      if (.not. error) then
        call findsubtimer (subname40, i, j, error)
        if (.not. error) then
#ifdef show_tictocs
          call show_tictoc('toc', i, j, s)
#endif
          if (active(i) .eq. j) then
            elapsedtime(i,j) = elapsedtime(i,j) + s - starttime(i,j)
            elapsedtime(i,1) = elapsedtime(i,1) + s - starttime(i,j)
            call tictoc_nesting ('toc',i,j,s)
          else 
            jj = active(i)
            if (print_timer_warnings) then
              write (stdout,'(a,a17,tr1,a32)')
     &           'WARNING:    improperly nested toc:  ', timer(i)
     &,           subtimer(i,j)
              if (jj .ne. 0) then
                write (stdout,'(tr12,a,a17,tr1,a32)')
     &           'active tic is:          ', timer(i), subtimer(i,jj)
              else
                write (stdout,'(tr12,a,a17)')
     &           'no tic active in group: ', timer(i)
              end if
            end if
            call inc_timer_warn
c
c           see if subtimer(i,j) is on stack for group i
c
            itop = istackptr(i)
            iptr = istackptr(i)
            keep_going = .true.
            do while (iptr .ge. 1 .and. keep_going)
              iptr = iptr -1
              if (jstack(i,iptr) .eq. j) then
c
c               toc the active subtimer in this group
c
                j = active(i)
                write (stdout,'(tr12,a,a17,tr1,a32)')
     &              'toc-ing timer:          ', timer(i), subtimer(i,j)
                elapsedtime(i,j) = elapsedtime(i,j) +s - starttime(i,j)
                elapsedtime(i,1) = elapsedtime(i,1) +s - starttime(i,j)
                call tictoc_nesting ('toc',i,j,s)
c
c               pop stopped subtimers off the stack
c
                do iptr2=itop,iptr,-1
                  j = active(i)
                  write (stdout,'(tr12,a,a17,tr1,a32)')
     &              'toc-ing timer:          ', timer(i), subtimer(i,j)
                  call tictoc_nesting ('toc',i,j,s)
                end do
                keep_going = .false.
              end if
            end do
            if (keep_going) then
c           
c           subtimer(i,j) is not on stack(i)
c
              if (print_timer_warnings) then
                write (stdout,'(tr12,a,a)')
     &             'toc ignored:  '
     &,            'one cannot toc a timer that is not currently tic-ed'
              end if
              call inc_timer_warn
            endif
          end if
        else
          group_overflow(i) = .true.
        end if
      else
        ngroup_overflow = .true.
      end if
#endif
      return
      end


      subroutine ticr (name, subname)
#ifdef timing
c
c=======================================================================
c     restart a timer which has been stopped without incrementing the
c     number of times it has been started.
c
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
c
# include "timer.h"
      character*(*) name, subname
      character*40 name40, subname40
      logical error

      name40 = name
      subname40 = subname
      call findtimer (name40, i, error)
      if (.not. error) then
        call findsubtimer (subname40, i, j, error)
        if (.not. error) then
#ifdef show_tictocs
          call show_tictoc('ticr', i, j, 0.0)
c         show_tictoc will print a time that is slightly too early, but
c         this printout must precede warnings based on this tic.
#endif
          call tictoc_nesting ('ticr',i,j,s)
c
c         tictoc nesting will return a starttime s for this tic
c
          nresumes(i,j) = nresumes(i,j) + 1
          nresumes(i,1) = nresumes(i,1) + 1
          starttime(i,j) = s
        else
          group_overflow(i) = .true.
        end if
      else
        ngroup_overflow = .true.
      end if
#endif
      return
      end


      subroutine showtimers
#ifdef timing
c
c=======================================================================
c     show a summary of timings.
c
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
c
# include "stdunits.h"
# include "timer.h"
      parameter (epsilon=1.0e-8)
      character*40 header(10)

c     final timing calibration for tic/toc calls
c     time sets of 1000 null tic/toc pairs

      print_timer_warnings = .false.
      print_tictocs = .false.
c
      call tic ('tic/toc2', '1000 calls + 1000 overhead2')
      do i=1,1000
        call tic ('tic/toc2', '1000 empty calls/pp')
        call toc ('tic/toc2', '1000 empty calls/pp')
      end do
      call toc ('tic/toc2', '1000 calls + 1000 overhead2')

      call tic ('tic/toc2a', '1000 calls + 1000 overheads')
      do i=1,1000
        call tic ('tic/toc2', '1000 empty calls')
        call toc ('tic/toc2', '1000 empty calls')
      end do
      call toc ('tic/toc2a', '1000 calls + 1000 overheads')

c     find total time

      s = second()
      elapsedtime(0,1) = s - starttime(0,1)

c     print summary

      header(1) = '  group'
      header(2) = '  sub group'
      header(3) = '   calls'
      header(4) = '   time(sec)'
      header(5) = '  time/call'
      header(6) = ' % group'
      header(7) = ' % total'
c# ifdef debug_timer
      write (stdout,'(a17,a32,a8,a13,a11,2a8,/)')
     &       (header(i),i=1,7)
c# endif
      do i=0,ntimers
        if (timer(i) .ne. 'unused' .and.
     &      timer(i) .ne. 'don''t use') then
          if (timer(i) .eq. 'tic/toc') ioverhead1 = i
          if (timer(i) .eq. 'tic/toc2') ioverhead2 = i
          do j=1,nsubtimers
            if (subtimer(i,j) .ne. 'unused'  .and.
     &          subtimer(i,j) .ne. 'don''t use') then
              if (abs(elapsedtime(i,1)) .le. epsilon)
     &          elapsedtime(i,1) = epsilon
c# ifdef debug_timer
              if (i .eq. 0 .or.
     &            .not.(i .eq. ioverhead1  .and. j .eq. 1) .and.
     &            .not.(i .eq. ioverhead2  .and. j .eq. 1) .and.
     &            .not.(j .eq. 1 .and. subtimer(i,3) .eq. 'unused'))then
                write (stdout,'(a17,a32,i8,f13.6,f11.6,2f8.3)')
     &              timer(i), subtimer(i,j),
     &              ncalls(i,j), elapsedtime(i,j),
     &              elapsedtime(i,j)/ncalls(i,j),
     &              100*elapsedtime(i,j)/elapsedtime(i,1),
     &              100*elapsedtime(i,j)/elapsedtime(0,1)
              endif
c# endif
            end if
          end do
# ifdef debug_timer
          write (stdout,*) ' '
# endif
        end if
      end do
c
c     print group and subgroup overflow warnings
c
      if (ngroup_overflow) then
        write (stdout,'(//,a,/,a,i3)')
     &       'WARNING:  Too many timer groups.  Some were lost',
     &       '          Increase ntimers=',ntimers
      end if
      do i=1,ntimers
        if (group_overflow(i)) then
          write (stdout,'(//,a,a,a,/,a,i3)')
     &    'WARNING:  Too many subgroups in group ',timer(i),
     &    '. Some were lost',
     &       '          Increase nsubtimers=',nsubtimers
        end if
      end do

c     adjust group and subgroup times
c
c     subtract nested tic/toc times from individual and group totaltimes
c     and compensate for increased name search times for later groups
c     and subgroups.

      totaltime = elapsedtime(0,1)

c     internal times

      time1pushpop = elapsedtime(ioverhead1,3) / 1000
      time2pushpop = elapsedtime(ioverhead2,3) / 1000
      time1 = elapsedtime(ioverhead1,4) / 1000
      time2 = elapsedtime(ioverhead2,4) / 1000

c     external times

      tictoctime1pushpop = (elapsedtime(ioverhead1,2) +
     &           elapsedtime(ioverhead1,3) - time1pushpop) / 1000
      tictoctime2pushpop = (elapsedtime(ioverhead2,2) +
     &           elapsedtime(ioverhead2,3) - time2pushpop) / 1000
      outertime1pushpop = (elapsedtime(ioverhead1,2)
     &          - time1pushpop) / 1000
      outertime2pushpop = (elapsedtime(ioverhead2,2)
     &          - time2pushpop) / 1000
      tictoctime3 = (elapsedtime(2,2) - time1) / 1000
      tictoctime4 = (elapsedtime(3,2) - time1) / 1000
      tictoctime5 = (elapsedtime(4,2) - time1) / 1000
      tictoctimelast = (elapsedtime(ioverhead2+1,2) - time2) / 1000

      alpha1 = tictoctime4 - tictoctime3
      alpha2 = tictoctime5 - tictoctime4
      alphapushpop = (tictoctime2pushpop - tictoctime1pushpop) /
     &                        (ioverhead2 - ioverhead1)
      alpha = (tictoctimelast - tictoctime3) /
     &             (ioverhead2 - ioverhead1)

c     tictoctime3 is at (1,4) and tictoctime1 at (1,3)
      tictoctime1 = tictoctime3 - alpha1
      pushpopextra1 = outertime1pushpop - tictoctime1
      pushpopinnerextra1 = time1pushpop - time1
# ifdef debug_timer
      write (stdout,'(a20,f13.9)') 
     & 'tictoctime1pushpop=',tictoctime1pushpop
      write (stdout,'(a20,f13.9)') 
     & 'tictoctime2pushpop=',tictoctime2pushpop
      write (stdout,'(a20,f13.9)') 'tictoctime1=',tictoctime1
      write (stdout,'(a20,f13.9)') 'tictoctime3=',tictoctime3
      write (stdout,'(a20,f13.9)') 'tictoctime4=',tictoctime4
      write (stdout,'(a20,f13.9)') 'tictoctime5=',tictoctime5
      write (stdout,'(a20,f13.9)') 'tictoctimelast=',tictoctimelast
      write (stdout,'(a20,f13.9)') 'time1pushpop=',time1pushpop
      write (stdout,'(a20,f13.9)') 'time2pushpop=',time2pushpop
      write (stdout,'(a20,f13.9)') 'time1=',time1
      write (stdout,'(a20,f13.9)') 'time2=',time2
      write (stdout,'(a20,f13.9)') 'alpha1=',alpha1
      write (stdout,'(a20,f13.9)') 'alpha2=',alpha2
      write (stdout,'(a20,f13.9)') 'alphapushpop=',alphapushpop
      write (stdout,'(a20,f13.9)') 'alpha=',alpha
      write (stdout,'(a20,f13.9)') 'pushpopextra1=',pushpopextra1
      write (stdout,'(a20,f13.9)') 
     & 'pushpopinnerextra1=',pushpopinnerextra1
# endif
      do i=0,ntimers
        if (timer(i) .ne. 'unused' .and.
     &      timer(i) .ne. 'don''t use') then
          do j=1,nsubtimers
            if (subtimer(i,j) .ne. 'unused'  .and.
     &          subtimer(i,j) .ne. 'don''t use') then
# ifdef debug_timer
              write (stdout,'(2i4,5i10,6f12.6)')
     &           i,j,ncalls(i,j)+nresumes(i,j)
     &           , nested_calls(i,j), weighted_calls(i,j)
     &           , npops(i,j),npushes(i,j)
     &           , elapsedtime(i,j)
     &           , - tictoctime1 * nested_calls(i,j)
     &           , - alpha*weighted_calls(i,j)
     &           , - pushpopextra1 * npops(i,j)
     &           , - pushpopinnerextra1 * npushes(i,j)
     &           ,  - ((ioverhead2-i)*time1 + (i-ioverhead1)*time2) /
     &             (ioverhead2-ioverhead1) * (ncalls(i,j)+nresumes(i,j))
# endif
c
c             correct elapsed time estimate
c
#ifdef real_seconds
              elapsedtime(i,j) = elapsedtime(i,j)
     &                           - tictoctime1 * nested_calls(i,j)
     &                           - alpha*weighted_calls(i,j)
     &                           - pushpopextra1 * npops(i,j)
     &                           - pushpopinnerextra1 * npushes(i,j)
     &       - ((ioverhead2-i)*time1 + (i-ioverhead1)*time2) /
     &             (ioverhead2-ioverhead1) * (ncalls(i,j)+nresumes(i,j))
#endif
# ifdef debug_timer
              write (stdout,'(2i4,5i10,6f12.6)')
     &           i,j,ncalls(i,j)+nresumes(i,j)
     &           , nested_calls(i,j), weighted_calls(i,j)
     &           , npops(i,j), npushes(i,j)
     &           , elapsedtime(i,j)
# endif
            end if
          end do
        end if
      end do

      subtimer(0,1) = 'computation time'
      subtimer(0,2) = 'calibration tic/toc calls'
      subtimer(0,3) = 'timing tic/toc calls'
      ncalls(0,2) = 8008
      ncalls(0,3) = nested_calls(0,1) - 8008
      tictoc_time = (totaltime - elapsedtime(0,1)) / nested_calls(0,1)
      elapsedtime(0,2) = 8008 * tictoc_time
      elapsedtime(0,3)=totaltime - elapsedtime(0,1) - elapsedtime(0,2)

# ifdef real_seconds
      write (stdout,'(///,a,tr20,a,/)') 'ADJUSTED TIMES',
     &      '(adjusted for overhead in calls to tic and toc)'
# else
      write (stdout,'(///,a,tr20,a,/)') 'UNADJUSTED TIMES',
     &      '(times include overhead in calls to tic and toc)'
# endif
      write (stdout,'(a17,a32,a8,a13,a11,2a8,/)')
     &       (header(i),i=1,7)
      do i=0,ntimers
        if (timer(i) .ne. 'unused' .and.
     &      timer(i) .ne. 'don''t use') then
          if (timer(i)(1:7) .ne. 'tic/toc') then
            do j=1,nsubtimers
              if (subtimer(i,j) .ne. 'unused'  .and.
     &          subtimer(i,j) .ne. 'don''t use') then
                if (abs(elapsedtime(i,1)) .le. epsilon) then
                  elapsedtime(i,1) = epsilon
                endif
                if (i .eq. 0 .or.
     &            .not.(j .eq. 1 .and. subtimer(i,3) .eq. 'unused'))then
                  write (stdout,'(a17,a32,i8,f13.6,f11.6,2f8.3)')
     &              timer(i), subtimer(i,j),
     &              ncalls(i,j), elapsedtime(i,j),
     &              elapsedtime(i,j)/ncalls(i,j),
     &              100*elapsedtime(i,j)/elapsedtime(i,1),
     &              100*elapsedtime(i,j)/elapsedtime(0,1)
                  if (istackptr(i) .ne. 0) then
c                    count times subtimer(i,j) appears on stack(i)
                     nstacked = 0
                     do jj=istackptr(i),1,-1
                       if (jstack(i,jj) .eq. j) nstacked = nstacked + 1
                     end do
                     if (nstacked .gt. 0) then
                       write (stdout,*)
     &                  'WARNING:  ', nstacked, ' tics still pending'
                     end if
                  end if
                end if
              end if
            end do
            write (stdout,*) ' '
          end if
        end if
      end do
c
# ifndef real_seconds
      write (stdout,'(///,tr10,a)')
     &  'Clock times were read in integer seconds'
      write (stdout,'(tr10,a)')
     & 'Expect errors in timing of 1 second/call'
# endif
      write (stdout,'(///,tr10,a,i10,a,/)') 'There were a total of ',
     &    nested_calls(0,1), ' pairs of calls to tic/toc in this run.'
      write (stdout,'(tr10,a,f11.9)') 
     &    'The estimated average time per tic/toc pair is ',
     &    (totaltime - elapsedtime(0,1)) / nested_calls(0,1)
# ifdef debug_timer
      write (stdout,'(tr10,a,f11.9)') 
     &    'The estimated average time per tic/toc pair is ',
     &    tictoctime1 + alpha*weighted_calls(0,1) / nested_calls(0,1)
     &                + (ioverhead2*time1 -ioverhead1*time2) /
     &                  (ioverhead2-ioverhead1) * ncalls(0,1)
c
      write (stdout,'(tr10,a,f11.9)') 
     &    'The estimated average time per tic/toc pair is ',
     &                      - (  - tictoctime1 * nested_calls(0,1)
     &                           - alpha*weighted_calls(0,1)
     &                           - pushpopextra1 * npops(0,1)
     &                           - pushpopinnerextra1 * npushes(0,1)
     &       - ((ioverhead2-i)*time1 + (i-ioverhead1)*time2) /
     &             (ioverhead2-ioverhead1) * (ncalls(0,1)+nresumes(0,1))
     &                       ) / nested_calls(0,1)
# endif
c
c     overflow warnings
c
      if (ngroup_overflow) then
        write (stdout,'(//,a,/,a,i3)')
     &       'WARNING:  Too many timer groups.  Some were lost',
     &       '          Increase ntimers=',ntimers
      end if
      do i=1,ntimers
        if (group_overflow(i)) then
          write (stdout,'(//,a,a,a,/,a,i3)')
     &    'WARNING:  Too many subgroups in group ',timer(i),
     &    '. Some were lost',
     &       '          Increase nsubtimers=',nsubtimers
        end if
      end do
#endif
c
      return
      end



#ifdef timing
      subroutine findtimer (name, i, error)
c=======================================================================
c     look up/insert a group timer name
c
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
# include "timer.h"
      character*(*) name
      character*40 name40
      logical error

      name40 = name
      do i=1,ntimers
        if (timer(i) .eq. name40) then
          error = .false.
          return
        else if (timer(i) .eq. 'unused') then
          timer(i) = name
          call findsubtimer ('total', i, j, error)
          error = .false.
          return
        else if (timer(i) .eq. 'don''t use') then
          error = .true.
          return
        end if
      end do
      error = .true.
      return
      end


      subroutine findsubtimer (subname, i, j, error)
c=======================================================================
c     look up/insert a subgroup timer name
c
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
# include "timer.h"
      character*(*) subname
      character*40 subname40
      logical error

      subname40 = subname
      do j=1,nsubtimers
        if (subtimer(i,j) .eq. subname40) then
          error = .false.
          return
        else if (subtimer(i,j) .eq. 'unused') then
          subtimer(i,j) = subname
          error = .false.
          return
        else if (subtimer(i,j) .eq. 'don''t use') then
          error = .true.
          return
        end if
      end do
      error = .true.
      return
      end


      subroutine show_tictoc(variety, i, j, s)
c=======================================================================
c     show requests to start or stop a timer with approximate time.
c     this routine is also used to show implied suspends and resumes.
c
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
# include "stdunits.h"
# include "timer.h"
      character*(*) variety
# ifdef show_tictocs
      if (print_tictocs) then
         if (s .eq. 0.0) then
           ss = second()
         else
           ss = s
         end if
        write (stdout,'(a12,tr1,a17,tr1,a32,a,f10.6)') variety,timer(i)
     &,       subtimer(i,j), ' at ', ss - starttime(0,1)
      end if
# endif
      return
      end 

      subroutine tictoc_nesting (variety,i,j,s)
c=======================================================================
c     manage stacks of suspended timers for each group.
c     keep statistics used in time adjustments.
c
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
# include "timer.h"
      character*(*) variety
      if (variety .eq. 'toc') then
c
c       if when a toc (i,j) happens, there is a timer(i,jj) in suspended
c       execution on the stack,resume it after completing the toc (i,j).
c
        call pop_timer (i,jj,s)
        active(i) = jj
c
c       update nested_calls for all active timers and subtimers
c       since tics must be completed by tocs, count only tocs, so
c       nested_calls(ii,jj) = 1 means one nested pair of tic/toc calls.
c
        nested_calls(0,1) = nested_calls(0,1) + 1
        weighted_calls(0,1) = weighted_calls(0,1) + i+j-4
        if (active(i) .gt. 0) then
          npops(0,1) = npops(0,1) + 1
        end if
        do ii=1,ntimers
          if (active(ii) .gt. 0) then
            nested_calls(ii,1) = nested_calls(ii,1) + 1
            weighted_calls(ii,1) = weighted_calls(ii,1) + i+j-4
            jj = active(ii)
            nested_calls(ii,jj) = nested_calls(ii,jj) + 1
            weighted_calls(ii,jj) = weighted_calls(ii,jj) + i+j-4
            if (ii .eq. i) then
              npops(ii,jj) = npops(ii,jj) + 1
              npops(ii,1) = npops(ii,1) + 1
              npushes(i,j) = npushes(i,j) + 1
              npushes(i,1) = npushes(i,1) + 1
            end if
          end if
        end do
      end if
c
c     if a tic (i,j) happens and another timer (i,jj) in the same group
c     is currently active, then toc timer(i,jj) and push jj on stack(i).
c     tic is counting on check_nesting to return a value for s, its
c     starttime, either from the clock or the toc time of pushed
c     timer (i,jj)
c
      if (variety .ne. 'toc') then
        if (active(i) .gt. 0) then
          jj = active(i)
          active(i) =  j
          call push_timer (i,j,jj,s)
c         if stack is full, push_timer may reject push and 
c         restore active(i)=jj
        else
          active(i) =  j
          s = second()
        end if
      end if
c
      return
      end



      subroutine inc_timer_warn
c=======================================================================
c     increment number of warnings and turn off printing at 1000
c
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
# include "stdunits.h"
# include "timer.h"
      ntictoc_errors = ntictoc_errors + 1
      if (ntictoc_errors .gt. max_timer_warnings .and. 
     &                              print_timer_warnings) then
        write (stdout,'(//,a,i5,a)') '*****TIC/TOC WARNINGS OFF:  ',
     &       max_timer_warnings, 'WARNINGS issued about tic/toc nesting'
      print_timer_warnings = .false.
      end if
      return
      end



      subroutine push_timer (i, j, jj, s)
c=======================================================================
c     toc timer(i,jj) and push it on stack(i)
c
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
# include "stdunits.h"
# include "timer.h"
      character*40 t40, st40
c
c     if there is room on the stack,
c       toc timer (i,jj) and push it on stack(i).
c       the stoptime s is returned to use as a starttime for the
c       tic (i,j) that causes the push.
c     else
c       refuse the push and set the current subtimer j to the one jj
c       that was supposed to be pushed on the stack.
c     end if

      if (istackptr(i) .lt. maxstack) then
        istackptr(i) = istackptr(i) + 1
        jstack(i,istackptr(i)) = jj
        s = second()
        elapsedtime(i,jj) = elapsedtime(i,jj) + s - starttime(i,jj)
        elapsedtime(i,1) = elapsedtime(i,1) + s - starttime(i,jj)
#ifdef show_tictocs
        call show_tictoc('suspend',i,jj,s)
#endif
      else
        write (stdout,'(a)') 'WARNING:    Tic/tocs nested too deeply'
        t40 = timer(i)
        st40 = subtimer(i,j)
        write (stdout,'(tr12,a,a,a,a,a)')
     &           'timer (',t40(1:len_trim(t40)),', ',
     &           st40(1:len_trim(st40)),') refused'
        write (stdout,'(tr12,a)') 'Increase maxstack above ', maxstack
        call inc_timer_warn
        j = jj
        active(i) = j
      end if
c
      return
      end


      subroutine pop_timer (i, jj, s)
c=======================================================================
c     pop a timer(i,jj) off stack(i) and restart it.
c
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
# include "timer.h"
c
c     pop a timer (i,jj) off stack(i) and restart it.
c     the starttime s originated as the stoptime of the
c     toc (i,j) that initiated the pop.
c
      if (istackptr(i) .gt. 0) then
        jj = jstack(i,istackptr(i))
        istackptr(i) = istackptr(i) - 1
        starttime(i,jj) = s
#ifdef show_tictocs
        call show_tictoc('resume',i,jj,s)
#endif
      else
        jj = 0
      end if
c
      return
      end
#endif


#ifdef cray_t3e_disabled
!not required: Balaji (vb@gfdl.gov)
      function second()
c
c=======================================================================
c     read a system clock and return a time in real seconds
c     author: Thomas Neumann  e-mail=> Thomas.Neumann@io-warnemuende.de
c=======================================================================
c
# include "timer.h"
      dimension ihms(3)
c
      call secondr(second)
      if (second .lt. prev_second) then
        second = second + 86400.0
      end if
      prev_second = second
      return
      end
#endif

#ifdef sgi_mipspro
      function second()
c
c=======================================================================
c     read a system clock and return a time in real seconds
c
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
c
# include "timer.h"
      dimension ihms(3)
c
      call itime (ihms)
      second = 3600.0 * ihms(1) + 60.0 * ihms(2) + ihms(3)
      if (second .lt. prev_second) then
        second = second + 86400.0
      end if
      prev_second = second
      return
      end
#endif

#ifdef rs6000
      function second()
c
c=======================================================================
c     read a system clock and return a time in real seconds             
c
c     author: Martin Schmidt  e-mail=> Martin.Schmidt@io-warnemuende.de
c=======================================================================
c
      type iar
         sequence
            integer*4 ihr
            integer*4 imin
            integer*4 isec
         end type
       type (iar) ihms
# include "stdunits.h"
       call itime_ (ihms)
       second = 3600.0 * ihms%ihr + 60.0 * ihms%imin + ihms%isec
       if (second .lt. prev_second) then
         second = second + 86400.0
       end if
       prev_second = second
       return
       end
#endif
#if !defined SGICRAY && !defined rs6000
      function second()
c
c=======================================================================
c     read a system clock and return a time in real seconds
c
c     author:     supply your system clock function here
c=======================================================================
c
# include "stdunits.h"
      write (stdout,*) ' '
      write (stdout,*)
     &   ' Please supply your system clock function in "timer.F"'
      write (stdout,*) ' under "function second()"'
      write (stdout,*) ' '
      second = 0.0
      return
      end
#endif


