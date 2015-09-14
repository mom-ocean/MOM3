#ifdef test_vmix
      program vmix1d
c
c=======================================================================
c
c
c     To test various  vertical mixing schemes in a 
c     simple one dimensional model (at one latitude and longitude
c     within the MOM grid): 
c     
c        1) setup the grid. (see grids.F)
c
c        2) compile and run this module using the script "run_vmix"
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use memory_window
      use arrays_2d
      use time_manager_mod

      logical error, vmixset, mixing
      parameter (ifdmax=100)
      character*30 cifdef(ifdmax)
# include "param.h"
# include "accel.h"
# include "coord.h"
# include "grdvar.h"
# include "iounit.h"
# include "mw.h"
# include "rinum.h"
# include "scalar.h"
# include "state.h"
# include "switch.h"
# include "tmngr.h"
# include "vmixc.h"
      dimension tt(km,2)
# include "fdifm.h"
# include "fdift.h"
# include "taskrows.h"
c
      allocate (kmt(imt,jmt))
      allocate (cori(jmt,2))
#ifdef kppvmix
      jstask=1
      jetask=jmt
      allocate (hblt(imt,jstask:jetask))
      allocate (hmxl(imt,jstask:jetask))
#endif

      stdout = 6
c
      write (stdout,'(//,25x,a/)')
     &' T E S T I N G   A   1-D   V E R S I O N   O F   V M I X'
# ifdef rot_grid
c
      write (stdout,*) '  "rot_grid" is incompatible with this driver' 
      if (imt .eq. imt) call abort()
      end
# endif
c
c     initialize physical constants
c
      radius   = 6371.0e5
      grav     = 980.6
      rho0     = 1.035
c
c     initialize time step accelerators
c
      do k=1,km
        dtxcel(k) = 1.0
      enddo
c
c-----------------------------------------------------------------------
c     set up the grids in x (longitude), y (latitude), and z (depth)
c     corresponding to Arakawa "b" gird system
c-----------------------------------------------------------------------
c
      call grids
c
c-----------------------------------------------------------------------
c     compute density coefficients based on depth of grid points
c-----------------------------------------------------------------------
c
      call eqstate (zt, km, ro0, to, so, c, tmink, tmaxk, smink, smaxk)
c
c-----------------------------------------------------------------------
c     initialize clock and calendar
c-----------------------------------------------------------------------
c
      year0  = 1994
      month0 = 1
      day0   = 1
c
      hour0  = 0
      min0   = 0
      sec0   = 0
      eqyear = .false.
      eqmon  = .false.
      refrun  = .false.
      refinit = .true.
      refuser = .false.
      if (refuser) then
        ryear  = 1900
        rmonth = 1
        rday   = 1
        rhour  = 0
        rmin   = 0
        rsec   = 0
      end if
c      runlen   = 120.0
      runlen   = 20.0
      rununits = 'days'
      idayrestart  = 0
      msrestart = 0
c      dtts    = 3600.0
c      dtuv    = 3600.0
      dtts    = 300.0
      dtuv    = 300.0
      segtim  = dtts/86400.0
c
c     note: nmix is set but mixing timesteps are ignored if we
c     do a robert time filter every time step
c
      nmix    = 11
c#define robert_time_filter
#ifdef robert_time_filter
        write (stdout,*) ' Note: robert time filter is applied every ts'
#else
      write(stdout,*)" Note: a forward ts is taken every ",nmix," steps"
#endif
      itt=1
      call tmngri (year0, month0, day0, hour0, min0, sec0
     &,            ryear, rmonth, rday, rhour, rmin, rsec
     &,            idayrestart, msrestart
     &,            runlen, rununits, rundays, dtts)
c
c-----------------------------------------------------------------------
c     prescribe some initial stratification
c-----------------------------------------------------------------------
c
c      z0 = 30.0e2
c      hh = 80.0e2
      z0 = 30.0e2
      hh = 200.0e2
      zm = zt(km)
      t0 = 7.5
      t1 = 10.0
      write (stdout,'(/,10x, a/)') 'Initial conditions'
      do k=1,km
        tt(k,1) = t0*(1.0 - tanh((zt(k)-hh)/z0)) + t1*(1.0-zt(k)/zm)
	tt(k,2) = 0.0349 - 0.035
#ifdef ncar_kpp_test_case_pure_convection
        if (zt(k).lt.25000.0) then
          tt(k,1) = 22.0 - 0.1*(zt(k)-zt(1))/100.0
          tt(k,2) = 0.0349 - 0.035
        else
          tt(k,1) = tt(k-1,1)
          tt(k,2) = 0.0349 - 0.035
        end if
#endif
#ifdef ncar_kpp_test_case_wind_mixing
        if (k .le. 10) then
          if (zt(k) .lt. 2000.0) then
            tt(k,1) = 13.0
            depth_mixed_layer=zt(k)
          else 
            tt(k,1) = 13.0 - (13.0-7.0)
     &           *(zt(k)-depth_mixed_layer)
     &           /(zt(10) - depth_mixed_layer)
          end if
        else
            tt(k,1) = 7.0
        end if
	tt(k,2) = 0.0349 - 0.035
#endif
c
	write (stdout, '("k=",i4," zt(k)=",g12.4
     &       ," T=",g12.4," S=",g12.4)') k,zt(k),tt(k,1),tt(k,2)

      enddo
c
c     set time levels
c
      taum1 = -1
      tau   =  0
      taup1 = +1
c
c     set I.C. for u,v,t,s  and set the land/sea masks to all sea (1.0)
c
      do j=1,jmw
        do i=1,imt
          kmt(i,j) = km
	  do k=1,km
	    umask(i,k,j) = 1.0
	    tmask(i,k,j) = 1.0
	    do n=1,2
	      u(i,k,j,n,tau)   = 0.0
	      u(i,k,j,n,taum1) = 0.0
	      u(i,k,j,n,taup1) = 0.0
	      t(i,k,j,n,tau)   = tt(k,n)
	      t(i,k,j,n,taum1) = tt(k,n)
	      t(i,k,j,n,taup1) = tt(k,n)
	    enddo
	  enddo
	enddo
      enddo
c
c-----------------------------------------------------------------------
c     pick a latitude for the one dimensional model
c-----------------------------------------------------------------------
c
      alat = 45.0
      alat = 25.0
      jrow = indp (alat, yu, jmt)
c
c-----------------------------------------------------------------------
c     set the coriolis factors for this latitude and choose
c     a centered coriolis term (explicit coriolis)
c-----------------------------------------------------------------------
c
      omega        = pi/43082.0
      do j=1,jmt
        cori(j,1) = c2*omega*sine(jrow)
        cori(j,2) = -c2*omega*sine(jrow)
      enddo
      gcor         = 1.0
c
c-----------------------------------------------------------------------
c     initialize the vertical mixing scheme
c-----------------------------------------------------------------------
c
      error   = .false.
      nifdef  = 0
      vmixset = .false.
      call ioinit
c
#ifdef ppvmix
      call ppmixi (error, cifdef, ifdmax, nifdef, vmixset)
#endif
c
#ifdef kppvmix
      call kppmixi (error, cifdef, ifdmax, nifdef, vmixset)
#endif
      if (error) call abort()
c
c-----------------------------------------------------------------------
c     integrate equations for all k at one point (i,j)
c     set "joff=1" so j=1 is not the southern boundary (land point)
c     "j" is the row in the MW and "jrow" is the latitude row on disk
c      j=2 corresponds to the latitude of jrow
c-----------------------------------------------------------------------
c
      is   = imt/2
      ie   = imt/2
      joff = 1
c
      write(stdout, '(/a,g14.7/)') '  1-D model latitude = ', yu(jrow)
c
      do itt=1,100000
c
c       set all switches to control when things happens 
c
        call tmngr (dtts)
	call set_time_switches
	mixing = .not. leapfrog
c
c       broadcast the central point to cover the memory window
c
        icent = is
        jcent = 2
        do n=1,2
	  do j=1,jmw
	    do i=is-1,ie+1
	      do k=1,km
	        u(i,k,j,n,tau) = u(icent,k,jcent,n,tau)
		u(i,k,j,n,taup1) = u(icent,k,jcent,n,taup1)
	        t(i,k,j,n,tau) = t(icent,k,jcent,n,tau)
		t(i,k,j,n,taup1)   = t(icent,k,jcent,n,taup1)
	      enddo
	    enddo
	  enddo
	enddo
#ifdef robert_time_filter
c
c       robert time filter for velocity
c
        smooth=.01
        do n=1,2
	  do j=1,jmw
	    do i=is-1,ie+1
	      do k=1,km
		u(i,k,j,n,tau)   = u(i,k,j,n,tau)
     &          + smooth*(0.5*(u(i,k,j,n,taup1) + u(i,k,j,n,taum1))
     &          - u(i,k,j,n,tau)) 
		t(i,k,j,n,tau)   = t(i,k,j,n,tau)
     &          + smooth*(0.5*(t(i,k,j,n,taup1) + t(i,k,j,n,taum1))
     &          - t(i,k,j,n,tau)) 
	      enddo
	    enddo
	  enddo
	enddo
#endif
c
c       move "tau"   variables ==> "tau-1" positions
c            "tau+1" variables ==> "tau"   positions
c
        do n=1,2
	  do j=1,jmw
	    do i=is-1,ie+1
	      do k=1,km
	        u(i,k,j,n,taum1) = u(i,k,j,n,tau)
		u(i,k,j,n,tau)   = u(i,k,j,n,taup1)
	        t(i,k,j,n,taum1) = t(i,k,j,n,tau)
		t(i,k,j,n,tau)   = t(i,k,j,n,taup1)
	      enddo
	    enddo
	  enddo
	enddo
c
#ifdef robert_time_filter
        dtu = 2.0*dtuv
	dtt = 2.0*dtts
#else
        if (mixing) then
c
c       move "tau"   variables ==> "tau-1" positions and set timestep
c
          do n=1,2
	    do j=1,jmw
	      do i=is-1,ie+1
	        do k=1,km
	          u(i,k,j,n,taum1) = u(i,k,j,n,tau)
	          t(i,k,j,n,taum1) = t(i,k,j,n,tau)
	        enddo
	      enddo
	    enddo
	  enddo
c
          dtu = dtuv
	  dtt = dtts
        else
          dtu = 2.0*dtuv
	  dtt = 2.0*dtts
        endif
#endif
c
c       prescribe some surface and bottom fluxes
c       wind components are cgs
c
	windx = -7.0e2
	windy = 3.0e2
	ro    = 1.0e-3
	cd    = 1.2e-3
	wind  = sqrt(windx**2 + windy**2)
	k     = 1
	airt  = tt(k,1) - 1.0
	damp  = 1.0/(3.0*86400.0)
        do j=1,jmw
          do i=is-1,ie+1
            stf(i,j,1) = damp*(airt - t(i,k,j,1,taum1))*dzt(1)
	    stf(i,j,2) = 0.0
	    btf(i,j,1) = 0.0
	    btf(i,j,2) = 0.0
            smf(i,j,1) = ro*cd*wind*(windx - u(i,k,j,1,taum1))
	    smf(i,j,2) = ro*cd*wind*(windy - u(i,k,j,2,taum1))
	    bmf(i,j,1) = 0.
	    bmf(i,j,2) = 0.
#ifdef ncar_kpp_test_case_pure_convection
            stf(i,j,1) = -100.0/41860.0
	    stf(i,j,2) = 0.0
	    btf(i,j,1) = 0.0
	    btf(i,j,2) = 0.0
            smf(i,j,1) = 0.0
	    smf(i,j,2) = 0.0
	    bmf(i,j,1) = 0.
	    bmf(i,j,2) = 0.
#endif
#ifdef ncar_kpp_test_case_wind_mixing
            stf(i,j,1) = 0.0
	    stf(i,j,2) = 0.0
	    btf(i,j,1) = 0.0
	    btf(i,j,2) = 0.0
            smf(i,j,1) = ro*cd*wind*(windx - u(i,k,j,1,taum1))
	    smf(i,j,2) = ro*cd*wind*(windy - u(i,k,j,2,taum1))
	    bmf(i,j,1) = 0.
	    bmf(i,j,2) = 0.
#endif
	  enddo
        enddo
c
c-----------------------------------------------------------------------
c       calculate mixing coefficients for i=is..ie and j=2..jmw-1
c       (variables are loaded for j=1,jmw. calculations from j=2,jmw-1.
c       limits is-1,ie+1 are needed to insure that data used for 
c       visc_cbu calculation is correct.)
c        output is:
c       "diff_cbt" = diffusive coeff at bottom of "t" cell
c       "visc_cbu" = diffusive coeff at bottom of "u" cell
c-----------------------------------------------------------------------
c
#ifdef ppvmix
        call ppmix (joff, 1, jmw, is, ie)
#endif
#ifdef kppvmix
        call kppmix (joff, 1, jmw, is, ie)
#endif
c
c-----------------------------------------------------------------------
c                   SOLVE THE 1-D TRACER EQUATIONS
c       construct diffusive flux for i=is..ie and and solve for j=2
c       "diff_fb" = diffusive flux at bottom of T cell 
c-----------------------------------------------------------------------
c
# ifdef kppvmix
c
c--------------------------------------------------------------------
c     add nonlocal transport term ghats to source
c--------------------------------------------------------------------
c
	do n=1,2
          nt2 = min(n,nsmix)
          do j=2,2
	    k=1
	    do i=is-1,ie+1
	      source(i,k,j) =  stf(i,j,n) * dztr(k)
     &                      * (- diff_cbt(i,k,j,nt2) * ghats(i,k,j))
	    enddo
            do i=is-1,ie+1
              do k=2,kmm1
	        source(i,k,j) = stf(i,j,n) * dztr(k)
     &                       * (diff_cbt(i,k-1,j,nt2) * ghats(i,k-1,j)
     &                        - diff_cbt(i,k  ,j,nt2) * ghats(i,k  ,j))
	      enddo
              source(i,km,j)=0.0
            enddo
          enddo
        enddo
# endif


	do n=1,2
	  nt2 = min(n,nsmix)
          do j=2,2
            do k=1,km-1
              do i=is,ie
                diff_fb(i,k,j) = diff_cbt(i,k,j,nt2)*(
     &                    t(i,k,j,n,taum1) - t(i,k+1,j,n,taum1))*dzwr(k)
              enddo
            enddo
            do i=is,ie
              diff_fb(i,0,j)  = stf(i,j,n)
              diff_fb(i,km,j) = btf(i,j,n)
            enddo
c
c           solve tracer eqns for each depth
c
	    i = is
	    do k=1,km
              t(i,k,j,n,taup1) = t(i,k,j,n,taum1)
     &             + dtt*(
     &             (diff_fb(i,k-1,j)-diff_fb(i,k,j))*dztr(k)
# ifdef source_term
     &             + source(i,k,j)
# endif
     &             )
c
c             set lateral boundaries
c
              do jj=1,jmw
                do ii=is-1,ie+1
                  t(ii,k,jj,n,taup1) = t(i,k,j,n,taup1)
	        enddo
	      enddo
            enddo                   
          enddo                   
        enddo
c
c-----------------------------------------------------------------------
c                   SOLVE THE 1-D MOMENTUM EQUATIONS
c       construct diffusive flux for i=is..ie and and solve for j=2
c       "diff_fb" = diffusive flux at bottom of U cell 
c-----------------------------------------------------------------------
c
	do n=1,2
          do j=2,2
            do k=1,km-1
              do i=is,ie
                diff_fb(i,k,j) = visc_cbu(i,k,j)*
     &                          (u(i,k,j,n,taum1) - u(i,k+1,j,n,taum1))
     &                          *dzwr(k)
              enddo
            enddo
            do i=is,ie
              diff_fb(i,0,j)  = smf(i,j,n)
              diff_fb(i,km,j) = bmf(i,j,n)
            enddo
c
c           solve momentum eqns for each depth
c
	    i = is
	    do k=1,km
              u(i,k,j,n,taup1) = u(i,k,j,n,taum1) + dtu*( 
     &                           CORIOLIS(i,k,j,n)
     &           +  (diff_fb(i,k-1,j) - diff_fb(i,k,j))*dztr(k)
     &             )
c
c             set lateral boundaries
c
              do jj=1,jmw
                do ii=is-1,ie+1
                  u(ii,k,jj,n,taup1) = u(i,k,j,n,taup1)
	        enddo
	      enddo
            enddo                   
          enddo
	enddo
c
c       show some results
c
        if (eoday) then
          i = is
	  j = 2
          write(stdout,'(/a,i6,a,a)')
     &    ' End of day at itt=',itt, '  ',stamp
          write(stdout,'(5(a,g11.4))')
     &    ' taux=',smf(i,j,1), ' tauy=',smf(i,j,2), ' hflx=',stf(i,j,1)
#ifdef kppvmix
     &,   ' hmxl=',hmxl(i,j+joff), ' hblt=',hblt(i,j+joff)
#endif
          do k=1,km
	    write(stdout,'(a,i2,15(a,g11.4,1x))') 
     &     'k=',k,' u=',u(i,k,j,1,taup1),' v=',u(i,k,j,2,taup1), ' t='
     &,     t(i,k,j,1,taup1), ' diff_cbt=',diff_cbt(i,k,j,1)
     &,     ' visc_cbu=', visc_cbu(i,k,j), 'riu=',riu(i,k,j)
     &,     'rit=',rit(i,k,j)
#ifdef kppvmix
     &,     'ghats=',ghats(i,k,j),'diff_cbt=',diff_cbt(i,k,j,1)
     &,     'ghats*diff_cbt=',ghats(i,k,j)*diff_cbt(i,k,j,1)
#endif
	  enddo
        endif
        if (eorun) stop
      enddo
      end

# ifdef timing
      subroutine tic (a, b)
      character*(*) a, b
      entry toc (a, b)
      return
      entry ticr (a, b)
      return
      end
# endif
#else
      subroutine vmix1d
      return
      end
#endif
