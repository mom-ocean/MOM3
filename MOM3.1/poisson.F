#ifdef test_poisson
# define driver_only
# include "congrad.F"
# include "grids.F"
# include "iomngr.F"
# include "size_check.F"
# include "timer.F"
# include "topog.F"
# include "tropic.F"
# include "util.F"
      program driver
c
c=======================================================================
c
c                    P O I S S O N      M O D U L E
c
c     This module solves TWO TYPES of elliptic equations.
c
c***********************************************************************
c
c     TYPE I is a STREAM FUNCTION equation of the form:
c
c     grad{(1/h)*grad(dpsi)} - 2dt*acor*{grad(f/h) X grad(dpsi)} = forc
c
c     solved for "dpsi" with "no slip"  boundary conditions where...
c
c     h    = h(x,y) = topography on "u,v" points
c     forc = forc(x,y) = forcing on "t" points
c     f    = f(y) = coriolis factor 2*omega*sin(phi) on "u,v" points
c     dpsi = psi(tau+1) - psi(tau-1) where psi is the stream function
c            defined on "t" points. tau refers to the time level.
c     2dt  = twice the time step
c     acor = (0..1) for controlling the implicit treatment of the
c            coriolis term
c
c
c***********************************************************************
c
c     TYPE II is a PRESSURE EQUATION of the form:
c
c                    grad{h*grad(surfpress)}  = forc
c
c          with "free slip" boundary conditions where:
c
c     h    = h(x,y) = topography on "u,v" points
c     forc = forc(x,y) = forcing on "t" points
c
c
c     To run a stand-alone test of the Poisson solvers, use the
c     "run_poisson" script
c
c
c     author:             c.h.goldberg      e-mail=> chg@gfdl.gov
c     based on code by:   r.c.pacanowski    e-mail=> rcp@gfdl.gov
c
c=======================================================================
c
#include "size.h"
#include "stdunits.h"
#include "coord.h"
c
      dimension dxt(imt), dyt(jmt), dxu(imt), dyu(jmt), csu(jmt)
      dimension cst(jmt)
      dimension phi(jmt), phit(jmt)
      dimension f(jmt)     
      dimension h(imt,jmt), hr(imt,jmt)
      dimension cf(imt,jmt,-1:1,-1:1)
      dimension forc(imt,jmt), dpsi(imt,jmt)
      dimension res(imt,jmt), kmt(imt,jmt), guess(imt,jmt)
      dimension iofs(mnisle), iperm(maxipp), jperm(maxipp)
      dimension nippts(mnisle)
      dimension surfpres(imt,jmt)
      dimension seaht(imt,jmt)
      dimension solution(imt,jmt)
#if defined partial_cell
      dimension ht(imt,jmt)
#endif

      dimension du(imt,jmt), dv(imt,jmt)
c
      dimension zu(imt,jmt,2)
      dimension zumod(imt,jmt,2)
      dimension tauu(imt,jmt), tauv(imt,jmt)
      dimension psx(imt,jmt), psy(imt,jmt)
c
      dimension map (imt,jmt)
      logical imask (-mnisle:mnisle)
c
      parameter (maxtest=100)
      character*16 testlist(maxtest)
      character*32 atest
c
      logical converged, solver, noslip
      logical show_sea_heights
      character*16 variable
c
      stdout = 6
c
# ifdef cyclic
      write (stdout,*) ' option "cyclic" enabled'
# endif
# ifdef symmetry
      write (stdout,*)
     &  ' option "symmetry" enabled for equator at j=',jmt-1
# endif
# ifdef rot_grid
      write (stdout,*)
     & ' option "rot_grid" is incompatible with this driver'
      if (imt .eq. imt) call abort()
      end
# endif
c
c-----------------------------------------------------------------------
c     set global domain
c-----------------------------------------------------------------------
c
      js = 1
      je = jmt
c
c-----------------------------------------------------------------------
c     initialize individual and sor optimizaion test lists
c-----------------------------------------------------------------------
c
      ntest  = 0
c
c-----------------------------------------------------------------------
c     ==> USER INPUT SECTION #1:  choose solvers to test
c-----------------------------------------------------------------------
c
c     test the stream function formulation (no slip assumed)
c 
      call select ('stream function', testlist, ntest)
c
c     turn on 5 point numerics
c
      call select ('5 point', testlist, ntest)
      call select ('congrad', testlist, ntest)
c
c     try 9 point solvers with 5 point coefficients (zero corner terms)
c
c     call select ('59 test', testlist, ntest)
c     call select ('congrad', testlist, ntest)
c
c     turn on 9 point numerics
c
c     call select ('9 point', testlist, ntest)
c     call select ('congrad', testlist, ntest)
c
c     test the surface pressure formulation using 9 point numerics
c     (surface pressure automatically turns on 9 point numerics)
c     (5 point numerics do not conserve energy)
c
      call select ('surface pressure', testlist, ntest)
      call select ('congrad', testlist, ntest)
c
c-----------------------------------------------------------------------
c     choose individual solver test parameters
c     c2dtsf = 2*dt
c     acor = (0..1) for implicit treatment of coriolis term
c     mxscan = maximum allowable iterations
c     precision = minimum solver step / expected_dpsimax
c     crit   = estimated maximum error in dpsi -- all solvers
c-----------------------------------------------------------------------
c
      c2dtsf = 2*3600.0
      acor   = 0.0
      mxscan = 6000
      precision            = 1.0e-7
      expected_dpsimax     = 1.0e12
      tolrsf              = precision * expected_dpsimax
      expected_surfpresmax = 1.0e3
      tolrsp              = precision * expected_surfpresmax
c
c-----------------------------------------------------------------------
c     ==> END USER INPUT SECTION #1
c-----------------------------------------------------------------------
c
      call inittimers
c
c-----------------------------------------------------------------------
c     calculate coordinates for "t" and "u" grid cells.
c-----------------------------------------------------------------------
c
      maxlen = max(imt, jmt, km)
      call gcoord (maxlen, imt2, jmt2, km2
     &,            dxtdeg, dytdeg, dxudeg, dyudeg
     &,            dzt, dzw, xt, xu, yt, yu, zt, zw)
c
c     check that returned grid sizes match those in file "size.h"
c
      call size_check (imt2, jmt2, km2, 'driver for poisson', 'stop')
c
c-----------------------------------------------------------------------
c     define some constants
c-----------------------------------------------------------------------
c
      pi     = 4.0*atan(1.0)
      omega  = pi/43082.0
      crad   = pi/180.0
      degtcm = pi/180.0
c
c-----------------------------------------------------------------------
c     convert grid resolution to cm
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        dyt(jrow) = dytdeg(jrow)*degtcm
        dyu(jrow) = dyudeg(jrow)*degtcm
      enddo
c
      do i=1,imt
        dxt(i) = dxtdeg(i)*degtcm
        dxu(i) = dxudeg(i)*degtcm
      enddo
#ifdef cyclic
      dxt(1)   = dxt(imt-1)
      dxt(imt) = dxt(2)
      dxu(1)   = dxu(imt-1)
      dxu(imt) = dxu(2)
#endif
c
c-----------------------------------------------------------------------
c     build latitude arrays for grid
c-----------------------------------------------------------------------
c
      do j=1,jmt
        if (yu(j) .eq. 90.0) yu(j) = yu(j) - 1.e-4
        phi(j) = yu(j)*crad
	csu(j) = cos(phi(j))
	phit(j) = yt(j)*crad
	cst(j) = cos(phit(j))
      enddo
c
c-----------------------------------------------------------------------
c       build coriolis term
c-----------------------------------------------------------------------
c
        do j=1,jmt
	  f(j) = 2.0*omega*sin(phi(j))
	enddo
c
c-----------------------------------------------------------------------
c     generate "kmt" = number of vertical levels on "t" cells
c-----------------------------------------------------------------------
c
      call topog (kmt, map, xt, yt, zt, xu, yu, zw, imt, jmt, km, ht)
c
      write (stdout,*) ' '
      write (stdout,*) ' Land/sea mask for testing POISSON MODULE:'
      write (stdout,*) ' '
c
      write (stdout,*) ' '
      do n=1,ntest
        write (stdout,'(a,a)')  '  ==> Selecting ', testlist(n)
      enddo
      write (stdout,*) ' '
c
c-----------------------------------------------------------------------
c     calculate island perimeters
c-----------------------------------------------------------------------
c
      call isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle
     &,                    imt, jmt, km, mnisle, maxipp)
c
c-----------------------------------------------------------------------
c     count ocean points
c-----------------------------------------------------------------------
c
      nocean = 0
      do j=2,jmt-1
        do i=2,imt-1
          if (map(i,j) .le. 0) nocean = nocean + 1
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     set mask for island perimeters on which to perform calculations
c     imask(-n) = .false.  [no equations ever on dry land mass n]
c     imask(0)  = .true.   [equations at all mid ocean points]
c     imask(n)  = .true./.false [controls whether there will be
c                                equations on the ocean perimeter of
c                                land mass n]
c     note: land mass 1 is the northwest-most land mass
c           usually includes the "north pole", and at low resolutions,
c           the "main continent"
c     for the numbering of the other landmasses, see generated map(i,j)
c-----------------------------------------------------------------------
c
      do isle=-mnisle,mnisle
        if (isle .ge. 0 .and. isle .le. nisle) then
          imask(isle) = .true.
        else
          imask(isle) = .false.
        endif
      enddo
c
c-----------------------------------------------------------------------
c     ==> USER INPUT SECTION #2: user-specified changes to island mask 
c-----------------------------------------------------------------------
c    
c     there are problems if imask is set .true. for a nonexistent
c     island. it is recommended that no changes be made in this section
c
c      imask(1) = .false.
c      imask(2) = .false.
c
c-----------------------------------------------------------------------
c     ==> END USER INPUT SECTION #2:
c-----------------------------------------------------------------------
c    
c     print diagnostic information
c
      write (stdout,*) ' => Number of ocean points = ', nocean
      do isle=-mnisle,mnisle
        if (imask(isle)) then
	  if (isle .eq. 0) then
            write (stdout,'(a)') 
     &        '=> calculations enabled for mid ocean points'
	  else
            write (stdout,'(2a,i3)') 
     &        '=> calculations enabled for ocean ',
     &                      'perimeter of land mass',isle
	  endif
	endif
      enddo
      do isle=0,nisle
        if (.not. imask(isle)) then
            write (stdout,'(2a,i3)') 
     &       '=> calculations disabled for ocean ',
     &                      'perimeter of land mass',isle
	endif
      enddo
c
c-----------------------------------------------------------------------
c     ==> USER INPUT SECTION #3: "dpsi" normalization 
c-----------------------------------------------------------------------
c    
c     imain is the land mass on which dpsi is normalized to 0
c     if imain is 0, then dpsi is not normalized.
c
      imain = min(2,nisle)
      do isle=1,nisle
        if (.not.imask(isle)) imain = isle
      enddo
c
c-----------------------------------------------------------------------
c     ==> END USER INPUT SECTION #3:
c-----------------------------------------------------------------------
c    
      if (imain .gt. 0 .and. imain .le. nisle) then
        write (stdout,'(a,i4)')
     &   'dpsi normalized to zero on land mass',imain
      else if (imain .eq. 0) then
        write (stdout,*) 'no normalization on dpsi'
      else
        write (stdout,*)
     &          'ERROR: illegal value for choice of normalization ',
     &           'land mass, imain =', imain
      endif
      write (stdout,*) 
     &         ' (user may set "imain" to any valid land mass number)'
c
c-----------------------------------------------------------------------
c     calculate h and reciprocal h at u/v points
c-----------------------------------------------------------------------
c
      do j=1,jmt-1
        do i=1,imt-1
	  kmu = min(kmt(i,j), kmt(i+1,j), kmt(i+1,j+1), kmt(i,j+1))
	  if (kmu .ne. 0) then
	    h(i,j)  = zw(kmu)
	    hr(i,j) = 1.0/h(i,j)
	  else
	    h(i,j)  = 0.0
	    hr(i,j) = 0.0
	  endif
	enddo
      enddo
      call border (h, 'u even', js, je)
      call border (hr, 'u even', js, je)
c
c-----------------------------------------------------------------------
c     show the value of acor
c-----------------------------------------------------------------------
c
      write (stdout,*) ' '
      if (acor .ne. 0.0) then
        write (stdout,*)
     &      ' => Note: implicit coriolis term ... acor=',acor
      else
        write (stdout,*) ' => No implicit coriolis term'
      endif
      write (stdout,*) ' '
c
c-----------------------------------------------------------------------
c     construct a simplified zonal wind stress (dynes/cm**2)
c-----------------------------------------------------------------------
c
# define hellerman_winds
#ifdef hellerman_winds
c
c     Zonally averaged annual mean Hellerman & Rosenstein wind stress
c
      do j=1,jmt
        do i=1,imt
          if (h(i,j) .ne. 0) then
            ulat = phi(j) / crad
            call hellerman (ulat, taux, tauy)
            tauu(i,j) = taux 
            tauv(i,j) = tauy
          else
            tauu(i,j) = 0.0
            tauv(i,j) = 0.0
          endif
        enddo
      enddo
#else
c
c     constant west-east winds
c
      do j=1,jmt
        do i=1,imt
          if (j.ge.2 .and. j.le.jmt-2) then
            tauu(i,j) = 1.0
            tauv(i,j) = 0.0
          endif
        enddo
      enddo
#endif
c
c-----------------------------------------------------------------------
c     calculate vertically averaged forcing 
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do i=1,imt
          zu(i,j,1) = tauu(i,j)
          zu(i,j,2) = tauv(i,j)
        enddo
      enddo
      do j=1,jmt
        do i=1,imt
          zumod(i,j,1) = zu(i,j,1)*c2dtsf
          zumod(i,j,2) = zu(i,j,2)*c2dtsf
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     set initial guess/estimate
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do i=1,imt
          guess(i,j) = 0.0
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     initialize
c-----------------------------------------------------------------------
c
      show_sea_heights = .true.
      nisle_sp = 0
      nisle_sf = nisle
c
c=======================================================================
c     test individual solvers
c=======================================================================
c
      do n=1,ntest
        write (stdout,'(/,i2,a,a)')  n, '. Testing ', testlist(n)
        solver = .false.
c
        write (atest,'(i2,a,a)') n,'.',testlist(n)
#ifdef timing
        call tic ('solver test', atest)
#endif
c
        call do_test (testlist(n)
     &,                map, imask 
     &,                zu, hr, dxu, dyu, phi, csu, f
     &,                zumod 
     &,                forc
     &,                acor, omega, c2dtsf
     &,                cf
     &,                h 
     &,                mxscan, mscan, tolrsf, tolrsp
     &,                iperm, jperm, iofs, nisle, nippts
     &,                crit, estimated_error
     &,                mode, variable
     &,                converged
     &,                guess, dpsi, res
     &,                solution, surfpres
     &,                npt
     &,                noslip, nisle_sp, nisle_sf
     &,                solver)
c
        if (solver) then
c
c         estimate convergence rate and theoretical maximum change in
c         dpsi if iterated to infinity, assuming geometric convergence.
c
          write (stdout,'(a,a,a,e15.8)')
     &           'Estimated maximum error in ', variable, ' = ' 
     &,          estimated_error
c
          if (mode .eq. 1) then
c
c           stream function: zero dpsi on land mass imain
c
            if (imain .gt. 0) then
              dpsi1 = dpsi(iperm(iofs(imain)+1), jperm(iofs(imain)+1))
              call con_adjust (dpsi, dpsi1, map)
              call con_adjust (solution, dpsi1, map)
            endif
c
c           stream function: fill land areas with constant dpsi
c
            call fill_land (dpsi, map, noslip 
     &,                     nisle, iperm, jperm, iofs, nippts)
            call fill_land (solution, map, noslip 
     &,                     nisle, iperm, jperm, iofs, nippts)
c
c           construct delta u velocity and delta v velocity  
c
            call ddyu (dpsi, du, dyu)
	    do j=2,jmt-1
	      do i=2,imt-1
	        du(i,j) = -du(i,j)*hr(i,j)
	      enddo
	    enddo
            call border (du, 'u even', js, je)
c
            call ddxu (dpsi, dv, dxu, csu)
	    do j=2,jmt-1
	      do i=2,imt-1
	        dv(i,j) = dv(i,j)*hr(i,j)
	      enddo
	    enddo
            call border (dv, 'u odd', js, je)
          endif
c
          if (mode .eq. 2) then
c
c           make average surfpres = zero
c
	    call zero_level (surfpres, variable, map, dxt, dyt, cst)
c
c           calculate grad (surfpres)
c
            call ddxu (surfpres, psx, dxu, csu)
            call border (psx, 'u even', js, je)
	    call ddyu (surfpres, psy, dyu)
            call border (psy, 'u odd', js, je)
	  endif
c
c-----------------------------------------------------------------------
c         show a slice through solution
c-----------------------------------------------------------------------
c
          call slice (solution, variable, res, imt/2)
c
          if (mode .eq. 2) then
c
c            write (stdout,'(/a/)') '   psx'
c            call slice (psx, 'psx', res, imt/2)
c
c            write (stdout,'(/a/)') '   psy'
c            call slice (psy, 'psy', res, imt/2)
c
c           surface pressure solution is subject to a checkerboard
c           null space (especially at low resolutions)
c           remove the null space. however, global checkerboard
c           removal may not remove local checkerboard patterns.
c
	    call checkerboard (surfpres, map, js, je)
c
c           correct average surfpres to zero since checkerboard
c           correction can alter mean
c
	    call zero_level (surfpres, variable, map, dxt, dyt, cst)
c
c           show a slice through solution
c
            write (stdout,'(/a/)')
     &        ' surf pressure uncheckerboarded with zero mean'
            call slice (surfpres, variable, res, imt/2)
c
            if (show_sea_heights) then
c
c             convert surface pressure to virtual sea heights (cm)
c
              grav=980.6
              do j=1,jmt
                do i=1,imt
                  seaht(i,j) = surfpres(i,j) / grav
                  res(i,j)   = res(i,j) / grav
                enddo
              enddo
c
              write (stdout,'(/a/)') '  EQUIVALENT SEA HEIGHTS (cm) '
              call slice (seaht, 'sea height', res, imt/2)
            endif
          endif
        endif
c
#ifdef timing
        call toc ('solver test', atest)
#endif
      enddo
c
      call showtimers
      call abort()
      end







      subroutine select (test, testlist, ntest)
c
c     builds a list of solver tests
c
#include "stdunits.h"
      parameter (maxtest=100)
      character*(*) testlist(maxtest), test
c
      ntest = ntest + 1
      if (ntest .gt. maxtest) then
        write (stdout,*) ' ==> ERROR: too many tests of Poisson solvers'
        write (stdout,*) '     Increase maxtest above ',maxtest
        write (stdout,*) '     in main program and subroutine select'
        call abort()
      else
        testlist(ntest) = test
      endif
      return
      end


        subroutine do_test (test
     &,                map, imask 
     &,                zu, hr, dxu, dyu, phi, csu, f
     &,                zumod 
     &,                forc
     &,                acor, omega, c2dtsf
     &,                cf
     &,                h 
     &,                mxscan, mscan, tolrsf, tolrsp
     &,                iperm, jperm, iofs, nisle, nippts
     &,                crit, estimated_error
     &,                mode, variable
     &,                converged
     &,                guess, dpsi, res
     &,                solution, surfpres
     &,                npt
     &,                noslip, nisle_sp, nisle_sf
     &,                solver)
c
#include "size.h"
#include "stdunits.h"
      dimension dxu(imt), dyu(jmt), phi(jmt), f(jmt), csu(jmt)
      dimension h(imt,jmt), hr(imt,jmt)
      dimension cf(imt,jmt,-1:1,-1:1)
      dimension forc(imt,jmt), dpsi(imt,jmt)
      dimension res(imt,jmt), kmt(imt,jmt), guess(imt,jmt)
      dimension iofs(mnisle), iperm(maxipp), jperm(maxipp)
      dimension nippts(mnisle)
      dimension surfpres(imt,jmt)
      dimension solution(imt,jmt)
c
      dimension zu(imt,jmt,2)
      dimension zumod(imt,jmt,2)
c
      dimension map (imt,jmt)
      logical imask (-mnisle:mnisle)
c
      parameter (maxtest=32)
      character*16 test, prev_test
c
      logical converged, solver, noslip
      character*11 conv
      character*16 variable
      character*8 bc_symm
      data prev_test /'none'/
c
      js = 1
      je = jmt
c
        if (test .eq. 'none') then
c
          write (stdout,*)  'no individual solver tests requested'
          solver = .false.
	  noslip = .true.
c
        else if (test .eq. 'stream function') then
c
c-----------------------------------------------------------------------
c         define stream function forcing on t grid boxes
c-----------------------------------------------------------------------
c
          call sfforc (zu, dxu, dyu, csu, forc)
c
          solver  = .false.
          mode    = 1
          variable = 'dpsi'
	  bc_symm  = 't odd'
	  noslip = .true.
	  nisle = nisle_sf
	  crit = tolrsf
c
        else if (test .eq. '5 point') then
c
c-----------------------------------------------------------------------
c         initialize the coefficients for the elliptic solvers using
c         stream function 5 point numerics  (all 5 pt solvers)
c         ne, nw, se, and sw coefficients are set to zero
c         to permit testing 9 pt solvers on 5 pt coefficients.
c-----------------------------------------------------------------------
c
          call sfc5pt (acor, f, c2dtsf, dxu, dyu, csu, hr, cf)
c
          npt    = 5
          solver = .false.
          mode     = 1
          variable = 'dpsi'
	  bc_symm  = 't odd'
c
          call subset (forc, cf(1,js,0,0), js, je, nerror)
	  if (test .ne. prev_test) then
              write (stdout,*) 
     &       'forc nonzero at ',nerror,' pts with cf(i,j,0,0)=0'
	  endif
c
        else if (test .eq. '9 point') then
c
          if (mode .eq. 1) then
c
c           initialize the coefficients for the elliptic solvers using
c           stream function 9 point numerics  (all 9 pt solvers)
c
            call sfc9pt (acor, f, c2dtsf, dxu, dyu, csu, hr, cf)
c
            npt = 9
            solver = .false.
c
            call subset (forc, cf(1,js,0,0), js, je, nerror)
	    if (test .ne. prev_test) then
              write (stdout,*) 
     &       'forc nonzero at ',nerror,' pts with cf(i,j,0,0)=0'
	    endif
          else if (mode .eq. 2) then
            npt = 9
            solver = .false.
          endif
        else if (test .eq. 'surface pressure') then
c
c-----------------------------------------------------------------------
c         initialize the coefficients for the surface pressure equation
c         9 point coefficients
c-----------------------------------------------------------------------
c
          call spc9pt (dxu, dyu, csu, h, cf, js, je)
c
c         calculate "forcing" (i.e. right side) for surface pressure
c
          call spforc (zumod, dxu, dyu, csu, h, forc, js, je)
c
          npt = 9
          solver = .false.
          mode = 2
          variable = 'surfpres'
	  bc_symm  = 't even'
	  noslip = .false.
	  nisle = nisle_sp
	  crit = tolrsp
c
          call subset (forc, cf(1,js,0,0), js, je, nerror)
	  if (test .ne. prev_test) then
            write (stdout,*) 
     &     'forc nonzero at ',nerror,' pts with cf(i,j,0,0)=0'
	  endif
        else if (test .eq. '59 test') then
c
c-----------------------------------------------------------------------
c         switch to testing 9 pt solvers on (5 pt) coefficients
c-----------------------------------------------------------------------
c
          npt = 9
          solver = .false.
        else if (test .eq. 'congrad') then
c
c-----------------------------------------------------------------------
c         "conjugate gradient method"
c-----------------------------------------------------------------------
c
#ifdef timing
          call tic ('solvers','congrad')
#endif
          call congr  (npt, variable, bc_symm
     &,                      guess, solution, forc, res
     &,                      cf, js, je
     &,                      mxscan        , mscan     , crit
     &,                      imask, iperm, jperm, iofs, nisle, nippts
     &,                      converged
     &,                      estimated_error
     &                       )
#ifdef timing
          call toc ('solvers','congrad')
#endif
c
          write (stdout,*)  ' '
          write (stdout,'(a26,tr1,a11,tr1,a,i6,a,e10.3)')
     &              '=> congrad:',conv(converged) 
     &,              'scans =',mscan, ' estimated err=',estimated_error
          write (stdout,'(tr20,a,i1,a)') 'uses ', npt, ' point operator'
          solver = .true.
        endif
c
c-----------------------------------------------------------------------
c       end of test selection if block
c-----------------------------------------------------------------------
c
        if (solver) then
c
c         copy solution to appropriate variable
c
	  do j=1,jmt
	    do i=1,imt
	      if      (mode .eq. 1) then
		dpsi(i,j)     = solution(i,j)
	      else if (mode .eq. 2) then
		surfpres(i,j) = solution(i,j)
	      endif
	    enddo
	  enddo
c
c         print options selected
c
          if      (mode .eq. 1) then
	    write (stdout,'(a)') 'solving for stream function'
	  else if (mode .eq. 2) then
	    write (stdout,'(a)') 'solving for surface pressure'
	  endif
	  if (noslip) then
	    write (stdout,'(a)')'solution is no slip on island boundary'
	  else
	    write (stdout,'(a)')
     &              'solution is free slip on island boundary'
	  endif
	endif
c
	prev_test = test
c
      return
      end





      subroutine slice (dpsi, alpha_dpsi, res, islice)
c
c-----------------------------------------------------------------------
c     show a slice through "dpsi"
c-----------------------------------------------------------------------
c
#include "size.h"
#include "stdunits.h"
      dimension dpsi(imt,jmt)
      dimension res (imt,jmt)
      character*(*) alpha_dpsi
c
      write (stdout,'(/,a,/)')
     & ' => Here is a north-south slice through the solution...'
      do j=jmt,1,-1
        write (stdout,
     &    '(a,i3,a,a10,a,i3,a,i3,a,e14.7,a,i3,a,i3,a,e14.7)' )
     &    " j=", j,", ",alpha_dpsi(1:LEN_TRIM(alpha_dpsi)),
     &    "(",islice,",",j,") =",
     &    dpsi(islice,j),
     &    ", res(",islice,",",j,") =", res(islice,j) 
      enddo
      return
      end



      subroutine hellerman (ulat, taux, tauy)
c
c=======================================================================
c     this routine linearly interpolates global, zonal mean values of
c     ocean model surface boundary conditions (SST, salinity, WSX, WSY)
c     having 4.5 deg spacing, to the requested latitude.
c     
c     author:      k. dixon      e-mail=> kd@gfdl.gov
c=======================================================================
c
      parameter (nolat = 40, nolatp = nolat+1)
      parameter (c1 = 1.0, c0 = 0.0, p5=0.5)
c
      dimension olatv(nolatp)
      dimension wsxobs(nolatp), wsyobs(nolatp)
c
c     taux and tauy are in units of dynes per square centimeter
c     olatv  = latitude points for observed data
c     dolat  = latitude spacing for observed data
c     ylatt  = latitude where t,s boundary conditions are desired
c     ylatv  = latitude where windstress boundary conditions are desired
c
c     "observed" windstress data are based on global, annual mean,
c     zonally averaged values from Hellerman and Rosenstein (1981).
c     some smoothing was done.
c
c     references:
c       Hellerman, S, and M. Rosenstein, normal monthly wind stress
c     over the world ocean with error estimates, J. Phys, Oceanogr., 13,
c     1093-1104,1983.
c
      data wsxobs /  0.00,
     &               0.00,  0.00, -0.02,  0.15,  0.31,
     &               0.50,  0.82,  1.08,  1.23,  1.16,
     &               0.84,  0.41,  0.02, -0.35, -0.55,
     &              -0.67, -0.64, -0.46, -0.29, -0.19,
     &              -0.16, -0.33, -0.52, -0.59, -0.55,
     &              -0.32,  0.09,  0.42,  0.56,  0.76,
     &               0.81,  0.65,  0.29,  0.06, -0.10,
     &              -0.05, -0.03,  0.05,  0.10,  0.01/
c
      data wsyobs /  .000,
     &               .000,  .009,  .032,  .005, -.023,
     &              -.075, -.155, -.202, -.230, -.179,
     &              -.049,  .093,  .214,  .294,  .344,
     &               .383,  .364,  .269,  .189,  .178,
     &               .125, -.122, -.213, -.251, -.259,
     &              -.202, -.189, -.179, -.183, -.009,
     &               .023,  .053, -.048, -.185, -.225,
     &              -.097, -.050, -.023, -.006,  .000/
c
c---------------------------------------------------------------------
c     set latitudes of sst and salinity observations
c     and set latitudes of windstress observations
c---------------------------------------------------------------------
c
      dolat = 180.0/nolat
      rdolat = c1/dolat
      do n=1,nolat
        olatv(n) = -90.0 + (n-1.0)*dolat
      enddo
      olatv(nolatp) = -90.0 + (nolat)*dolat
c
c---------------------------------------------------------------------
c   use linear interpolation to produce the estimated surface boundary
c   condition values for wind stress components at u,v row j
c---------------------------------------------------------------------
c
      ylatv = ulat
c
      if (ylatv .le. olatv(1)) then
        nn = 1
        d = c0
      elseif (ylatv .ge. olatv(nolatp)) then
        nn = nolatp - 1
        d = dolat
      else
        do i=2,nolatp
          if (ylatv .le. olatv(i)) then
            nn = i - 1
            d  = ylatv - olatv(nn)
            goto 301
          endif
        enddo
      endif
301   continue
      taux = (wsxobs(nn)*(dolat - d) + wsxobs(nn+1)*d) *rdolat
      tauy = (wsyobs(nn)*(dolat - d) + wsyobs(nn+1)*d) *rdolat
c
      return
      end


#endif






      subroutine border (v, sym, js, je)
c
c-----------------------------------------------------------------------
c     adjust borders of an array for cyclic and symmetry settings
c
c     symmetry conditions are:
c     't odd':  "t" grid variable: asymmetric reflection at north
c     't even':  "t" grid variable: symmetric reflection at north
c     'u odd':  "u" grid variable: asymmetric reflection at north
c     'u even':  "u" grid variable: symmetric reflection at north
c-----------------------------------------------------------------------
c
      character*(*) sym
#include "size.h"
#include "stdunits.h"
c
      dimension v(imt,js:je)
c
c     set southern border
c
      if (js .eq. 1) then
        do i=2,imt-1
          v(i,1) = 0.0
        enddo
      endif
c
      if (je .eq. jmt) then
#ifdef symmetry
        if     (sym .eq. 't even') then
          do i=2,imt-1
            v(i,jmt) = v(i,jmt-1)
          enddo
        elseif (sym .eq. 't odd') then
          do i=2,imt-1
            v(i,jmt) = -v(i,jmt-1)
          enddo
        elseif (sym .eq. 'u even') then
          do i=2,imt-1
            v(i,jmt)   = v(i,jmt-2)
          enddo
        elseif (sym .eq. 'u odd') then
          do i=2,imt-1
            v(i,jmt)   = -v(i,jmt-2)
	    v(i,jmt-1) = 0.0
          enddo
        else
          write (stdout,*) '=> Error: illegal sym =',sym
	  call abort()
        endif
#else
        do i=2,imt-1
          v(i,jmt) = 0.0
        enddo
#endif
      endif
c
#ifdef cyclic
      do j=js,je
        v(1,j) = v(imt-1,j)
        v(imt,j) = v(2,j)
      enddo
#else
      do j=js,je
        v(1,j) = 0.0
        v(imt,j) = 0.0
      enddo
#endif
      return
      end



      subroutine checkerboard (solution, map, js, je)
c
c-----------------------------------------------------------------------
c     removes "checkboard" null space from an array "solution"
c
c     author              c.h. goldberg    e-mail => chg@gfdl.gov
c-----------------------------------------------------------------------
c
#include "size.h"
#include "stdunits.h"
c
      dimension solution(imt,js:je)
      dimension map(imt,js:je)
      dimension sum(0:1,0:1)
      dimension noceansum(0:1,0:1)
      dimension correction(0:1,0:1)
c
      do j1=0,1
        do i1=0,1
	  sum(i1,j1) = 0.0
	  noceansum(i1,j1) = 0
        enddo
      enddo
c
      do j1=0,1
        do i1=0,1
          do j=js+1+j1,je-1,2
            do i=2+i1,imt-1,2
              sum(i1,j1) = sum(i1,j1) + solution(i,j)
            enddo
          enddo
        enddo
      enddo
c
      do j1=0,1
        do i1=0,1
          do j=js+1+j1,je-1,2
            do i=2+i1,imt-1,2
	      if (map(i,j) .le. 0) then
		noceansum(i1,j1) = noceansum(i1,j1) + 1
              endif
            enddo
          enddo
        enddo
      enddo
c
      redsum   = sum(0,0) + sum(1,1)
      blacksum = sum(1,0) + sum(0,1)
      nred     = noceansum(0,0) + noceansum(1,1)
      nblack   = noceansum(1,0) + noceansum(0,1)
      diff = redsum/nred - blacksum/nblack
      c    = diff / 2.0
c
      write (stdout,*) ' '
      write (stdout,'(a,i6,a,i6,a,e14.7)')
     &          '=> checkerboard: nred = ',nred, ', nblack = ',nblack
     &,         ', removing a checkerboard correction of ', c
c
      correction (0,0) = -c
      correction (1,1) = -c
      correction (1,0) =  c
      correction (0,1) =  c
c
      do j1=0,1
        do i1=0,1
          do j=js+1+j1,je-1,2
            do i=2+i1,imt-1,2
              if (map(i,j) .le. 0) then
		solution(i,j) = solution(i,j) + correction(i1,j1)
              endif
            enddo
          enddo
        enddo
      enddo
c
      return
      end





      subroutine fill_land (solution, map, noslip,
     &                      nisle, iperm, jperm, iofs, nippts)
#include "size.h"
c
      dimension solution(imt,jmt)
      dimension map(imt,jmt)
      dimension iperm(maxipp)
      dimension jperm(maxipp)
      dimension nippts(mnisle)
      dimension iofs(mnisle)
      logical noslip
c
      if (noslip) then
	call fill_land1 (solution, map,
     &                   nisle, iperm, jperm, iofs, nippts)
      else
	call fill_land2 (solution, map,
     &                   nisle, iperm, jperm, iofs, nippts)
      endif

      return
      end


      subroutine fill_land1 (solution, map, 
     &                      nisle, iperm, jperm, iofs, nippts)
c
c=======================================================================
c
c     fills each land area with the [presumed constant] value
c     of solution along its ocean perimeter.
c
c     author              c.h. goldberg    e-mail => chg@gfdl.gov
c=======================================================================
c
#include "size.h"
c
      dimension solution(imt,jmt)
      dimension map(imt,jmt)
      dimension iperm(maxipp)
      dimension jperm(maxipp)
      dimension nippts(mnisle)
      dimension iofs(mnisle)
c
      do j=2,jmt-1
        do i=2,imt-1
	  if (map(i,j) .gt. 0) then
	    isle = map(i,j)
	    fill = solution(iperm(iofs(isle)+1),jperm(iofs(isle)+1))
	    solution(i,j) = fill
	  endif
	enddo
      enddo
c
      call mirror_adjust (solution)
c
      return
      end


      subroutine fill_land2 (solution, map, 
     &                      nisle, iperm, jperm, iofs, nippts)
c
c=======================================================================
c
c     fills the boundary cells of each land area with the value
c     of solution at the adjacent ocean perimeter point.
c     only NSEW directions are searched [no diagonal directions]
c     in case of multiple ocean perimeter points, their average
c     is used.
c
c     author              c.h. goldberg    e-mail => chg@gfdl.gov
c=======================================================================
c
#include "size.h"
c
      logical last_pass
      dimension solution(imt,jmt)
      dimension map(imt,jmt)
      dimension iperm(maxipp)
      dimension jperm(maxipp)
      dimension nippts(mnisle)
      dimension iofs(mnisle)
c     
      last_pass = .false.
      if (.not. last_pass) then
        last_pass = .true.

        do j=2,jmt-1
          do i=2,imt-1
            if (map(i,j) .gt. 0) then
              isle = map(i,j)
              nbrs = 0
              sum = 0.0
              if (map(i,j+1) .eq. -isle) then
                nbrs = nbrs + 1
                sum = sum + solution(i,j+1)
              endif
	      if (map(i+1,j) .eq. -isle) then
	        nbrs = nbrs + 1
	        sum = sum + solution(i+1,j)
	      endif
	      if (map(i,j-1) .eq. -isle) then
	        nbrs = nbrs + 1
	        sum = sum + solution(i,j-1)
	      endif
	      if (map(i-1,j) .eq. -isle) then
	        nbrs = nbrs + 1
	        sum = sum + solution(i-1,j)
	      endif
	      if (nbrs .gt. 0) then
	        solution(i,j) = sum / nbrs
c               last_pass = .false.
	      endif
	    endif
	  enddo
        enddo
c 
        call mirror_adjust (solution)
c
      endif 
c
      return
      end



      subroutine mirror_adjust (solution)
c
c=======================================================================
c
c     fills each border cell with the value adjacent to it
c
c     author              c.h. goldberg    e-mail => chg@gfdl.gov
c
c
c=======================================================================
c
#include "size.h"
c
      dimension solution(imt,jmt)

      call border (solution, 't odd', js, je)
      do i=1,imt
	solution(i,1)   = solution(i,2)
	solution(i,jmt) = solution(i,jmt-1)
      enddo
#ifndef cyclic
      do j=2,jmt-1
	solution(1,j)   = solution(2,j)
	solution(imt,j) = solution(imt-1,j)
      enddo
      solution(  1,  1) = solution(    2,    2)
      solution(  1,jmt) = solution(    2,jmt-1)
      solution(imt,  1) = solution(imt-1,    2)
      solution(imt,jmt) = solution(imt-1,jmt-1)
#endif
      return
      end





      subroutine zero_level (surfpres, variable, map, dxt, dyt, cst)
#include "size.h"
#include "stdunits.h"
      dimension surfpres(imt,jmt)
      dimension map(imt,jmt)
      dimension dxt(imt), dyt(jmt), cst(jmt)
      character*(*) variable
c
c     this does not correctly handle multiple basins
c
      sum = 0.0
      area_ocean = 0.0
      do j=2,jmt-1
        do i=2,imt-1
	  if (map(i,j) .le. 0) then
            area = dxt(i)*cst(j)*dyt(j)
            sum = sum + surfpres(i,j)*area
            area_ocean = area_ocean + area
	  endif
        enddo
      enddo
      surfpres0 = sum / area_ocean
      call con_adjust (surfpres, surfpres0, map)
      write (stdout,'(a,e14.7,a,a/)')
     &  '=> zero_level: removing a mean of ', surfpres0, ' from '
     &, variable
      return
      end


      subroutine ddxu (tquant, uquant, dxu, cosu)
c
c=======================================================================
c
c     Calculates x partial derivative of field tquant
c     Answer is centered at u/v points
c
c     author           c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "size.h"
      dimension tquant(imt,jmt), uquant(imt,jmt)
      dimension dxu(imt), cosu(jmt)
c
c     calculate partial derivative = ddx (tquant)
c
      call diffdxu (tquant, uquant)
      do j=1,jmt-1
        do i=1,imt-1
	  uquant(i,j) = uquant(i,j) / (dxu(i)*cosu(j))
	enddo
      enddo
c
      return
      end




      subroutine diffdxu (tquant, uquant)
c
c=======================================================================
c
c     Calculates x partial difference of field tquant
c     Answer is centered at u/v points
c
c     author           c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "size.h"
      dimension tquant(imt,jmt), uquant(imt,jmt)
      dimension cddxu(0:1,0:1), cddyu(0:1,0:1)
      dimension cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)
c
c-----------------------------------------------------------------------
c     set locally needed constants
c-----------------------------------------------------------------------
c
      c0    = 0.0
      p5    = 0.5
c
c-----------------------------------------------------------------------
c     construct coefficients for partial differences. a partial
c     difference in "x" is defined as an "x" difference of a quantity 
c     which is averaged in "y". (and symmetrically for "y" differences).
c     Note that this is an x difference and NOT an x derivative.
c     partial differences of quantities on the "t" grid are defined on
c     the "u" grid and visa versa.
c     therefore partial differences at: 
c     u/v points (i,j), involve nearby t/s points with subscripts: 
c        (i  ,j+1)    (i+1,j+1)
c        (i  ,j  )    (i+1,j  )
c     t/s points (i,j), involve nearby u/v points with subscripts:
c        (i-1,j  )    (i  ,j  )
c        (i-1,j-1)    (i  ,j-1)
c     thus if qu(i,j) is defined on u/v points, its partial
c     difference ddxqt = ddxt(qu) is defined on t/s points and has the
c     value 
c     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
c                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
c-----------------------------------------------------------------------
c
      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5
c
      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5
c
      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5
c
      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5
c
c-----------------------------------------------------------------------
c     calculate partial difference 
c        diffdx (tquant) = deltax (bary (tquant))
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do i=1,imt
	  uquant(i,j) = 0.0
	enddo
      enddo
c
      do j1=0,1
        do i1=0,1
          do j=1,jmt-1
            do i=1,imt-1
	      uquant(i,j) = uquant(i,j) + cddxu(i1,j1)*tquant(i+i1,j+j1)
            enddo
          enddo
        enddo
      enddo
c
      return
      end



      subroutine ddyu (tquant, uquant, dyu)
c
c=======================================================================
c
c     Calculates y partial derivative of field tquant
c     Answer is centered at u/v points
c
c     author           c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "size.h"
      dimension tquant(imt,jmt), uquant(imt,jmt)
      dimension dyu(jmt)
c
      call diffdyu (tquant, uquant)
      do j=1,jmt-1
        do i=1,imt-1
	  uquant(i,j) = uquant(i,j) / dyu(j)
	enddo
      enddo
c
      return
      end



      subroutine diffdyu (tquant, uquant)
c
c=======================================================================
c
c     Calculates y partial difference of field tquant
c     Answer is centered at u/v points
c
c     author           c. h. goldberg       chg@gfdl.gov
c=======================================================================
c
#include "size.h"
      dimension tquant(imt,jmt), uquant(imt,jmt)
      dimension cddxu(0:1,0:1), cddyu(0:1,0:1)
      dimension cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)
c
c-----------------------------------------------------------------------
c     set locally needed constants
c-----------------------------------------------------------------------
c
      c0    = 0.0
      p5    = 0.5
c
c-----------------------------------------------------------------------
c     construct coefficients for partial differences. a partial
c     difference in "x" is defined as an "x" difference of a quantity 
c     which is averaged in "y". (and symmetrically for "y" differences).
c     Note that this is an x difference and NOT an x derivative.
c     partial differences of quantities on the "t" grid are defined on
c     the "u" grid and visa versa.
c     therefore partial differences at: 
c     u/v points (i,j), involve nearby t/s points with subscripts: 
c        (i  ,j+1)    (i+1,j+1)
c        (i  ,j  )    (i+1,j  )
c     t/s points (i,j), involve nearby u/v points with subscripts:
c        (i-1,j  )    (i  ,j  )
c        (i-1,j-1)    (i  ,j-1)
c     thus if qu(i,j) is defined on u/v points, its partial
c     difference ddxqt = ddxt(qu) is defined on t/s points and has the
c     value 
c     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
c                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
c-----------------------------------------------------------------------
c
      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5
c
      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5
c
      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5
c
      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5
c
c-----------------------------------------------------------------------
c     calculate partial difference
c        diffdy (tquant) = deltay (barx (tquant))
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do i=1,imt
	  uquant(i,j) = 0.0
	enddo
      enddo
c
      do j1=0,1
        do i1=0,1
          do j=1,jmt-1
            do i=1,imt-1
	      uquant(i,j) = uquant(i,j) + cddyu(i1,j1)*tquant(i+i1,j+j1)
            enddo
          enddo
        enddo
      enddo
c
      return
      end




      subroutine con_adjust (dpsi, dpsi1, map)
c
c-----------------------------------------------------------------------
c     the constant dpsi1 is subtracted from dpsi(i,j) at all
c     ocean points (i.e., where map(i,j) .le. 0)
c-----------------------------------------------------------------------
c
#include "size.h"
      dimension dpsi(imt,jmt)
      dimension map(imt,jmt)
c
      do j=1,jmt
        do i=1,imt
          if (map(i,j) .le. 0) then
            dpsi(i,j) = dpsi(i,j) - dpsi1
          endif
        enddo
      enddo
      return
      end


      function conv (converged)
c  
c-----------------------------------------------------------------------
c     converts logical to character form for printing
c-----------------------------------------------------------------------
c
      character*(*) conv
      logical converged
c
      if (converged) then
        conv = '[converged]'
      else
        conv = '[diverged] '
      endif
c
      return
      end

      function relerr (x, y)  
      if (x .eq. y) then
	relerr = 0.0
      else
        xymin = min(abs(x), abs(y))
        xymax = max(abs(x), abs(y))
	if (xymin .gt. 1.0e-23 * xymax) then
	  relerr = abs(x-y)/xymin
        else
	  relerr = xymax
	endif
      endif

      return
      end

