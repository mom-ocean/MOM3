#ifdef test_convect
# include "util.F"
# include "grids.F"
# include "denscoef.F"
# include "size_check.F"
# include "state.F"
#ifdef timing
# include "timer.F"
#endif

      program driver
c
c=======================================================================
c
c                        CONVECTION MODULE
c
c     To test various convection schemes in a simple one dimensional 
c     model
c     
c        1) compile and run this module using the "run_convect" script
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use memory_window
# include "param.h"
# include "accel.h"
# include "coord.h"
# include "grdvar.h"
# include "mw.h"
# include "scalar.h"
# include "state.h"
# include "switch.h"
# include "tmngr.h"
      dimension kmt(imt,jmt), tt(km,2)
c
      stdout = 6
c
      write (stdout,'(//,25x,a/)')
     &' T E S T I N G   A   1-D   M O D E L   O F   C O N V E C T I O N'
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
c     prescribe some initial stratification
c-----------------------------------------------------------------------
c
      z0 = 30.0e2
      hh = 80.0e2
      zm = zt(km)
      t0 = 7.5
      t1 = 10.0
      write (stdout,'(/,10x, a/)') 'Initial conditions: two bubbles'
      do k=1,km
        if (k .eq. km/2) then
          tt(k,1) = tt(1,1)
          tt(k,2) = 1.1*tt(1,2)
        else if (k .eq. km) then
          tt(k,1) = tt(4,1)
          tt(k,2) = 0.9*tt(4,2)
        else
          tt(k,1) = t0*(1.0 - tanh((zt(k)-hh)/z0)) + t1*(1.0-zt(k)/zm)
          tt(k,2) = 0.0349 - 0.035
        endif
         write (stdout,*) 'k=',k,' zt(k)=',zt(k), ' temp=',tt(k,1)
     &, ' salt=',tt(k,2)
      enddo
c
c     set I.C. for t,s
c
      taum1 = -1
      taup1 = 0
      do j=1,jmw
        do i=1,imt
          do k=1,km
            do n=1,2
              t(i,k,j,n,taum1) = tt(k,n)
              t(i,k,j,n,taup1) = tt(k,n)
            enddo
          enddo
        enddo
      enddo
      do j=1,jmt
        do i=1,imt
          kmt(i,j) = km
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     integrate equations for all k at one point (i,j)
c     "j" is the row in the MW
c-----------------------------------------------------------------------
c
      i         = imt/2
      j         = 2
      is        = i
      ie        = i
      js        = j
      je        = j
      joff      = 0
      ncon      = 3
c
      write (stdout,'(/a,i2)')
     &     '  CASE 1: standard convection with ncon =', ncon
      write (stdout,'(/a)') '  CASE 2: full convection'
c
      do itt=1,8
c
c       do old style convection
c
        call convct (t(1,1,1,1,taum1), ncon, joff, js, je, is, ie, kmt)
c
c       do full convection
c
        call convct2 (t(1,1,1,1,taup1), joff, js, je, is, ie, kmt)
c
c       show some results
c
        write (stdout,'(/32x,a,i6,/14x,a,i2,a,27x,a/)')
     & ' After time step =',itt, 'CASE 1 (ncon=',ncon,')', 'CASE 2'
        do k=1,km
          write (stdout,'(a,i2,3x, 2(a,g12.5),10x,2(a,g12.5))') 
     &    'k=',k,' t =',t(i,k,j,1,taum1),' s =',t(i,k,j,2,taum1)
     &,   ' t =',t(i,k,j,1,taup1), ' s =',t(i,k,j,2,taup1)
        enddo
      enddo
      stop
      end
#endif



      subroutine convct (ts, ncon, joff, js, je, istrt, iend, kmt)
#if !defined implicitvmix || defined isoneutralmix
# ifdef partial_cell
      use memory_window
# endif
# if defined sponges && !defined sponges_old
      use xsponge_mod
# endif
# include "param.h"
      parameter (is=2, ie=imt-1)
# include "accel.h"
# ifdef partial_cell
#  include "mw.h"
# endif
c
c-----------------------------------------------------------------------
c     standard explicit convection scheme
c     convectively adjust water column if gravitationally unstable
c
c     inputs:
c
c     ncon  = number of passes through convection routine
c     joff  = offset between "j" in MW and "jrow" latitude on disk
c     js    = starting row in MW
c     je    = ending row in MW
c     is    = starting longitude index
c     ie    = ending longitude index
c
c     Note: istrt,iend are currently bypassed. instead, is and ie are
c           set as parameters to optimize performance
c     kmt   = number of ocean T-cells in the vertical
c     ts    = temperature and salinity before convection
c
c     outputs:
c
c     ts    = tracers after convection
c
c-----------------------------------------------------------------------
c
      dimension ts(imt,km,jmw,nt), temp(imt,km,jmw)
      dimension kmt(imt,js+joff:je+joff)
c
# ifdef timing
      call tic ('tracer', 'convection: convct')
# endif
c
c     ks=1: compare lev. 1 to 2; 3 to 4; etc.
c     ks=2: compare lev. 2 to 3; 4 to 5; etc.
c
      do nn=1,ncon
        do ks=1,2
c
c         find density for rows
c
          call statec (ts(1,1,1,1), ts(1,1,1,2), temp(1,1,1)
     &,                max(js,2), je, is, ie, ks)
c
c         set "heavy water" in land to stop convection
c
          dense = 1.e15
          do j=js,je
            jrow = j + joff
# if defined sponges && !defined sponges_old
            if (spng(jrow) .ne. 0.0) cycle
# endif
            do i=is,ie
              k = kmt(i,jrow) + 1
              if (k .le. km) then
                temp(i,k,j) = dense
              endif
            enddo
          enddo
c
c         if unstable,  mix tracers on adjoining levels
c
          do n=1,nt
            do j=js,je
# if defined sponges && !defined sponges_old
              jrow = j + joff
              if (spng(jrow) .ne. 0.0) cycle
# endif
              do k=ks,kmm1,2
                do i=is,ie
                  if (temp(i,k,j) .gt. temp(i,k+1,j)) then
# ifdef partial_cell
                    dhtxcl = dht(i,k,j)/dtxcel(k)
                    dhtxcl2 = dht(i,k+1,j)/dtxcel(k+1)
                    ts(i,k,j,n)   = (dhtxcl*ts(i,k,j,n) +
     &                               dhtxcl2*ts(i,k+1,j,n))
     &                             /(dhtxcl + dhtxcl2)
                    ts(i,k+1,j,n) = ts(i,k,j,n)
# else
                    ts(i,k,j,n)   = (dztxcl(k)*ts(i,k,j,n) +
     &                          dztxcl(k+1)*ts(i,k+1,j,n))*dzwxcl(k)
                    ts(i,k+1,j,n) = ts(i,k,j,n)
# endif
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c
      do n=1,nt
        do j=js,je
          call setbcx (ts(1,1,j,n), imt, km)
        enddo
      enddo
c
# ifdef timing
      call toc ('tracer', 'convection: convct')
# endif
#endif
      return
      end


      subroutine convct2 (ts, joff, js, je, istrt, iend, kmt)

!=======================================================================
!     The following convection scheme is an alternative to the standard 
!     scheme. In contrast to the standard scheme, it totally removes 
!     all gravitational instability in the water column. It does that 
!     in one pass, so the parameter ncon becomes irrelevant if this 
!     option is selected. Since most convection propagates downward the
!     scheme looks downward first and follows any instability (upward or
!     downward) before checking the other direction. The routine mixes 
!     passive tracers only after the entire instability is found. The 
!     scheme is similar to that used by Rahmstorf (jgr 96,6951-6963) and 
!     by Marotzke (jpo 21,903-907). It is discussed in a note to Ocean 
!     Modelling (101). It uses as much cpu time as 1-3 passes of the 
!     standard scheme, depending on the amount of static instability 
!     found in the model, and is much faster than using "implicitvmix".
!
!     based on code by s.rahmstorf   e-mail: rahmstorf@pik-potsdam.de
!     Recoded by m.eby in June, 2000 e-mail: eby@uvic.ca
!     To recover the exact same numerical values as the original
!     Rahmstorf code, look for the two "Rahmstorf" comments below.      
!      
!     inputs:
!
!     joff   = offset between "j" in MW and "jrow" latitude on disk
!     js     = starting row in MW
!     je     = ending row in MW
!     is     = starting longitude index
!     ie     = ending longitude index
!     kmt    = number of ocean "t" boxes in the vertical
!
!     Note: istrt,iend are currently bypassed. instead, is and ie are
!           set as parameters to optimize performance
!     ts     = tracers before convection
!
!     outputs:
!
!     ts     = tracers after convection
!
!     other previously undefined variables:
!
!     chk_la = logical flag to check level above kt
!     chk_lb = logical flag to check level below kb
!     kb     = bottom level of (potential) instability
!     kbo    = bottom level of ocean
!     kt     = top level of (potential) instability
!     ktot   = total number of levels convecting in the column
!     kven   = number of levels that ventilated in the column
!     la     = test level above kt
!     lb     = test level below kb
!     rl     = lower level density referenced to lower level
!     ru     = upper level density referenced to lower level
!     tmx    = mixed tracer (1=temp, 2=salt, 3=other)
!     tsm    = sum of tracers (weighted by thickness) in the instability
!     zsm    = total thickness of the instability
!
!=======================================================================

#ifdef partial_cell
      use memory_window
#endif
#if defined sponges && !defined sponges_old
      use xsponge_mod
#endif
#include "param.h"
      integer is, ie
      parameter (is=2, ie=imtm1)
#include "accel.h"
#ifdef partial_cell
# include "mw.h"
#endif
#include "state.h"

      integer i, iend, istrt, j, je, joff, jrow, js, k, kb, kbo
      integer kmt(imt,js+joff:je+joff), kt, ktot, kven, l, la, lb, n
      logical chk_la, chk_lb
      real    rl, ru, tmx(3), ts(imt,km,jmw,nt), tsm(3), zsm

#include "dens.h"
#ifdef timing

      call tic ('tracer', 'convection: convct2')
#endif

      do j=js,je
        jrow = j + joff
#if defined sponges && !defined sponges_old
        if (spng(jrow) .ne. c0) cycle
#endif
        do i=is,ie
!         ktot = 0
          kbo = kmt(i,jrow)

!         search for unstable regions starting from the top
          kt = 1
          kb = 2
          do while (kt .lt. kbo)
            ru = dens (ts(i,kt,j,1)-to(kb),ts(i,kt,j,2)-so(kb),kb)
            rl = dens (ts(i,kb,j,1)-to(kb),ts(i,kb,j,2)-so(kb),kb)

!           sum the first pair found in an unstable region
            if (ru .gt. rl) then
              chk_la = .true.
              chk_lb = .true.
#ifdef partial_cell
              dztxcl(kt) = dht(i,kt,j)/dtxcel(kt)
              dztxcl(kb) = dht(i,kb,j)/dtxcel(kb)
#endif
              zsm = dztxcl(kt) + dztxcl(kb)
              tsm(1) = ts(i,kt,j,1)*dztxcl(kt) + ts(i,kb,j,1)*dztxcl(kb)
              tmx(1) = tsm(1)/zsm
              tsm(2) = ts(i,kt,j,2)*dztxcl(kt) + ts(i,kb,j,2)*dztxcl(kb)
              tmx(2) = tsm(2)/zsm

              do while (chk_lb .or. chk_la)

!               check for an unstable level (lb) below kb
                if (kb .ge. kbo) chk_lb = .false.
                do while (chk_lb)
                  chk_lb = .false.
                  lb = kb + 1 
                  ru = dens (tmx(1)-to(lb),      tmx(2)-so(lb),      lb)
                  rl = dens (ts(i,lb,j,1)-to(lb),ts(i,lb,j,2)-so(lb),lb)
                  if (ru .gt. rl) then
!                   add new level to sums
                    kb = lb
#ifdef partial_cell
                    dztxcl(kb) = dht(i,kb,j)/dtxcel(kb)
#endif
                    zsm = zsm + dztxcl(kb)
                    tsm(1) = tsm(1) + ts(i,kb,j,1)*dztxcl(kb)
                    tmx(1) = tsm(1)/zsm
                    tsm(2) = tsm(2) + ts(i,kb,j,2)*dztxcl(kb)
                    tmx(2) = tsm(2)/zsm
                    chk_la = .true.
                    if (kb .lt. kbo) chk_lb = .true.
                  endif
                enddo

!               check for an unstable level (la) above kt
! to get equivalent of original Rahmstorf code, uncomment the next line
!                chk_la = .true.
                if (kt .le. 1) chk_la = .false.
                do while (chk_la)
                  chk_la = .false.
                  la = kt - 1
                  ru = dens (ts(i,la,j,1)-to(kt),ts(i,la,j,2)-so(kt),kt)
                  rl = dens (tmx(1)-to(kt),      tmx(2)-so(kt),      kt)
                  if (ru .gt. rl) then
!                   add new level to sums
                    kt = la
#ifdef partial_cell
                    dztxcl(kt) = dht(i,kt,j)/dtxcel(kt)
#endif
                    zsm = zsm + dztxcl(kt)
                    tsm(1) = tsm(1) + ts(i,kt,j,1)*dztxcl(kt)
                    tmx(1) = tsm(1)/zsm
                    tsm(2) = tsm(2) + ts(i,kt,j,2)*dztxcl(kt)
                    tmx(2) = tsm(2)/zsm
                    chk_lb = .true.
! to get equivalent of original Rahmstorf code, comment out the next line
                    if (kt .gt. 1) chk_la = .true.
                  endif
                enddo
              enddo

!             mix all tracers from kt to kb
              do k=kt,kb
                ts(i,k,j,1) = tmx(1)
                ts(i,k,j,2) = tmx(2)
              enddo
              do n=3,nt
                tsm(3) = c0
                do k=kt,kb
                  tsm(3) = tsm(3) + ts(i,k,j,n)*dztxcl(k)
                enddo
                tmx(3) = tsm(3)/zsm 
                do k=kt,kb
                  ts(i,k,j,n) = tmx(3)
                enddo
              enddo

!             some possible diagnostics
!              ktot = ktot + kb - kt + 1
!              if (kt .eq. 1) kven = kb

              kt = kb + 1
            else
              kt = kb
            endif

!           continue the search for other unstable regions
            kb = kt + 1
          enddo
        enddo

        do n=1,nt
          call setbcx (ts(1,1,j,n), imt, km)
	enddo
      enddo
#ifdef timing

      call toc ('tracer', 'convection: convct2')
#endif

      return
      end
