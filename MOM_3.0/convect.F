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
c
c-----------------------------------------------------------------------
c     The following convection scheme is an alternative to the
c     standard scheme. In contrast to the standard scheme,
c     it totally removes all gravitational instability in the
c     water column. It does that in one pass, so the parameter
c     ncon becomes irrelevant if this option is selected.
c     The scheme is equivalent to those used by Rahmstorf 
c     (jgr 96,6951-6963) and by Marotzke (jpo 21,903-907).
c     It is discussed in a note to Ocean Modelling (101). It uses
c     as much cpu time as 1-3 passes of the standard scheme, 
c     depending on the amount of static instability found in the
c     model, and is much faster than using "implicitvmix".
c
c     Written by Stefan Rahmstorf while at Institut fuer Meereskunde,
c     Kiel, Germany.            Comments welcome:
c                         rahmstorf@pik-potsdam.de
c
c     inputs:
c
c     kmt   = number of ocean T-cells in the vertical
c     joff  = offset between "j" in MW and "jrow" latitude on disk
c     js    = starting row in MW
c     je    = ending row in MW
c     is    = starting longitude index
c     ie    = ending longitude index
c
c     Note: istrt,iend are currently bypassed. instead, is and ie are
c           set as parameters to optimize performance
c     ts    = temperature and salinity before convection
c
c     outputs:
c
c     ts    = tracers after convection
c
c     definition of internal variables:
c
c     kcon = maximum number of levels at this location
c     lcon = counts levels down
c     lcona = upper layer of a convective part of water column
c     lconb = lower layer of a convective part of water column
c     rhoup = density referenced to same level
c     rholo = density referenced to level below
c                (note that densities are not absolute!)
c     dztsum = sum of layer thicknesses
c     trasum = sum of layer tracer values
c     tramix = mixed tracer value after convection
c     lctot = total of number of levels involved in convection
c     lcven = number of levels ventilated (convection to surface)
c     note: lctot can in rare cases count some levels twice, if they
c           get involved in two originally separate, but then
c           overlapping convection areas in the water column! It
c           is a control parameter; the sensible parameter to plot
c           is lcven. Lcven is 0 on land, 1 on ocean points with no
c           convection, and anything up to kbot on convecting points. 
c
c     author: Stefan Rahmstorf rahmstorf@pik-potsdam.de
c-----------------------------------------------------------------------
c
#ifdef partial_cell
      use memory_window
#endif
#if defined sponges && !defined sponges_old
      use xsponge_mod
#endif
#include "param.h"
      parameter (is=2, ie=imt-1)
#include "accel.h"
#include "mw.h"
#include "state.h"
      dimension kmt(imt,js+joff:je+joff), ts(imt,km,jmw,nt)
      dimension trasum(nt), rhoup(imt,km), rholo(imt,km)
      dimension lctot(imt), lcven(imt)
#include "dens.h"
c
c     check each row column by column; note that 'goto 1310'
c     finishes a particular column and moves to the next one.
c
#ifdef timing
      call tic ('tracer', 'convection: convct2')
#endif
c
      do j=js,je
c
c       find density of entire row for stability determination
c
        jrow = j + joff
# if defined sponges && !defined sponges_old
        if (spng(jrow) .ne. 0.0) cycle
# endif
        do l=1,kbot-1
          do i=is,ie
            l1=l+1
            tup = ts(i,l1,j,1) - to(l1)
            sup = ts(i,l1,j,2) - so(l1)
            tlo = ts(i, l,j,1) - to(l1)
            slo = ts(i, l,j,2) - so(l1)
            rhoup(i,l1) = dens (tup, sup, l1)
            rholo(i,l)  = dens (tlo ,slo ,l1)
          enddo
        enddo
c
        do i=is,ie
          kcon     = kmt(i,jrow)
          lctot(i) = 0
          lcven(i) = 0
          if (kcon .eq. 0) goto 1310
          lcven(i) = 1
          lcon       = 0
c
c         1. initial search for uppermost unstable pair; if none is
c            found, move on to next column
c
          do k=kcon-1,1,-1
            if (rholo(i,k) .gt. rhoup(i,k+1)) lcon = k
          enddo
c
          if (lcon .eq. 0) goto 1310
c
1319      lcona = lcon
          lconb = lcon + 1
c
c         2. mix the first two unstable layers
c
#ifdef partial_cell
          dhtxcla = dht(i,lcona,j)/dtxcel(lcona)
          dhtxclb = dht(i,lconb,j)/dtxcel(lconb)
          dztsum  = dhtxcla + dhtxclb
#else
          dhtxcla = dztxcl(lcona)
          dhtxclb = dztxcl(lconb)
          dztsum = dztxcl(lcona) + dztxcl(lconb)
#endif
          do n=1,nt
            trasum(n)        = ts(i,lcona,j,n)*dhtxcla + 
     &                         ts(i,lconb,j,n)*dhtxclb
            tramix           = trasum(n) / dztsum
            ts(i,lcona,j,n) = tramix
            ts(i,lconb,j,n) = tramix
          enddo
c
c         3. test layer below lconb
c
1306      continue
          if (lconb .eq. kcon) goto 1308
c
          l1 = lconb + 1
          rholo(i,lconb) = dens (ts(i,lconb,j,1)-to(l1)
     &,                          ts(i,lconb,j,2)-so(l1), l1)
c
          if (rholo(i,lconb) .gt. rhoup(i,l1)) then
            lconb = lconb+1
#ifdef partial_cell
            dhtxclb = dht(i,lconb,j)/dtxcel(lconb)
            dztsum  = dztsum + dhtxclb
#else
            dhtxclb = dztxcl(lconb)
            dztsum = dztsum + dztxcl(lconb)
#endif
            do n=1,nt
              trasum(n) = trasum(n) + ts(i,lconb,j,n)*dhtxclb
              tramix = trasum(n) / dztsum
              do lmix=lcona,lconb
                ts(i,lmix,j,n) = tramix
              enddo
            enddo
            goto 1306
          end if
c
c         4. test layer above lcona
c
1308      continue
          if (lcona .gt. 1) then
            l1 = lcona-1
            rholo(i,l1) = dens(ts(i,l1,j,1)-to(lcona) 
     &,                        ts(i,l1,j,2)-so(lcona),lcona)
            rhoup(i,lcona) = dens(ts(i,lcona,j,1)-to(lcona)
     &,                           ts(i,lcona,j,2)-so(lcona),lcona)
            if (rholo(i,lcona-1) .gt. rhoup(i,lcona)) then
              lcona = lcona-1
#ifdef partial_cell
              dhtxcla = dht(i,lcona,j)/dtxcel(lcona)
              dztsum  = dztsum + dhtxcla
#else
              dhtxcla = dztxcl(lcona)
              dztsum = dztsum + dhtxcla
#endif
              do n=1,nt
                trasum(n) = trasum(n) + ts(i,lcona,j,n)*dhtxcla
                tramix = trasum(n) / dztsum 
                do lmix=lcona,lconb
                  ts(i,lmix,j,n) = tramix
                enddo
              enddo
              goto 1306
            end if
          end if
c
c         5. remember the total number of levels mixed by convection
c            in this water column, as well as the ventilated column
c
          lctot(i) = lctot(i) + lconb - lcona + 1
          if (lcona .eq. 1) lcven(i) = lconb - lcona + 1
c
c         6. resume search if step 3. and 4. have been passed and this
c            unstable part of the water column has thus been removed,
c            i.e. find further unstable areas further down the column
c
          if (lconb .eq. kcon) goto 1310
          lcon = lconb
c
1302      continue
          lcon = lcon + 1
c
          if (lcon .eq. kcon) goto 1310
c
          if (rholo(i,lcon) .le. rhoup(i,lcon+1)) goto 1302
c
          goto 1319

1310      continue
        enddo
c
        do n=1,nt
          call setbcx (ts(1,1,j,n), imt, km)
        enddo
      enddo
c
#ifdef timing
      call toc ('tracer', 'convection: convct2')
#endif
      return
      end
