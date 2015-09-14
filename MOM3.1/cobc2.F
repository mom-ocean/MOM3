#if defined obc_west || defined obc_east
      subroutine cobc2 (joff, js, je, i, n)
c
c=======================================================================
c     compute tracers at "tau+1" for i= 2 and/or imt-1
c     using open boundary conditions similar to Stevens (1990).
c
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       i    = boundary longitude (2 or imt-1)
c       n    = tracer component (this routine is included in the
c              "n" loop of tracer.F)
c
c
c     similar to "cobc.F" but at meridional boundaries
c
c     author:   Arne Biastoch   e-mail: abiastoch@ifm.uni-kiel.de
c=======================================================================
c
      use arrays_2d
      use memory_window
      use sbc_info
      use sbc_ocean
# include "param.h"
# include "accel.h"
# include "cobc.h"
# include "coord.h"
# include "emode.h"
# include "grdvar.h"
# include "hmixc.h"
# ifdef isoneutralmix
#  include "isopyc.h"
# endif
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# include "vmixc.h"
      dimension twodt(km)
      dimension ADV_Tx_obc(km,jmw)
# include "fdift.h"
c
# ifdef timing
      call tic ('ocean', 'cobc2')
# endif
c
c-----------------------------------------------------------------------
c
# ifdef obc_west
c
      if (i .eq. 2) then
c
c-----------------------------------------------------------------------
c       1) compute the advective velocity "uad" 
c          at the east face of the T cells
c
c       2) compute phase velocity at the western boundary: c1w
c-----------------------------------------------------------------------
c
c
        do k=1,km
          do j=js,je
            jrow = j + joff
            var = -dxu(i+1)*csu(jrow)/dtts
c	  
            uad(k,j) = (u(i,k,j,1,tau)*dyt(j)+u(i,k,j-1,1,tau)
     &                 *dyt(j-1))/(dyt(j)+dyt(j-1))*tmask(i,k,j)
            if (uad(k,j) .gt. c0) uad(k,j) = c0
c	     
	    var2(k,j,n) = t(i+2,k,j,n,taum1)-t(i+1,k,j,n,taum1)
            if (var2(k,j,n) .eq. c0) then
              c1w(k,j,n) = var
	    else
              c1w(k,j,n) = var*(t(i+1,k,j,n,tau)-t(i+1,k,j,n,taum1))
     &                          /var2(k,j,n)*tmask(i,k,j)
              if (c1w(k,j,n) .gt. c0)  c1w(k,j,n) = c0
              if (c1w(k,j,n) .lt. var) c1w(k,j,n) = var
	    endif
c	     
            if (jrow .eq. 2)     uad(k,j)=c0
            if (jrow .eq. jmtm1) uad(k,j)=c0
          enddo
        enddo
      endif
# endif
c
# ifdef obc_east
c
      if (i .eq. imtm1) then
c
c-----------------------------------------------------------------------
c       1) compute the advective velocity "uad" 
c          at the west face of the T cells
c
c       2) compute phase velocity at the western boundary: c1e
c-----------------------------------------------------------------------
c
c
        do k=1,km
          do j=js,je
            jrow = j + joff
            var = dxu(i-2)*csu(jrow)/dtts 
c	  
            uad(k,j) = (u(i-1,k,j,1,tau)*dyt(j)+u(i-1,k,j-1,1,tau)
     &                 *dyt(j-1))/(dyt(j)+dyt(j-1))*tmask(i,k,j)
            if (uad(k,j) .lt. c0) uad(k,j) = c0
c	     
	    var2(k,j,n) = t(i-1,k,j,n,taum1)-t(i-2,k,j,n,taum1)
	    if (var2(k,j,n) .eq. c0) then
              c1e(k,j,n) = var
	    else
	      c1e(k,j,n) = -var*(t(i-1,k,j,n,tau)
     &                -t(i-1,k,j,n,taum1))/var2(k,j,n)*tmask(i,k,j)
              if (c1e(k,j,n) .lt. c0)  c1e(k,j,n) = c0
              if (c1e(k,j,n) .gt. var) c1e(k,j,n) = var
	    endif
c	     
            if (jrow .eq. 2)     uad(k,j) = c0
            if (jrow .eq. jmtm1) uad(k,j) = c0
          enddo
        enddo
      endif
# endif
c
c-----------------------------------------------------------------------
c       diffusive flux across northern face of T cells
c       (use constant horizontal diffusion)
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=1,km
	    do ii=i-1,i
c
c             only del**2 diffusion 
c
              jrow = j + joff
              diff_fe(ii,k,j)  = ah*cstr(jrow)*dxur(ii)*
     &                           (t(ii+1,k,j,n,tau) - t(ii,k,j,n,tau))
            enddo
          enddo
        enddo
        do j=js,je
          do k=1,km
            diff_fn(i,k,j)  = ah*csu(jrow)*dyur(jrow)*
     &                         (t(i,k,j+1,n,tau) - t(i,k,j,n,tau))
            diff_fn(i,k,j-1)= ah*csu(jrow-1)*dyur(jrow-1)*
     &                         (t(i,k,j,n,tau) - t(i,k,j-1,n,tau))
          enddo
        enddo
c
c-----------------------------------------------------------------------
c       diffusive flux across bottom face of T cells
c       use constant vertical diffusion
c-----------------------------------------------------------------------
c
        nt2 = min(n,nsmix)
        do j=js,je
          do k=1,km-1
            diff_fb(i,k,j) = diff_cbt(i,k,j,nt2)*dz_wtr(i,k,j)*
     &                       (t(i,k,j,n,tau) - t(i,k+1,j,n,tau))
          enddo
        enddo
c
c       for vertical diffusion of tracers:
c       set k=0 elements to reflect surface tracer flux.
c       set bottom level "kz" to reflect insulation condition.
c       for vertical advection of tracers:
c       set k=0 elements of vertical advection array to the rigid
c       lid condition (w(0) = 0). also assume w(km)=0. If "kz" <
c       "km" and there is a bottom slope, w(kz) may be non zero.
c
        do j=js,je
          jrow            = j + joff
          kb              = kmt(i,jrow)
          diff_fb(i,0,j)  = stf(i,j,n)
          diff_fb(i,kb,j) = btf(i,j,n)
        enddo
c
c
c-----------------------------------------------------------------------
c       calculate lateral diffusion and advection of tracers
c-----------------------------------------------------------------------
c
# ifdef obc_west
c
        if (i .eq. 2) then
c
c         radiation condition at the western wall
c         and diffusion (del**2) of tracer
c
          do j=js,je
            do k=1,km
              jrow = j + joff
              ADV_Tx_obc(k,j) = -(c1w(k,j,n)+uad(k,j))
     &         *(t(i+1,k,j,n,tau)-t(i,k,j,n,tau))*dxur(i)*csur(jrow)
            enddo
          enddo
        endif
# endif
# ifdef obc_east
c
        if (i .eq. imtm1) then
c
c         radiation condition at the eastern wall
c         and diffusion (del**2) of tracer
c
          do j=js,je
            do k=1,km
              jrow = j + joff
              ADV_Tx_obc(k,j) = -(c1e(k,j,n)+uad(k,j))
     &         *(t(i,k,j,n,tau)-t(i-1,k,j,n,tau))*dxur(i-1)*csur(jrow)
            enddo
          enddo
        endif
# endif
c
c-----------------------------------------------------------------------
c       construct tracer source terms here
c       force n/s wall to observed values
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=1,km
            source(i,k,j) = c0
          enddo
        enddo
c
# ifndef orlanski
        call obcsponge2 (joff, js, je, i, n, t(1,1,1,1,tau),
     &                  uad, source)
# endif
c
c-----------------------------------------------------------------------
c       calculate the new tracer quantities allowing for implicit
c       treatment of vertical diffusion
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=1,km
            twodt(k) = dtts*dtxcel(k)
            t(i,k,j,n,taup1) = t(i,k,j,n,tau) + twodt(k)*
     &                        (DIFF_Tx(i,k,j) + DIFF_Tz(i,k,j) 
     &                       + DIFF_Ty(i,k,j)
     &                       + ADV_Tx_obc(k,j) 
     &                       + source(i,k,j)
     &                        )*tmask(i,k,j)
          enddo
        enddo
c
# ifdef trace_indices
        if (n.eq.1) write (stdout,'(2x,6(a,i4))')
     &    "=> In  cobc2: js=",js," je=",je," joff=",joff
     &,   " jrows=",js+joff," to ",je+joff, " at i= ",i
# endif
c
#if defined implicitvmix || defined isoneutralmix
c
c-----------------------------------------------------------------------
c       add dT/dt component due to implicit vertical diffusion
c-----------------------------------------------------------------------
c
        call ivdift (joff, js, je, i, i, n, twodt)
#endif
c
c-----------------------------------------------------------------------
c       construct diagnostics associated with tracer "n"
c-----------------------------------------------------------------------
c
        call diagt1 (joff, js, je, i, i, n, twodt)
c
c-----------------------------------------------------------------------
c     construct diagnostics (for total dT/dt)
c-----------------------------------------------------------------------
c
c     only for second call of cobc2 (n=2)
      if (n.eq.2) then
c
      idiag = 1
      call diagt2 (joff, js, je, i, i, idiag)
c
#if !defined simple_sbc && !defined smhi_sbc
c
c-----------------------------------------------------------------------
c     if needed, construct the Atmos S.B.C.(surface boundary conditions)
c     averaged over this segment
c     eg: SST and possibly SSS
c-----------------------------------------------------------------------
c
      isst  = index_of_sbc ('sst')
      isss  = index_of_sbc ('sss')
      if (isst .ne. 0 .or. isss .ne. 0) then
        call asbct (joff, js, je, i, i, isst, isss)
      endif
#endif
      endif
c
# ifdef timing
      call toc ('ocean', 'cobc2')
# endif
      return
      end

# ifndef orlanski
      subroutine obcsponge2 ( joff, js, je, i, n, tm1, uad, source )
c
c=======================================================================
c     newtonian damping variables for obc regions. damp to prescribed
c     values if there are inflow conditions. data must be prepared
c     using the "mkobc"routines included in the programs for working
c     with the MOM dataset.
c
c     input:
c       i    = (2,imt-1) open boundary row
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       n    = (1,2) = (T,S) tracer component
c       tm1  = tracer at "tau"
c       uad  = advective vel. for tracers at open walls
c            ==> restore to boundary values for inflow conditions 
c
c     output:
c       source = newtonian damping term
c
c     author:   Arne Biastoch   e-mail: abiastoch@ifm.uni-kiel.de
c=======================================================================
c
      use xobc_mod
# include "param.h"
# include "iounit.h"
# include "obc_data.h"
# include "switch.h"
# include "tmngr.h"
c
      dimension source(imt,km,2:jmw-1)
      dimension tm1(imt,km,jmw,nt)
      dimension uad(km,jmt)
c
# ifdef timing
        call tic ('cobc2', 'obc sponge2')
# endif
      if ((joff + js .eq. 2) .and. n .eq. 1) then
c      
c-----------------------------------------------------------------------
c       decide whether to read sponge data or not
c-----------------------------------------------------------------------
c
        begtim = position_within_data (model_time, obc_start_time
     &,                              obc_end_time, .true.)
        methodobc = 3
        call timeinterp (begtim, indxob2, tobc2, obc2dpm, 12, .true.
     &,           methodobc, inextdobc2, iprevdobc2, wprevobc2,
     &            readob2, inextobc2, iprevobc2)
c
c       read in the next data record from disk when needed
c
        if (readob2) then
          call getunit (ionew6, 'obc2.ud', opt_obc2)
          read (ionew6, rec=inextdobc2) obctnext, spdpmn, im, kk, jm
# ifdef obc_west	  
     &,  obcw,(((obbuf_west (j,k,nn,inextobc2),j=1,jmt),k=1,km),nn=1,2)
# endif     
# ifdef obc_east
     &,  obce,(((obbuf_east (j,k,nn,inextobc2),j=1,jmt),k=1,km),nn=1,2)
# endif     
          write (stdout,'(/a,i3,a,i2,a,i2,a,g14.7,1x,a/)')
     &    '=> read obc2 record =',inextdobc2,' into buffer =',inextobc2
     &,   ' method #',methodobc,' at day =', begtim, stamp
     
        call relunit (ionew6)
        endif
      endif
      if (n .le. 2) then
c
c-----------------------------------------------------------------------
c       construct newtonian damping term using obc data
c-----------------------------------------------------------------------
c
# ifdef obc_west
	  if (i .eq. 2 .and. obcw .ne. c0) then
	    tnext = c1-wprevobc2
	    do j=js,je
              jrow = joff + j
	      do k=1,km
                if (uad(k,j) .eq. 0.) then
	           data = tnext    *obbuf_west(jrow,k,n,inextobc2) 
     &                  + wprevobc2*obbuf_west(jrow,k,n,iprevobc2)
                   source(i,k,j) = source(i,k,j) - 
     &                               obcw*(tm1(i,k,j,n) - data)
                endif
              enddo
	    enddo
	  endif
# endif
c
# ifdef obc_east
	  if (i .eq. imtm1 .and. obce .ne. c0) then
	    tnext = c1-wprevobc2
	    do j=js,je
              jrow = joff + j
	      do k=1,km
                if (uad(k,j) .eq. 0.) then
	           data = tnext    *obbuf_east(jrow,k,n,inextobc2) 
     &                  + wprevobc2*obbuf_east(jrow,k,n,iprevobc2)
                   source(i,k,j) = source(i,k,j) - 
     &                               obce*(tm1(i,k,j,n) - data)
                endif
              enddo
	    enddo
	  endif
# endif
      endif
# ifdef timing
        call toc ('cobc2', 'obc sponge2')
# endif
      return
      end

# endif
c
#else
      subroutine cobc2 (joff, js, je, i, n)
      return
      end
#endif
