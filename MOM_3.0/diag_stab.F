      subroutine stabi
c
c-----------------------------------------------------------------------
c     initialization for stability testing  
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
# include "diag_stab.h"
      tdig     = 1.e-4
      cflcrt   = 1.5
      maxcfl   = 3
      cflons   = 0.0
      cflone   = 360.0
      cflats   = -90.0
      cflate   = 90.0
      cfldps   = 0.0
      cfldpe   = 6000.0e2
      return
      end

      subroutine stab (j, jrow)
c
c-----------------------------------------------------------------------
c     test for various measures of stability  
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
#ifdef stability_tests
      use arrays_2d
      use memory_window
# include "param.h"
# include "accel.h"
# include "coord.h"
# include "diag_stab.h"
# include "docnam.h"
# include "grdvar.h"
# include "hmixc.h"
# include "iounit.h"
# ifdef isoneutralmix
#  include "isopyc.h"
# endif 
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# include "state.h"
# include "tmngr.h"
# include "vmixc.h"
c
# ifdef timing
      call tic ('diagnostic', 'stability tests')
# endif
c
c-----------------------------------------------------------------------
c     perform stability tests only within specified lat,lon and depth
c     region
c-----------------------------------------------------------------------
c
      if (jrow .ge. jscfl .and. jrow .le. jecfl) then
c
c-----------------------------------------------------------------------
c     scan for CFL violations. save locations of closest approach to
c     local CFL limit.
c-----------------------------------------------------------------------
c
        cl = cflcrt * p5
        cosur = max(csur(jrow),epsln)
        do k=kscfl,kecfl
          dtmax = max(dtuv, dtts*dtxcel(k))
          f1    = dtmax*dyur(jrow)
          f2    = dtmax*dzwr(k)
	  f3    = dtmax*cosur
          do i=iscfl,iecfl
            cflu  = abs(f3*dxur(i)*u(i,k,j,1,tau))
            cflv  = abs(f1*u(i,k,j,2,tau))
            cflwu = abs(f2*adv_vbu(i,k,j))*umask(i,k,j)
            cflwt = abs(f2*adv_vbt(i,k,j))*tmask(i,k,j)
            if (cflu .ge. cl .or. cflv .ge. cl .or. cflwu .ge. cl .or.
     &        cflwt .ge. cl) then
              write (stdout,'(/,a33,i4,a1,i3,a1,i3,a13,f6.3,/)')
     &         ' ==> CFL exceeded at coordinate (i,j,k) = (',i
     &,        ',',jrow,',',k,') by factor =',cflcrt
              umax  = p5*csu(jrow)*dxu(i)/dtmax
              pcflu = abs(100.0*u(i,k,j,1,tau)/umax)
              vmax  = p5*dyu(jrow)/dtmax
              pcflv = abs(100.0*u(i,k,j,2,tau)/vmax)
              wmax  = p5*dzw(k)/dtmax
              pcflwu= abs(100.0*adv_vbu(i,k,j)/wmax)
              pcflwt= abs(100.0*adv_vbt(i,k,j)/wmax)
              write (stdout,'(a,f8.2,a,g15.8,a)')
     &         ' u reached   ', pcflu,' % of the CFL limit (',umax,')'
              write (stdout,'(a,f8.2,a,g15.8,a)')
     &         ' v reached   ', pcflv,' % of the CFL limit (',vmax,')'
              write (stdout,'(a,f8.2,a,g15.8,a)')
     &       ' adv_vbu reached', pcflwu,' % of the CFL limit (',wmax,')'
              write (stdout,'(a,f8.2,a,g15.8,a)')
     &       ' adv_vbt reached', pcflwt,' % of the CFL limit (',wmax,')'
              is = max(1,i-3)
              ie = min(imt,i+3)
              ks = max(1,k-3)
              ke = min(km,k+3)
              scl = c0
              fx  = 1.0e-2
              write (stdout,9100)'u velocity', itt
     &,        jrow, yu(jrow), xu(is), xu(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (u(1,1,j,1,tau), imt, is, ie, ks, ke, scl)
c  
              write (stdout,9100) 'v velocity', itt
     &,        jrow, yu(jrow), xu(is), xu(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (u(1,1,j,2,tau), imt, is, ie, ks, ke, scl)
c  
              write (stdout,9100)  'adv_vbu ', itt
     &,        jrow, yu(jrow), xu(is), xu(ie), fx*zw(ks), fx*zw(ke), scl
              call matrix (adv_vbu(1,1,j), imt, is, ie, ks, ke, scl)
c
              write (stdout,9100)  'adv_vbt ', itt
     &,        jrow, yt(jrow), xt(is), xt(ie), fx*zw(ks), fx*zw(ke), scl
              call matrix (adv_vbt(1,1,j), imt, is, ie, ks, ke, scl)
c
              do m=1,nt
                write (stdout,9100) trname(m), itt
     &,        jrow, yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
                call matrix (t(1,1,j,m,tau), imt, is, ie, ks, ke, scl)
              enddo
              numcfl = numcfl + 1
c
c           turn off checking on this time step when max is exceeded
c
              if (numcfl .gt. maxcfl) stabts = .false.
            endif
          enddo
        enddo
        do k=kscfl,kecfl
          dtmax = max(dtuv, dtts*dtxcel(k))
          vmax  = p5*dyu(jrow)/dtmax
          wmax  = p5*dzw(k)/dtmax
          do i=iscfl,iecfl
            umax  = p5*csu(jrow)*dxu(i)/dtmax
            if (abs(100.0*u(i,k,j,1,tau)/umax) .gt. cflup) then
              cflup = abs(100.0*u(i,k,j,1,tau)/umax)
              cflum = u(i,k,j,1,tau)
              icflu = i
              jcflu = jrow
              kcflu = k
            endif
            if (abs(100.0*u(i,k,j,2,tau)/vmax) .gt. cflvp) then
              cflvp = abs(100.0*u(i,k,j,2,tau)/vmax)
              cflvm = u(i,k,j,2,tau)
              icflv = i
              jcflv = jrow
              kcflv = k
            endif
            if (abs(100.0*umask(i,k,j)*adv_vbu(i,k,j)/wmax) .gt. cflwup)
     &        then
              cflwup = abs(100.0*adv_vbu(i,k,j)/wmax)
              cflwum = adv_vbu(i,k,j)
              icflwu = i
              jcflwu = jrow
              kcflwu = k
            endif
            if (abs(100.0*tmask(i,k,j)*adv_vbt(i,k,j)/wmax) .gt. cflwtp)
     &        then
              cflwtp = abs(100.0*adv_vbt(i,k,j)/wmax)
              cflwtm = adv_vbt(i,k,j)
              icflwt = i
              jcflwt = jrow
              kcflwt = k
            endif
          enddo
        enddo
c
c-----------------------------------------------------------------------
c     look for max peclet numbers using velocities at "tau".
c     look for max reynolds numbers using velocities at "tau".
c     (note: if the prandtl number is different than 1, only the
c            vertical diffusve coeff for Temperature is used.) 
c-----------------------------------------------------------------------
c
        do k=kscfl,min(kecfl,kbot)
          do i=iscfl,iecfl
	    ramb = c1/(visc_cbu(i,k,j) + epsln)
	    rahb = c1/(diff_cbt(i,k,j,1) + epsln)
c
#if defined velocity_horz_mix_const
            rame = c1/(visc_ceu + epsln)
            ramn = c1/(visc_cnu + epsln)
#endif
#ifdef velocity_horz_mix_var
            rame = c1/(visc_ceu(i,k,j) + epsln)
            ramn = c1/(visc_cnu(i,k,j) + epsln)
#endif
#ifdef bryan_lewis_horizontal
            rahe = c1/(ahh(k) + epsln)
            rahn = c1/(ahh(k) + epsln)
#endif
#ifdef tracer_horz_mix_const
            rahe = c1/(diff_cet + epsln)
            rahn = c1/(diff_cnt + epsln)
#endif
#ifdef tracer_horz_mix_var
            rahe = c1/(diff_cet(i,k,j) + epsln)
            rahn = c1/(diff_cnt(i,k,j) + epsln)
#endif
#ifdef isoneutralmix
            rahe = c1/(ah+ahisop*fzisop(k))
            rahn = c1/(ah+ahisop*fzisop(k))
#endif
c
            reyx = abs(u(i,k,j,1,tau)*csu(jrow)*dxu(i))*rame
            if (reyx .gt. reynx) then
              ireynx = i
              jreynx = jrow
              kreynx = k
              reynx  = reyx
              reynu  = u(i,k,j,1,tau)
              reynmu = c1/rame
            endif
            reyy = abs(u(i,k,j,2,tau)*dyu(jrow))*ramn
            if (reyy .gt. reyny) then
              ireyny = i
              jreyny = jrow
              kreyny = k
              reyny  = reyy
              reynv  = u(i,k,j,2,tau)
              reynmv = c1/ramn
            endif
	    kk = min(k+1,km)
	    if (k .ge. kmu(i,jrow)) then
	      reyz = c0
	    else 
              reyz =umask(i,kk,j)*abs(adv_vbu(i,k,j)*dzw(k))*ramb
            endif
	    if (reyz .gt. reynz) then
              ireynz = i
              jreynz = jrow
              kreynz = k
              reynz  = reyz
              reynw  = adv_vbu(i,k,j)
              reynmw = c1/ramb
            endif
            pecx = abs(u(i,k,j,1,tau)*csu(jrow)*dxu(i))*rahe
            if (pecx .gt. peclx) then
              ipeclx = i
              jpeclx = jrow
              kpeclx = k
              peclx  = pecx
              peclu  = u(i,k,j,1,tau)
              peclmu = c1/rahe
            endif
            pecy = abs(u(i,k,j,2,tau)*dyu(jrow))*rahn
            if (pecy .gt. pecly) then
              ipecly = i
              jpecly = jrow
              kpecly = k
              pecly  = pecy
              peclv  = u(i,k,j,2,tau)
              peclmv = c1/rahn
            endif
	    kk = min(k+1,km)
	    if (k .ge. kmt(i,jrow)) then
	      pecz = 0.0
	    else 
              pecz =tmask(i,kk,j)*abs(adv_vbt(i,k,j)*dzw(k))*rahb
            endif
	    if (pecz .gt. peclz) then
              ipeclz = i
              jpeclz = jrow
              kpeclz = k
              peclz  = pecz
              peclw  = adv_vbt(i,k,j)
              peclmw = c1/rahb
            endif
          enddo
        enddo
c
c-----------------------------------------------------------------------
c     look for ficticious creation of local extremum for tracers
c     by finding local min and max tracer at "tau" and comparing to
c     tracer at "tau+1"
c-----------------------------------------------------------------------
c
        call getunit (iostab, 'iostab'
     &,             'formatted sequential append')
        ks = max(2,kscfl)
        ke = min(km-1,kecfl)
        is = max(2,iscfl)
        ie = min(imt-1,iecfl)
        do n=1,nt
          do k=ks,ke
            do i=is,ie
	      if (tmask(i,k,j) .ne. c0) then
	        tbig = max(t(i,k,j,n,tau),t(i,k,j,n,taum1))
	        tsml = min(t(i,k,j,n,tau),t(i,k,j,n,taum1))
	        do jj=j-1,j+1,2
	          if (tmask(i,k,jj) .ne. c0) then
	  	    tbig = max(tbig,t(i,k,jj,n,tau),t(i,k,jj,n,taum1))
	  	    tsml = min(tsml,t(i,k,jj,n,tau),t(i,k,jj,n,taum1))
	  	  endif
	        enddo
c
	        do ii=i-1,i+1,2
	          if (tmask(ii,k,j) .ne. c0) then
		    tbig = max(tbig,t(ii,k,j,n,tau),t(ii,k,j,n,taum1))
		    tsml = min(tsml,t(ii,k,j,n,tau),t(ii,k,j,n,taum1))
		  endif
	        enddo
c
	        do kk=k-1,k+1,2
	          if (tmask(i,kk,j) .ne. c0) then
		    tbig = max(tbig,t(i,kk,j,n,tau),t(i,kk,j,n,taum1))
		    tsml = min(tsml,t(i,kk,j,n,tau),t(i,kk,j,n,taum1))
		  endif
	        enddo
c
                tcrit = tdig*abs(t(i,k,j,n,taup1))
                if (tmask(i,k,j) .ne. c0 .and. 
     &            ((t(i,k,j,n,taup1) .gt. tbig + tcrit)
     &	          .or. (t(i,k,j,n,taup1) .lt. tsml - tcrit))) then
		  write (iostab,'(i4, i4, i4, i2, 3g14.7)')
     &             i, k, jrow, n, t(i,k,j,n,taup1), tsml, tbig
                endif
	      endif
	    enddo
	  enddo
        enddo
        call relunit (iostab)
c      
c-----------------------------------------------------------------------
c     look for temperatures and salinities outside ranges used for
c     specifying density coefficients
c-----------------------------------------------------------------------
c
        call getunit (iobadt, 'iobadt'
     &,             'formatted sequential append')
        call getunit (iobads, 'iobads'
     &,             'formatted sequential append')
        ks = max(1,kscfl)
        ke = min(km,kecfl)
        is = max(2,iscfl)
        ie = min(imt-1,iecfl)
        do k=1,km
          do i=is,ie
	    if (tmask(i,k,j) .ne. c0) then
              tlocal = t(i,k,j,1,taup1)
# if defined salinity_psu 
              slocal = t(i,k,j,2,taup1)
# else
              slocal = t(i,k,j,2,taup1)*1000.0 + 35.0
# endif
	      n = 1
              if (tlocal .lt. tmink(k) .or. tlocal .gt. tmaxk(k)) then
	  	write (iobadt,'(i4, i4, i4, i2, 3g14.7)')
     &          i, k, jrow, n, tlocal, tmink(k), tmaxk(k)
              endif
	      n = 2
              if (slocal .lt. smink(k) .or. slocal .gt. smaxk(k)) then
		write (iobads,'(i4, i4, i4, i2, 3g14.7)')
     &             i, k, jrow, n, slocal, smink(k), smaxk(k)
              endif
	    endif
	  enddo
        enddo
        call relunit (iobadt)
        call relunit (iobads)
c      
# ifdef timing
        call toc ('diagnostic', 'stability tests')
# endif
9100    format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     &  ,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     &  ,', scaling=',1pg10.3)

      endif

#endif
      return
      end

