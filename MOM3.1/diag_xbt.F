      subroutine xbti
c
c=======================================================================
c     initialize all  XBT positions and arrays
c
c     inputs:
c
c     xbtlon = longitude of XBT
c     xbtlat = latitude of XBT
c     xbtdpt = depth of XBT
c
c     outputs:
c
c     numxbt = number of XBTs
c     xname  = names of quantities measured by the XBTs
c     nxbtts = time step counter for averaging XBTs
c     ixbt   = nearest model grid point to xbtlon
c     jxbt   = nearest model grid point to xbtlat
c     kxbt   = nearest model grid point to xbtdpt
c     nsxbt  = starting XBT number for jrow
c     nexbt  = ending XBT number for jrow
c     axbt   = space for accumulating XBT data for averaging
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c=======================================================================
c
#ifdef save_xbts
      logical errors
# include "param.h"
# include "coord.h"
# include "iounit.h"
# include "diag_xbt.h"
      dimension numb(maxxbt), numrow(jmt)
c
c-----------------------------------------------------------------------
c     initialize positions of XBT stations. there may be up to
c     "maxxbt" stations specified.
c
c      example: deploy 15 XBT`s at latitude 10S every 10 degs in the
c               pacific starting at 135E and drop them down to 4000m.
c               then deploy 15 more at the equator.
c               
c     10S pacific deployment
c
c      n = 0
c      do  i=1,15
c        n = n + 1
c        xbtlon(n) = 135.0 + (i-1)*10.0
c	 xbtlat(n) = -10.0
c	 xbtdpt(n) = 4000.0e2
c      enddo
c
c     equatorial pacific deployment
c
c      do  i=1,15
c        n = n + 1
c        xbtlon(n) = 135.0 + (i-1)*10.0
c	 xbtlat(n) = 0.0
c	 xbtdpt(n) = 4000.0e2
c      enddo
c-----------------------------------------------------------------------
c
# ifdef timing
      call tic ('diagnostic', 'xbt')
# endif
c
      write (stdout,8900)
c
c     initialize all XBTs to "not set" condition
c
      do n=1,maxxbt
	xbtlon(n) = epsln
	xbtlat(n) = epsln
	xbtdpt(n) = epsln
      enddo
c
c-----------------------------------------------------------------------
c     begin user deployment of XBTs
c-----------------------------------------------------------------------
c
c
c     test deployment of 2 XBTs near equator to depth of near 40m
c
      n = 0
      do i=1,2
        n = n + 1
        xbtlon(n) = 180.0 + (i-1)*10.0
	xbtlat(n) = 0.0 - 3*i
	xbtdpt(n) = 40.0e2
      enddo
c
c-----------------------------------------------------------------------
c     end user deployment of XBTs
c-----------------------------------------------------------------------
c
c
c     verify that the number of XBTs doesn`t exceed the "maxxbt" limit
c
      if (n .gt. maxxbt) then
        write (stdout,*) ' '
        write (stdout,*) ' Error: ==> number of XBT stations (',n
     &,                  ') exceeds parameter "maxxbt"'
	call abort()
      else if (n .lt. maxxbt) then
        write (stdout,*) ' '
        write (stdout,*) ' Warning: ==> save space by reducing'
     &,                  ' parameter "maxxbt" to ',n,' in "diag_xbt.h"'
      endif
      if (n .eq. 0) then
        write (stdout,*) ' '
        write (stdout,*) ' Error: ==> no XBT locations were specified.'
     &,                  ' set them in routine "xbti"'
	call abort()
      endif
c
c-----------------------------------------------------------------------
c     specify names of quantities measured by the XBT
c-----------------------------------------------------------------------
c
      xnameu(1,1) = ' dU/dt      '
      xnameu(1,2) = ' dV/dt      '
      xnameu(2,1) = ' -(P)x      '
      xnameu(2,2) = ' -(P)y      '
      xnameu(3,1) = ' -(UU)x     '
      xnameu(3,2) = ' -(UV)x     '
      xnameu(4,1) = ' -(VU)y     '
      xnameu(4,2) = ' -(VV)y     '
      xnameu(5,1) = ' -(WU)z     '
      xnameu(5,2) = ' -(WV)z     '
      xnameu(6,1) = ' DIFF_Ux    '
      xnameu(6,2) = ' DIFF_Vx    '
      xnameu(7,1) = ' DIFF_Uy    '
      xnameu(7,2) = ' DIFF_Vy    '
      xnameu(8,1) = ' DIFF_Uz    '
      xnameu(8,2) = ' DIFF_Vz    '
      xnameu(9,1) = ' DIFF_Umet  '
      xnameu(9,2) = ' DIFF_Vmet'
      xnameu(10,1) = ' fV         '
      xnameu(10,2) = ' -fU        '
      xnameu(11,1) = '   source   '
      xnameu(11,2) = '   source   '
      xnameu(12,1) = ' -(surf P)x '
      xnameu(12,2) = ' -(surf P)y '
      xnameu(13,1) = 'ADV_Umetric '
      xnameu(13,2) = '-ADV_Vmetric'
      xnameu(14,1) = ' -U(U)x     '
      xnameu(14,2) = ' -U(V)x     '
      xnameu(15,1) = ' -V(U)y     '
      xnameu(15,2) = ' -V(V)y     '
      xnameu(16,1) = ' -W(U)z     '
      xnameu(16,2) = ' -W(V)z     '
      xnameu(17,1) = ' average U  '
      xnameu(17,2) = ' average V  '
c
      xnamex(1) = ' Surf T flux'
      xnamex(2) = ' Taux       '
      xnamex(3) = ' Tauy       '
      xnamex(4) = ' average W  '
c
      xnamet(1) = ' dT/dt      '
      xnamet(2) = ' -(UT)x     '
      xnamet(3) = ' -(VT)y     '
      xnamet(4) = ' -(WT)z     '
      xnamet(5) = ' DIFF_Tx    '
      xnamet(6) = ' DIFF_Ty    '
      xnamet(7) = ' DIFF_Tz    '
      xnamet(8) = ' T source   '
      xnamet(9) = ' T convect  '
      xnamet(10) = ' T filter   '
      xnamet(11) = ' -U(T)x     '
      xnamet(12) = ' -V(T)y     '
      xnamet(13) = ' -W(T)z     '
      xnamet(14) = ' chg var  T '
      xnamet(15) = ' average T  '
c
c-----------------------------------------------------------------------
c     initialize counter for the number of time steps per average
c-----------------------------------------------------------------------
c
      nxbtts = 0
c
c-----------------------------------------------------------------------
c     convert XBT positions to nearest model grid points and set
c     unused XBTs to arbitrary value. don`t allow more than "kmxbt"
c     levels to be sampled at any position.
c-----------------------------------------------------------------------
c
      errors = .false.
      maxkm = 1
      do n=1,maxxbt
        if (xbtlon(n) .ne. epsln) then
	  numxbt  = n
	  ixbt(n) = indp(xbtlon(n), xt(2), imt-2) + 1
	  jxbt(n) = indp(xbtlat(n), yt(2), jmt-2) + 1
	  kxbt(n) = indp(xbtdpt(n), zt, km)
	  maxkm   = max(kxbt(n),maxkm)
	  numb(n) = (jxbt(n)-1)*imt + ixbt(n)
	  if (kxbt(n) .gt. kmxbt) then
	    errors = .true.
	    write (stdout,'(/a,i4,a,f10.2,a/a,i3,a,i3,a,i3/)') 
     &      '=> Error: XBT station #',n, ' is set for a depth of '
     &,                xbtdpt(n),' cm in xbti.F'
     &,     '          The nearest model level is k=',kxbt(n)
     &,                '. Increase "kmxbt" in diag_xbt.h from ', kmxbt 
     &,                ' to ',kxbt(n)
	  endif
        else
	  ixbt(n) = 0
	  jxbt(n) = 0
	  kxbt(n) = 0
	  numb(n) = 0
	endif
      enddo
      if (maxkm .lt. kmxbt) then
        write (stdout,*) ' '
        write (stdout,*) ' Warning: ==> save space by reducing'
     &,                  ' parameter "kmxbt" to ',maxkm,' "diag_xbt.h"'
      endif
      if (errors) then
        call abort()
      endif
c
c-----------------------------------------------------------------------
c     sort the XBTs from south to north and west to east
c-----------------------------------------------------------------------
c
      do nsort=1,1000
        isort = 0
	do n=2,numxbt
	  if (numb(n) .lt. numb(n-1)) then
	    isort = 1
	    call iswapx (ixbt(n), ixbt(n-1))
	    call iswapx (jxbt(n), jxbt(n-1))
	    call iswapx (kxbt(n), kxbt(n-1))
	    call iswapx (numb(n), numb(n-1))
	  endif
	enddo
	if (isort .eq. 0) go to 10
      enddo
10    continue
      do n=1,numxbt
	write (stdout,9000) n, yt(jxbt(n)), xt(ixbt(n))
     &,                       kxbt(n), zt(kxbt(n))*0.01
      enddo
c
c-----------------------------------------------------------------------
c     count the number of XBTs on each model latitude. 
c     nsxbt is the starting XBT number
c     nexbt is the ending XBT number 
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        numrow(jrow) = 0
      enddo
      do n=1,numxbt
        numrow(jxbt(n)) = numrow(jxbt(n)) + 1
      enddo
      n = 0
      do jrow=1,jmt
        if (numrow(jrow) .ne. 0) then
	  nsxbt(jrow) = n + 1
	  nexbt(jrow) = nsxbt(jrow) + numrow(jrow) - 1
	  n           = n + numrow(jrow)
	else
	  nsxbt(jrow) = 0
	  nexbt(jrow) = 0
	endif
      enddo
c
c-----------------------------------------------------------------------
c     initialize all time accumulators to zero
c-----------------------------------------------------------------------
c
      do num=1,numxbt
	do k=1,kmxbt
          do n=1,nt
	    do m=1,ntxbt
	      txbt(k,m,n,num) = c0
            enddo
          enddo
          do n=1,2
	    do m=1,nuxbt
	      uxbt(k,m,n,num) = c0
            enddo
          enddo
	  xbtw(k,num) = c0
        enddo
	uxbtsf(1,num) = c0
	uxbtsf(2,num) = c0
	do n=1,nt
	  txbtsf(n,num) = c0
	enddo
      enddo
c
# ifdef timing
      call toc ('diagnostic', 'xbt')
# endif
c
      return
8900  format (//,20x,' X B T    S T A T I O N    L O C A T I O N S'/)
9000  format (1x, ' XBT station #',i4, ' is at lat =',f6.2,', lon ='
     &,      f6.2, ', for ',i3,' levels down to a depth of ',f6.1,' m')
      end
      subroutine iswapx (i, j)
      itemp = i
      i     = j
      j     = itemp
      return
      end


      subroutine txbt1 (joff, js, je, n)
c
c=======================================================================
c     input:
c       joff  = offset relating "j" in the MW to latitude "jrow"
c       js    = starting row in the MW
c       je    = ending row in the MW
c       n     = tracer component
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use memory_window
# include "param.h"
# include "accel.h"
# include "grdvar.h"
# include "hmixc.h"
# ifdef isoneutralmix
#  include "isopyc.h"
# endif
# include "mw.h"
# include "scalar.h"
# include "vmixc.h"
# include "diag_xbt.h"
# include "fdift.h"
c
# ifdef timing
      call tic ('diagnostic', 'xbt')
# endif
c
      if (n .gt. nt) then
        write (stdout,*) '=> Error: n=',n,' in txbt1.F'
        call abort()
      endif
c
      do j=js,je
        jrow = j + joff
        if (nsxbt(jrow) .ne. 0) then
c
c-----------------------------------------------------------------------
c         accumulate data for the nth XBT for one time step
c-----------------------------------------------------------------------
c
          do nth=nsxbt(jrow),nexbt(jrow)
            i   = ixbt(nth)
	    lev = kxbt(nth)
	    do k=1,lev
c
c-----------------------------------------------------------------------
c             tracer
c-----------------------------------------------------------------------
c
              term = tmask(i,k,j)*t(i,k,j,n,tau)
	      txbt(k,15,n,nth) = txbt(k,15,n,nth) + term 
c
c-----------------------------------------------------------------------
c             d(tracer)/dt
c-----------------------------------------------------------------------
c
	      r2dt = c1/(c2dtts*dtxcel(k))
              term = tmask(i,k,j)*(t(i,k,j,n,taup1) -
     &                             t(i,k,j,n,taum1))*r2dt
	      txbt(k,9,n,nth) = txbt(k,9,n,nth) + term 
c
c-----------------------------------------------------------------------
c             zonal advection (flux form) of tracer
c-----------------------------------------------------------------------
c
              term = -tmask(i,k,j)*ADV_Tx(i,k,j)
	      txbt(k,2,n,nth) = txbt(k,2,n,nth) + term 
c
c-----------------------------------------------------------------------
c             pure zonal advection of tracer
c-----------------------------------------------------------------------
c
c             - U(T)x = T(U)x - (UT)x
c
              dudx = (adv_vet(i,k,j)-adv_vet(i-1,k,j))*dxtr(i)
     &               *cstr(jrow)
# if defined gm_advect
     &              +(adv_vetiso(i,k,j)-adv_vetiso(i-1,k,j))*dxtr(i)
     &               *cstr(jrow)
# endif
              term = tmask(i,k,j)*(t(i,k,j,n,tau)*dudx - ADV_Tx(i,k,j))
	      txbt(k,11,n,nth) = txbt(k,11,n,nth) + term 
c
c-----------------------------------------------------------------------
c             meridional advection (flux form) of tracer
c-----------------------------------------------------------------------
c
              term = -tmask(i,k,j)*ADV_Ty(i,k,j)
	      txbt(k,3,n,nth) = txbt(k,3,n,nth) + term 
c
c-----------------------------------------------------------------------
c             pure meridional advection of tracer
c-----------------------------------------------------------------------
c
c             - V(T)y = T(V)y - (VT)y
c
              dvdy = (adv_vnt(i,k,j)-adv_vnt(i,k,j-1))*dytr(jrow)
     &               *cstr(jrow)
# if defined gm_advect
     &             + (adv_vntiso(i,k,j)-adv_vntiso(i,k,j-1))*dytr(jrow)
     &               *cstr(jrow)
# endif
              term = tmask(i,k,j)*(t(i,k,j,n,tau)*dvdy - ADV_Ty(i,k,j))
	      txbt(k,12,n,nth) = txbt(k,12,n,nth) + term 
c
c-----------------------------------------------------------------------
c             vertical advection (flux form) of tracer
c-----------------------------------------------------------------------
c
              term = -tmask(i,k,j)*ADV_Tz(i,k,j)
	      txbt(k,4,n,nth) = txbt(k,4,n,nth) + term 
c
c-----------------------------------------------------------------------
c             pure vertical advection of tracer
c-----------------------------------------------------------------------
c
c             - W(T)z = T(W)z - (WT)z
c
              dwdz = (adv_vbt(i,k-1,j)-adv_vbt(i,k,j))
# ifdef partial_cell
     &                                                /dht(i,k,j)
# else
     &                                                *dztr(k)
# endif
#  if defined gm_advect
     &              + (adv_vbtiso(i,k-1,j)-adv_vbtiso(i,k,j))
#   ifdef partial_cell
     &                                                       /dht(i,k,j)
#   else
     &                                                       *dztr(k)
#   endif
#  endif
              term = tmask(i,k,j)*(t(i,k,j,n,tau)*dwdz - ADV_Tz(i,k,j))
	      txbt(k,13,n,nth) = txbt(k,13,n,nth) + term 
c
c-----------------------------------------------------------------------
c             zonal diffusion of tracer
c-----------------------------------------------------------------------
c
              term = tmask(i,k,j)*DIFF_Tx(i,k,j)
	      txbt(k,5,n,nth) = txbt(k,5,n,nth) + term 
c
c-----------------------------------------------------------------------
c             meridional diffusion of tracer
c-----------------------------------------------------------------------
c
              term = tmask(i,k,j)*DIFF_Ty(i,k,j)
	      txbt(k,6,n,nth) = txbt(k,6,n,nth) + term 
c
c-----------------------------------------------------------------------
c             vertical diffusion of tracer
c-----------------------------------------------------------------------
c
              term = tmask(i,k,j)*DIFF_Tz(i,k,j)
# if defined implicitvmix || defined isoneutralmix
     &               + tmask(i,k,j)*zzi(i,k,j)
# endif
	      txbt(k,7,n,nth) = txbt(k,7,n,nth) + term 
c
c-----------------------------------------------------------------------
c             tracer source term
c-----------------------------------------------------------------------
c
# ifdef source_term
              term = tmask(i,k,j)*source(i,k,j)
	      txbt(k,8,n,nth) = txbt(k,8,n,nth) + term 
# else
	      txbt(k,8,n,nth) = 0.0 
# endif
c
            enddo
c
c-----------------------------------------------------------------------
c           surface tracer flux
c-----------------------------------------------------------------------
c
            k    = 1
            term = tmask(i,k,j)*stf(i,j,n)
	    txbtsf(n,nth) = txbtsf(n,nth) + term
          enddo
        endif
      enddo
c
# ifdef timing
      call toc ('diagnostic', 'xbt')
# endif
c
      return
      end


      subroutine txbt2 (joff, js, je, iterm)
c
c=======================================================================
c     Accumulate parts of d(tracer)/dt and change in tracer variance
c
c     input:
c       joff  = offset relating "j" in the MW to latitude "jrow"
c       js    = starting row in the MW
c       je    = ending row in the MW
c       iterm = 1  => total change
c       iterm = 10 => change due to filtering      
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use memory_window
# include "param.h"
# include "accel.h"
# include "mw.h"
# include "scalar.h"
# include "diag_xbt.h"
c
# ifdef timing
      call tic ('diagnostic', 'xbt')
# endif
c
      if (iterm .ne. 1 .and. iterm .ne. 10) then
        write (stdout,*) '=>Error: iterm=',iterm,' in txbt2'
	call abort()
      endif
c
      do j=js,je
        jrow  = j + joff
        if (nsxbt(jrow) .ne. 0) then
	  do nth=nsxbt(jrow),nexbt(jrow)
	    i   = ixbt(nth)
	    lev = kxbt(nth)
	    do n=1,nt
              do k=1,lev
c
c-----------------------------------------------------------------------
c               d/dt(tracer)
c-----------------------------------------------------------------------
c
                r2dt = c1/(c2dtts*dtxcel(k))
                term = tmask(i,k,j)*(t(i,k,j,n,taup1) -
     &                               t(i,k,j,n,taum1))*r2dt
	        txbt(k,iterm,n,nth) = txbt(k,iterm,n,nth) + term 
c
c-----------------------------------------------------------------------
c               change in variance of tracer
c-----------------------------------------------------------------------
c
                if (iterm .eq. 1) then
                  term = tmask(i,k,j)*(t(i,k,j,n,taup1)**2-
     &                                 t(i,k,j,n,taum1)**2)
	          txbt(k,14,n,nth) = txbt(k,14,n,nth) + term 
                endif
              enddo
            enddo
          enddo
	endif
      enddo
c
# ifdef timing
      call toc ('diagnostic', 'xbt')
# endif
c
      return
      end


      subroutine uxbt1 (joff, js, je, n)
c
c=======================================================================
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       n    = (1,2) for (u,v) component of velocity
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use arrays_2d
      use memory_window
# include "param.h"
# include "grdvar.h"
# include "hmixc.h"
# include "mw.h"
# include "scalar.h"
# include "vmixc.h"
# include "diag_xbt.h"
# include "fdifm.h"
c
# ifdef timing
      call tic ('diagnostic', 'xbt')
# endif
      do j=js,je
        jrow = j + joff
        if (nsxbt(jrow) .ne. 0) then
c
c
c         only allow n = 1 (u component) or 2 (v component)
c
          if (n .gt. 2) then
            write (stdout,*) '=> Error: n=',n,' in uxbt1.F'
            call abort()
          endif
c
c-----------------------------------------------------------------------
c         accumulate data for the nth XBT for one time step
c-----------------------------------------------------------------------
c
          do nth=nsxbt(jrow),nexbt(jrow)
            i   = ixbt(nth)
	    lev = kxbt(nth)
	    do k=1,lev
c
c-----------------------------------------------------------------------
c             pressure term
c-----------------------------------------------------------------------
c
              term = -umask(i,k,j)*grad_p(i,k,j,n)
	      uxbt(k,2,n,nth) = uxbt(k,2,n,nth) + term 
c
c-----------------------------------------------------------------------
c             zonal advection (flux form) of momentum
c-----------------------------------------------------------------------
c
              term = -umask(i,k,j)*ADV_Ux(i,k,j)
	      uxbt(k,3,n,nth) = uxbt(k,3,n,nth) + term 
c
c-----------------------------------------------------------------------
c             pure zonal advection of momentum
c-----------------------------------------------------------------------
c
c             - U(U)x = U(U)x - (UU)x (when n=1)
c             - U(V)x = V(U)x - (UV)x (when n=2)
c
              dudx = (adv_veu(i,k,j)-adv_veu(i-1,k,j))*dxur(i)
     &               *csur(jrow)
              term = umask(i,k,j)*(u(i,k,j,n,tau)*dudx - ADV_Ux(i,k,j))
	      uxbt(k,14,n,nth) = uxbt(k,14,n,nth) + term 
c
c-----------------------------------------------------------------------
c             advective metric term
c-----------------------------------------------------------------------
c
              term = ADV_metric(i,k,j,n)
	      uxbt(k,13,n,nth) = uxbt(k,13,n,nth) + term 
c
c-----------------------------------------------------------------------
c             meridional advection (flux form) of momentum
c-----------------------------------------------------------------------
c
              term = -umask(i,k,j)*ADV_Uy(i,k,j)
	      uxbt(k,4,n,nth) = uxbt(k,4,n,nth) + term 
c
c-----------------------------------------------------------------------
c             pure meridional advection of momentum
c-----------------------------------------------------------------------
c
c             - V(U)y = U(V)y - (VU)y (when n=1)
c             - V(V)y = V(V)y - (VV)y (when n=2)
c
              dvdy = (adv_vnu(i,k,j)-adv_vnu(i,k,j-1))*dyur(jrow)
     &               *csur(jrow)
              term = umask(i,k,j)*(u(i,k,j,n,tau)*dvdy - ADV_Uy(i,k,j))
	      uxbt(k,15,n,nth) = uxbt(k,15,n,nth) + term 
c
c-----------------------------------------------------------------------
c             vertical advection (flux form) of momentum
c-----------------------------------------------------------------------
c
              term = -umask(i,k,j)*ADV_Uz(i,k,j)
	      uxbt(k,5,n,nth) = uxbt(k,5,n,nth) + term 
c
c-----------------------------------------------------------------------
c             pure vertical advection of momentum
c-----------------------------------------------------------------------
c
c             - W(U)z = U(W)z - (WU)z (when n=1)
c             - W(V)z = V(W)z - (WV)z (when n=2)
c
              dwdz = (adv_vbu(i,k-1,j)-adv_vbu(i,k,j))
# ifdef partial_cell
     &                                                /dhu(i,k,j)
# else
     &                                                *dztr(k)
# endif
              term = umask(i,k,j)*(u(i,k,j,n,tau)*dwdz - ADV_Uz(i,k,j))
	      uxbt(k,16,n,nth) = uxbt(k,16,n,nth) + term 
c
c-----------------------------------------------------------------------
c             zonal diffusion of momentum
c-----------------------------------------------------------------------
c
              term = umask(i,k,j)*DIFF_Ux(i,k,j)
	      uxbt(k,6,n,nth) = uxbt(k,6,n,nth) + term 
c
c-----------------------------------------------------------------------
c             meridional diffusion of momentum
c-----------------------------------------------------------------------
c
              term = umask(i,k,j)*DIFF_Uy(i,k,j)
	      uxbt(k,7,n,nth) = uxbt(k,7,n,nth) + term 
c
c-----------------------------------------------------------------------
c             diffusive metric term
c-----------------------------------------------------------------------
c
              term = umask(i,k,j)*DIFF_metric(i,k,j,n)
	      uxbt(k,9,n,nth) = uxbt(k,9,n,nth) + term 
c
c-----------------------------------------------------------------------
c             vertical diffusion of momentum
c-----------------------------------------------------------------------
c
              term = umask(i,k,j)*DIFF_Uz(i,k,j)
# ifdef implicitvmix
     &              +umask(i,k,j)*zzi(i,k,j)
# endif
	      uxbt(k,8,n,nth) = uxbt(k,8,n,nth) + term 
c
c-----------------------------------------------------------------------
c             coriolis term
c-----------------------------------------------------------------------
c
              term = umask(i,k,j)*CORIOLIS(i,k,j,n)
	      uxbt(k,10,n,nth) = uxbt(k,10,n,nth) + term 
c
c-----------------------------------------------------------------------
c             accumulate the source term
c-----------------------------------------------------------------------
c
# ifdef source_term
              term = umask(i,k,j)*source(i,k,j)
	      uxbt(k,11,n,nth) = uxbt(k,11,n,nth) + term 
# else
	      uxbt(k,11,n,nth) = 0.0 
# endif
c
c-----------------------------------------------------------------------
c             accumulate u, v, and w
c-----------------------------------------------------------------------
c
              term = umask(i,k,j)*u(i,k,j,n,tau)
	      uxbt(k,17,n,nth) = uxbt(k,17,n,nth) + term 
c
              if (n .eq. 2) then
                term = p5*(adv_vbu(i,k,j)+adv_vbu(i,k-1,j))*umask(i,k,j)
	        xbtw(k,nth) = xbtw(k,nth) + term 
              endif
            enddo
c
c-----------------------------------------------------------------------
c           accumulate the surface momentum flux
c-----------------------------------------------------------------------
c
            k    = 1
            term = umask(i,k,j)*smf(i,j,n)
            uxbtsf(n,nth) = uxbtsf(n,nth) + term
c
c-----------------------------------------------------------------------
c           update accumulation counter once per time step
c-----------------------------------------------------------------------
c
            if (nth .eq. numxbt .and. (n .eq. 2))  nxbtts = nxbtts + 1
          enddo
        endif
      enddo
c
# ifdef timing
      call toc ('diagnostic', 'xbt')
# endif
c
      return
      end




      subroutine uxbt2 (joff, js, je, c2dtuv, acor)
c
c=======================================================================
c     Accumulate  d/dt and the implicit coriolis terms in the 
c     momentum equations for XBTs
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       c2dtuv = (2*dtuv,dtuv) on (lpfrod,mixing) time steps
c       acor   = implicit factor 
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
c
      use arrays_2d
      use memory_window
# include "param.h"
# include "mw.h"
# include "diag_xbt.h"
c
# ifdef timing
      call tic ('diagnostic', 'xbt')
# endif
c
c-----------------------------------------------------------------------
c     local constants
c-----------------------------------------------------------------------
c
      r2dt = c1/c2dtuv
c
c-----------------------------------------------------------------------
c     d/dt of velocity (external mode part of tau+1 will be added 
c     later when the external mode is solved)
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow  = j + joff
        if (nsxbt(jrow) .ne. 0) then
	  do nth=nsxbt(jrow),nexbt(jrow)
	    i   = ixbt(nth)
	    lev = kxbt(nth)
	    do n=1,2
              do k=1,lev
c
c-----------------------------------------------------------------------
c               d/dt of velocity (external mode part of tau+1 will be  
c               added later when the external mode is solved)
c-----------------------------------------------------------------------
c
                term = umask(i,k,j)*(u(i,k,j,n,taup1) - 
     &                               u(i,k,j,n,taum1))*r2dt
	        uxbt(k,1,n,nth) = uxbt(k,1,n,nth) + term 
c
c
c-----------------------------------------------------------------------
c               implicit coriolis term (external mode part will be 
c               added later when external mode is solved)
c-----------------------------------------------------------------------
c
                if (acor .ne. c0) then
#ifdef rot_grid
                  term = umask(i,k,j)*acor*cori(i,jrow,n)*
#else
                  term = umask(i,k,j)*acor*cori(jrow,n)*
#endif
     &                     (u(i,k,j,3-n,taup1) - u(i,k,j,3-n,taum1))
	          uxbt(k,10,n,nth) = uxbt(k,10,n,nth) + term 
	        endif
              enddo
            enddo
	  enddo
	endif
      enddo
c
# ifdef timing
      call toc ('diagnostic', 'xbt')
# endif
c
      return
      end


      subroutine uxbt3
c
c=======================================================================
c     Accumulate surface pressure gradients, external mode
c     part of d/dt and external moe part of implicit coriolis term
c     for XBTs
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c=======================================================================
c
      use arrays_2d
      use memory_window
# include "param.h"
# include "emode.h"
# include "grdvar.h"
# include "mw.h"
# include "scalar.h"
# include "diag_xbt.h"
      parameter (is=1, ie=1, js=1, je=1)
      dimension psgrad(is:ie,js:je,2)
c
# ifdef timing
      call tic ('diagnostic', 'xbt')
# endif
c
      r2dtuv = c1/c2dtuv
      do n=1,numxbt
        i      = ixbt(n)
	jrow   = jxbt(n)
	lev    = kxbt(n)
        kz     = kmu(i,jrow)
        if (kz .ne. 0) then
c
c         construct the surface pressure gradients for pt (i,jrow)
c         note: it will be stored in psgrad(is,js,:)
c
          call calc_psgrad (psgrad, uext, vext, jrow, jrow, i, i)
#ifdef rot_grid
          atosp = acor*cori(i,jrow,1)
#else
          atosp = acor*cori(jrow,1)
#endif
          kl = min(kz,lev)
          do k=1,kl
            uxbt(k,1,1,n)  = uxbt(k,1,1,n)  + uext*r2dtuv
            uxbt(k,1,2,n)  = uxbt(k,1,2,n)  + vext*r2dtuv
            uxbt(k,12,1,n) = uxbt(k,12,1,n) - psgrad(is,js,1)
            uxbt(k,12,2,n) = uxbt(k,12,2,n) - psgrad(is,js,2)
            uxbt(k,10,1,n) = uxbt(k,10,1,n) + vext*atosp
            uxbt(k,10,2,n) = uxbt(k,10,2,n) - uext*atosp
          enddo
        endif
      enddo
c
# ifdef timing
      call toc ('diagnostic', 'xbt')
# endif
c
      return
      end


      subroutine xbto
c
c=======================================================================
c     XBT output: average and save all XBT data
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      character*128 file_name, file_stamp
# include "param.h"
# include "coord.h"
# include "grdvar.h"
# include "iounit.h"
# include "scalar.h"
# include "switch.h"
# include "tmngr.h"
# include "diag_xbt.h"
      dimension errt(nt)
c
# ifdef timing
      call tic ('diagnostic', 'xbt')
# endif
c
c     reltim defines the end of the averaged data in years
c
      reltim = prelyr
      period = dtts*nxbtts
c
      if (ioxbt .eq. stdout .or. ioxbt .lt. 0) then
	write (stdout,'(//,20x,a,/)')
     &  'A V E R A G E D    X B T   S T A T I O N S'
	write (stdout,8800) numxbt, itt, period*secday, stamp  
      endif
c
      do n=1,numxbt
        i     = ixbt(n)
	jrow  = jxbt(n)
	lev   = kxbt(n)
c
        if (ioxbt .eq. stdout .or. ioxbt .lt. 0) then
	  write (stdout,8900)  n, yt(jrow), xt(i)
	endif
c
c-----------------------------------------------------------------------
c       construct time mean quantities
c-----------------------------------------------------------------------
c
        rnavg = c1/nxbtts
c
c-----------------------------------------------------------------------
c       average the data for xbt number n, write it out, then zero
c       the accumulators
c-----------------------------------------------------------------------
c
	do i=1,nt
	  do m=1,ntxbt
	    do k=1,lev
	      txbt(k,m,i,n) = rnavg*txbt(k,m,i,n)
            enddo
          enddo
	  txbtsf(i,n) = rnavg*txbtsf(i,n)
        enddo
	do i=1,2
	  do m=1,nuxbt
	    do k=1,lev
	      uxbt(k,m,i,n) = rnavg*uxbt(k,m,i,n)
            enddo
          enddo
	  uxbtsf(i,n) = rnavg*uxbtsf(i,n)
        enddo
	do k=1,lev
	  xbtw(k,n) = rnavg*xbtw(k,n)
	enddo
c
c       construct change due to convection and filtering
c
        do i=1,nt
          do k=1,lev
	    tconv = txbt(k,10,i,n) - txbt(k,9,i,n)
	    tfilt = txbt(k,1,i,n) - txbt(k,10,i,n)
	    txbt(k,9,i,n)  = tconv
	    txbt(k,10,i,n) = tfilt
	  enddo  
	enddo  
c
c-----------------------------------------------------------------------
c       write out the results
c-----------------------------------------------------------------------
c
        if (ioxbt .eq. stdout .or. ioxbt .lt. 0) then
	  do k=1,lev
	    write (stdout,8700) k, zt(k)*0.01
            do m=1,nuxbt
	      write (stdout,9000) xnameu(m,1), uxbt(k,m,1,n)
     &,                           xnameu(m,2), uxbt(k,m,2,n)
            enddo
            erru = c0
            errv = c0
            do lll=2,13
              erru = erru + uxbt(k,lll,1,n)
              errv = errv + uxbt(k,lll,2,n)
	    enddo
	    erru = uxbt(k,1,1,n) - erru
	    errv = uxbt(k,1,2,n) - errv
	    write (stdout,9000) ' error      ', erru
     &,                         ' error      ', errv 
	    if (k .eq. 1) then
              write (stdout,9000) xnamex(2),uxbtsf(1,n)
     &,                           xnamex(3),uxbtsf(2,n)
	    endif
	    write (stdout,*) ' '
	    write (stdout,*) '           ',xnamex(4),' = ',xbtw(k,n)
c
	    write (stdout,9001) (i,i=1,nt)
            do m=1,ntxbt
	      write (stdout,9002) xnamet(m), (txbt(k,m,i,n),i=1,nt)
            enddo
	    do i=1,nt
	      errt(i) = c0
	      do lll=2,10
	        errt(i) = errt(i) + txbt(k,lll,i,n)
	      enddo
	      errt(i) = errt(i) - txbt(k,1,i,n)
	    enddo
	    write (stdout,9002) '  error     ', (errt(i),i=1,nt)
	    if (k .eq. 1) then
              write (stdout,9002) xnamex(1), (txbtsf(i,n),i=1,nt)
	    endif
          enddo
        endif
c
      enddo

      file_name = file_stamp ('xbt', xstamp, '.dta')
      call xbt_out (trim(file_name), reltim, period*secday,stamp)
c
c-----------------------------------------------------------------------
c       zero the accumulators
c-----------------------------------------------------------------------
c
      do n=1,numxbt
	do i=1,nt
	  do m=1,ntxbt
	    do k=1,lev
	      txbt(k,m,i,n) = c0
            enddo
          enddo
	  txbtsf(i,n) = c0
        enddo
	do i=1,2
	  do m=1,nuxbt
	    do k=1,lev
	      uxbt(k,m,i,n) = c0
            enddo
          enddo
	  uxbtsf(i,n) = c0
        enddo
	do k=1,lev
	  xbtw(k,n) = c0
	enddo
      enddo
c
c-----------------------------------------------------------------------
c     zero the "averaging" counter for the next averaging period
c-----------------------------------------------------------------------
c
      nxbtts = 0
c
# ifdef timing
      call toc ('diagnostic', 'xbt')
# endif
c
      return
8700  format (/,1x,' at model level ',i3,': depth = '
     &,f8.2,'m'/) 
8800  format (//1x,'===> XBT diagnostic: ',i4, ' stations at ts = ',i7
     &,', averaged over ', f8.2,' days ending on ',a32)
8900  format (//1x,' XBT station #',i4,' location:  lat =',f6.2
     &,       ', lon =',f6.2)
9000  format (1x, a12,' = ',1pe14.7,2x, a12,' = ',1pe14.7,2x
     &,           a12,' = ',1pe14.7,2x, a12,' = ',1pe14.7,2x)
9001  format (/1x,14x, 3x,'tracer # ',i2,2x, 3x,'tracer # ',i2,2x
     &,                3x,'tracer # ',i2,2x, 3x,'tracer # ',i2,2x)
9002  format (1x,a12, 2x,1pe14.7, 2x,1e14.7, 2x,1e14.7, 2x,1e14.7)
9999  format (/1x
     & ,'===> ',i4,' XBT"s written to file xbt.dta on ts = ',i10
     &,        ', averaged over = ',f7.2, 'days ', a,/)
      end


#else
      return
      end
#endif
