      subroutine adv_flux_lin (joff, js, je, is, ie, n)
c
c=======================================================================
c     Linearized advective tracer flux
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c=======================================================================
c
#ifdef linearized_advection
      use memory_window
# include "param.h"
# include "mw.h"
      parameter (istrt=2, iend=imt-1)
c
c-----------------------------------------------------------------------
c     limit the longitude indices based on those from the argument list
c     Note: this is currently bypassed. istrt and iend are set as
c           parameters to optimize performance
c-----------------------------------------------------------------------
c
c     istrt = max(2,is)
c     iend  = min(imt-1,ie)
c
c-----------------------------------------------------------------------
c     advective flux across eastern and northern face of T cells
c     is zero due to linearization about state of no motion
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c     calculate 2*advective flux across bottom face of T cells.
c     (It`s done this way for performance issues)
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=1,km-1
          do i=istrt,iend
	    adv_fb(i,k,j)  = adv_vbt(i,k,j)*(tbarz(k,n) + tbarz(k+1,n))
          enddo
        enddo
      enddo
#endif
      return
      end



      subroutine adv_flux_2nd (joff, js, je, is, ie, n)
c
c=======================================================================
c     2nd order advective tracer flux
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c        n   = tracer 
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c=======================================================================
c
      use memory_window
#include "param.h"
#ifdef second_order_tracer_advection
# if defined gm_advect
#  include "isopyc.h" 
# endif
# include "mw.h"
      parameter (istrt=2, iend=imt-1)
c
c-----------------------------------------------------------------------
c     calculate 2*advective flux across east, north & bottom face of
c      T cells. (It`s done this way for performance issues)
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=1,km
          do i=istrt-1,iend
            totvel = adv_vet(i,k,j)
# if defined gm_advect 
     &              + adv_vetiso(i,k,j)
# endif
	    adv_fe(i,k,j) = totvel*(t(i,  k,j,n,tau) + 
     &                              t(i+1,k,j,n,tau))
          enddo
        enddo
      enddo
      do j=js-1,je
        do k=1,km
          do i=istrt,iend
            totvel = adv_vnt(i,k,j)
# if defined gm_advect 
     &              + adv_vntiso(i,k,j)
# endif
	    adv_fn(i,k,j) = totvel*(t(i,k,j,n,tau) + 
     &                              t(i,k,j+1,n,tau))
          enddo
        enddo
      enddo
      do j=js,je
        do k=1,km-1
          do i=istrt,iend
            totvel = adv_vbt(i,k,j)
# if defined gm_advect 
     &              + adv_vbtiso(i,k,j)
# endif
	    adv_fb(i,k,j)  = totvel*(t(i,k,  j,n,tau) +
     &                               t(i,k+1,j,n,tau))
          enddo
        enddo
      enddo
#endif
      return
      end



      subroutine adv_flux_bottom_upwind (joff, js, je, is, ie, n)
c
c=======================================================================
c     upwind advective tracer flux
c     Use only over the bottom grid cells
c     Must be used in conjunction with another advection 
c     scheme for the interior cells
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c        n   = tracer 
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c     author:   s.m.griffies         e-mail smg@gfdl.gov
c=======================================================================
c
      use arrays_2d
#include "param.h"
#ifdef bottom_upwind
      parameter (istrt=2, iend=imt-1)
# include "grdvar.h"
# if defined gm_advect
#  include "isopyc.h" 
# endif
# include "mw.h"
      limit(j) = min(j+joff,jmt) - joff
c
      lag = taum1
      kupwind = 1
c
c-----------------------------------------------------------------------
c     calculate 2*advective flux across eastern face of T cells.
c     (It`s done this way for performance issues)
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        do i=istrt,iend
          do kk=1,kupwind
            k   = max(1,kmt(i,jrow)-kk+1)
            totvel = adv_vet(i,k,j)
# if defined gm_advect
     &              + adv_vetiso(i,k,j)
# endif
            upos = (totvel + abs(totvel))
     &             *tmask(i,k,j)*tmask(i+1,k,j)
            uneg = (totvel - abs(totvel))
     &             *tmask(i+1,k,j)*tmask(i,k,j)
c
	    adv_fe(i,k,j) =  upos*t(i,  k,j,n,lag) 
     &                     + uneg*t(i+1,k,j,n,lag)
          enddo
        enddo
        call setbcx (adv_fe(1,1,j), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     calculate 2*advective flux across northern face of T cells.
c     (It`s done this way for performance issues)
c-----------------------------------------------------------------------
c
      do j=js-1,je
        jrow = j + joff
        do i=istrt,iend
          do kk=1,kupwind
            k   = max(1,kmt(i,jrow)-kk+1)
            totvel = adv_vnt(i,k,j)
# if defined gm_advect
     &              + adv_vntiso(i,k,j)
# endif
            vpos = (totvel + abs(totvel))
     &             *tmask(i,k,j)*tmask(i,k,j+1)
            vneg = (totvel - abs(totvel))
     &             *tmask(i,k,j+1)*tmask(i,k,j)
c
	    adv_fn(i,k,j) =    vpos*t(i,k,j  ,n,lag) 
     &                       + vneg*t(i,k,j+1,n,lag)
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     calculate 2*advective flux across bottom face of T cells.
c     (It`s done this way for performance issues)
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        do i=istrt,iend
          do kk=1,kupwind
            k   = max(1,kmt(i,jrow)-kk)
            kp1 = min(k+1,km)
            totvel = adv_vbt(i,k,j)
# if defined gm_advect
     &              + adv_vbtiso(i,k,j)
# endif
            wpos = (totvel + abs(totvel))
     &             *tmask(i,kp1,j)*tmask(i,k,j)
            wneg = (totvel - abs(totvel))
     &             *tmask(i,k,j)*tmask(i,kp1,j)
c
	    adv_fb(i,k,j)  =   wneg*t(i,k  ,j,n,lag) 
     &                       + wpos*t(i,kp1,j,n,lag)
          enddo
        enddo
      enddo 
#endif
      return
      end



      subroutine adv_flux_quick (joff, js, je, is, ie, n)
c
c=======================================================================
c     3rd order advective tracer flux
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c        n   = tracer 
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c=======================================================================
c
#ifdef quicker
      use memory_window
# include "param.h"
      parameter (istrt=2, iend=imt-1)
# include "grdvar.h"
# if defined gm_advect
#  include "isopyc.h" 
# endif
# include "mw.h"
# include "taskrows.h"
      limit(j) = min(j+joff,jmt) - joff
c
# ifdef ncar_upwind3
      lag = tau
# else
      lag = taum1
# endif
c
c-----------------------------------------------------------------------
c     calculate 2*advective flux across eastern face of T cells.  (It`s
c     done this way for performance issues).  In grid points near land
c     points and near bottom and surface, project the curvature coeff of
c     the land point to the center point.
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=1,kbot
          kl = k
          do i=istrt,iend-1
	    ip2 = i + 2
            totvel = adv_vet(i,k,j)
# if defined gm_advect 
     &              + adv_vetiso(i,k,j)
# endif
            upos = p5*(totvel + abs(totvel))
     &                             *tmask(i,kl,j)*tmask(i+1,kl,j)
            uneg = p5*(totvel - abs(totvel))
     &                             *tmask(i+1,kl,j)*tmask(i,kl,j)
c
c           eastmsk is one only in grid points just east of land:
            eastmsk=tmask(i,kl,j)*(1-tmask(i-1,kl,j))
c           westmsk is one only in grid points just west of land:
            westmsk=tmask(i+1,kl,j)*(1-tmask(ip2,kl,j))
	    adv_fe(i,k,j) =
     &                  totvel*(
     &                          quick_x(i,1)*t(i,  k,j,n,tau) 
     &                        + quick_x(i,2)*t(i+1,k,j,n,tau))
     &                  - upos*(curv_xp(i,1)*t(i+1,k,j,n,lag)
     &                         +curv_xp(i,2)*t(i  ,k,j,n,lag)
     &                         +curv_xp(i,3)
     &        *(t(i-1,k,j,n,lag)*(1-eastmsk)+t(i  ,k,j,n,lag)*eastmsk))
     &                  - uneg*(curv_xn(i,1)
     &        *(t(ip2,k,j,n,lag)*(1-westmsk)+t(i+1,k,j,n,lag)*westmsk)
     &                         +curv_xn(i,2)*t(i+1,k,j,n,lag)
     &                         +curv_xn(i,3)*t(i  ,k,j,n,lag))
          enddo
        enddo
c
        do k=1,kbot
          kl = k
          i=iend
	  ip2 = 3
          totvel = adv_vet(i,k,j)
# if defined gm_advect
     &              + adv_vetiso(i,k,j)
# endif
          upos = p5*(totvel + abs(totvel))
     &                             *tmask(i,kl,j)*tmask(i+1,kl,j)
          uneg = p5*(totvel - abs(totvel))
     &                             *tmask(i+1,kl,j)*tmask(i,kl,j)
c
c         eastmsk is one only in grid points just east of land:
          eastmsk=tmask(i,kl,j)*(1-tmask(i-1,kl,j))
c         westmsk is one only in grid points just west of land:
          westmsk=tmask(i+1,kl,j)*(1-tmask(ip2,kl,j))
	  adv_fe(i,k,j) =
     &                totvel*(
     &                          quick_x(i,1)*t(i,  k,j,n,tau)
     &                        + quick_x(i,2)*t(i+1,k,j,n,tau))
     &                  - upos*(curv_xp(i,1)*t(i+1,k,j,n,lag)
     &                         +curv_xp(i,2)*t(i  ,k,j,n,lag)
     &                         +curv_xp(i,3)
     &        *(t(i-1,k,j,n,lag)*(1-eastmsk)+t(i  ,k,j,n,lag)*eastmsk))
     &                  - uneg*(curv_xn(i,1)
     &        *(t(ip2,k,j,n,lag)*(1-westmsk)+t(i+1,k,j,n,lag)*westmsk)
     &                         +curv_xn(i,2)*t(i+1,k,j,n,lag)
     &                         +curv_xn(i,3)*t(i  ,k,j,n,lag))
        enddo
        call setbcx (adv_fe(1,1,j), imt, kbot)
      enddo
c
c-----------------------------------------------------------------------
c     calculate 2*advective flux across northern face of T cells.
c     (It`s done this way for performance issues)
c-----------------------------------------------------------------------
c
      if (js+joff .eq. jscomp) then
        if (joff .le. 0) then
          do j=1,js-1
	    do k=1,kbot
	      do i=2,imt-1
	        adv_f4n(i,k,j,n) = c0
	      enddo
	    enddo
	  enddo
          jsfn = js
        else
          jsfn = js-1
	endif
      else
        jsfn = js
      endif
      do j=jsfn,je
        jrow = j + joff
        jp2 = limit(j+2)
        do k=1,kbot
          kl = k
          do i=istrt,iend
            totvel = adv_vnt(i,k,j)
# if defined gm_advect
     &              + adv_vntiso(i,k,j)
# endif
            vpos = p5*(totvel + abs(totvel))
     &                             *tmask(i,kl,j)*tmask(i,kl,j+1)
            vneg = p5*(totvel - abs(totvel))
     &                             *tmask(i,kl,j+1)*tmask(i,kl,j)
c
c           rnormsk is one only in grid points just north of land:
            rnormsk=tmask(i,kl,j)*(1-tmask(i,kl,j-1))
c           soutmsk is one only in grid points just south of land:
            soutmsk=tmask(i,kl,j+1)*(1-tmask(i,kl,jp2))
	    adv_f4n(i,k,j,n) = totvel*(
     &                          quick_y(jrow,1)*t(i,k,j  ,n,tau)
     &                        + quick_y(jrow,2)*t(i,k,j+1,n,tau))
     &                  - vpos*(curv_yp(jrow,1)*t(i,k,j+1,n,lag)
     &                         +curv_yp(jrow,2)*t(i,k,j  ,n,lag)
     &                         +curv_yp(jrow,3)*
     &         (t(i,k,j-1,n,lag)*(1-rnormsk)+t(i,k,j  ,n,lag)*rnormsk))
     &                  - vneg*(curv_yn(jrow,1)*
     &         (t(i,k,jp2,n,lag)*(1-soutmsk) +t(i,k,j+1,n,lag)*soutmsk)
     &                         +curv_yn(jrow,2)*t(i,k,j+1,n,lag)
     &                         +curv_yn(jrow,3)*t(i,k,j  ,n,lag))
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     calculate 2*advective flux across bottom face of T cells.
c     (It`s done this way for performance issues)
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=2,kbot-2
          do i=istrt,iend
            totvel = adv_vbt(i,k,j)
# if defined gm_advect
     &              + adv_vbtiso(i,k,j)
# endif
            wpos = p5*(totvel + abs(totvel))
     &             *tmask(i,k+2,j)*tmask(i,k+1,j)             
            wneg = p5*(totvel - abs(totvel))
     &                            *tmask(i,k,j)*tmask(i,k+1,j)
c
c           upmsk is one only in grid points just above land:
            upmsk=tmask(i,kl,j)*(1-tmask(i,kl-1,j))
	    adv_fb(i,k,j)  = totvel*(
     &                          quick_z(k,1)*t(i,k  ,j,n,tau) 
     &                        + quick_z(k,2)*t(i,k+1,j,n,tau))
     &                  - wneg*(curv_zp(k,1)*t(i,k+1,j,n,lag)
     &                         +curv_zp(k,2)*t(i,k  ,j,n,lag)
     &                         +curv_zp(k,3)*t(i,k-1,j,n,lag))
     &                  - wpos*(curv_zn(k,1)*t(i,k+2,j,n,lag)
     &                         +curv_zn(k,2)*t(i,k+1,j,n,lag)
     &                         +curv_zn(k,3)*
     &           (t(i,k  ,j,n,lag)*(1-upmsk)+t(i,k+1,j,n,lag)*upmsk))
          enddo
        enddo
	k=1
        do i=istrt,iend
          totvel = adv_vbt(i,k,j)
# if defined gm_advect
     &            + adv_vbtiso(i,k,j)
# endif
          wpos = p5*(totvel + abs(totvel))
     &         *tmask(i,k+2,j)*tmask(i,k+1,j)*tmask(i,k,j)

          wneg = p5*(totvel - abs(totvel))
     &         *tmask(i,k,j)*tmask(i,k+1,j)
c
	  adv_fb(i,k,j)  = totvel*(
     &                        quick_z(k,1)*t(i,k  ,j,n,tau)
     &                      + quick_z(k,2)*t(i,k+1,j,n,tau))
     &                - wneg*(curv_zp(k,1)*t(i,k+1,j,n,lag)
     &                         +curv_zp(k,2)*t(i,k  ,j,n,lag)
     &                         +curv_zp(k,3)*t(i,k  ,j,n,lag))
     &                - wpos*(curv_zn(k,1)*t(i,k+2,j,n,lag)
     &                       +curv_zn(k,2)*t(i,k+1,j,n,lag)
     &                       +curv_zn(k,3)*t(i,k  ,j,n,lag))
        enddo
	k=kbot-1
        do i=istrt,iend
          totvel = adv_vbt(i,k,j)
# if defined gm_advect
     &            + adv_vbtiso(i,k,j)
# endif
          wneg = p5*(totvel - abs(totvel))
     &             *tmask(i,k-1,j)*tmask(i,k,j)*tmask(i,k+1,j)
          wpos = p5*(totvel + abs(totvel))
     &                            *tmask(i,k+1,j)*tmask(i,k,j)
c
	  adv_fb(i,k,j)  = totvel*(
     &                        quick_z(k,1)*t(i,k  ,j,n,tau) 
     &                      + quick_z(k,2)*t(i,k+1,j,n,tau))
     &                - wneg*(curv_zp(k,1)*t(i,k+1,j,n,lag)
     &                       +curv_zp(k,2)*t(i,k  ,j,n,lag)
     &                       +curv_zp(k,3)*t(i,k-1,j,n,lag))
     &                - wpos*(curv_zn(k,1)*t(i,k+1,j,n,lag)
     &                       +curv_zn(k,2)*t(i,k+1,j,n,lag)
     &                       +curv_zn(k,3)*t(i,k  ,j,n,lag))
        enddo
      enddo
#endif
      return
      end



      subroutine adv_flux_4th (joff, js, je, is, ie, n)
c
c=======================================================================
c     4th order advective tracer flux
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c        n   = tracer 
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c=======================================================================
c
#if defined fourth_order_tracer_advection
      use memory_window
# include "param.h"
      parameter (istrt=2, iend=imt-1)
# if defined gm_advect
#  include "isopyc.h" 
# endif
# include "mw.h"
# include "taskrows.h"
      limit(j) = min(j+joff,jmt) - joff
c
c-----------------------------------------------------------------------
c     limit the longitude indices based on those from the argument list
c     Note: this is currently bypassed. istrt and iend are set as
c           parameters to optimize performance
c-----------------------------------------------------------------------
c
c     istrt = max(2,is)
c     iend  = min(imt-1,ie)
c
c-----------------------------------------------------------------------
c     calculate 2*advective flux across eastern face of T cells.
c     (It`s done this way for performance issues)
c-----------------------------------------------------------------------
c
      a2nd = 1.0
      b2nd = 0.0
      a4th = 7.0/6.0
      b4th = -1.0/6.0
c
      do j=js,je
        do k=1,km
          do i=istrt,iend-1
            mask = tmask(i-1,k,j)*tmask(i+2,k,j)
            totvel = adv_vet(i,k,j)
# if defined gm_advect 
     &              + adv_vetiso(i,k,j)
# endif
            adv_fe(i,k,j) = totvel*(
     &       (a2nd*(1.0-mask) + a4th*mask)*(t(i,  k,j,n,tau) + 
     &                                      t(i+1,k,j,n,tau))+
     &       (b2nd*(1.0-mask) + b4th*mask)*(t(i-1,  k,j,n,tau) + 
     &                                      t(i+2,k,j,n,tau)))
          enddo
          i = iend
# ifdef cyclic
	  ip2 = 3
# else
	  ip2 = imt
# endif
          mask = tmask(i-1,k,j)*tmask(ip2,k,j)
          totvel = adv_vet(i,k,j)
# if defined gm_advect 
     &            + adv_vetiso(i,k,j)
# endif
          adv_fe(i,k,j) = totvel*(
     &       (a2nd*(1.0-mask) + a4th*mask)*(t(i,  k,j,n,tau) + 
     &                                      t(i+1,k,j,n,tau))+
     &       (b2nd*(1.0-mask) + b4th*mask)*(t(i-1,  k,j,n,tau) + 
     &                                      t(ip2,k,j,n,tau)))
          adv_fe(1,k,j) = adv_fe(imt-1,k,j)
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     2*advective flux across northern face of T cells is built
c     into ADV_Ty. (It`s done this way for performance issues)
c-----------------------------------------------------------------------
c
      if (js+joff .eq. jscomp) then
        if (joff .le. 0) then
          do j=1,js-1
	    do k=1,km
	      do i=2,imt-1
	        adv_f4n(i,k,j,n) = c0
	      enddo
	    enddo
	  enddo
          jsfn = js
        else
          jsfn = js-1
	endif
      else
        jsfn = js
      endif
      do j=jsfn,je
        jp2 = limit(j+2)
        do k=1,km
          do i=istrt,iend
            mask = tmask(i,k,j-1)*tmask(i,k,jp2)
            totvel = adv_vnt(i,k,j)
# if defined gm_advect 
     &              + adv_vntiso(i,k,j)
# endif
            adv_f4n(i,k,j,n) = totvel*(
     &       (a2nd*(1.0-mask) + a4th*mask)*(t(i,  k,j,n,tau) + 
     &                                      t(i,k,j+1,n,tau))+
     &       (b2nd*(1.0-mask) + b4th*mask)*(t(i,  k,j-1,n,tau) + 
     &                                      t(i,k,jp2,n,tau)))
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     calculate 2*advective flux across bottom face of T cells.
c     (It`s done this way for performance issues)
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=2,km-2
            do i=istrt,iend
              mask = tmask(i,k-1,j)*tmask(i,k+2,j)
              totvel = adv_vbt(i,k,j)
# if defined gm_advect 
     &                + adv_vbtiso(i,k,j)
# endif
              adv_fb(i,k,j) = totvel*(
     &         (a2nd*(1.0-mask) + a4th*mask)*(t(i,  k,j,n,tau) + 
     &                                        t(i,k+1,j,n,tau))+
     &         (b2nd*(1.0-mask) + b4th*mask)*(t(i,k-1,j,n,tau) + 
     &                                        t(i,k+2,j,n,tau)))
            enddo
          enddo
	  k = 1
	  m = km-1
          do i=istrt,iend
            totvel = adv_vbt(i,k,j)
# if defined gm_advect 
     &              + adv_vbtiso(i,k,j)
# endif
            adv_fb(i,k,j) = totvel*(t(i,k  ,j,n,tau) + 
     &                              t(i,k+1,j,n,tau))
c
            totvel = adv_vbt(i,m,j)
# if defined gm_advect 
     &              + adv_vbtiso(i,m,j)
# endif
            adv_fb(i,m,j) = totvel*(t(i,m  ,j,n,tau) + 
     &                              t(i,m+1,j,n,tau))
          enddo
        enddo
#endif
      return
      end      

      subroutine adv_flux_fct (joff, js, je, is, ie, n)
c
c=======================================================================
c     computes advective fluxes using a flux-corrected transport scheme
c
c        for reference see
c        Gerdes, R., C. Koeberle and J. Willebrand, 1991
c        the influence of numerical advection schemes on the results of 
c        ocean general circulation models. Clim Dynamics 5, 211-226 
c        and
c        Zalesak, S. T., 1979: Fully multidimensional flux-corrected 
c        transport algorithms for fluids. J. Comp. Phys. 31, 335-362. 
c
c     input:
c       joff  = offset relating "j" in the MW to latitude "jrow"
c       js    = starting row in the MW
c       je    = ending row in the MW
c       jstrt = starting row in the MW for fct
c       jend  = ending row in the MW for fct
c       is    = starting longitude index in the MW
c       ie    = ending longitude index in the MW
c       istrt = max(2,starting longitude index in the MW)
c       iend  = min(imt-1,ending longitude index in the MW)
c       n     = tracer index
c
c     output:
c       adv_fn = 2*advective flux across northern face of T-cell
c       adv_fe = 2*advective flux across eastern face of T-cell
c       adv_fb = 2*advective flux across bottom face of T-cell
c                ready for use in statement functions as defined 
c                in "fdift.h"
c
c     author:   C. Koeberle       e-mail ckoeberl@awi-bremerhaven.de
c               thanks to  M. Eby who did extensive debugging of a
c               very rough beta version and who also added the 
c               Gent_McWilliams portion and diagnostics to show
c               differences between fct & 2nd order centered advection
c=======================================================================
c               
#ifdef tracer_advection_fct
      use memory_window
      use mpp_io_mod
      character*128 file_name, file_stamp
# if defined fct_dlm1 || !defined fct_dlm2
#  define dlimiter_1
# endif
# if defined fct_dlm2 && !defined fct_dlm1
#  define dlimiter_2
# endif
# include "param.h"
# include "accel.h"
# include "coord.h"
# include "diag.h"
# include "grdvar.h"
# ifdef gm_advect
#  include "isopyc.h"
# endif
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# include "tmngr.h"
c
c---------------------------------------------------------------------
c     dimension local data
c---------------------------------------------------------------------
c
      dimension twodt(km), dcf(imt), Trmin(imt), Trmax(imt)
     &,         Cpos(imt), Cneg(imt), flxlft(imt), flxrgt(imt)
     &,         Rpl(imt,km), Rmn(imt,km), tmaski(imt,km,jmw)
# ifdef fct_3d
     &,         tmpext(imt,km,2)
# endif
# ifdef dlimiter_1
     &,         t_lo(imt,km)
# endif
# include "fdift.h"
c
c-----------------------------------------------------------------------
c      FCT statement functions
c
c      calculate some variables needed as one-dimensional delimiters to 
c      raw anti-diffusive fluxes 
c
c      1) calculate P+/P-: maximal/minimal possible growth/decay of T
c      2) calculate Q+/Q-: maximal/minimal allowed growth/decay of T 
c           the appropriate choice of this "allowed" is the key
c           for the success of this scheme...
c      3) calculate ration R+- of Q+- to P+-, that is maximal/minimal 
c         possible change of T if no limit would be active, 
c         must be at least 1
c
c-----------------------------------------------------------------------
c
c     here the calculations of the flux delimiters are done as statement
c     functions for simplicity always loop along one row of POINTS
c 
c     maximal/minimal possible change at grid point
c     this variable is related to a point
c
      Pplus(i)  = c2dtts*dcf(i)*(max(c0,flxlft(i))-min(c0,flxrgt(i)))
      Pminus(i) = c2dtts*dcf(i)*(max(c0,flxrgt(i))-min(c0,flxlft(i)))
c
c     maximal/minimal feasible change at grid point
c     this variable is related to a point
c
# ifdef dlimiter_1
      Qplus(i)  = Trmax(i) - t_lo(i,k)
      Qminus(i) = t_lo(i,k) - Trmin(i)
# endif
# ifdef dlimiter_2
      Qplus(i)  = Trmax(i) - t_lo(i,k,j+1,n)
      Qminus(i) = t_lo(i,k,j+1,n) - Trmin(i)
# endif
c
c     ratio of Q and P to get delimiter
c     this variable is related to a point
c
      Rpls(i)  = min(1.,tmask(i,k,j+1)*Qplus(i)/(Pplus(i)+epsln))
      Rmns(i)  = min(1.,tmask(i,k,j+1)*Qminus(i)/(Pminus(i)+epsln))
c
c     j_limit1(j) limits the maximum MW row to "latitude" row "jmt-1"
c     j_limit2(j) limits the maximum MW row to "latitude" row "jmt"
c
      j_limit1(j) = min(j+joff,jmt-1) - joff
      j_limit2(j) = min(j+joff,jmt)   - joff
c
c-----------------------------------------------------------------------
c     limit the indices based on those from the argument list
c-----------------------------------------------------------------------
c
      istrt   = max(2,is)
      iend    = min(imt-1,ie)
      istrtm1 = istrt - 1
      iendp1  = iend + 1
      jstrt   = js
      jend    = min(je,jmt-1-joff)
c
c-----------------------------------------------------------------------
c     initialization when calculating jrow 2 
c-----------------------------------------------------------------------
c
      if (joff + js .eq. 2) then
        jstrt = js - 1
        do k=1,km
          do i=istrt-1,iend
	    adv_fn(i,k,jstrt) = c0
	    anti_fn(i,k,jstrt,n) = c0
            R_plusY(i,k,jstrt,n) = c0
            R_minusY(i,k,jstrt,n) = c0
# ifdef fct_3d
            R_plus3(i,k,jstrt,n) = c0
            R_minus3(i,k,jstrt,n) = c0
#  endif
          enddo
        enddo
      endif
c
c-----------------------------------------------------------------------
c     Initialize some diagnostics fields
c-----------------------------------------------------------------------
c
#  ifdef save_fct
      if (fctts .and. .not. euler2 .and. n .le. 2) then
        do j=js,jend
          do i=1,imt
            do k=1,km
              fct_minus_ctr(i,k,j)  = c0
              flux_diffX(i,k,j)     = c0
              flux_diffY(i,k,j)     = c0
              flux_diffZ(i,k,j)     = c0
            enddo
          enddo
        enddo
      endif
#  endif
c
c-----------------------------------------------------------------------
c     create an inverse land mask
c-----------------------------------------------------------------------
c
      do j=1,jmw
        do k=1,km
          do i=1,imt
            tmaski(i,k,j) = c1 - tmask(i,k,j)
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     calculate 2*advective (low order scheme) flux across northern, 
c     eastern and bottom faces of T cells
c-----------------------------------------------------------------------
c
      jlast = j_limit1(jend+1)
      do j=js-1,jlast
        do k=1,km
          do i=istrt,iend
            totadv = adv_vnt(i,k,j)
# if defined gm_advect
     &             + adv_vntiso(i,k,j)
# endif
            adv_fn(i,k,j) = totadv*
     &                       (t(i,k,j,n,taum1) + t(i,k,j+1,n,taum1))
     &                      + abs(totadv)*
     &                       (t(i,k,j,n,taum1) - t(i,k,j+1,n,taum1))
          enddo
        enddo
      enddo
c
      jlast = j_limit1(jend+1)
      do j=js,jlast
        do k=1,km
          do i=istrtm1,iend
            totadv = adv_vet(i,k,j)
# if defined gm_advect
     &             + adv_vetiso(i,k,j)
# endif
            adv_fe(i,k,j) = totadv*
     &                       (t(i,k,j,n,taum1) + t(i+1,k,j,n,taum1))
     &                      + abs(totadv)*
     &                       (t(i,k,j,n,taum1) - t(i+1,k,j,n,taum1))
          enddo
        enddo
c
        do k=1,kmm1
          do i=istrt,iend
            totadv = adv_vbt(i,k,j)
# if defined gm_advect
     &             + adv_vbtiso(i,k,j)
# endif
            adv_fb(i,k,j) = totadv*
     &                      (t(i,k+1,j,n,taum1) + t(i,k,j,n,taum1))
     &                    + abs(totadv)*
     &                      (t(i,k+1,j,n,taum1) - t(i,k,j,n,taum1))
          enddo
        enddo
        do i=istrt,iend
          adv_fb(i,0,j) = adv_vbt(i,0,j)*c2*t(i,1,j,n,taum1)
          adv_fb(i,km,j) = c0
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     main j loop
c-----------------------------------------------------------------------
c
      do j=jstrt,jend
        jrow = (j+1) + joff
        jp2  = j_limit2(j+2)
        jp1  = j_limit1(j+1)
c
c-----------------------------------------------------------------------
c       solve for "tau+1" tracer at center of T cells in row j+1
c       - low order solution -
c-----------------------------------------------------------------------
c
        do k=1,km
	  twodt(k) = c2dtts*dtxcel(k)
          do i=istrt,iend
# ifdef dlimiter_1
	    t_lo(i,k) = (t(i,k,j+1,n,taum1) - twodt(k)
# endif
# ifdef dlimiter_2
	    t_lo(i,k,j+1,n) = (t(i,k,j+1,n,taum1) - twodt(k)
# endif
     &                      *(ADV_Tx(i,k,jp1) + ADV_Ty(i,k,jp1) + 
     &                        ADV_Tz(i,k,jp1))*tmask(i,k,j+1))
          enddo
        enddo
# ifdef dlimiter_1
        call setbcx (t_lo(1,1), imt, km)
# endif
# ifdef dlimiter_2
        call setbcx (t_lo(1,1,j+1,n), imt, km)
# endif
c
c-----------------------------------------------------------------------
c       next calculate raw antidiffusive fluxes, that is high order 
c       scheme flux (leap frog) minus the low order (upstream) 
c-----------------------------------------------------------------------
c
        do k=1,km
          do i=istrtm1,iend
            totadv = adv_vet(i,k,jp1)
# if defined gm_advect
     &             + adv_vetiso(i,k,jp1)
# endif
            anti_fe(i,k,j+1,n) = totadv*(t(i,k,j+1,n,tau) +
     &                           t(i+1,k,j+1,n,tau)) - adv_fe(i,k,jp1)
          enddo
          do i=istrt,iend
            totadv = adv_vnt(i,k,jp1)
# if defined gm_advect
     &             + adv_vntiso(i,k,jp1)
# endif
            anti_fn(i,k,j+1,n) = totadv*(t(i,k,j+1,n,tau) + 
     &                           t(i,k,jp2,n,tau)) - adv_fn(i,k,jp1)
          enddo
        enddo
c
        do k=1,kmm1
          do i=istrt,iend
            totadv = adv_vbt(i,k,jp1)
# if defined gm_advect
     &             + adv_vbtiso(i,k,jp1)
# endif
            anti_fb(i,k,j+1,n) = totadv*(t(i,k,j+1,n,tau) +
     &                           t(i,k+1,j+1,n,tau)) - adv_fb(i,k,jp1)
     &                           *tmask(i,k,j+1)
          enddo
        enddo
        do i=istrt,iend
          anti_fb(i,0,j+1,n) = adv_vbt(i,0,jp1)*c2*t(i,1,j+1,n,taum1)
          anti_fb(i,km,j+1,n) = c0
        enddo
c
c-----------------------------------------------------------------------
c       now calculate and apply one-dimensional delimiters to these 
c       raw antidiffusive fluxes 
c
c       1) calculate T*, that are all halfway neighbors of T
c       2) calculate ratio R+- of Q+- to P+-, that is maximal/minimal 
c          possible change of T if no limit would be active, 
c          must be at least 1
c       3) choose correct ratio depending on direction of flow as a 
c          delimiter
c       4) apply this delimiter to raw antidiffusive flux
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c       delimit x-direction
c-----------------------------------------------------------------------
c
        do k=1,km
c
c         prepare some data for use in statement function   
# ifdef dlimiter_1
c
c         running mean of two adjacent points
c
          do i=istrt,iendp1
            Trmax(i) = p5*(t(i-1,k,j+1,n,tau) + t(i,k,j+1,n,tau))
          enddo
# endif
c
c         extremum of low order solution central point and adjacent 
c         halfway neighbours; check for land
c
          do i=istrt,iend
# ifdef dlimiter_1
            fxa = tmask(i-1,k,j+1)*Trmax(i) +
     &            tmaski(i-1,k,j+1)*t_lo(i,k)
            fxb = tmask(i+1,k,j+1)*Trmax(i+1) +
     &            tmaski(i+1,k,j+1)*t_lo(i,k)
            Trmax(i) = max(fxa,fxb,t_lo(i,k))
            Trmin(i) = min(fxa,fxb,t_lo(i,k))
# endif
# ifdef dlimiter_2
            fxa = tmask(i-1,k,j+1)*t_lo(i-1,k,j+1,n) +
     &            tmaski(i-1,k,j+1)*t_lo(i,k,j+1,n)
            fxb = tmask(i+1,k,j+1)*t_lo(i+1,k,j+1,n) +
     &            tmaski(i+1,k,j+1)*t_lo(i,k,j+1,n)
            Trmax(i) = max(fxa,fxb,t_lo(i,k,j+1,n))
            Trmin(i) = min(fxa,fxb,t_lo(i,k,j+1,n))
# endif
# ifdef fct_3d
            tmpext(i,k,1) = Trmax(i)
            tmpext(i,k,2) = Trmin(i)
# endif
            dcf(i) = cstdxt2r(i,j+1)
            flxlft(i) = anti_fe(i-1,k,j+1,n)
            flxrgt(i) = anti_fe(i,k,j+1,n)
          enddo
c
c         calculate ratio R           
c
          do i=istrt,iend
            Rpl(i,k) = Rpls(i)
            Rmn(i,k) = Rmns(i)
          enddo
          call setbcx (Rpl(1,k), imt, 1)
          call setbcx (Rmn(1,k), imt, 1)
c
c         calculate delimiter using ratio at adjacent points         
c
          do i=istrt,iendp1
            Cpos(i-1) = min(Rpl(i,k),Rmn(i-1,k))
            Cneg(i-1) = min(Rpl(i-1,k),Rmn(i,k))
          enddo
# ifndef tst_fct_his
c
c         finally apply appropriate delimiter to flux
c
          do i=istrtm1,iend
            anti_fe(i,k,j+1,n) = p5*((Cpos(i) + Cneg(i))
     &                               *anti_fe(i,k,j+1,n) +
     &                               (Cpos(i) - Cneg(i))
     &                               *abs(anti_fe(i,k,j+1,n)))
          enddo
# endif
        enddo
c
c-----------------------------------------------------------------------
c       delimit y-direction
c-----------------------------------------------------------------------
c
        do k=1,km
c
c         prepare some data for use in statement function       
c
          do i=istrt,iend
# ifdef dlimiter_1
            fxa = p5*tmask(i,k,j)*(t(i,k,j,n,tau) + 
     &            t(i,k,j+1,n,tau)) +
     &            tmaski(i,k,j)*t_lo(i,k)
            fxb = p5*tmask(i,k,jp2)*(t(i,k,j+1,n,tau) +
     &            t(i,k,jp2,n,tau)) +
     &            tmaski(i,k,jp2)*t_lo(i,k)
            Trmax(i) = max(fxa,fxb,t_lo(i,k))
            Trmin(i) = min(fxa,fxb,t_lo(i,k))
# endif
# ifdef dlimiter_2
            fxa = tmask(i,k,j)*t_lo(i,k,j,n) +
     &            tmaski(i,k,j)*t_lo(i,k,j+1,n)
            fxb = p5*tmask(i,k,jp2)*(t(i,k,j+1,n,tau) +
     &            t(i,k,jp2,n,tau)) +
     &            tmaski(i,k,jp2)*t_lo(i,k,j+1,n)
            Trmax(i) = max(fxa,fxb,t_lo(i,k,j+1,n))
            Trmin(i) = min(fxa,fxb,t_lo(i,k,j+1,n))
# endif
# ifdef fct_3d
            tmpext(i,k,1) = max(Trmax(i),tmpext(i,k,1))
            tmpext(i,k,2) = min(Trmin(i),tmpext(i,k,2))
# endif
            dcf(i) = cstdyt2r(jrow)
            flxlft(i) = anti_fn(i,k,j,n)
            flxrgt(i) = anti_fn(i,k,j+1,n)
          enddo
c
c         calculate ratio R, related to a point          
c
          do i=istrt,iend
            R_plusY (i,k,j+1,n) = Rpls(i)
            R_minusY(i,k,j+1,n) = Rmns(i)
          enddo
c
c         calculate delimiter using ratio at adjacent points
c
          do i=istrt,iend
            Cpos(i) = min(R_plusY(i,k,j+1,n),R_minusY(i,k,j,n))
            Cneg(i) = min(R_plusY(i,k,j,n),R_minusY(i,k,j+1,n))
          enddo
# ifndef tst_fct_his
c
c         finally get delimiter c dependent on direction of flux and 
c         apply it to raw antidiffusive flux
c
          do i=istrt,iend
            anti_fn(i,k,j,n) = p5*((Cpos(i) + Cneg(i))
     &                             *anti_fn(i,k,j,n) +
     &                             (Cpos(i) - Cneg(i))
     &                             *abs(anti_fn(i,k,j,n)))
          enddo
# endif
        enddo
c
c-----------------------------------------------------------------------
c       delimit z-direction
c-----------------------------------------------------------------------
c
        do k=1,km
c
c         prepare some data for use in statement function
c
          do i=istrt,iend
            dcf(i) = 
# ifdef partial_cell
     &               p5/dht(i,k,j)
# else
     &               dzt2r(k)
# endif
            flxlft(i) = anti_fb(i,k,j+1,n)
            flxrgt(i) = anti_fb(i,k-1,j+1,n)
            if (k .gt. 1)then
# ifdef dlimiter_1
              fxa = p5*tmask(i,k-1,j+1)*
     &              (t(i,k-1,j+1,n,tau) + t(i,k,j+1,n,tau)) + 
     &              tmaski(i,k-1,j+1)*t_lo(i,k)
# endif
# ifdef dlimiter_2
              fxa = tmask(i,k-1,j+1)*t_lo(i,k-1,j+1,n) +
     &              tmaski(i,k-1,j+1)*t_lo(i,k,j+1,n)
# endif
            else
# ifdef dlimiter_1
              fxa = t_lo(i,k)
# endif
# ifdef dlimiter_2
              fxa = t_lo(i,k,j+1,n)
# endif
            endif
# ifdef dlimiter_1
            if (k .lt. km) then
              fxb = p5*tmask(i,k+1,j+1)*
     &              (t(i,k,j+1,n,tau)+t(i,k+1,j+1,n,tau)) +
     &              tmaski(i,k+1,j+1)*t_lo(i,k)
# endif
# ifdef dlimiter_2
            if (k .lt. km) then
              fxb = tmask(i,k+1,j+1)*t_lo(i,k+1,j+1,n) +
     &              tmaski(i,k+1,j+1)*t_lo(i,k,j+1,n)
# endif
            else
# ifdef dlimiter_1
              fxb = t_lo(i,k)
# endif
# ifdef dlimiter_2
              fxb = t_lo(i,k,j+1,n)
# endif
            endif
# ifdef dlimiter_1
            Trmax(i) = max(fxa,fxb,t_lo(i,k))
            Trmin(i) = min(fxa,fxb,t_lo(i,k))
# endif
# ifdef dlimiter_2
            Trmax(i) = max(fxa,fxb,t_lo(i,k,j+1,n))
            Trmin(i) = min(fxa,fxb,t_lo(i,k,j+1,n))
# endif
# ifdef fct_3d
            tmpext(i,k,1) = max(Trmax(i),tmpext(i,k,1))
            tmpext(i,k,2) = min(Trmin(i),tmpext(i,k,2))
# endif
          enddo
c
c         calculate delimiter using ratio at adjacent points
c         this variable is related to an arc (between two points,
c         the same way as fluxes are defined.)
c
          do i=istrt,iend
            Rpl(i,k) = Rpls(i)
            Rmn(i,k) = Rmns(i)
          enddo
c
        enddo
c
        do k=1,kmm1
c
c         calculate delimiter using ratio at adjacent points
c         this variable is related to an arc (between two points,
c         the same way as fluxes are defined.)
c
          do i=istrt,iend
            Cneg(i) = min(Rpl(i,k+1),Rmn(i,k))
            Cpos(i) = min(Rpl(i,k),Rmn(i,k+1))
          enddo
# ifndef tst_fct_his
c
c         finally get delimiter c dependent on direction of flux and 
c         apply it to raw antidiffusive flux
c
          do i=istrt,iend
            anti_fb(i,k,j+1,n) = p5*((Cpos(i)+Cneg(i))
     &                               *anti_fb(i,k,j+1,n) +
     &                               (Cpos(i)-Cneg(i))
     &                               *abs(anti_fb(i,k,j+1,n)))
          enddo
# endif
        enddo
        do i=istrt,iend
          anti_fb(i,0,j+1,n)  = c0
          anti_fb(i,km,j+1,n) = c0
        enddo
# ifdef fct_3d
c
c-----------------------------------------------------------------------
c     then calculate and apply 3-d delimiter to just delimited 
c     antidiffusive fluxes   
c-----------------------------------------------------------------------
c
        do k=1,km
c
c         prepare some data for use in statement function
c
          do i=istrt,iend
            Trmax(i) = tmpext(i,k,1)
            Trmin(i) = tmpext(i,k,2)
          enddo
c
          do i=istrt,iend
            R_plus3(i,k,j+1,n) = min(1.,tmask(i,k,j+1)*Qplus(i)/
     &                           (epsln+c2dtts*(
     &                           cstdxt2r(i,j+1)
     &                             *(max(c0,anti_fe(i-1,k,j+1,n)) - 
     &                               min(c0,anti_fe(i,k,j+1,n))) +
     &                           cstdyt2r(jrow)
     &                             *(max(c0,anti_fn(i,k,j,n)) -
     &                               min(c0,anti_fn(i,k,j+1,n))) +
#  ifdef partial_cell
     &                          (p5/dht(i,k,j))
#  else
     &                           dzt2r(k)
#  endif
     &                             *(max(c0,anti_fb(i,k,j+1,n)) -
     &                               min(c0,anti_fb(i,k-1,j+1,n)))
     &                           )))
c
            R_minus3(i,k,j+1,n) = min(1.,tmask(i,k,j+1)*Qminus(i)/
     &                            (epsln+c2dtts*(
     &                            cstdxt2r(i,j+1)
     &                              *(max(c0,anti_fe(i,k,j+1,n)) -
     &                                min(c0,anti_fe(i-1,k,j+1,n))) +
     &                            cstdyt2r(jrow)
     &                              *(max(c0,anti_fn(i,k,j+1,n)) -
     &                                min(c0,anti_fn(i,k,j,n))) +
#  ifdef partial_cell
     &                          (p5/dht(i,k,j))
#  else
     &                           dzt2r(k)
#  endif
     &                              *(max(c0,anti_fb(i,k-1,j+1,n)) -
     &                                min(c0,anti_fb(i,k,j+1,n)))
     &                            )))
c
          enddo
        enddo
        call setbcx (R_plus3(1,1,j+1,n), imt, km)
        call setbcx (R_minus3(1,1,j+1,n), imt, km)
#  ifndef tst_fct_his
c
c       finally apply 3-d delimiters to precorrected fluxes
c
        do k=1,km
          do i=istrt,iendp1
            Cpos(i-1) = min(R_plus3(i,k,j+1,n),R_minus3(i-1,k,j+1,n))
            Cneg(i-1) = min(R_plus3(i-1,k,j+1,n),R_minus3(i,k,j+1,n))
          enddo
          do i=istrtm1,iend
            anti_fe(i,k,j+1,n) = p5*((Cpos(i) + Cneg(i))
     &                               *anti_fe(i,k,j+1,n) +
     &                               (Cpos(i) - Cneg(i))
     &                               *abs(anti_fe(i,k,j+1,n)))
          enddo
c
          do i=istrt,iend
            Cpos(i) = min(R_plus3(i,k,j+1,n),R_minus3(i,k,j,n))
            Cneg(i) = min(R_plus3(i,k,j,n),R_minus3(i,k,j+1,n))
          enddo
          do i=istrt,iend
            anti_fn(i,k,j,n) = p5*((Cpos(i) + Cneg(i))
     &                             *anti_fn(i,k,j,n) +
     &                             (Cpos(i) - Cneg(i))
     &                             *abs(anti_fn(i,k,j,n)))
          enddo
        enddo
c
        do k=1,kmm1
          do i=istrt,iend
            Cneg(i) = min(R_plus3(i,k+1,j+1,n),R_minus3(i,k,j+1,n))
            Cpos(i) = min(R_plus3(i,k,j+1,n),R_minus3(i,k+1,j+1,n))
          enddo
          do i=istrt,iend
            anti_fb(i,k,j+1,n) = p5*((Cpos(i) + Cneg(i))
     &                               *anti_fb(i,k,j+1,n) +
     &                               (Cpos(i) - Cneg(i))
     &                               *abs(anti_fb(i,k,j+1,n)))
          enddo
        enddo
c
#  endif
# endif
c
c-----------------------------------------------------------------------
c       complete advective fluxes by adding low order fluxes to 
c       delimited antidiffusive fluxes    
c-----------------------------------------------------------------------
c
        do k=1,km
          do i=istrtm1,iend
# ifdef tst_fct_los
            anti_fe(i,k,j+1,n) = adv_fe(i,k,jp1)
# else
            anti_fe(i,k,j+1,n) = anti_fe(i,k,j+1,n) + adv_fe(i,k,jp1)
# endif
          enddo
          do i=istrt,iend
# ifdef tst_fct_los
            anti_fn(i,k,j,n) = adv_fn(i,k,j)
# else
            anti_fn(i,k,j,n) = (anti_fn(i,k,j,n) + adv_fn(i,k,j))
# endif
     &                         *tmask(i,k,j)
# ifdef tst_fct_los
            anti_fb(i,k,j+1,n) = adv_fb(i,k,jp1)
# else
            anti_fb(i,k,j+1,n) = (anti_fb(i,k,j+1,n) + adv_fb(i,k,jp1))
# endif
     &                           *tmask(i,k,j+1)
          enddo
        enddo
c
      enddo
c
c-----------------------------------------------------------------------
c     set 2*corrected advective fluxes across northern, eastern and 
c     bottom faces of T cells
c-----------------------------------------------------------------------
c
c
      do j=js-1,jend
        do k=1,km
          do i=istrt,iend
            adv_fn(i,k,j) = anti_fn(i,k,j,n)
          enddo
        enddo
      enddo
c
      do j=js,jend
        do k=1,km
          do i=istrtm1,iend
	    adv_fe(i,k,j) = anti_fe(i,k,j,n)
          enddo
        enddo
c
        do k=1,kmm1
          do i=istrt,iend
	    adv_fb(i,k,j) = anti_fb(i,k,j,n)
          enddo
        enddo
      enddo
c
# if defined save_fct
      if (fctts .and. .not. euler2 .and. n .le. 2) then
c
c       calculate the change in tracer per time step due to fct 
c       above that of centered, to see where fct is important.
c       Do so for the total contribution to the updated tracer 
c       (data placed in array tmaski),
c       as well as look at differences for each of the three 
c       flux components (flux_diffX, flux_diffY, flux_diffZ).
c
        do j=js,jend
          jrow = j + joff
          jp1 = j + 1
          do k=1,km
	    twodt(k) = c2dtts*dtxcel(k)
            do i=istrt,iend
#  if defined gm_advect
              fct_minus_ctr(i,k,j) = twodt(k)*cstdxt2r(i,j)*(
     &                   (adv_fe(i-1,k,j) - (adv_vet(i-1,k,j) +
     &                   adv_vetiso(i-1,k,j))*
     &                   (t(i-1,k,j,n,tau) + t(i,k,j,n,tau))) -
     &                   (adv_fe(i,k,j) - (adv_vet(i,k,j) +
     &                   adv_vetiso(i,k,j))*
     &                   (t(i,k,j,n,tau) + t(i+1,k,j,n,tau))))
              fct_minus_ctr(i,k,j) = fct_minus_ctr(i,k,j)
     &                   + twodt(k)*cstdyt2r(jrow)*(
     &                   (adv_fn(i,k,j-1) - (adv_vnt(i,k,j-1) +
     &                   adv_vntiso(i,k,j-1))*
     &                   (t(i,k,j-1,n,tau) + t(i,k,j,n,tau))) -
     &                   (adv_fn(i,k,j) - (adv_vnt(i,k,j) +
     &                   adv_vntiso(i,k,j))*
     &                   (t(i,k,j,n,tau) + t(i,k,j+1,n,tau))))
              fct_minus_ctr(i,k,j) = fct_minus_ctr(i,k,j) + twodt(k)*
#   ifdef partial_cell
     &                                            (p5/dht(i,k,j))*(
#   else
     &                                             dzt2r(k)*(
#   endif
     &                   (anti_fb(i,k,j,n) - (adv_vbt(i,k,j) +
     &                   adv_vbtiso(i,k,j))*
     &                   (t(i,k,j,n,tau) + t(i,k+1,j,n,tau))) -
     &                   (adv_fb(i,k-1,j) - (adv_vbt(i,k-1,j) +
     &                   adv_vbtiso(i,k-1,j))*
     &                   (t(i,k-1,j,n,tau) + t(i,k,j,n,tau))))
c
              flux_diffX(i,k,j) =
     &             p5*tmask(i,k,j)*tmask(i+1,k,j)
     &             *(
     &             adv_fe(i,k,j)
     &           - (adv_vetiso(i,k,j)+adv_vet(i,k,j))
     &             *(t(i,k,j,n,tau) + t(i+1,k,j,n,tau))
     &             )
              flux_diffY(i,k,j) =
     &             p5*tmask(i,k,j)*tmask(i,k,j+1)
     &             *(
     &             adv_fn(i,k,j)
     &           - (adv_vnt(i,k,j)+adv_vntiso(i,k,j))
     &             *(t(i,k,j,n,tau) + t(i,k,j+1,n,tau))
     &             *tmask(i,k,j)*tmask(i,k,j+1)
     &             )
              flux_diffZ(i,k,j) =
     &             p5*tmask(i,k,j)*tmask(i,k+1,j)
     &             *(
     &             adv_fb(i,k,j)
     &           - (adv_vbt(i,k,j)+adv_vbtiso(i,k,j))
     &             *(t(i,k,j,n,tau) + t(i,k+1,j,n,tau))
     &             )
c
#  else
c
              fct_minus_ctr(i,k,j) = twodt(k)*(
     &                   ((adv_fe(i-1,k,j) - adv_vet(i-1,k,j)*
     &                   (t(i-1,k,j,n,tau) + t(i,k,j,n,tau))) -
     &                   (adv_fe(i,k,j) - adv_vet(i,k,j)*
     &                   (t(i,k,j,n,tau) + t(i+1,k,j,n,tau))))*
     &                   cstdxt2r(i,j) +
     &                   ((adv_fn(i,k,j-1) - adv_vnt(i,k,j-1)*
     &                   (t(i,k,j-1,n,tau) + t(i,k,j,n,tau))) -
     &                   (adv_fn(i,k,j) - adv_vnt(i,k,j)*
     &                   (t(i,k,j,n,tau) + t(i,k,j+1,n,tau))))*
     &                   cstdyt2r(jrow) +
     &                   ((adv_fb(i,k,j) - adv_vbt(i,k,j)*
     &                   (t(i,k,j,n,tau) + t(i,k+1,j,n,tau))) -
     &                   (adv_fb(i,k-1,j) - adv_vbt(i,k-1,j)*
     &                   (t(i,k-1,j,n,tau) + t(i,k,j,n,tau))))*
#   ifdef partial_cell
     &                   (p5/dht(i,k,j))
#   else
     &                   dzt2r(k)
#   endif
     &                                   )
              flux_diffX(i,k,j) =
     &             p5*(
     &             adv_fe(i,k,j)
     &           - adv_vet(i,k,j)*(t(i,k,j,n,tau) + t(i+1,k,j,n,tau))
     &             )
              flux_diffY(i,k,j) =
     &             p5*(
     &             adv_fn(i,k,j)
     &           - adv_vnt(i,k,j)*(t(i,k,j,n,tau) + t(i,k,j+1,n,tau))
     &             )
              flux_diffZ(i,k,j) =
     &             p5*(
     &             adv_fb(i,k,j)
     &           - adv_vbt(i,k,j)*(t(i,k,j,n,tau) + t(i,k+1,j,n,tau))
     &             )
#  endif
            enddo
          enddo
	  call setbcx (fct_minus_ctr(1,1,j), imt, km)
	  call setbcx (flux_diffX(1,1,j), imt, km)
	  call setbcx (flux_diffY(1,1,j), imt, km)
	  call setbcx (flux_diffZ(1,1,j), imt, km)
        enddo
        call fct_diag (joff, js, jend, n, .true.)
      endif
# endif
      return
      end

#else
      return
      end
#endif
