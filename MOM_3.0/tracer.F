      subroutine tracer (joff, js, je, is, ie)
c
c=======================================================================
c     compute tracers at "tau+1" for rows js through je in the MW.
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
      use arrays_2d
      use memory_window
      use sbc_info
      use sbc_ocean
#if defined bbl_ag && defined sponges && !defined sponges_old
      use xsponge_mod
#endif
#include "param.h"
      parameter (istrt=2, iend=imt-1)
#ifdef bbl_ag
# include "bbl.h"
#endif
#include "accel.h"
#if defined obc_north
# include "cobc.h"
#endif
#include "coord.h"
#include "diag_region.h"
#include "emode.h"
#include "grdvar.h"
#include "hmixc.h"
#ifdef isoneutralmix
# include "isopyc.h"
#endif
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#ifdef bbl_ag
# include "state.h"
#endif
#include "vmixc.h"
      dimension twodt(km)
      integer tlev
#if defined save_density_terms || defined bbl_ag
# include "dens.h"
#endif
#include "fdift.h"
c
c-----------------------------------------------------------------------
c     bail out if starting row exceeds ending row
c-----------------------------------------------------------------------
c
      if (js .le. je) then
c
# ifdef timing
        call tic ('ocean', 'tracer')
        call tic ('tracer', 'misc')
# endif
c
        do k=1,km
	  twodt(k) = c2dtts*dtxcel(k)
        enddo
c
c-----------------------------------------------------------------------
c     limit the longitude indices based on those from the argument list
c     Note: this is currently bypassed. istrt and iend are set as
c           parameters to optimize performance
c-----------------------------------------------------------------------
c
c      istrt = max(2,is)
c      iend  = min(imt-1,ie)
c
c-----------------------------------------------------------------------
c     solve for one tracer at a time
c     n = 1 => temperature
c     n = 2 => salinity
c     n > 2 => passive tracers (if applicable)
c-----------------------------------------------------------------------
c
        do n=1,nt
c
#ifdef linearized_advection
c
c-----------------------------------------------------------------------
c       calculate 2*linear advective tracer flux
c-----------------------------------------------------------------------
c
          call adv_flux_lin (joff, js, je, is, ie, n)
#endif
#ifdef second_order_tracer_advection
c
c-----------------------------------------------------------------------
c       calculate 2* 2nd order advective tracer flux
c-----------------------------------------------------------------------
c
          call adv_flux_2nd (joff, js, je, is, ie, n)
#endif
#ifdef tracer_advection_fct
c
c-----------------------------------------------------------------------
c       calculate 2* FCT tracer flux
c-----------------------------------------------------------------------
c
          call adv_flux_fct(joff, js, je, is, ie, n)
#endif
#ifdef fourth_order_tracer_advection
c
c-----------------------------------------------------------------------
c       calculate 2* 4th order advective tracer flux
c-----------------------------------------------------------------------
c
          call adv_flux_4th (joff, js, je, is, ie, n)
#endif
#ifdef quicker
c
c-----------------------------------------------------------------------
c       calculate 2* 3rd order advective tracer flux
c-----------------------------------------------------------------------
c
          call adv_flux_quick (joff, js, je, is, ie, n)
#endif
#ifdef bottom_upwind
c
c-----------------------------------------------------------------------
c       calculate 2* upwind advective tracer flux at ocean bottom
c-----------------------------------------------------------------------
c
          call adv_flux_bottom_upwind (joff, js, je, is, ie, n)
#endif
c
c-----------------------------------------------------------------------
c       calculate diffusive flux across eastern and northern faces 
c       of T cells due to various parameterizations for diffusion.
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c       Horizontal tracer flux across the eastern and northern faces
c       of T cells. If partial cells enabled, then produce thickness
c       weighted fluxes. 
c-----------------------------------------------------------------------
c
c         east face
c
          do j=js,je
            do k=1,kbot
              do i=istrt-1,iend              
#if defined tracer_horz_laplacian
                diff_fe(i,k,j) = (T_i(i,k,j,n,1)
     &                           -T_i(i,k,j,n,0))
#endif
#if defined tracer_horz_biharmonic
                diff_fe(i,k,j) = (del2_tracer(i+1,k,j,n)
     &                           -del2_tracer(i,k,j,n))
#endif
#ifdef tracer_horz_mix_const
     &                           *ah_cstdxur(i,j)
#endif
#ifdef tracer_horz_mix_var
     &                           *diff_cet(i,k,j)*cstdxur(i,j)
#endif
#ifdef bryan_lewis_horizontal
     &                           *ahh(k)*cstdxur(i,j)
#endif
#if defined partial_cell
     &                           *dhte(i,k,j)
#endif
              enddo
            enddo
          enddo
c
c        north face 
c
          do j=js-1,je
            jrow = j + joff
	    jp1 = min(jrow+1,jmt-1)-joff
#ifdef tracer_horz_mix_const
            ah_csu_dyur =  diff_cnt*csu_dyur(jrow)
#endif
            do k=1,kbot
              do i=istrt,iend
#if defined tracer_horz_laplacian
                diff_fn(i,k,j) = (T_j(i,k,j,n,1)
     &                           -T_j(i,k,j,n,0))
#endif
#if defined tracer_horz_biharmonic
                diff_fn(i,k,j) = (del2_tracer(i,k,jp1,n)
     &                           -del2_tracer(i,k,j,n))
#endif
#ifdef tracer_horz_mix_const
     &                           *ah_csu_dyur
#endif
#ifdef tracer_horz_mix_var
     &                           *diff_cnt(i,k,j)*csu_dyur(jrow)
#endif
#ifdef bryan_lewis_horizontal
     &                           *ahh(k)*csu_dyur(jrow)
#endif
#if defined partial_cell
     &                           *dhtn(i,k,j)
#endif
              enddo
            enddo
#ifdef bbl_ag
            do i=istrt,iend
              diff_fn(i,km,j)=c0
            enddo
#endif            
          enddo
c 
#if defined obc_south || defined obc_north
          do j=js-1,je
# if defined obc_south
            if (j+joff .eq. 2) then
              do k=1,km
                do i=istrt,iend
                  diff_fn(i,k,j-1) = c0
                enddo
              enddo
            endif
# endif
# if defined obc_north
            if (j+joff .eq. jmtm1) then
              do k=1,km
                do i=istrt,iend
                  diff_fn(i,k,j) = c0
                enddo
              enddo
            endif
          enddo
# endif
c
#endif
c
c-----------------------------------------------------------------------
c       calculate diffusive flux across bottom face of T cells
c-----------------------------------------------------------------------
c
          nt2 = min(n,nsmix)
          do j=js,je
            do k=1,kbot-1
              do i=istrt,iend
                diff_fb(i,k,j) = diff_cbt(i,k,j,nt2)*dz_wtr(i,k,j)*
     &                         (t(i,k,j,n,taum1) - t(i,k+1,j,n,taum1))
              enddo
            enddo
          enddo
c
#ifdef isoneutralmix
c
c-----------------------------------------------------------------------
c       compute isopycnal diffusive flux through east, north, and 
c       bottom faces of T cells.  If partial cells enabled, then 
c       produce thickness weighted fluxes. 
c-----------------------------------------------------------------------
c
          call isoflux (joff, js, je, is, ie, n)
#endif
c
#ifdef source_term
c
c-----------------------------------------------------------------------
c       set source term for T cells
c-----------------------------------------------------------------------
c
          do j=js,je
            do k=1,km
              do i=istrt,iend
	        source(i,k,j) = c0
              enddo
            enddo
          enddo
c
# if defined sponges
c
c-----------------------------------------------------------------------
c       add newtonian damping term for sponge layers to source
c-----------------------------------------------------------------------
c
          call sponge1 (joff, js, je, istrt, iend, n, t(1,1,1,1,taum1)
     &,               source)
# endif
c
# if defined shortwave && !defined simple_sbc
c
c-----------------------------------------------------------------------
c       incorporate short wave penetration into source
c-----------------------------------------------------------------------
c
          if (n .eq. 1) then
	    call swflux0 (joff, js, je, istrt, iend, source)
	  endif
# endif
c
# ifdef xlandmix
c
c-----------------------------------------------------------------------
c       add crossland mixing of tracers into source
c-----------------------------------------------------------------------
c
          call xland (joff, js, je, is, ie, n)
# endif
c
# ifdef kppvmix
c
c--------------------------------------------------------------------
c     add nonlocal transport term ghats to source
c--------------------------------------------------------------------
c
          nt2 = min(n,nsmix)
          do j=js,je
	    k=1
	    do i=istrt,iend
#  if !defined kpp_nonlocal_trans_no_match 
	      source(i,k,j) = source(i,k,j) + stf(i,j,n) * dztr(k)
     &                      * (- diff_cbt(i,k,j,nt2) * ghats(i,k,j))
#  else
	      source(i,k,j) = source(i,k,j) + stf(i,j,n) * dztr(k)
     &                      * (-                       ghats(i,k,j))
#  endif
	    enddo
            do k=2,kbot-1
	      do i=istrt,iend
#  if !defined kpp_nonlocal_trans_no_match
	        source(i,k,j) = source(i,k,j) + stf(i,j,n) * dztr(k)
     &                       * (diff_cbt(i,k-1,j,nt2) * ghats(i,k-1,j)
     &                        - diff_cbt(i,k  ,j,nt2) * ghats(i,k  ,j))
#  else
	        source(i,k,j) = source(i,k,j) + stf(i,j,n) * dztr(k)
     &                       * (                        ghats(i,k-1,j)
     &                        -                         ghats(i,k  ,j))
#  endif
	      enddo
            enddo
          enddo
# endif
          
#endif source_term
c
#if defined bbl_ag
c
c         Entrainment at top face of bbl (at k=km)
c
          if (n .eq. 1) then
            do j=js-1,je
              jrow = j+joff
	      do i=istrt-1,iend
	        kb        = max(kmt(i,jrow),1)
	        rhoz(i,j) = rho(i,kb,j) - rho(i,km,j)
                uz2(i,j)  = (u(i,km,j,1,tau) - u(i,kb,j,1,tau))**2
     &                   +(u(i,km,j,2,tau) - u(i,kb,j,2,tau))**2
	      enddo
            enddo
c
            do j=js,je
	      do i=istrt,iend
	        Kbbl(i,j) = entrain_kbbl_constant
	        do i1=-1,0
	          do j1=-1,0
       	            Rit = -grav*rhoz(i,j)*dzt(km)/(uz2(i+i1,j+j1)+epsln)
	            Kbbl(i,j) = Kbbl(i,j)+0.25*entrain_kbbl/(1.0 + Rit)
                  enddo
                enddo
              enddo
            enddo
          endif
c
c         entrainment across top face of BBL (at k=km)
c
	  do j=js,je
	    jrow = j+joff
	    do i=istrt,iend
	      kb = max(kmt(i,jrow),1)
              entrain(i,j) =(Kbbl(i,j)*(t(i,kb,j,n,tau)-t(i,km,j,n,tau))
     &                    *dztr(km))*tmask(i,km,j)
            enddo
          enddo
#endif
c
c-----------------------------------------------------------------------
c       set surface and bottom vert b.c. on T cells for diffusion
c       and advection. for isopycnal diffusion, set adiabatic boundary
c       conditions.
c       note: the b.c. at adv_fb(i,k=bottom,j) is set by the above code.
c             However, it is not set when k=km so it is set below.
c             adv_fb(i,km,j) is always zero (to within roundoff).
c-----------------------------------------------------------------------
c
          do j=js,je
            jrow   = j + joff
            do i=istrt,iend
	      kb              = kmt(i,jrow)
c
c             vertical diffusion fluxes              
c              
              diff_fb(i,0,j)  = stf(i,j,n)
              diff_fb(i,kb,j) = btf(i,j,n)
#if defined isoneutralmix
              diff_fbiso(i,kb,j) = c0
#endif
c
c             vertical advection fluxes              
c              
	      adv_fb(i,0,j)   = adv_vbt(i,0,j)*(t(i,1,j,n,tau) +
     &                                          t(i,1,j,n,tau)) 
#if defined explicit_fresh_water_flux
              if (n .eq. 2) then 
  	        adv_fb(i,0,j)   = (adv_vbt(i,0,j)+sfft(i,jrow))
     &                      *(t(i,1,j,n,tau) + t(i,1,j,n,tau))
              endif
#endif
	      adv_fb(i,km,j)  = adv_vbt(i,km,j)*(t(i,km,j,n,tau) +
     &                                          t(i,km,j,n,tau))  
c
c             bbl modifications
c              
#if defined bbl_ag
              if (kb .gt. 0
# ifdef bbl_mask
     &           .and.mask_t_bbl(i,jrow).eq.1
# endif
     &             ) then
                diff_fb(i,kb,j)   = entrain(i,j)
                diff_fb(i,km-1,j) = entrain(i,j)
	        adv_fb(i,kb,j)  = (adv_vbt(i,kb,j)+abs(adv_vbt(i,kb,j)))
     &            	        *t(i,km,j,n,taum1)+
     &                            (adv_vbt(i,kb,j)-abs(adv_vbt(i,kb,j)))
     &                          *t(i,kb,j,n,taum1)
                adv_fb(i,km-1,j) = adv_fb(i,kb,j)
	      endif
#endif
            enddo
          enddo
c
c         
c-----------------------------------------------------------------------
c       solve for "tau+1" tracer using statement functions to represent
c       each component of the calculation 
c-----------------------------------------------------------------------
c
#if defined obc_west || defined obc_east || defined obc_north
# include "tracer_obc.inc"
#else
c
c       1st: solve using all components which are treated explicitly
c
          do j=js,je
            jrow   = j + joff
            do k=1,kbot
              do i=istrt,iend
	        t(i,k,j,n,taup1) = t(i,k,j,n,taum1) + twodt(k)*(
     &              DIFF_Tx(i,k,j) + DIFF_Ty(i,k,j) + DIFF_Tz(i,k,j)
     &             - ADV_Tx(i,k,j) -  ADV_Ty(i,k,j) -  ADV_Tz(i,k,j)
# ifdef source_term
     &             + source(i,k,j)
# endif
     &                            )*tmask(i,k,j)
              enddo
            enddo
          enddo
c
# if defined bbl_ag
c
c-----------------------------------------------------------------------
c         solve for "tau+1" tracer in bbl (k=km)
c-----------------------------------------------------------------------
c
c         horizontal advective fluxes
c
          do j=js,je
            jrow = j+joff
	    do i=istrt-1,iend
	      adv_fe(i,km,j) =
     &	                 ((adv_vet(i,km,j)+abs(adv_vet(i,km,j)))
     &                       *t(i,km,j,n,taum1)
     &                     +(adv_vet(i,km,j)-abs(adv_vet(i,km,j)))
     &                       *t(i+1,km,j,n,taum1))
     &                      *tmask(i+1,km,j)*tmask(i,km,j)
#  ifdef detrain
     &                      *oppsignx(i,j)
#  endif
	    enddo
	    call setbcx (adv_fe(1,km,j), imt, 1)
          enddo
#  ifdef detrain
          do j=js,je
            jrow = j+joff
	    do i=istrt-1,iend
	      kb     = kmt(i,jrow)
	      kb1    = kmt(i+1,jrow)
	      deeper = (1.0+sign(1,kb1-kb))*p5
	      kbhor  = min(kb,kb1)
              adv_fdete(i,j,n)=
     &	                   ((adv_vet(i,km,j)+abs(adv_vet(i,km,j)))
     &                       *t(i,km,j,n,taum1)
     &                     +(adv_vet(i,km,j)-abs(adv_vet(i,km,j)))
     &                       *t(i+1,kbhor+1,j,n,taum1))
     &                      *tmask(i+1,km,j)*tmask(i,km,j)*
     &                      (1.0-oppsignx(i,j))*deeper
              adv_fdetw(i+1,j,n)=
     &	                   ((adv_vet(i,km,j)+abs(adv_vet(i,km,j)))
     &                       *t(i,kbhor+1,j,n,taum1)
     &                     +(adv_vet(i,km,j)-abs(adv_vet(i,km,j)))
     &                       *t(i+1,km,j,n,taum1))
     &                      *tmask(i+1,km,j)*tmask(i,km,j)*
     &                      (1.0-oppsignx(i,j))*(1.0-deeper)
	    enddo
	    call setbcx (adv_fdete(1,j,n), imt, 1)
	    call setbcx (adv_fdetw(1,j,n), imt, 1)
          enddo
#  endif
c
	  do j=js-1,je
	    jrow = j+joff
	    do i=istrt,iend
	      kb     = kmt(i,jrow)
	      kb1    = kmt(i,jrow+1)
	      deeper = (1.0+sign(1,kb1-kb))*p5
	      kbhor  = min(kb,kb1)
       	      adv_fn(i,km,j) =
     &        	       ((adv_vnt(i,km,j)+abs(adv_vnt(i,km,j)))
     &                       *t(i,km,j,n,taum1)
     &                       +(adv_vnt(i,km,j)-abs(adv_vnt(i,km,j)))
     &                       *t(i,km,j+1,n,taum1))
     &                      *tmask(i,km,j+1)*tmask(i,km,j)
#  ifdef detrain
     &                      *oppsigny(i,j)
#  endif
	    enddo
          enddo
#  ifdef detrain
	  do j=js-1,je 
	    jrow = j+joff
	    do i=istrt,iend
	      kb     = kmt(i,jrow)
	      kb1    = kmt(i,jrow+1)
	      deeper = (1.0+sign(1,kb1-kb))*p5
	      kbhor  = min(kb,kb1)
       	      adv_fdetn(i,j,n) =
     &			     ((adv_vnt(i,km,j)+abs(adv_vnt(i,km,j)))
     &                       *t(i,km,j,n,taum1)
     &                       +(adv_vnt(i,km,j)-abs(adv_vnt(i,km,j)))
     &                       *t(i,kbhor+1,j+1,n,taum1))
     &                      *tmask(i,km,j+1)*tmask(i,km,j)
     &                      *(1-oppsigny(i,j))*deeper
       	      adv_fdets(i,j+1,n) =
     &	                     ((adv_vnt(i,km,j)+abs(adv_vnt(i,km,j)))
     &                       *t(i,kbhor+1,j,n,taum1)
     &                       +(adv_vnt(i,km,j)-abs(adv_vnt(i,km,j)))
     &                       *t(i,km,j+1,n,taum1))
     &                      *tmask(i,km,j+1)*tmask(i,km,j)
     &                      *(1-oppsigny(i,j))*(1-deeper)
	    enddo
          enddo
#  endif
c
#  ifdef quicker
c         set adv_f4n in the bbl to be used for diagnostics
c         purposes (in gyre.F):
	  do j=js-1,je 
	    jrow = j+joff
	    do i=istrt,iend
              adv_f4n(i,km,j,n)=adv_fn(i,km,j)
	    enddo
          enddo
#  endif
c
c         update tracers
c 
          do j=js,je
	    jrow = j+joff
	    do i=istrt,iend
	      t(i,km,j,n,taup1) = (t(i,km,j,n,taum1) + c2dtts*(
     &                        -(adv_fe(i,km,j)-adv_fe(i-1,km,j)
#  ifdef detrain
     &                        +adv_fdete(i,j,n)-adv_fdetw(i,j,n)
#  endif
     &                            )*cstr(jrow)*dxt2r(i)
#  if defined partial_cell
     &                             /dht(i,km,j)
#  endif
     &                        -(adv_fn(i,km,j)-adv_fn(i,km,j-1)
#  ifdef detrain
     &                        +adv_fdetn(i,j,n)-adv_fdets(i,j,n)
#  endif
     &                           )*cstdyt2r(jrow)
#  if defined partial_cell
     &                             /dht(i,km,j)
#  endif
     &                      -(adv_fb(i,km-1,j))*dzt2r(km)
     &                      + entrain(i,j)*dztr(km)
     &                         ))*tmask(i,km,j)
	    enddo
	    call setbcx (t(1,km,j,n,taup1), imt, 1)
          enddo
c
#  if defined detrain
c
c---------------------------------------------------------------------
c     put detrainment from bottom boundary layer in proper location
c---------------------------------------------------------------------
c
          do j=js,je
	    jrow = j+joff
	    do i=istrt,iend
	      kb = kmt(i,jrow)
	      kbe = kmt(i+1,jrow)
	      kbw = kmt(i-1,jrow)
	      kbs = kmt(i,jrow-1)
	      kbn = kmt(i,jrow+1)
	      kbhore = min(kb,kbe)
	      kbhorw = min(kb,kbw)
	      kbhorn = min(kb,kbn)
	      kbhors = min(kb,kbs)
	      t(i,kbhore+1,j,n,taup1) = t(i,kbhore+1,j,n,taup1)-
     &          twodt(kbhore+1)*adv_fdetw(i+1,j,n)*dxt2r(i)*cstr(jrow)
#   if defined partial_cell
     &         /dht(i,kbhore+1,j)
#   else
     &         *dzt(km)/dzt(kbhore+1)
#   endif
	      t(i,kbhorw+1,j,n,taup1) = t(i,kbhorw+1,j,n,taup1)+
     &          twodt(kbhorw+1)*adv_fdete(i-1,j,n)*dxt2r(i)*cstr(jrow)
#   if defined partial_cell
     &         /dht(i,kbhorw+1,j)
#   else
     &         *dzt(km)/dzt(kbhorw+1)
#   endif
	      t(i,kbhorn+1,j,n,taup1) = t(i,kbhorn+1,j,n,taup1)-
     &          twodt(kbhorn+1)*adv_fdets(i,j+1,n)*cstdyt2r(jrow)
#   if defined partial_cell
     &         /dht(i,kbhorn+1,j)
#   else
     &         *dzt(km)/dzt(kbhorn+1)
#   endif
	      t(i,kbhors+1,j,n,taup1) = t(i,kbhors+1,j,n,taup1)+
     &          twodt(kbhors+1)*adv_fdetn(i,j-1,n)*cstdyt2r(jrow)
#   if defined partial_cell
     &         /dht(i,kbhors+1,j)
#   else
     &         *dzt(km)/dzt(kbhors+1)
#   endif
	    enddo
          enddo
#  endif
# endif
c
# ifdef save_density_terms
c
c       Save contribution of advection to the potential density tendency
c
          if (rhots .and. .not. euler2 .and. n .le. 2) then
c
            if (n .eq. 1) then
              do j=js,je
                do k=1,km
                  do i=1,imt
                    advrho(i,k,j) = c0
                  enddo
                enddo
              enddo
            endif
c
            do j=js,je
              do k=1,km
                do i=istrt,iend
                  advrho(i,k,j) = advrho(i,k,j)
     &               -drho(i,k,j,n)*tmask(i,k,j)
     &               *(ADV_Tx(i,k,j)+ADV_Ty(i,k,j)+ADV_Tz(i,k,j))
                enddo
              enddo
            enddo
c
          endif
c
# endif


# if defined implicitvmix || defined redi_diffusion
c
c       2nd: add in portion of vertical diffusion handled implicitly  
c       
          call ivdift (joff, js, je, istrt, iend, n, twodt)
# endif
#endif
c
          do j=js,je
            call setbcx (t(1,1,j,n,taup1), imt, km)
	  enddo
c
c-----------------------------------------------------------------------
c       construct diagnostics associated with tracer "n" 
c-----------------------------------------------------------------------
c
          call diagt1 (joff, js, je, istrt, iend, n, twodt)
c
c-----------------------------------------------------------------------
c       end of tracer component "n" loop
c-----------------------------------------------------------------------
c
        enddo
c
#ifdef bbl_ag
c
c     Convectively adjust BBL with level above (1 pass only)
c
        do j=js,je
	  jrow = j+joff
# if defined sponges && !defined sponges_old
          if (spng(jrow) .ne. 0.0) cycle
# endif
	  do i=istrt,iend
	    kb = kmt(i,jrow)
	    if (kb .gt. 0
# ifdef bbl_mask
     &           .and.mask_t_bbl(i,jrow).eq.1
# endif
     &           ) then
	      tbot = t(i,km,j,1,taup1)-to(kb)
	      tbotm1 = t(i,kb,j,1,taup1)-to(kb)
	      sbot = t(i,km,j,2,taup1)-so(kb)
	      sbotm1 = t(i,kb,j,2,taup1)-so(kb)
	      dbot = dens(tbot,sbot,kb)
	      dbotm1 = dens(tbotm1,sbotm1,kb)
	      if (dbot .lt. dbotm1) then
# ifdef partial_cell
	        dztsum = (dzt(km)+dht(i,kb,j))
                d1     = dht(i,kb,j)
# else
	        dztsum = (dzt(km)+dzt(kb))
                d1     = dzt(kb)
# endif
	        do n=1,nt
                  tav = (t(i,km,j,n,taup1)*dzt(km)
     &	          +t(i,kb,j,n,taup1)*d1)/dztsum
 	          t(i,km,j,n,taup1) = tav
 	          t(i,kb,j,n,taup1) = tav
                enddo
              endif
            endif
          enddo
        enddo
#endif
#if (!defined implicitvmix || defined isoneutralmix) && !defined kppvmix
c
c-----------------------------------------------------------------------
c     explicit convection: adjust column if gravitationally unstable
c-----------------------------------------------------------------------
c
# if defined oldconvect
        call convct (t(1,1,1,1,taup1), ncon, joff, js, je, istrt, iend
     &,            kmt(1,js+joff))
# else
        call convct2 (t(1,1,1,1,taup1), joff, js, je, istrt, iend
     &,              kmt(1,js+joff))
# endif
#endif
c
c-----------------------------------------------------------------------
c     construct diagnostics after convection 
c-----------------------------------------------------------------------
c
        idiag = 10
        call diagt2 (joff, js, je, istrt, iend, idiag)
c
#if defined fourfil || defined firfil
c
c-----------------------------------------------------------------------
c     filter tracer tendency and reconstruct T(tau+1) at high latitudes
c-----------------------------------------------------------------------
c
        if (istrt .eq. 2 .and. iend .eq. imt-1) then
          call filt (joff, js, je)
        else
          write (stdout,'(a)') 
     &    'Error: filtering requires is=2 and ie=imt-1 in tracer'
          stop '=>tracer'
        endif
#endif
        do n=1,nt
          do j=js,je
            call setbcx (t(1,1,j,n,taup1), imt, km)
          enddo
        enddo
c
c-----------------------------------------------------------------------
c     construct diagnostics after filtering (for total dT/dt) 
c-----------------------------------------------------------------------
c
        idiag = 1
        call diagt2 (joff, js, je, istrt, iend, idiag)
c
#if !defined simple_sbc && !defined robert_time_filter
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
          tlev = tau
          if (euler2) tlev = taum1
        call asbct (joff, js, je, istrt, iend, isst, isss
     &,               t(1,1,1,1,tlev))
        endif
#endif
c
#if defined debug_tracer
        do j=js,je
	  jrow = j+joff
          write(stdout,*) 'In tracer: j=',j,' jrow=',jrow
        call print_checksum (t(1,1,j,1,taup1), imt, km,'ck T(taup1) = ')
        call print_checksum (t(1,1,j,2,taup1), imt, km,'ck S(taup1) = ')
        enddo
#endif
#if defined trace_indices && !defined obc_west && !defined obc_east
        write (stdout,'(2x,5(a,i4))')
     &   "=> In tracer: js=",js," je=",je," joff=",joff
     &  ,", calculating latitude rows: ",js+joff," to ",je+joff
#endif
#ifdef timing
        call toc ('tracer', 'misc')
        call toc ('ocean', 'tracer')
#endif

      endif

      return
      end



      subroutine diagt1 (joff, js, je, is, ie, n, twodt)
c
c-----------------------------------------------------------------------
c     construct diagnostics associated with tracer component "n" 
c
c     input:
c       joff  = offset relating "j" in the MW to latitude "jrow"
c       js    = starting row in the MW
c       je    = ending row in the MW
c       is    = starting longitude index in the MW
c       ie    = ending longitude index in the MW
c       n     = (1,2) = (u,v) velocity component
c       twodt = (2*dtts,dtts) on (leapfrog,mixing) time steps
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      use mpp_io_mod
      use arrays_2d
      use memory_window
#include "param.h"
#include "accel.h"
#include "coord.h"
#include "diag.h"
#include "diag_region.h"
#include "emode.h"
#include "grdvar.h"
#include "hmixc.h"
#ifdef isoneutralmix
# include "isopyc.h"
#endif
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "vmixc.h"
#ifdef meridional_tracer_budget
      dimension stor(imt,km), div(imt,km), sorc(imt,km), dif(imt,km)
      dimension t1(imt), t2(imt), t3(imt), t4(imt)
#endif
#ifdef time_step_monitor
      dimension temp1(imt,km), temp2(imt,km), temp3(imt,km)
#endif
      dimension twodt(km)
#include "fdift.h"
c
#ifdef timing
      if (n .eq. 1) then
	call tic ('tracer', 'diagnostics')
      else
	call ticr ('tracer', 'diagnostics')
      endif
#endif
c
#if defined save_mixing_coeff
c
c-----------------------------------------------------------------------
c     diagnostic: estimate mixing coefficients on east, north, and
c                 bottom face of T cells from the flux
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c               S. M. Griffies         e-mail smg@gfdl.gov
c-----------------------------------------------------------------------
c
      if (cmixts .and. n .eq. 1 .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'mixing coeff')
# endif
c
c       Lateral diffusivities
c 
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=2,imt-1
# if defined tracer_horz_mix_const
              ce(i,k,j,2) = diff_cet*tmask(i,k,j)*tmask(i+1,k,j)
              cn(i,k,j,2) = diff_cnt*tmask(i,k,j)*tmask(i,k,j+1)
# endif
# if defined tracer_horz_mix_var
              ce(i,k,j,2) = diff_cet(i,k,j)*tmask(i,k,j)*tmask(i+1,k,j)
              cn(i,k,j,2) = diff_cnt(i,k,j)*tmask(i,k,j)*tmask(i,k,j+1)
# endif
# if defined bryan_lewis_horizontal
              ce(i,k,j,2) = ahh(k)*tmask(i,k,j)*tmask(i+1,k,j)
              cn(i,k,j,2) = ahh(k)*tmask(i,k,j)*tmask(i,k,j+1)
# endif
# if defined isoneutralmix
              ce(i,k,j,2) = tmask(i+1,k,j)*agm(i,jrow)
              cn(i,k,j,2) = tmask(i,k,j+1)*agm(i,jrow)
# endif
            enddo
          enddo
        enddo
c
c       Vertical diffusivities  
c
        do j=js,je
          jrow = j + joff
          do k=1,km-1
            do i=2,imt-1
              cb(i,k,j,2)=diff_cbt(i,k,j,1)*tmask(i,k,j)*tmask(i,k+1,j)
# if defined redi_diffusion
              cb(i,k,j,2)=tmask(i,k+1,j)*(cb(i,k,j,2) - K33(i,k,j))
# endif
            enddo
          enddo
          do i=2,imt-1
            cb(i,km,j,2) = 0.0
          enddo
        enddo
c
        do j=js,je
          call setbcx (ce(1,1,j,2), imt, km)
          call setbcx (cn(1,1,j,2), imt, km)
          call setbcx (cb(1,1,j,2), imt, km)
        enddo
# ifdef timing
        call toc ('diagnostic', 'mixing coeff')
# endif
      endif
#endif
c
# if defined save_convection
c
c-----------------------------------------------------------------------
c     diagnostic: initialize temperature before convection
c
c     author:   S.M. Griffies       e-mail smg@gfdl.gov
c-----------------------------------------------------------------------
c
      if (exconvts .and. n .eq. 1 .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'explicit convection')
# endif
        do j=js,je
          do k=1,km
            do i=1,imt
              excnv0(i,k,j) = t(i,k,j,1,taup1)
            enddo
          enddo
        enddo
# ifdef timing
        call toc ('diagnostic', 'explicit convection')
# endif
      endif
c
#endif
c
#ifdef time_step_monitor
c
c-----------------------------------------------------------------------
c     diagnostic: integrate |d(tracer)/dt|  and tracer variance on "tau"
c                 globally
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c               (based on diagnostic by M. Cox)
c-----------------------------------------------------------------------
c
      if (tsits .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'time step monitor')
# endif

        do j=js,je
          jrow = j + joff
          r2dt    = c1/c2dtts
          cosdyt  = cst(jrow)*dyt(jrow)
          do k=1,km
            fx = r2dt/dtxcel(k)
            do i=is,ie
              dvolume = dxt(i)*cosdyt*tmask(i,k,j)  
#  ifdef partial_cell
     &                    *dht(i,k,j)
#  else
     &                    *dzt(k)
#  endif
c
c             Compute terms used for global tracer diagnostics 
c
#  if defined explicit_free_surface_nonlinear
              temp3(i,k) = t(i,k,j,n,taup1)*dvolume
#  else
              temp3(i,k) = t(i,k,j,n,tau)*dvolume
#  endif      
              temp1(i,k) = t(i,k,j,n,tau)**2*dvolume
              temp2(i,k) = abs(t(i,k,j,n,taup1)-t(i,k,j,n,taum1))*
     &                     dvolume*fx
            enddo
            do i=is,ie
              tbar(k,n,jrow)   = tbar(k,n,jrow) + temp3(i,k)
              travar(k,n,jrow) = travar(k,n,jrow) + temp1(i,k)
              dtabs(k,n,jrow)  = dtabs(k,n,jrow) + temp2(i,k)
            enddo
          enddo
        enddo

# if defined explicit_free_surface_nonlinear
c
c       save the surface tracer values for getting the
c       surface contribution to tbar inside diag_fs
c
        do j=js,je
          jrow=j+joff
          do i=is,ie
            surftracer(i,jrow,n,1) = t(i,1,j,n,taup1)
            surftracer(i,jrow,n,2) = t(i,1,j,n,tau)
          enddo
        enddo
# endif

# ifdef timing
        call toc ('diagnostic', 'time step monitor')
# endif
      endif
#endif
c
#ifdef tracer_averages
c
c-----------------------------------------------------------------------
c     diagnostic: accumulate tracers for averages under horizontal
c                 regions (use units of meters, rather than cm) 
c
c     author:   K. Dixon       e-mail kd@gfdl.gov
c-----------------------------------------------------------------------
c
      if (tavgts .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'tracer averages')
# endif
        do j=js,je
          jrow = j + joff
          do i=is,ie
            mask = mskhr(i,jrow)
            if (mask .ne. 0) then
              boxar = cst(jrow)*dxt(i)*dyt(jrow)*tmask(i,1,j)*0.0001
              sumbf(mask,n) = sumbf(mask,n) + stf(i,j,n)*boxar
              do k=1,km
                sumbk(mask,k,n) = sumbk(mask,k,n) + t(i,k,j,n,tau)
     &                             *boxar*tmask(i,k,j)*0.01
# ifdef partial_cell
     &                    *dht(i,k,j)
# else
     &                    *dzt(k)
# endif
             enddo
            endif
          enddo
        enddo
# ifdef timing
        call toc ('diagnostic', 'tracer averages')
# endif
      endif
#endif
c
#ifdef tracer_yz
c
c----------------------------------------------------------------------
c     diagnostic: integrate tracer equations in longitude
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c----------------------------------------------------------------------
c
      if (tyzts .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'tracer_yz')
# endif
        do j=js,je
          jrow = j + joff
          do k=1,km
            rtwodt = c1/twodt(k)
            sumdx = c0
	    do m=1,5
	      tyz(jrow,k,n,m) = c0
	    enddo
            do i=is,ie
              delx  = dxt(i)*tmask(i,k,j)
	      sumdx = sumdx + delx
              tyz(jrow,k,n,1) = tyz(jrow,k,n,1)+ delx*t(i,k,j,n,tau)
              tyz(jrow,k,n,2) = tyz(jrow,k,n,2) + rtwodt*delx*
     &                            (t(i,k,j,n,taup1) - t(i,k,j,n,taum1))
              tyz(jrow,k,n,3) =  tyz(jrow,k,n,3) 
     &                           - delx*(ADV_Ty(i,k,j) + ADV_Tz(i,k,j))
              tyz(jrow,k,n,4) =  tyz(jrow,k,n,4) 
     &                         + delx*(DIFF_Ty(i,k,j) + DIFF_Tz(i,k,j))
# ifdef source_term
              tyz(jrow,k,n,5) = tyz(jrow,k,n,5) + delx*source(i,k,j)
# endif
            enddo
	    sumdxr = c1/(sumdx + epsln)
	    do m=1,5
	      tyz(jrow,k,n,m) = tyz(jrow,k,n,m)*sumdxr
	    enddo
          enddo
        enddo
	if (js+joff .eq. 2) then
	  do m=1,5
	    do k=1,km
	      tyz(1,k,n,m)   = c0
	      tyz(jmt,k,n,m) = c0
	    enddo
	  enddo
	endif
# ifdef timing
        call toc ('diagnostic', 'tracer_yz')
# endif
      endif
#endif
c
#ifdef meridional_tracer_budget
c
c----------------------------------------------------------------------
c     diagnostic: integrate equations in depth, lon, and time for 
c                 each basin and jrow. basin # 0 is used to catch
c                 values in land areas
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c----------------------------------------------------------------------
c
      if (tmbperts .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'meridional tracer budget')
# endif
        do j=js,je
          jrow = j + joff
          if (jrow .eq. 2 .and. n .eq. 1) numtmb = numtmb + 1
          cosdyt = cst(jrow)*dyt(jrow)
	  do k=1,km
	    do i=is,ie
              stor(i,k) = 0.0
              div(i,k)  = 0.0
              sorc(i,k) = 0.0
              dif(i,k)  = 0.0
            enddo
	  enddo
          do k=1,km
            rtwodt = c1/twodt(k)
            do i=is,ie
              dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
              dxdydz          = dxdy
# ifdef partial_cell
     &                          *dht(i,k,j)
# else
     &                          *dzt(k)
# endif
              stor(i,k) = stor(i,k) + rtwodt*dxdydz*
     &                            (t(i,k,j,n,taup1) - t(i,k,j,n,taum1))
              div(i,k)  = div(i,k)  - dxdydz*ADV_Ty(i,k,j)
# ifdef source_term
              sorc(i,k) = sorc(i,k) + dxdydz*source(i,k,j)
# endif
              dif(i,k)  = dif(i,k)  + dxdydz*DIFF_Ty(i,k,j)
            enddo
          enddo
c
c         the accumulation is done this way for issues of speed
c
          do i=is,ie
            t1(i) = stor(i,1)
            t2(i) = div(i,1)
            t3(i) = sorc(i,1)
            t4(i) = dif(i,1)
          enddo
c
	  do k=2,km
	    do i=is,ie
              t1(i) = t1(i) + stor(i,k)
              t2(i) = t2(i)  + div(i,k)
              t3(i) = t3(i) + sorc(i,k)
              t4(i) = t4(i)  + dif(i,k)
            enddo
	  enddo
c
	  do i=is,ie
            m               = msktmb(i,jrow)
            tstor(jrow,n,m) = tstor(jrow,n,m) + t1(i)
            tdiv(jrow,n,m)  = tdiv(jrow,n,m)  + t2(i)
            tsorc(jrow,n,m) = tsorc(jrow,n,m) + t3(i)
            tdif(jrow,n,m)  = tdif(jrow,n,m)  + t4(i)
          enddo
c
	  k = 1
	  do i=is,ie
            m               = msktmb(i,jrow)
            dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
            tflux(jrow,n,m) = tflux(jrow,n,m) + dxdy*stf(i,j,n)
          enddo
          if (n .eq. 1) then
	    do k=1,km
	      do i=is,ie
                m               = msktmb(i,jrow)
                dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
                dxdydz          = dxdy
# ifdef partial_cell
     &                            *dht(i,k,j)
# else
     &                            *dzt(k)
# endif
                smdvol(jrow,m)  = smdvol(jrow,m)  + dxdydz
              enddo
	    enddo
          endif
        enddo
# ifdef timing
        call toc ('diagnostic', 'meridional tracer budget')
# endif
      endif
#endif
c
#ifdef gyre_components
c
c-----------------------------------------------------------------------
c     diagnostic: compute the northward transport components of
c                 each tracer
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c               (based on diagnostic by K. Bryan)
c-----------------------------------------------------------------------
c
      if (gyrets .and. eots)  call gyre (joff, js, je, is, ie, n)
#endif
c
#ifdef term_balances
c
c-----------------------------------------------------------------------
c     diagnostic: integrate r.h.s. terms in the tracer equations  
c                 over specified regional volumes.
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (trmbts .and. eots)  call ttb1 (joff, js, je, is, ie, n)
#endif
c
#ifdef save_xbts
c
c-----------------------------------------------------------------------
c     diagnostic: accumulate r.h.s. terms in the tracer equations
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (xbtperts .and. eots) call txbt1 (joff, js, je, n)
#endif
c
#ifdef timing
      call toc ('tracer', 'diagnostics')
#endif
      return
      end




      subroutine diagt2 (joff, js, je, is, ie, idiag)
c
c-----------------------------------------------------------------------
c     construct d(tracer)/dt diagnostics 
c
c     input:
c       joff  = offset relating "j" in the MW to latitude "jrow"
c       js    = starting row in the MW
c       je    = ending row in the MW
c       is    = starting longitude index in the MW
c       ie    = ending longitude index in the MW
c       idiag = 1  => total tracer change
c       idiag = 10 => change of tracer due to filtering(also convection)
c-----------------------------------------------------------------------
c
      use mpp_io_mod
      use memory_window
      character*128 file_name, file_stamp
#include "param.h"
#include "coord.h"
#include "diag.h"
#include "iounit.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"
c
#ifdef timing
      call ticr ('tracer', 'diagnostics')
#endif
c
#if defined save_convection
c
c-----------------------------------------------------------------------
c     diagnostic: save tracer time change due to explicit convection
c     idiag = 10 signifies diagnostics immediately after convection
c
c     author:   S.M. Griffies       e-mail smg@gfdl.gov
c-----------------------------------------------------------------------
c
      if (exconvts .and. idiag .eq. 10 .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'explicit convection')
# endif
        rdt = c1/c2dtts
c
c       excnv1 = epsln over land cells and d(convect)/dt over ocean
c
        do j=js,je
          do k=1,km
            do i=1,imt
              excnv1(i,k,j) = tmask(i,k,j)*
     &         (t(i,k,j,1,taup1)-excnv0(i,k,j))*rdt
     &                       + (c1-tmask(i,k,j))*epsln
            enddo
          enddo
        enddo
        call conv_out (joff, js, je, .true.)

# ifdef timing
        call toc ('diagnostic', 'explicit convection')
# endif
      endif
#endif
c
#ifdef term_balances
c
c-----------------------------------------------------------------------
c     diagnostic: integrate d/dt(tracer) over specified regional volumes
c                  after convection and filtering
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (trmbts .and. eots) call ttb2 (joff, js, je, is, ie, idiag)
#endif
c
#ifdef save_xbts
c
c-----------------------------------------------------------------------
c     diagnostic: accumulate d/dt(tracer) after convection
c                 and filtering
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (xbtperts .and. eots) call txbt2 (joff, js, je, idiag)
#endif
#ifdef timing
      call toc ('tracer', 'diagnostics')
#endif
      return
      end



#if !defined simple_sbc
      subroutine asbct (joff, js, je, is, ie, isst, isss, tt)
c
c-----------------------------------------------------------------------
c     construct the Atmos S.B.C. (surface boundary conditions) SST
c     and possibly SSS
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       isst   = index for sst 
c       isss   = index for sss 
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      use arrays_2d
      use memory_window
      use sbc_info
      use sbc_ocean
# include "param.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
      dimension tt(imt,km,jmw,nt)
# ifdef timing
      call tic ('tracer', 'constructing SBC')
# endif
c
c     initialize the Atmos S.B.C. at the start of each ocean segment
c     (do not alter values in land)
c
      if (eots .and. osegs) then
        do j=js,je
          jrow   = j + joff
          if (isst .ne. 0) then
            do i=is,ie
# if defined coupled
              if (kmt(i,jrow) .ne. 0) sbcocn(i,jrow,isst) = c0
# else
              sbcocn(i,jrow,isst) = c0
# endif
            enddo
          endif
          if (isss .ne. 0) then
            do i=is,ie
# if defined coupled
              if (kmt(i,jrow) .ne. 0) sbcocn(i,jrow,isss) = c0
# else
              sbcocn(i,jrow,isss) = c0
# endif
            enddo
          endif
        enddo
      endif
c
c     accumulate surface tracers for the Atmos S.B.C. every time step
c
      if (eots) then
        do j=js,je
          jrow   = j + joff
	  if (isst .ne. 0) then
            do i=is,ie
              sbcocn(i,jrow,isst) = sbcocn(i,jrow,isst)+tt(i,1,j,1)
     &                              *tmask(i,1,j)
            enddo
          endif
          if (isss .ne. 0) then
            do i=is,ie
              sbcocn(i,jrow,isss) = sbcocn(i,jrow,isss)+tt(i,1,j,2)
     &                              *tmask(i,1,j)
            enddo
	  endif
        enddo
      endif
c
c     average the surface tracers for the Atmos S.B.C. at the end of
c     each ocean segment. (do not alter values in land)
c
      if (eots .and. osege) then
        rts = c1/ntspos
        do j=js,je
          jrow   = j + joff
          if (isst .ne. 0) then
            do i=is,ie
# if defined coupled
              if (kmt(i,jrow) .ne. 0) sbcocn(i,jrow,isst) = 
     &                                rts*sbcocn(i,jrow,isst)
# else
              sbcocn(i,jrow,isst) = rts*sbcocn(i,jrow,isst)
# endif
            enddo
          endif
          if (isss .ne. 0) then
            do i=is,ie
# if defined coupled
              if (kmt(i,jrow) .ne. 0) sbcocn(i,jrow,isss) = 
     &	                              rts*sbcocn(i,jrow,isss)
# else
              sbcocn(i,jrow,isss) = rts*sbcocn(i,jrow,isss)
# endif
            enddo
          endif
        enddo
      endif
# ifdef timing
      call toc ('tracer', 'constructing SBC')
# endif
      return
      end
#endif


#if defined implicitvmix || defined redi_diffusion
      subroutine ivdift (joff, js, je, is, ie, n, twodt)
c
c-----------------------------------------------------------------------
c     solve vertical diffusion of tracers implicitly
c
c     input:
c       joff  = offset relating "j" in the MW to latitude "jrow"
c       js    = starting row in the MW
c       je    = ending row in the MW
c       is    = starting longitude index in the MW
c       ie    = ending longitude index in the MW
c       n     = tracer component
c       twodt = (2*dtts, dtts) on (leapfrog, mixing) time steps
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      use arrays_2d
      use memory_window
# include "param.h"
# include "mw.h"
# include "switch.h"
# include "vmixc.h"
      dimension twodt(km)
      dimension dcb(imt,km,2:jmw-1)
c
c     store terms to compute implicit vertical mixing on
c     diagnostic time steps
c
# ifdef timing
      call tic ('tracer', 'implicit vert diff')
# endif
c
# ifdef save_xbts
      if (eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = t(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif
# else
#  ifdef term_balances
      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = t(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif
#  endif
# endif
c
      nt2 = min(n,nsmix)
      do j=js,je
        do k=1,km
          do i=1,imt
            dcb(i,k,j) = diff_cbt(i,k,j,nt2)
          enddo
        enddo
      enddo
# ifdef tcvmix
      call invtri (tp1(1,1,1,n), stf(1,1,n), btf(1,1,n), vdca(???)
     &,              twodt, kmt, tmask(1,1,1), is, ie, joff, js, je)
      if (imt .gt. 1) then
        write (stdout,*) 
     &     'Error: ivdif.F is not converted for ifvef tcvmix'
	call abort()
      endif
# else
      call invtri (t(1,1,1,n,taup1), stf(1,1,n), btf(1,1,n)
     &, dcb(1,1,2), twodt, kmt, tmask(1,1,1), is, ie
     &, joff, js, je, 'tracer')
# endif
c
c     compute residual implicit vertical mixing
c
# ifdef save_xbts
      if (eots) then
        do j=js,je
          do k=1,km
            rc2dt = c1/twodt(k)
            do i=is,ie
              zzi(i,k,j) = rc2dt*(t(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
	  enddo
	enddo
      endif
# else
#  ifdef term_balances
      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            rc2dt = c1/twodt(k)
            do i=is,ie
              zzi(i,k,j) = rc2dt*(t(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
	  enddo
	enddo
      endif
#  endif
# endif
# ifdef timing
      call toc ('tracer', 'implicit vert diff')
# endif
      return
      end
#endif





#ifdef sponges
      subroutine sponge1 (joff, js, je, is, ie, n, tm1, tsource)
c
c=======================================================================
c     newtonian damping variables for sponge regions adjacent to
c     artificial southern and northern boundaries in limited domain
c     basins. data must be prepared using the run_sponge script 
c     included in PREP_DATA.
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       n    = (1,2) = tracer component
c       tm1  = tracer at "tau-1"
c
c     output:
c       tsource = newtonian damping term
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use xsponge_mod
      use memory_window
# include "param.h"
# include "iounit.h"
# if defined equatorial_thermocline
#  include "coord.h"
# endif
# include "mw.h"
# include "switch.h"
# include "tmngr.h"
c
      dimension tsource(imt,km,2:jmw-1)
      dimension tm1(imt,km,jmw,nt)
c
# ifdef timing
        call tic ('tracer', 'sponge boundary')
# endif
c
# if defined equatorial_thermocline
c
c       only appropriate for idealized equatorial studies 
c
c       USER INPUT:
c       spng_width = width (deg) of sponges
c       spng_damp  = reciprocal of damping period (1/sec)
c
        spng_width = 3.0
        spng_damp  = 1.0/(5.0*86400.0)
c
        do j=js,je
          jrow = joff + j
          if ((yt(jrow) - yu(1)) .le. spng_width
     &      .or. (yu(jmt-1) - yt(jrow)) .le. spng_width) then
            do k=1,km
              do i=is,ie
	        tsource(i,k,j) = tsource(i,k,j) - 
     &                  spng_damp*(tm1(i,k,j,n) - tbarz(k,n))
              enddo
            enddo
          endif
        enddo
# else
c
      if ((joff + js .eq. 2) .and. n .eq. 1 .and. is .eq. 2) then
c
c-----------------------------------------------------------------------
c       decide whether to read sponge data or not
c-----------------------------------------------------------------------
c
        begtim = position_within_data (model_time, sp_start_time
     &,                                sp_end_time, .true.)
        method = 3
        call timeinterp (begtim, indxsp, tspng, spgdpm, 12, .true.
     &,           method, inextd, iprevd, wprev, readsp, inext, iprev)
c
c       read in the next data record from disk when needed
c
        if (readsp) then
# ifdef ssread_sswrite
          lrec = max_rows*2*imt*km
          call oget (ionew2, lrec, inextd, bufmb(1,1,1,1,inext))
# else
          call getunit (ionew2, 'sponges', opt_sponge) 
          read (ionew2, rec=inextd) bufmb1 
          do in=1,max_rows
            do k=1,km
              do i=1,imt
                bufmb(i,k,in,1,inext) = bufmb1(i,k,in,1)
                bufmb(i,k,in,2,inext) = bufmb1(i,k,in,2)
              enddo
           enddo
          enddo
# endif
          write (stdout,'(/a,i3,a,i2,a,i2,a,a/)')
     &    '=> read sponge record =',inextd,' into buffer =', inext
     &,   ' method #',method,' at ', stamp
          write (stdout,*)
     &   ' sponge cksum =', checksum (bufmb(1,1,1,1,inext)
     &,                               imt*km, max_rows*2)
# ifndef ssread_sswrite
          call relunit (ionew2)
# endif
        endif
      endif
      if (n .le. 2) then
c
c-----------------------------------------------------------------------
c       construct newtonian damping term using sponge data
c-----------------------------------------------------------------------
c
	do j=js,je
	  jrow = joff + j
	  if (spng(jrow) .ne. c0) then
	    j_sponge = sponge_row(jrow)
	    tnext = c1-wprev
	    do k=1,km
	      do i=is,ie
	        data1 = tnext*bufmb(i,k,j_sponge,n,inext) 
     &                + wprev*bufmb(i,k,j_sponge,n,iprev)
                tsource(i,k,j) = tsource(i,k,j) - 
     &                           spng(jrow)*(tm1(i,k,j,n) - data1)
              enddo
	    enddo
	  endif
        enddo
      endif
# endif
# ifdef timing
        call toc ('tracer', 'sponge boundary')
# endif
      return
      end
#endif

#if defined shortwave && !defined simple_sbc
      subroutine swflux0 (joff, js, je, is, ie, source)
c
c=======================================================================
c       incorporate short wave penetration via the "source" 
c       term. note that the divergence of shortwave for the first
c       level "divpen(1)" is compensating for the effect of having
c       the shortwave component already included in the total
c       surface tracer flux "stf(i,j,1)"
c
c       incorporating the penetrative shortwave radiative flux into
c       the vertical diffusive flux term directly is not correct when
c       vertical diffusion is solved implicitly. The tri-diagonal
c       solver takes diffusion coefficients at tops and bottoms of
c       cells but not the fluxes.
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c
c     output:
c       source = source term penetrative short wave
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c=======================================================================
c
      use sbc_info
      use sbc_ocean
#include "param.h"
#include "cshort.h"      
c
      dimension source(imt,km,2:jmw-1)
# ifdef timing
        call tic ('tracer', 'shortwave penetration')
# endif
c
      isw = index_of_sbc ('short wave')
      if (isw .ne. 0) then
        do j=js,je
          jrow   = j + joff
	  do k=1,km
	    do i=is,ie
              source(i,k,j) = source(i,k,j)
     &                       + sbcocn(i,jrow,isw)*divpen(k) 
            enddo
          enddo
        enddo
      endif
# ifdef timing
        call toc ('tracer', 'shortwave penetration')
# endif
      return
      end
#endif
