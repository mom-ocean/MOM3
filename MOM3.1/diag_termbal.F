      subroutine utb1 (joff, js, je, is, ie, n)
c
c=======================================================================
c     accumulate terms in the momentum equations over the
c     volume of the specified regions
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       n    = velocity component
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#ifdef term_balances
      use mpp_io_mod
      use arrays_2d
      use memory_window
# include "param.h"
# include "taskrows.h"
# include "coord.h"
# include "diag_region.h"
# include "diag.h"
# include "grdvar.h"
# include "hmixc.h"
# include "mw.h"
# include "scalar.h"
# include "vmixc.h"
# include "fdifm.h"
c
# ifdef timing
      call tic ('diagnostic', 'term balances')
# endif
c
c-----------------------------------------------------------------------
c     set local constants
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        fx = csu(jrow)*dyu(jrow)
# ifdef symmetry
        if (jrow .eq. jmtm1) fx = fx*p5
# endif
c
c-----------------------------------------------------------------------
c       accumulate terms for all regions within the current jrow 
c-----------------------------------------------------------------------
c
        do k=1,km
          do i=is,ie
            nreg = nhreg*(mskvr(k)-1) + mskhr(i,jrow)
            if (nreg .gt. 0 .and. mskhr(i,jrow) .gt. 0) then
              boxvol = fx*dxu(i)
# ifdef partial_cell
     &                          *dhu(i,k,j)
# else
     &                          *dzt(k)
# endif
c
c-----------------------------------------------------------------------
c             pressure term
c-----------------------------------------------------------------------
c
              term = -umask(i,k,j)*grad_p(i,k,j,n)
              call addto (termbm(k,2,n,nreg), term*boxvol)
c
c-----------------------------------------------------------------------
c             zonal advection (flux form) of momentum
c-----------------------------------------------------------------------
c
              term = -umask(i,k,j)*ADV_Ux(i,k,j)
              call addto (termbm(k,3,n,nreg), term*boxvol)
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
              call addto (termbm(k,14,n,nreg), term*boxvol)
c
c-----------------------------------------------------------------------
c             advective metric term
c-----------------------------------------------------------------------
c
              term = ADV_metric(i,k,j,n)
              call addto (termbm(k,13,n,nreg), term*boxvol)
c
c-----------------------------------------------------------------------
c             meridional advection (flux form) of momentum
c-----------------------------------------------------------------------
c
              term = -umask(i,k,j)*ADV_Uy(i,k,j)
              call addto (termbm(k,4,n,nreg), term*boxvol)
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
              call addto (termbm(k,15,n,nreg), term*boxvol)
c
c-----------------------------------------------------------------------
c             vertical advection (flux form) of momentum
c-----------------------------------------------------------------------
c
              term = -umask(i,k,j)*ADV_Uz(i,k,j)
              call addto (termbm(k,5,n,nreg), term*boxvol)
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
              call addto (termbm(k,16,n,nreg), term*boxvol)
c
c-----------------------------------------------------------------------
c             zonal diffusion of momentum
c-----------------------------------------------------------------------
c
              term = umask(i,k,j)*DIFF_Ux(i,k,j)
              call addto (termbm(k,6,n,nreg), term*boxvol)
c
c-----------------------------------------------------------------------
c             meridional diffusion of momentum
c-----------------------------------------------------------------------
c
              term = umask(i,k,j)*DIFF_Uy(i,k,j)
              call addto (termbm(k,7,n,nreg), term*boxvol)
c
c-----------------------------------------------------------------------
c             diffusive metric term
c-----------------------------------------------------------------------
c
              term = umask(i,k,j)*DIFF_metric(i,k,j,n)
              call addto (termbm(k,9,n,nreg), term*boxvol)
c
c-----------------------------------------------------------------------
c             vertical diffusion of momentum
c-----------------------------------------------------------------------
c
              term = umask(i,k,j)*DIFF_Uz(i,k,j)
# ifdef implicitvmix
     &              +umask(i,k,j)*zzi(i,k,j)
# endif
              call addto (termbm(k,8,n,nreg), term*boxvol)
c
c-----------------------------------------------------------------------
c             coriolis term
c-----------------------------------------------------------------------
c
              term = umask(i,k,j)*CORIOLIS(i,k,j,n)
              call addto (termbm(k,10,n,nreg), term*boxvol)
c
c-----------------------------------------------------------------------
c             accumulate the source term
c-----------------------------------------------------------------------
c
# ifdef source_term
              term = umask(i,k,j)*source(i,k,j)
              call addto (termbm(k,11,n,nreg), term*boxvol)
# else
              termbm(k,11,n,nreg) = 0.0
# endif
c
c-----------------------------------------------------------------------
c             accumulate u, v, and w
c-----------------------------------------------------------------------
c
              term = umask(i,k,j)*u(i,k,j,n,tau)
              call addto (termbm(k,17,n,nreg), term*boxvol)
c
              if (n .eq. 2) then
                term = p5*(adv_vbu(i,k,j)+adv_vbu(i,k-1,j))*umask(i,k,j)
	        call addto (avgw(nreg), term*boxvol)
              endif
c
c-----------------------------------------------------------------------
c             accumulate the surface momentum flux
c-----------------------------------------------------------------------
c
              if (k .eq. 1) then
                term = umask(i,k,j)*smf(i,j,n)
                call addto (smflx(n,nreg), term*fx*dxu(i))
	      endif
            endif
          enddo        
        enddo
      enddo
c
# ifdef timing
      call toc ('diagnostic', 'term balances')
# endif
c
      return
      end


      subroutine utb2 (joff, js, je, is, ie, c2dtuv, acor)
c
c=======================================================================
c     accumulate external mode parts of d/dt and the implicit coriolis
c     term in the momentum equations over the volume in the specified
c     regions
c
c     input:
c       joff   = offset relating "j" in the MW to latitude "jrow"
c       js     = starting row in the MW
c       je     = ending row in the MW
c       is     = starting longitude index in the MW
c       ie     = ending longitude index in the MW
c       c2dtuv = (2*dtuv,dtuv) on (lpfrod,mixing) time steps
c       acor   = implicit factor 
c
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use mpp_io_mod
      use arrays_2d
      use memory_window
# include "param.h"
# include "taskrows.h"
# include "coord.h"
# include "diag_region.h"
# include "diag.h"
# include "grdvar.h"
# include "mw.h"
# include "fdifm.h"
c
# ifdef timing
      call tic ('diagnostic', 'term balances')
# endif
c
c
c-----------------------------------------------------------------------
c       local constants
c-----------------------------------------------------------------------
c
      r2dt = c1/c2dtuv
c
      do j=js,je
        jrow = j + joff
        fx   = csu(jrow)*dyu(jrow)
# ifdef symmetry
        if (jrow .eq. jmtm1) fx = fx*p5
# endif
        do n=1,2
          do k=1,km
            do i=is,ie
              nreg   = nhreg*(mskvr(k)-1) + mskhr(i,jrow)
              if (nreg .gt. 0 .and. mskhr(i,jrow) .gt. 0) then
                boxvol = fx*dxu(i)
# ifdef partial_cell
     &                            *dhu(i,k,j)
# else
     &                            *dzt(k)
# endif
c
c-----------------------------------------------------------------------
c               d/dt of velocity (external mode part of tau+1 will be  
c               added later when the external mode is solved)
c-----------------------------------------------------------------------
c
                term = umask(i,k,j)*(u(i,k,j,n,taup1) - 
     &                               u(i,k,j,n,taum1))*r2dt
                call addto (termbm(k,1,n,nreg), term*boxvol)
c
c-----------------------------------------------------------------------
c               implicit coriolis term (external mode part will be added
c               later when external mode is solved)
c-----------------------------------------------------------------------
c
                if (acor .ne. c0) then
# ifdef rot_grid
                  term = umask(i,k,j)*acor*cori(i,jrow,n)*
# else
                  term = umask(i,k,j)*acor*cori(jrow,n)*
# endif
     &                       (u(i,k,j,3-n,taup1) - u(i,k,j,3-n,taum1))
                  call addto (termbm(k,10,n,nreg), term*boxvol)
	        endif
	      endif
            enddo        
          enddo        
        enddo
      enddo
c
# ifdef timing
      call toc ('diagnostic', 'term balances')
# endif
      return
      end



      subroutine utb3
c
c=======================================================================
c     accumulate external mode parts of d/dt, the implicit coriolis
c     term and the surface pressure gradient sover the volume in the 
c     specified regions.
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use mpp_io_mod
      use arrays_2d
      use memory_window
# include "param.h"
# include "coord.h"
# include "diag_region.h"
# include "grdvar.h"
# include "mw.h"
# include "scalar.h"
# include "diag_termbal.h"
      parameter (is=1, ie=1, js=1, je=1)
      dimension psgrad(is:ie,js:je,2)
c
# ifdef timing
      call tic ('diagnostic', 'term balances')
# endif
c
      do jrow=1,jmt-1
        fddt  = csu(jrow)*dyu(jrow)/c2dtuv
        fspr  = csu(jrow)*dyu(jrow)
# ifdef rot_grid
        do i=2,imtm1
          atosp = acor*cori(i,jrow,1)
          f1    = atosp*csu(jrow)*dyu(jrow)
# else
        atosp = acor*cori(jrow,1)
        f1    = atosp*csu(jrow)*dyu(jrow)
        do i=2,imtm1
# endif
          kz = kmu(i,jrow)
          if (kz .ne. 0) then
            do k=1,kz
              n = nhreg*(mskvr(k)-1) + mskhr(i,jrow)
              if (n .gt. 0 .and. mskhr(i,jrow) .gt. 0) then
c
c               construct the surface pressure gradients for pt (i,jrow)
c
		if (k .eq. 1) then
		  call calc_psgrad(psgrad, uext, vext, jrow, jrow, i, i)
                endif                                       
                boxfac = fddt*dxu(i)
# ifdef partial_cell
                j = jrow
                boxfac = fddt*dxu(i)*dhu(i,k,j)
# else
     &                              *dzt(k)
# endif
                boxspr = fspr*dxu(i)
# ifdef partial_cell
     &                              *dhu(i,k,j)
# else
     &                              *dzt(k)
# endif
                termbm(k,1,1,n)  = termbm(k,1,1,n)  + uext*boxfac
                termbm(k,1,2,n)  = termbm(k,1,2,n)  + vext*boxfac
                termbm(k,12,1,n) = termbm(k,12,1,n) - 
     &                              psgrad(is,js,1)*boxspr
                termbm(k,12,2,n) = termbm(k,12,2,n) - 
     &                              psgrad(is,js,2)*boxspr
                boxacr = f1*dxu(i)
# ifdef partial_cell
     &                            *dhu(i,k,j)
# else
     &                            *dzt(k)
# endif
                termbm(k,10,1,n) = termbm(k,10,1,n) + vext*boxacr
                termbm(k,10,2,n) = termbm(k,10,2,n) - uext*boxacr
              endif
	    enddo
          endif
	enddo
      enddo
c
# ifdef timing
      call toc ('diagnostic', 'term balances')
# endif
      return
      end




      subroutine ttb1 (joff, js, je, is, ie, n)
c
c=======================================================================
c     accumulate terms in the tracer equations over the volume in the
c     specified regions
c
c     input:
c       joff  = offset relating "j" in the MW to latitude "jrow"
c       js    = starting row in the MW
c       je    = ending row in the MW
c       is    = starting longitude index in the MW
c       ie    = ending longitude index in the MW
c       n     = tracer component
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use mpp_io_mod
      use arrays_2d
      use memory_window
# include "param.h"
# include "taskrows.h"
# include "accel.h"
# include "coord.h"
# include "diag_region.h"
# include "diag.h"
# include "grdvar.h"
# include "hmixc.h"
# ifdef isoneutralmix
#  include "isopyc.h"
# endif
# include "mw.h"
# include "scalar.h"
# include "vmixc.h"
# include "fdift.h"
c
# ifdef timing
      call tic ('diagnostic', 'term balances')
# endif
c
c-----------------------------------------------------------------------
c     limit the longitude indices
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
      do j=js,je
        jrow = j + joff
c
c-----------------------------------------------------------------------
c       set local constants
c-----------------------------------------------------------------------
c
        fx   = cst(jrow)*dyt(jrow)
        do k=1,km
          do i=istrt,iend
            nreg = nhreg*(mskvr(k)-1) + mskhr(i,jrow)
            if (nreg .gt. 0 .and. mskhr(i,jrow) .gt. 0) then
	      area   = fx*dxt(i)
              boxvol = area
# ifdef partial_cell
     &                     *dht(i,k,j)
# else
     &                     *dzt(k)
# endif
c
c-----------------------------------------------------------------------
c             tracer
c-----------------------------------------------------------------------
c
              term = tmask(i,k,j)*t(i,k,j,n,tau)
              call addto (termbt(k,15,n,nreg), term*boxvol)
              call addto (termbt(k,15,n,0),    term*boxvol)
c
c-----------------------------------------------------------------------
c             d(tracer)/dt
c-----------------------------------------------------------------------
c
	      r2dt = c1/(c2dtts*dtxcel(k))
              term = tmask(i,k,j)*(t(i,k,j,n,taup1) -
     &                             t(i,k,j,n,taum1))*r2dt
              call addto (termbt(k,9,n,nreg), term*boxvol)
              call addto (termbt(k,9,n,0),    term*boxvol)
c
c-----------------------------------------------------------------------
c             zonal advection (flux form) of tracer
c-----------------------------------------------------------------------
c
              term = -tmask(i,k,j)*ADV_Tx(i,k,j)
              call addto (termbt(k,2,n,nreg), term*boxvol)
              call addto (termbt(k,2,n,0),    term*boxvol)
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
              call addto (termbt(k,11,n,nreg), term*boxvol)
              call addto (termbt(k,11,n,0),    term*boxvol)
c
c-----------------------------------------------------------------------
c             meridional advection (flux form) of tracer
c-----------------------------------------------------------------------
c
              term = -tmask(i,k,j)*ADV_Ty(i,k,j)
              call addto (termbt(k,3,n,nreg), term*boxvol)
              call addto (termbt(k,3,n,0),    term*boxvol)
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
              call addto (termbt(k,12,n,nreg), term*boxvol)
              call addto (termbt(k,12,n,0),    term*boxvol)
c
c-----------------------------------------------------------------------
c             vertical advection (flux form) of tracer
c-----------------------------------------------------------------------
c
              term = -tmask(i,k,j)*ADV_Tz(i,k,j)
              call addto (termbt(k,4,n,nreg), term*boxvol)
              call addto (termbt(k,4,n,0),    term*boxvol)
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
              call addto (termbt(k,13,n,nreg), term*boxvol)
              call addto (termbt(k,13,n,0),    term*boxvol)
c
c-----------------------------------------------------------------------
c             zonal diffusion of tracer
c-----------------------------------------------------------------------
c
              term = tmask(i,k,j)*DIFF_Tx(i,k,j)
              call addto (termbt(k,5,n,nreg), term*boxvol)
              call addto (termbt(k,5,n,0),    term*boxvol)
c
c-----------------------------------------------------------------------
c             meridional diffusion of tracer
c-----------------------------------------------------------------------
c
              term = tmask(i,k,j)*DIFF_Ty(i,k,j)
              call addto (termbt(k,6,n,nreg), term*boxvol)
              call addto (termbt(k,6,n,0),    term*boxvol)
c
c-----------------------------------------------------------------------
c             vertical diffusion of tracer
c-----------------------------------------------------------------------
c
              term = tmask(i,k,j)*DIFF_Tz(i,k,j)
# if defined implicitvmix || defined redi_diffusion
     &               + tmask(i,k,j)*zzi(i,k,j)
# endif
              call addto (termbt(k,7,n,nreg), term*boxvol)
              call addto (termbt(k,7,n,0),    term*boxvol)
c
c-----------------------------------------------------------------------
c             tracer source term
c-----------------------------------------------------------------------
c
# ifdef source_term
              term = tmask(i,k,j)*source(i,k,j)
              call addto (termbt(k,8,n,nreg), term*boxvol)
              call addto (termbt(k,8,n,0),    term*boxvol)
# else
              termbt(k,8,n,nreg) = 0.0
              termbt(k,8,n,0)    = 0.0
# endif
c
              if (k .eq. 1) then
c
c-----------------------------------------------------------------------
c               surface tracer
c-----------------------------------------------------------------------
c
                term = tmask(i,k,j)*t(i,k,j,n,tau)
                call addto (asst(n,nreg), term*area)
                call addto (asst(n,0),    term*area)
c
c-----------------------------------------------------------------------
c               surface tracer flux
c-----------------------------------------------------------------------
c
                term = tmask(i,k,j)*stf(i,j,n)
                call addto (stflx(n,nreg), term*area)
                call addto (stflx(n,0), term*area)
              endif
            endif
          enddo
        enddo
      enddo
# ifdef timing
      call toc ('diagnostic', 'term balances')
# endif
      return
      end



      subroutine ttb2 (joff, js, je, is, ie, iterm)
c
c=======================================================================
c     accumulate d/dt and change in tracer in the tracer equations over
c     the volume in the specified regions
c
c     input:
c       joff  = offset relating "j" in the MW to latitude "jrow"
c       js    = starting row in the MW
c       je    = ending row in the MW
c       is    = starting longitude index in the MW
c       ie    = ending longitude index in the MW
c
c     iterm = 1  => total change
c     iterm = 10 => change from filtering      
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use mpp_io_mod
      use arrays_2d
      use memory_window
# include "param.h"
# include "taskrows.h"
# include "accel.h"
# include "coord.h"
# include "diag_region.h"
# include "diag.h"
# include "grdvar.h"
# include "mw.h"
# include "scalar.h"
c
# ifdef timing
      call tic ('diagnostic', 'term balances')
# endif
c
      if (iterm .ne. 1 .and. iterm .ne. 10) then
        write (stdout,*) '=>Error: iterm=',iterm,' in ttb2'
	call abort()
      endif
c
      do j=js,je
        jrow = j + joff
        fx   = cst(jrow)*dyt(jrow)
        do n=1,nt
          do k=1,km
	    r2dt = c1/(c2dtts*dtxcel(k))
            do i=is,ie
              nreg   = nhreg*(mskvr(k)-1) + mskhr(i,jrow)
              if (nreg .gt. 0 .and. mskhr(i,jrow) .gt. 0) then
	        area   = fx*dxt(i)
                boxvol = area
# ifdef partial_cell
     &                       *dht(i,k,j)
# else
     &                       *dzt(k)
# endif
c
c-----------------------------------------------------------------------
c               d/dt(tracer)
c-----------------------------------------------------------------------
c
                term = tmask(i,k,j)*(t(i,k,j,n,taup1) -
     &                               t(i,k,j,n,taum1))*r2dt
                call addto (termbt(k,iterm,n,nreg), term*boxvol)
                call addto (termbt(k,iterm,n,0),    term*boxvol)
c
c-----------------------------------------------------------------------
c               change in variance of tracer
c-----------------------------------------------------------------------
c
                if (iterm .eq. 1) then
                  term = tmask(i,k,j)*(t(i,k,j,n,taup1)**2-
     &                                 t(i,k,j,n,taum1)**2)
                  call addto (termbt(k,14,n,nreg), term*boxvol)
                  call addto (termbt(k,14,n,0),    term*boxvol)
                endif
              endif
	    enddo
	  enddo
        enddo
      enddo
c
# ifdef timing
      call toc ('diagnostic', 'term balances')
# endif
      return
      end

#else
      return
      end
#endif


