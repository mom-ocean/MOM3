      subroutine gyre (joff, js, je, is, ie, n)
c
c-----------------------------------------------------------------------
c     compute the northward transport components of each tracer
c
c     input:
c       joff  = offset relating "j" in the MW to latitude "jrow"
c       js    = starting row in the MW
c       je    = ending row in the MW
c       is    = starting longitude index in the MW
c       ie    = ending longitude index in the MW
c       n     = tracer component
c
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c     based on work by M.Cox and K. Bryan
c     modifications for basins: gokhan danabasoglu
c-----------------------------------------------------------------------
c
#ifdef gyre_components
      use mpp_io_mod
      use arrays_2d
      use memory_window
# include "param.h"
# include "coord.h"
# include "diag_region.h"
# include "diag.h"
# include "grdvar.h"
# include "hmixc.h"
# include "mw.h"
# include "scalar.h"
c
# ifdef timing
      call tic ('diagnostic', 'gyre components')
# endif
c
      do j=js,je
        jrow = j + joff
        if (jrow .lt. jmtm1) then
	  small = 1.e-10
          do k=1,km
            totdxn = small
            totdxs = small
            vbr    = c0
            tbrs   = c0
            tbrn   = c0
            do i=is,ie
              totdxn = totdxn + dxt(i)*tmask(i,k,j+1)
# ifdef partial_cell
     &                 *dht(i,k,j+1)
# endif
              totdxs = totdxs + dxt(i)*tmask(i,k,j)
# ifdef partial_cell
     &                 *dht(i,k,j)
# endif
              vbr    = vbr  + u(i,k,j,2,tau)*dxu(i)*csu(jrow)
# ifdef partial_cell
     &                 *dhu(i,k,j)
# endif
              tbrn   = tbrn + t(i,k,j+1,n,tau)*tmask(i,k,j+1)*dxt(i)
# ifdef partial_cell
     &                 *dht(i,k,j+1)
# endif
              tbrs   = tbrs + t(i,k,j,n,tau)*tmask(i,k,j)*dxt(i)
# ifdef partial_cell
     &                 *dht(i,k,j)
# endif
	    enddo
            tbrn          = tbrn/totdxn
            tbrs          = tbrs/totdxs
            ttn(1,jrow,n) = ttn(1,jrow,n) + vbr*p5*(tbrn+tbrs)
# if !defined partial_cell
     &                                                       *dzt(k)
# endif
            do i=is,ie
	      tempdiff_fn = diff_fn(i,k,j)*
     &                           tmask(i,k,j+1)*tmask(i,k,j)*dxt(i)
# if !defined partial_cell
     &                                 *dzt(k)
# endif
# if defined fourth_order_tracer_advection || defined quicker
	      tempadv_fn = p5*adv_f4n(i,k,j,n)*dxt(i)
# else
	      tempadv_fn = p5*adv_fn(i,k,j)*dxt(i)
# endif
# if !defined partial_cell
     &                                     *dzt(k)
# endif
              ttn(6,jrow,n)    = ttn(6,jrow,n) + tempadv_fn
              ttn(7,jrow,n)    = ttn(7,jrow,n) - tempdiff_fn
              ttn2(6,jrow,n,0) = ttn2(6,jrow,n,0) + tempadv_fn
              ttn2(7,jrow,n,0) = ttn2(7,jrow,n,0) - tempdiff_fn
              if (mskhr(i,jrow) .ne. 0) then
                ttn2(6,jrow,n,mskhr(i,jrow)) = 
     &                   ttn2(6,jrow,n,mskhr(i,jrow)) + tempadv_fn
                ttn2(7,jrow,n,mskhr(i,jrow)) = 
     &                   ttn2(7,jrow,n,mskhr(i,jrow)) - tempdiff_fn
              endif
            enddo
          enddo
c
#ifdef rot_grid
          do i=is,ie
            if (cori(i,jrow,1) .eq. c0 .and. jrow .gt. 1) then
              factor = c4*cori(i,jrow-1,1)
            else
              factor = c4*cori(i,jrow,1)
            endif
#else
          if (cori(jrow,1) .eq. c0 .and. jrow .gt. 1) then
            factor = c4*cori(jrow-1,1)
          else
            factor = c4*cori(jrow,1)
          endif
          do i=is,ie
#endif          
            totz = c0
            vbrz = c0
            tbrz = c0
            do k=1,km
	      mask = tmask(i,k,j)*tmask(i,k,j+1)
              vbrz = vbrz + adv_vnt(i,k,j)*dxt(i)
# if !defined partial_cell
     &                                           *dzt(k)
# endif
              tbrz = tbrz +mask*(t(i,k,j,n,tau)+t(i,k,j+1,n,tau))
# ifdef partial_cell
     &                           *min(dht(i,k,j),dht(i,k,j+1))
# else
     &                           *dzt(k)
# endif
              totz = totz + mask
# ifdef partial_cell
     &                          *min(dht(i,k,j),dht(i,k,j+1))
# else
     &                          *dzt(k)
# endif
            enddo
            if (totz .ne. c0) then
              tbrz = tbrz/totz
              ttn(3,jrow,n) = ttn(3,jrow,n) + vbrz*tbrz*p5
              ttn(5,jrow,n) = ttn(5,jrow,n) - (smf(i,j,1)*dxu(i) + 
     &                        smf(i-1,j,1)*dxu(i-1))*(t(i,1,j,n,tau)
     &                        +t(i,1,j+1,n,tau)-tbrz)
     &         *csu(jrow)/factor
            endif
          enddo
          ttn(2,jrow,n) = ttn(6,jrow,n)-ttn(1,jrow,n)
          ttn(4,jrow,n) = ttn(6,jrow,n)-ttn(3,jrow,n)-ttn(5,jrow,n)
          ttn(8,jrow,n) = ttn(6,jrow,n)+ttn(7,jrow,n)
          do ll=0,nhreg
            ttn2(8,jrow,n,ll) = ttn2(6,jrow,n,ll)+ttn2(7,jrow,n,ll)
          enddo
        endif
      enddo
# ifdef timing
      call toc ('diagnostic', 'gyre components')
# endif
#endif
      return
      end
