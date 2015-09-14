      subroutine rinum (joff, js, je, is, ie)
c
c=======================================================================
c
c     compute richardson number
c
c     inputs:
c
c      joff   = offset between rows in the MW and latitude rows
c      js     = starting row for loading variables to calculate
c               coefficients.
c      je     = ending row for loading variables to calculate
c               coefficients. calculations end at je-1
c      is     = starting index for calculating coefficients in the
c               longitude direction
c      ie     = ending index for calculating coefficients in the
c               longitude direction
c      km     = number of vertical levels
c      grav   = gravity (cm/sec**2)
c
c     outputs:
c
c      rit     = richardson number at bottom of T cells
c      riu     = richardson number at bottom of U cells
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use memory_window
#if defined smooth_Richardson_number
      use arrays_2d
#endif
      integer tlev
# include "param.h"
# include "coord.h"
# include "mw.h"
# include "rinum.h"
# include "scalar.h"
c
      dimension ro(imt,km,jmw)
c
c-----------------------------------------------------------------------
c     set local constants
c-----------------------------------------------------------------------
c
      fx    = -p25*grav
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
c-----------------------------------------------------------------------
c     set time level
c-----------------------------------------------------------------------
c
# ifdef implicitvmix
      tlev = tau
# else
      tlev = taum1
# endif
c
c-----------------------------------------------------------------------
c     compute density difference across bottom of T cells at tau-1
c-----------------------------------------------------------------------
c
      do ks=1,2
        call statec (t(1,1,1,1,tlev), t(1,1,1,2,tlev), ro(1,1,1)
     &,              js, je, istrt-1, iend+1, ks)
        do j=js,je
          do k=ks,kbot-1,2
            do i=istrt-1,iend+1
              rhom1z(i,k,j) = (ro(i,k,j) - ro(i,k+1,j))*tmask(i,k+1,j)
            enddo
          enddo
	enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute richardson numbers on bottom of U cells
c-----------------------------------------------------------------------
c
      jsriu = max(js,1+jbuf)-1
      do j=jsriu,je-1
        do k=1,kbot-1
          t1 = fx*dzw(k)
          do i=istrt-1,iend
# ifdef partial_cell
          t1 = fx*dhwt(i,k,j)
# endif
            riu(i,k,j) = t1*umask(i,k+1,j)*(
     &                   rhom1z(i,k,j+1) + rhom1z(i+1,k,j+1) +
     &                   rhom1z(i,k,j)   + rhom1z(i+1,k,j)) /
     &                   ((u(i,k,j,1,tlev) - u(i,k+1,j,1,tlev))**2 +
     &                    (u(i,k,j,2,tlev) - u(i,k+1,j,2,tlev))**2 + 
     &                    epsln)
          enddo
        enddo
	call setbcx (riu(1,1,j), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     compute richardson numbers on bottom of T cells as average
c     of four nearest richardson numbers on bottom of U cells.
c     (do not consider land cells in the average... only active ones) 
c-----------------------------------------------------------------------
c
      jstrt = max(js-1,1+jbuf)
      do j=jstrt,je-1
        do k=1,kbot-1
          do i=istrt,iend
	    active_cells = umask(i,k+1,j) + umask(i-1,k+1,j) +
     &                     umask(i,k+1,j-1) + umask(i-1,k+1,j-1) + epsln
            rit(i,k,j)   = (riu(i,k,j) + riu(i-1,k,j)
     &                    + riu(i,k,j-1) + riu(i-1,k,j-1))/active_cells
c
c           make sure no static instability exists which is not seen by
c           the Richardson number.  This may happen due to the
c           horizontal averaging used in calculating the Richardson
c           number:
c
            if (rit(i,k,j) .gt. 0.0 .and. rhom1z(i,k,j) .gt. 0.0) then
              rit(i,k,j)=-10.0
            endif
          enddo
        enddo
      enddo

#if defined smooth_Richardson_number
c-----------------------------------------------------------------------
c     smooth Richardson number in the vertical using a 121 filter:
c-----------------------------------------------------------------------
      num_v_smooth_Ri=1
      do mr = 1,num_v_smooth_Ri
        do j=jstrt,je-1
          jrow=j+joff
          do i=istrt,iend
            rit_prev    =  p25 * rit(i,1,j)
            kmax=kmt(i,jrow)
            if (kmax.gt.0) then
#ifdef bbl_ag
c             in bbl, this do loop still goes to kmax-2:
#endif
              do k=2,kmax-2
                tmp      =  rit(i,k,j)
                rit(i,k,j) = rit_prev + p5 * rit(i,k,j)
     &               + p25 * rit(i,k+1,j)
                rit_prev   =  p25 * tmp
              enddo
            endif
          enddo
        enddo
      enddo
#endif

c
      return
      end








