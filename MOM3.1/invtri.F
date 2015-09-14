      subroutine invtri (z, topbc, botbc, dcb, tdt, kmz, mask, is, ie
     &,                  joff, js, je, caller)
#if defined implicitvmix || defined redi_diffusion
c
c=======================================================================
c     solve the vertical diffusion equation implicitly using the
c     method of inverting a tridiagonal matrix as described in
c     Numerical Recipes in Fortran, The art of Scientific Computing,
c     Second Edition, Press, Teukolsky, Vetterling, Flannery, 1992
c     pages 42,43.
c     this routine assumes that the variables are defined at grid points
c     and the top and bottom b.c. are flux conditions.
c
c     inputs:
c     z         = right hand side terms
c     topbc     = top boundary condition
c     botbc     = bottom boundary condition
c     dcb       = vertical mixing coeff at base of cell
c     tdt       = 2 * timestep
c     kmz       = level indicator
c     mask      = land/sea mask
c     is        = index of starting longitude
c     ie        = index of ending longitude
c     js        = starting latitude row in MW
c     je        = ending latitude row in MW
c     joff      = offset between jrow on disk and j in the MW
c
c     outputs:
c     z         = returned solution
c
c     author:   Anand Gnanadesikan    a1g@gfdl.gov
c=======================================================================
c
# ifdef partial_cell
      use memory_window
# endif
# if defined sponges && !defined sponges_old
      use xsponge_mod
# endif
# include "param.h"
# include "grdvar.h"
# include "mw.h"
# include "taskrows.h"
# include "vmixc.h"
# ifdef partial_cell
      dimension dhw_ut(imt,0:km,2:jmw-1)
      dimension dhut(imt,km,2:jmw-1)
# endif
      character*(*) caller
      dimension z(imt,km,jmw)
      dimension topbc(imt,1:jmw), botbc(imt,1:jmw)
      dimension dcb(imt,km,2:jmw-1)
      dimension kmz(imt,jstask:jetask), tdt(km)
      real mask(imt,km,1:jmw)
c
      dimension a(imt,km,2:jmw-1), b(imt,km,2:jmw-1)
      dimension c(imt,0:km,2:jmw-1), d(imt,km,2:jmw-1)
      dimension f(imt,0:km,2:jmw-1), e(imt,km,2:jmw-1)
      dimension bet(imt,2:jmw-1)
c
      do j=js,je
        jrow = j + joff
# if defined sponges && !defined sponges_old
        if (spng(jrow) .ne. 0.0 .and. caller .eq. 'tracer') cycle
# endif
# ifdef partial_cell
        if (caller .eq. 'tracer') then
          do k=1,km
            do i=is,ie
	      dhut(i,k,j) = dht(i,k,j)
            enddo
	  enddo
          do k=0,km
            do i=is,ie
	      dhw_ut(i,k,j) = dhwt(i,k,j)
            enddo
	  enddo
        elseif (caller .eq. 'baroclinic') then
          do k=1,km
            do i=is,ie
	      dhut(i,k,j) = dhu(i,k,j)
            enddo
	  enddo
          do k=0,km
            do i=is,ie
	      dhw_ut(i,k,j) = min(dhwt(i,k,j),   dhwt(i+1,k,j)
     &,                           dhwt(i,k,j+1), dhwt(i+1,k,j+1))
            enddo
	  enddo
        else
          write (*,*) '=>Error: in invtri. caller = ',caller
          call abort()
        endif
# endif
        do k=1,km
          km1 = max(1,k-1)
	  kp1 = min(k+1,km)
          factu = dztur(k)*tdt(k)*aidif
	  factl = dztlr(k)*tdt(k)*aidif
          do i=is,ie
# ifdef partial_cell
            dhtur = c1/(dhw_ut(i,k-1,j)*dhut(i,k,j))
            dhtlr = c1/(dhw_ut(i,k,j)*dhut(i,k,j))
            factu = dhtur*tdt(k)*aidif
	    factl = dhtlr*tdt(k)*aidif
# endif
            a(i,k,j) = -dcb(i,km1,j)*factu*mask(i,k,j)
            c(i,k,j) = -dcb(i,k,j)*factl*mask(i,kp1,j)
            f(i,k,j) = z(i,k,j)*mask(i,k,j) 
            b(i,k,j) = c1 - a(i,k,j) - c(i,k,j)
          enddo
        enddo

        do i=is,ie
          a(i,1,j)  = c0
          c(i,km,j) = c0
          b(i,1,j)  = c1 - a(i,1,j) - c(i,1,j)
#if defined bbl_ag
          a(i,km,j) = c0
	  c(i,km-1,j) = c0
	  b(i,km-1,j) = c1 - a(i,km-1,j) - c(i,km-1,j)
#endif
          b(i,km,j) = c1 - a(i,km,j) - c(i,km,j)
c
c         top and bottom b.c.
c
          f(i,1,j)  = z(i,1,j) + topbc(i,j)*tdt(1)*aidif*mask(i,1,j)
# ifdef partial_cell
     &                          /dhut(i,1,j)
# else
     &                          *dztr(1)
# endif
          k=max(2,kmz(i,jrow))
          f(i,k,j)   = z(i,k,j) - botbc(i,j)*tdt(k)*aidif*mask(i,k,j)
# ifdef partial_cell
     &                          /dhut(i,k,j)
# else
     &                          *dztr(k)
# endif
	enddo
      enddo
c
c     decomposition and forward substitution
c
      eps = 1.e-30
      do j=js,je
# if defined sponges && !defined sponges_old
        jrow = j + joff
        if (spng(jrow) .ne. 0.0 .and. caller .eq. 'tracer') cycle
# endif
        do i=is,ie
	    bet(i,j) = mask(i,1,j)/(b(i,1,j) + eps)
	    z(i,1,j) = f(i,1,j)*bet(i,j)
        enddo
	do k=2,km
          do i=is,ie
	    e(i,k,j) = c(i,k-1,j)*bet(i,j)
	    bet(i,j) = mask(i,k,j)/(b(i,k,j) - a(i,k,j)*e(i,k,j) + eps)
	    z(i,k,j) = (f(i,k,j) - a(i,k,j)*z(i,k-1,j))*bet(i,j)
          enddo
        enddo
      enddo
c
c     back substitution
c
      do j=js,je
# if defined sponges && !defined sponges_old
        jrow = j + joff
        if (spng(jrow) .ne. 0.0 .and. caller .eq. 'tracer') cycle
# endif
        do k=km-1,1,-1
          do i=is,ie
	    z(i,k,j) = z(i,k,j) - e(i,k+1,j)*z(i,k+1,j)
          enddo
        enddo
      enddo
#endif
      return
      end
