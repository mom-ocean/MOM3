      subroutine statekpp (joff, js, je, is, ie)  
#ifdef kppvmix
c
c-----------------------------------------------------------------------
c     "statekpp" computes for one slab all necessary input arrays 
c     for the kpp mixing scheme based on the mom implementation of 
c     the equation of state ("z-model" version).
c
c     input:
c      t     = the input row of temperatures
c      s     = the input row of salinities  (model units:(ppt-35)/1000)
c
c     output:
c      rho   = potential density of surface layer           (g/cm^3)
c      dbloc = local delta buoyancy at interfaces           (cm/s^2)
c      dbsfc = delta buoyancy with respect to surface       (cm/s^2)
c      talpha= d(rho)/ d(potential temperature)           (g/cm^3/C)
c      sbeta = d(rho)/ d(salinity)                      (g/cm^3/PSU)
c         
c     related documentation: changes.density
c
c     written by: jan morzel
c     date:       february 10, 1995  (converted from "sigma" version)
c-----------------------------------------------------------------------
c
      use memory_window
# include "param.h"
# include "mw.h"
# include "scalar.h"
# include "state.h"
# include "vmixc.h"
# include "dens.h"
c
c input
c      real t(imt,km),        ! potential temperature     (C)
c     &     s(imt,km)         ! salinity                (model units)
c
c output
c      real rhosfc(imt,jmw),   ! potential density of sfc layer(g/cm^3)
c     &     talpha(imt,km,jmw),! d(rho)/ d(pot.temperature) (g/cm^3/C)
c     &     sbeta(imt,km,jmw), ! d(rho)/ d(salinity)    (g/cm^3/PSU)
c     &     dbloc(imt,km,jmw), ! local delta buoyancy at interfaces
c                                (cm/s^2)
c     &     dbsfc(imt,km,jmw)  ! delta buoyancy w/ respect to sfc
c                              (cm/s^2) note: 0. values at k=1 only  
c                              needed for mixed layer depth computation
c local
      real temp(imt,km,jmw),     ! temporary array for temperature (C)
     &     rhok,             ! density of t(k  ) & s(k  ) at depth k
     &     rhokm1,           ! density of t(k-1) & s(k-1) at depth k
     &     rho1              ! density of t(1  ) & s(1  ) at depth k
c
c-----------------------------------------------------------------------
c     set local constants
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
c-----------------------------------------------------------------------
c     because the Unesco equation of state is not valid below -2C, the
c     minimum temperature is set to -2C here for the first "ktchk" 
c     vertical levels 
c-----------------------------------------------------------------------
c
      do j=js,je
        ktchk = km
         do k=1,ktchk
           do i=istrt,iend
             if (t(i,k,j,1,tau).lt.-2.) then
               temp(i,k,j) = -2.
             else
               temp(i,k,j) = t(i,k,j,1,tau)
             endif
           enddo
         enddo
        if (ktchk .lt. km) then
          do k=ktchk+1,km
            do i=1,iend
              temp(i,k,j) = t(i,k,j,1,tau)
            enddo
          enddo
        endif
      enddo
c
c-----------------------------------------------------------------------
c     alpha and beta at centers of T cells
c-----------------------------------------------------------------------
c
      do j=js,je
	do k=1,km
	  do i=istrt,iend
	    tprime = temp(i,k,j)-to(k)
	    sprime = t(i,k,j,2,tau)-so(k)
	    talpha(i,k,j) = -drhodt(tprime, sprime, k)
	    sbeta(i,k,j)  =  drhods(tprime, sprime, k)
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     calculate density, in surface layer, and set dbsfc to zero
c-----------------------------------------------------------------------
c
      k=1
      do j=js,je
        do i=istrt,iend
	  tprime = temp(i,k,j)-to(k)
	  sprime = t(i,k,j,2,tau)-so(k)
          rhosfc(i,j)   = dens (tprime,sprime, k) + ro0(k) + 1.
          dbsfc(i,k,j)  = c0
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     calculate  gradients in interior layers
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=2,km
          do i=istrt,iend
            rhok   = dens (temp(i,k,j)-to(k),t(i,k,j,2,tau)-so(k),k)
     &                    + ro0(k) + 1.
            rhokm1 = dens (temp(i,k-1,j)-to(k),t(i,k-1,j,2,tau)-so(k),k)
     &                    + ro0(k) + 1.
            rho1   = dens (temp(i,1,j)-to(k),t(i,1,j,2,tau)-so(k),k)
     &                    + ro0(k) + 1.

            dbloc(i,k-1,j) = grav * (1.0 - rhokm1/rhok)
            dbsfc(i,k ,j ) = grav * (1.0 - rho1  /rhok)

          enddo
        enddo
      enddo
      do j=js,je
        do i=istrt,iend
            dbloc(i,km,j) = dbloc(i,km-1,j)
        enddo
      enddo
#endif
      return
      end
