      subroutine state (tt, s, ro, ro_full, js, je, is, ie, joff)

c
c=======================================================================
c
c     state computes normalized densities by using a 3rd
c     order polynomial fit to the knudsen formula, for each level
c     subtract normalizing constants from temperature and salinity
c     and compute polynomial approximation of knudsen density.
c
c     note.. for precision purposes, there is a depth dependent
c     constant subtracted from the density returned by this routine.
c     so... this routine should be used only for horizontal gradients
c     of density.
c
c     inputs:
c
c     t  = the input row of temperatures (potential deg C)
c     s  = the input row of salinities (units: (ppt-35)/1000)
c     js = starting row for computing density within the MW
c     je = ending row for computing density within the MW
c     is = starting longitude index for computing density within the MW
c     ie = ending longitude index for computing density within the MW
c
c     output:
c
c     rho = normalized densities
c     These densities are in cgs units(g/cm3) and represent 
c     the in situ density at a level minus a depth dependent
c     normalization. The complete in situ density is given by 
c     rho_complete(i,k,j) = dens (t(i,k,j)-to(k), s(i,k,j)-so(k), k) 
c                           + rho_norm(k)*10-3,
c     where rho_norm(k) are the depth dependent normalization densities
c     [in sigma units (density-1)*1000] 
c
c=======================================================================
c
      use arrays_2d, only : kmt
#include "param.h"
      dimension tt(imt,km,jmw), s(imt,km,jmw), ro(imt,km,jmw)
#if defined partial_cell || defined bbl_ag
c
c     compute density anomoly ro = UNESCO density - rho0 
c
# ifdef bbl_ag
#  include "bbl.h"
# endif
# include "coord.h"
# include "grdvar.h"
# include "scalar.h"
# include "switch.h"
# include "mw.h"
      dimension p(imt,km,jmw), psu(imt,km,jmw), ro_full(imt,km,jmw)
c
c-----------------------------------------------------------------------
c     construct pressure (dynes/cm^2) at each T cell grid point
c     for use in computing density in the UNESCO equation of state.
c-----------------------------------------------------------------------
c      
# if defined pressure_depth
c
c     ro must be the full density for this computation.
c
      call hydrostatic_pressure (p, ro_full, js, je, joff)

# else
      do j=js,je
        do k=1,km
	  do i=1,imt
	    p(i,k,j) = grav*rho0*zt(k)
	  enddo
	enddo
      enddo
#  if defined partial_cell
      do j=js,je
        jrow = j + joff
	do i=1,imt
	  kb = max(kmt(i,jrow),2)
	  p(i,kb,j) = grav*rho0*(zw(kb-1) + fracdz(kb,0)*dht(i,kb,j))
	enddo
      enddo
#  endif
#  if defined bbl_ag
      do j=js,je
        jrow = j + joff
	do i=1,imt
	  kb = max(kmt(i,jrow),1)
	  do k=kb+1,km
	    p(i,k,j) = grav*rho0*zt_bbl(i,j)
	  enddo
	enddo
      enddo
#  endif
# endif
c
c-----------------------------------------------------------------------
c     convert pressure from dynes/cm^2 to bars
c     convert salinity from model units to psu
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=1,km
	  do i=1,imt
# if defined salinity_psu
	    psu(i,k,j) = s(i,k,j)
# else
	    psu(i,k,j) = s(i,k,j)*1000.0 + 35.0
# endif
	    p(i,k,j)   = p(i,k,j)*1.e-6
	  enddo
	enddo
      enddo
c
c-----------------------------------------------------------------------
c     construct full density at T cell grid points
c-----------------------------------------------------------------------
c
      call unesco_density (tt, psu, p, ro, js, je)
c
c-----------------------------------------------------------------------
c     subtract mean density to produce density anomoly
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=1,km
	  do i=1,imt
	    ro(i,k,j) = ro(i,k,j) - rho0
	  enddo
	enddo
      enddo
c
# if defined bbl_ag
c
c     set density anomoly below ocean bottom to density anomoly in bbl
c
      do j=js,je
        jrow = j + joff
        do i=1,imt
	  kb = max(kmt(i,jrow),1)
	  do k=kb+1,km
	    ro(i,k,j) = ro(i,km,j)
	  enddo
	enddo
      enddo
# endif
#else
# include "state.h"
# include "dens.h"
c
# ifdef timing
      call tic ('state', 'state')
# endif
      do j=js,je
        do k=1,km
          do i=is,ie
            ro(i,k,j) = dens (tt(i,k,j)-to(k), s(i,k,j)-so(k), k)
          enddo
        enddo
      enddo
# ifdef timing
      call toc ('state', 'state')
# endif
#endif
      return
      end


      subroutine hydrostatic_pressure (p, ro, js, je, joff)
c
c-----------------------------------------------------------------------
c     compute hydrostatic pressure p(dynes/cm^2) at T cell grid points
c     if "ro" is an anomoly, p will be a hydrostatic pressure anomoly
c     if "ro" is a full density, p will be a hydrostatic pressure
c-----------------------------------------------------------------------
c
      use arrays_2d
      use memory_window
#include "param.h"
#ifdef bbl_ag
# include "bbl.h"
#endif
#include "coord.h"
#include "mw.h"
#include "scalar.h"
      dimension p(imt,km,jmw), ro(imt,km,jmw)
      
#ifdef pressure_gradient_average
c
c     construct density as in Brown and Campana (1978) (see manual)
c
      call state (t(1,1,1,1,taum1), t(1,1,1,2,taum1), rhotaum1(1,1,1)
     &,         ro,  js, je+1, istrt-1, iend+1)
      call state (t(1,1,1,1,taup1), t(1,1,1,2,taup1), rhotaup1(1,1,1)
     &,         ro,  js, je+1, istrt-1, iend+1)
      aprime = 0.25
      do j=js,je+1
        do k=1,km
          do i=1,imt
	    ro(i,k,j) = aprime*(rhotaum1(i,k,j) + rhotaup1(i,k,j))
     &       + (1.0-2.0*aprime)*ro(i,k,j)
          enddo
	enddo
      enddo
#endif
      p5grav = p5*grav
      do j=js,je
        jrow = j + joff
	do i=1,imt
          p(i,1,j) = ro(i,1,j)*grav*dzw(0)
	  kb = max(kmt(i,jrow),1)
          do k=2,kb
            p(i,k,j) = p(i,k-1,j) + (ro(i,k-1,j)+ro(i,k,j)) 
#ifdef partial_cell
     &                  *p5grav*dhwt(i,k-1,j)
#else
     &                  *p5grav*dzw(k-1)
#endif
          enddo
	  p_bottom = p(i,kb,j)
#ifdef bbl_ag
     &                + (ro(i,kb,j)+ro(i,km,j))*p5grav*dzw_bbl(i,j)
#endif
          do k=kb+1,km
            p(i,k,j) = p_bottom
          enddo
        enddo
	call setbcx (p(1,1,j), imt, km)
      enddo
      return
      end

#ifdef air_pressure
      subroutine add_air_pressure (js, je, is, ie)
      use memory_window
# include "size.h"
# include "mw.h"
      do j=js,je
	do i=is,ie
          do k=1,km
            pressure(i,k,j) = pressure(i,k,j) + smf(i,j,3)
          enddo
        enddo
	call setbcx (pressure(1,1,j), imt, km)
      enddo
      return
      end
#endif

      subroutine statec (t, s, rho, js, je, is, ie, ind)
c
c
c=======================================================================
c
c     statec computes, for one row, the normalized densities by using
c     a 3rd order polynomial fit to the knudsen formula. For
c     purposes of checking vertical stability between adjacent
c     levels, the reference depth for pressure dependence in
c     the knudsen formula must be held constant. that level is
c     determined by "ind".
c
c     inputs:
c
c     t   = the input row of temperatures (potential deg C)
#if defined salinity_psu
c     s   = the input row of salinities (units: psu)
#else
c     s   = the input row of salinities (units: (ppt-35)/1000)
#endif
c     js  = starting row for computing density within the MW
c     je  = ending row for computing density within the MW
c     is  = starting longitude index for computing density within the MW
c     ie  = ending longitude index for computing density within the MW
c     ind = 1 for comparing levels 1 to 2, 3 to 4, etc.
c           (coefficients for the lower of the 2 levels are used)
c           2 for comparing levels 2 to 3, 4 to 5, etc.
c           (coefficients for the lower of the 2 levels are used)
c
c     output:
c
c     rho = normalized densities
c     These densities are in cgs units(g/cm3) and represent 
c     the in situ density at a level minus a depth dependent
c     normalization. The complete in situ density is given by 
c     rho_complete(i,k,j) = dens (t(i,k,j)-to(k), s(i,k,j)-so(k), k) 
c                           + rho_norm(k)*10-3,
c     where rho_norm(k) are the depth dependent normalization densities
c     [in sigma units (density-1)*1000]
c
c=======================================================================
c
#include "param.h"
#include "state.h"
      dimension t(imt,km,jmw), s(imt,km,jmw), rho(imt,km,jmw)
#include "dens.h"
c
#ifdef timing
      call tic ('state', 'statec')
#endif
      if (ind .lt. 1 .or. ind .gt. 2) then
        write (stderr,99) ind
        call abort()
      endif
c
      do j=js,je
        do l=1,km,2
          if (ind .eq. 1) then
            k = min(l+1,km)
          else
            k = l
          endif
          do i=is,ie
            rho(i,l,j) = dens (t(i,l,j)-to(k), s(i,l,j)-so(k), k)
          enddo
        enddo
      enddo
c
      do j=js,je
        do l=2,km,2
          if (ind .eq. 1) then
            k = l
          else
            k = min(l+1,km)
          endif
          do i=is,ie
            rho(i,l,j) = dens (t(i,l,j)-to(k), s(i,l,j)-so(k), k)
          enddo
        enddo
      enddo
#ifdef timing
      call toc ('state', 'statec')
#endif
      return
   99 format(/' error => bad "ind" in statec: ind =',i10)
      end


# ifdef potential_density
      subroutine state_ref (t, s, rho, js, je, is, ie, kref)
c
c=======================================================================
c
c     Construct potential density referenced to the level kref.
c
c=======================================================================
c
#include "param.h"
#include "state.h"
      dimension t(imt,km,jmw), s(imt,km,jmw), rho(imt,km,jmw)
#include "dens.h"
c
#ifdef timing
      call tic ('state_ref', 'state_ref')
#endif
      do j=js,je
        do k=1,km
          do i=is,ie
            rho(i,k,j) = dens(t(i,k,j)-to(kref),s(i,k,j)-so(kref),kref)
          enddo
        enddo
      enddo
#ifdef timing
      call toc ('state_ref', 'state_ref')
#endif
c
      return
      end
# endif





