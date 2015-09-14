      subroutine setvbc (joff, js, je, is, ie)
c
c=======================================================================
c     set momentum and tracer vertical boundary conditions 
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c=======================================================================
c
      use arrays_2d
      use memory_window
      use sbc_info
      use sbc_ocean
#include "param.h"
#include "coord.h"
#include "grdvar.h"
#include "mw.h"
#include "scalar.h"
#ifdef simple_sbc
      dimension star(2+nt)
#endif
#ifdef river_inflow 
# include "river.h"
#endif
c
c-----------------------------------------------------------------------
c     bail out if starting row exceeds ending row
c-----------------------------------------------------------------------
c
      if (js .le. je) then
c
#ifdef timing
      call tic ('setvbc', 'setvbc')
#endif
c
c-----------------------------------------------------------------------
c     limit the longitude indices
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
c----------------------------------------------------------------------
c     set no flux condition for all tracers at surface & bottom.
c----------------------------------------------------------------------
c
      do n=1,nt
        do j=js,je
          do i=istrt,iend
            stf(i,j,n) = c0
            btf(i,j,n) = c0
          enddo
        enddo
      enddo
c
#if defined minimize_sbc_memory
c
c----------------------------------------------------------------------
c     construct "sbcocn" for latitudes corresponding to memory window
c     rows js through je.
c----------------------------------------------------------------------
c
      call atmos (joff, js, je)
#endif
c
c----------------------------------------------------------------------
c       apply surface tracer and momentum fluxes from the atmosphere
c----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
c
#ifdef simple_sbc
        do nn=1,2+nt
          star(nn) = c0
        enddo
        call bcest (yt(jrow), yu(jrow), star)
# if !defined salinity_psu
        star(4) = (star(4) - 35.0) * 0.001
# endif
c
c       star(1) is taux in dynes/cm**2
c       star(2) is tauy in dynes/cm**2
c       star(3) is SST in deg C
c       star(4) is SSS in parts per part
c       star(>4) is zero
c
#endif
#ifdef linear_tstar
        y_lat   = abs(yt(jrow))
        star(3) = 26.93 - 0.3618*y_lat
        star(4) = 34.72
# if !defined salinity_psu
        star(4) = (star(4) - 35.0) * 0.001
# endif
#endif
#ifdef equatorial_thermocline
        star(3) = tbarz(1,1)
        star(4) = tbarz(1,2)
#endif
c
c----------------------------------------------------------------------
c       surface tracer fluxes (in order from n=1..nt)
c----------------------------------------------------------------------
c
        it = index_of_sbc ('heatflux')
        do n=1,nt
	  itflx = it + (n-1)
          do i=istrt,iend
            stf(i,j,n) = c0
c
#if defined simple_sbc && defined restorst
c
c           surface fluxes are generated internally 
c
            stf(i,j,n) = damp_factor(n)*tmask(i,1,j)
     &                  *(star(n+2) - t(i,1,j,n,taum1))
# if defined explicit_free_surface && !defined explicit_free_surface_linear
     &                  *(dzt(1)+etat(i,jrow,1))/dzt(1)
# endif

# if defined explicit_fresh_water_flux
            if (n .eq. 2) then
              stf(i,j,n) = c0
            endif
# endif
#endif
c
#if !defined simple_sbc
c
c           surface fluxes are provided from an "atmos" model
c
            stf(i,j,n) = (sbcocn(i,jrow,itflx))*tmask(i,1,j)
#endif      
c
#if defined explicit_fresh_water_flux && defined river_inflow
            if (n .eq. 2) then
c             
c             Add tracer fluxes coming with the rivers
c             
              do nr=1,nriv 
                do m=1,nboxriv(nr)
                  ir = iriv(m,nr)
                  jr = jriv(m,nr)
                  if (jrow .eq. jr .and. i .eq. ir) then
                    stf(i,j,n) = rff(nr)*tmask(i,1,j)
     &                   *(rtf(nr,n)-t(i,1,j,n,taum1))
                  endif
                enddo 
              enddo
            endif
#endif      
          enddo
        enddo
c
c----------------------------------------------------------------------
c       surface stress
c----------------------------------------------------------------------
c
#ifdef simple_sbc
c
c       simple case where idealized winds are generated internally
c
# ifdef constant_taux
        star(1) = taux0
# endif
# ifdef constant_tauy
        star(2) = tauy0
# endif
        do i=istrt,iend
          smf(i,j,1) = star(1)*umask(i,1,j)
          smf(i,j,2) = star(2)*umask(i,1,j)
# ifdef analytic_zonal_winds
          y_lat = pi*abs(yu(jrow))/180.0
          smf(i,j,1) = (0.8*(-sin(6.0*y_lat)-1.0)
     &    +0.5*(tanh(5.0*pi-10.0*y_lat)+tanh(10.0*y_lat)))
     &    *umask(i,1,j)
          smf(i,j,2) = c0 
# endif
# ifdef air_pressure
          smf(i,j,3) = c0
# endif
        enddo
#endif
#if !defined simple_sbc
c
c       general case where winds come from an atmospheric model
c       
        do n=1,2
	  if (n .eq. 1) then
	    itau = index_of_sbc('taux')
	  else
	    itau = index_of_sbc('tauy')
	  endif
          if (itau .ne. 0) then
            do i=istrt,iend
              smf(i,j,n) = sbcocn(i,jrow,itau)*umask(i,1,j)
            enddo
          else
            do i=istrt,iend
              smf(i,j,n) = c0
            enddo
          endif
        enddo
# ifdef air_pressure
        ipres = index_of_sbc ('air press')
        if (ipres .ne. 0) then
          do i=istrt,iend
            smf(i,j,3) = sbcocn(i,jrow,ipres)*tmask(i,1,j)
          enddo
        else
          do i=istrt,iend
            smf(i,j,3) = c0
          enddo
        endif
# endif
#endif
      enddo
c
#ifdef obc
c
c----------------------------------------------------------------------
c       set all surface fluxes on the open boundaries to zero
c----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
# ifdef obc_south
        if (jrow .eq. 2) then
	  do i=1,imt
	    do n=1,nt
	      stf(i,j,n) = c0
	    enddo
	    do n=1,2
	      smf(i,j,n) = c0
	    enddo
	  enddo
	endif
# endif
# ifdef obc_north
        if (jrow .eq. jmtm1) then
	  do i=1,imt
	    do n=1,nt
	      stf(i,j,n) = c0
	    enddo
	  enddo
	  do i=1,imt
	    do n=1,2
	      smf(i,j,n) = c0
	    enddo
	  enddo
	endif
        if (jrow .eq. jmtm2) then
	  do i=1,imt
	    do n=1,2
	      smf(i,j,n) = c0
	    enddo
	  enddo
	endif
# endif
	do n=1,nt
# ifdef obc_west
	  stf(2,j,n)     = c0
# endif
# ifdef obc_east
	  stf(imtm1,j,n) = c0
# endif
	enddo
	do n=1,2
# ifdef obc_west
	  smf(2,j,n)     = c0
# endif	  
# ifdef obc_east
	  smf(imtm2,j,n) = c0
	  smf(imtm1,j,n) = c0
# endif	  
	enddo
      enddo
#endif
c
#if !defined bbl_ag
c
c----------------------------------------------------------------------
c       set bottom drag
c----------------------------------------------------------------------
c
      do n=1,2
        if (cdbot .eq. c0) then
          do j=js,je
            do i=istrt,iend
	      bmf(i,j,n) = c0
            enddo
          enddo
        else
          do j=js,je
            jrow = j + joff
            do i=istrt,iend
              kz = kmu(i,jrow)
              if (kz .ne. 0) then
                uvmag    = sqrt(u(i,kz,j,1,taum1)**2 + 
     &                          u(i,kz,j,2,taum1)**2)
                bmf(i,j,n) = cdbot*u(i,kz,j,n,taum1)*uvmag
              else
                bmf(i,j,n) = c0
              endif
            enddo
          enddo
        endif
      enddo
#endif
c
c----------------------------------------------------------------------
c     apply zonal boundary conditions
c----------------------------------------------------------------------
c
      do n=1,nt
        call setbcx (stf(1,js,n), imt, je-js+1)
        call setbcx (btf(1,js,n), imt, je-js+1)
      enddo
      do n=1,2
        call setbcx (smf(1,js,n), imt, je-js+1)
        call setbcx (bmf(1,js,n), imt, je-js+1)
      enddo
c
#ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     & "=> In setvbc: js=",js," je=",je," joff=",joff
     &," jrows=",js+joff," to ",je+joff
     & 
#endif
#ifdef timing
      call toc ('setvbc', 'setvbc')
#endif
      endif
c
      return
      end
