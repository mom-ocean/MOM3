      subroutine setsff (js, je)
#ifdef explicit_fresh_water_flux
c
c=======================================================================
c     set momentum and tracer vertical boundary conditions 
c
c     input:
c     js = global starting latitude index for the task
c     je = global ending latitude index for the task
c     is   = starting longitude index 
c     ie   = ending longitude index 
c=======================================================================
c
      use arrays_2d
      use sbc_info
      use sbc_ocean
# include "param.h"
# include "coord.h"
# include "grdvar.h"
# include "scalar.h"
# ifdef river_inflow
#  include "river.h"
# endif
c
c-----------------------------------------------------------------------
c     bail out if starting row exceeds ending row
c-----------------------------------------------------------------------
c
      if (js .gt. je) return
c
# ifdef timing
      call tic ('setsff', 'setsff')
# endif
c
c----------------------------------------------------------------------
c     set freshwater flux onto U and T points 
c----------------------------------------------------------------------
c
      sfft(:,:) = 0.0
      sffu(:,:) = 0.0
      
      do jrow=js,je-1
	dyn  = dun(jrow)*cst(jrow+1)
	dys  = dus(jrow)*cst(jrow)
	dyr  = dyur(jrow)*csur(jrow)
        do i=2,imt-1
	  asw = duw(i)*dys
	  anw = duw(i)*dyn
	  ase = due(i)*dys
	  ane = due(i)*dyn
# ifdef simple_sbc
          sfft(i,jrow)=3.e-6*min(kmt(i,jrow),1)
# else
          itflx = index_of_sbc ('fresh wtr')
          if (itflx .ne. 0) then
            sfft(i,jrow) = sbcocn(i,jrow,itflx)*emsk(i,jrow)
          endif
# endif   
          sffu(i,jrow) = dyr*dxur(i)*(
     &                    sfft(i,jrow)*asw + sfft(i+1,jrow)*ase
     &                  + sfft(i,jrow+1)*anw + sfft(i+1,jrow+1)*ane)
        enddo
      enddo
      
c
# ifdef river_inflow
c
c----------------------------------------------------------------------
c     add the river runoff to surface freshwater flux
c----------------------------------------------------------------------
c
      call getriver
      do nr=1,nriv 
        do m=1,nboxriv(nr)
          ir = iriv(m,nr)
          jr = jriv(m,nr)
          if (js.le.jr.and.je.ge.jr.and.1.lt.ir.and.imt.gt.ir) then
            sfft(ir,jr) = sfft(ir,jr) + rff(nr) 
#  ifdef coupled
c
c           add tracer fluxes due to river runoff
c
            it = index_of_sbc ('heatflux')
            do n=1,nt
	      itflx = it + (n-1)
	      if (itflx .ne. 0) then
		sbcocn(ir,jr,itflx) = sbcocn(ir,jr,itflx)
     &                                 + rtf(nr,n) * rff(nr)
	      endif
            enddo
#  endif             
          endif
        enddo
      enddo
# endif
# ifdef obc
c
c----------------------------------------------------------------------
c     set all surface fluxes on the open boundaries to zero
c----------------------------------------------------------------------
c
      do jrow=js,je
#  ifdef obc_south
        if (jrow .eq. 2) then
	  do i=1,imt
	    sfft(i,jrow) = c0
	  enddo
	endif
#  endif
#  ifdef obc_north
        if (jrow .eq. jmtm1) then
	  do i=1,imt
	    sfft(i,jrow) = c0
	  enddo
	endif
#  endif
#  ifdef obc_west
	sfft(2,jrow)     = c0
#  endif
#  ifdef obc_east
	stf(imtm1,jrow) = c0
#  endif
      enddo
# endif
c
c----------------------------------------------------------------------
c     apply zonal boundary conditions
c----------------------------------------------------------------------
c
# ifdef explicit_fresh_water_flux
      call setbcx (sfft(1,js), imt, je-js+1)
# endif
c
# ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     & "=> In setsff: js=",js," je=",je
     & 
# endif
# ifdef timing
      call toc ('setsff', 'setsff')
# endif
c
#endif
      return
      end
