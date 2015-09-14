#ifdef tracer_horz_biharmonic
      subroutine delsq_tracer (joff, js, je, is, ie)
c
c=======================================================================
c     compute del2 for tracers on MW rows "js" ... "je"
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c=======================================================================
c
      use memory_window
#include "param.h"
#include "grdvar.h"
#include "hmixc.h"
#include "mw.h"
#include "scalar.h"
#include "fdift.h"
c
c-----------------------------------------------------------------------
c     bail out if starting row exceeds ending row
c-----------------------------------------------------------------------
c
      if (js .gt. je) return
c
# ifdef timing
      call tic ('tracer', 'biharmonic tracer mixing')
# endif
c
c-----------------------------------------------------------------------
c     limit the MW row indices to what can be calculated
c-----------------------------------------------------------------------
c
      jstrt = max(js-1,2)
      jend  = je-1
c
c-----------------------------------------------------------------------
c     compute the laplacian tracer flux
c-----------------------------------------------------------------------
c
      do n=1,nt
c
c       diffusive flux across the eastern face of T cells
c       with unit diffusivity
c
        do j=jstrt,jend
          jrow = j + joff
          do k=1,km
            do i=is-1,ie              
              diff_fe(i,k,j) = (t_i(i,k,j,n,1) - t_i(i,k,j,n,0))
     &                        *cstdxur(i,j)
# if defined partial_cell
     &                        *dhte(i,k,j)
# endif
            enddo
          enddo
        enddo
c
c       diffusive flux across the northern face of T cells
c       with unit diffusivity
c
        do j=jstrt-1,jend
          jrow = max(j + joff,1)
          do k=1,km
            do i=is,ie
              diff_fn(i,k,j) = (t_j(i,k,j,n,1) - t_j(i,k,j,n,0))
     &	                       *csu_dyur(jrow)
# if defined partial_cell
     &                         *dhtn(i,k,j)
# endif
            enddo
          enddo
        enddo
c
c       Now compute del2, which is equivalent to 
c       minus the horizontal Laplacian diffusion operator.
c       Fluxes of del2 are computed in tracer.F
c
        do j=jstrt,jend
          jrow   = j + joff
          do k=1,km
            do i=is,ie
	      del2_tracer(i,k,j,n) = - DIFF_Tx(i,k,j) - DIFF_Ty(i,k,j)
            enddo
          enddo
	  call setbcx (del2_tracer(1,1,j,n), imt, km)
c
c         set del2 = 0 on southern wall
c
          if (jrow .eq. 2) then
            do k=1,km
              do i=is-1,ie+1
# ifdef obc_south
	        del2_tracer(i,k,j-1,n) = del2_tracer(i,k,j,n)
# else
	        del2_tracer(i,k,j-1,n) = c0
# endif
              enddo
            enddo
	  endif
c
c         set del2 = 0 on northern wall
c
          if (jrow-1 .eq. jmt-1) then
            do k=1,km
              do i=is-1,ie+1
# ifdef obc_north
	        del2_tracer(i,k,j,n) = del2_tracer(i,k,j-1,n)
# else
	        del2_tracer(i,k,j,n) = c0
# endif
              enddo
            enddo
          endif
c
        enddo
c
      enddo
# ifdef timing
      call toc ('tracer', 'biharmonic tracer mixing')
# endif

      return
      end
#else
      subroutine delsq_tracer_dummy
      return
      end
#endif


#ifdef velocity_horz_biharmonic
      subroutine delsq_velocity (joff, js, je, is, ie)
c
c=======================================================================
c     compute del2 of velocity on MW rows "js" ... "je"
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c=======================================================================
c
      use memory_window
#include "param.h"
#include "grdvar.h"
#include "hmixc.h"
#include "mw.h"
# ifdef partial_cell
      dimension pcsink(imt,km,jmw)
# endif
#include "scalar.h"
#include "fdifm.h"
c
c-----------------------------------------------------------------------
c     bail out if starting row exceeds ending row
c-----------------------------------------------------------------------
c
      if (js .gt. je) return
c
# ifdef timing
      call tic ('velocity', 'biharmonic velocity mixing')
# endif
c
c-----------------------------------------------------------------------
c     limit the MW row indices to what can be calculated
c-----------------------------------------------------------------------
c
      jstrt = max(js-1,2)
      jend  = je-1
c  
c-----------------------------------------------------------------------
c     compute del2 of velocities.  Use unit viscosity 
c-----------------------------------------------------------------------
c
      do n=1,2
c
c       diffusive flux across east face of U cells 
c       Use unit viscosity
c
        do j=jstrt,jend
          jrow = j + joff
          do k=1,km
            do i=is-1,ie              
              diff_fe(i,k,j) = (u(i+1,k,j,n,taum1) - u(i,k,j,n,taum1))
     &                         *csudxtr(i+1,j)
# ifdef partial_cell
     &                         *dhue(i,k,j)
# endif
            enddo
          enddo
        enddo
c
c       diffusive flux across north face of U cells
c       Use unit viscosity
c
        do j=jstrt-1,jend
          jrow = j + joff
          do k=1,km
            do i=is,ie
              diff_fn(i,k,j) = (u(i,k,j+1,n,taum1) - u(i,k,j,n,taum1))
     &                         *cst_dytr(jrow+1)
# ifdef partial_cell
     &                         *dhun(i,k,j)
# endif
            enddo
          enddo
        enddo
c
c       Sink due to partial cells
c       Use unit viscosity
c
# ifdef partial_cell
        if (n .eq. 1) then
          viscceu = c1
          visccnu = c1
          do j=jstrt,jend
            jrow = j + joff
            do k=1,km
              do i=is,ie      
                pcsink(i,k,j) = - viscceu*(csudxur(i,j)/dhu(i,k,j))*
     &          (csudxtr(i+1,j)*(dhu(i,k,j)-dhue(i,k,j))
     &         + csudxtr(i,j)*(dhu(i,k,j)-dhue(i-1,k,j)))
     &         - visccnu*(csudyur(jrow)/dhu(i,k,j))*
     &          (cst_dytr(jrow+1)*(dhu(i,k,j)-dhun(i,k,j))
     &         + cst_dytr(jrow)*(dhu(i,k,j)-dhun(i,k,j-1)))
              enddo
            enddo
          enddo
        endif
# endif
c
c       compute del2 
c       Use the friction metric from a constant unit viscosity
c
        do j=jstrt,jend
          jrow   = j + joff
          do k=1,km
            do i=is,ie
	      del2_vel(i,k,j,n) = -(  DIFF_Ux(i,k,j) + DIFF_Uy(i,k,j)
     &                              + am3(jrow)*u(i,k,j,n,taum1)
     &                              + am4(jrow,n)
     &          *(u(i+1,k,j,3-n,taum1) - u(i-1,k,j,3-n,taum1))*dxmetr(i)
# ifdef partial_cell
     &                                +pcsink(i,k,j)*u(i,k,j,n,taum1)
# endif
     &                             )*umask(i,k,j) 
            enddo
          enddo
	  call setbcx (del2_vel(1,1,j,n), imt, km)
c
c         set del2 = 0 on southern wall
c
          if (jrow .eq. 2) then
            do k=1,km
              do i=is-1,ie+1
# ifdef obc_south
	        del2_vel(i,k,j-1,n) = del2_vel(i,k,j,n)
# else
	        del2_vel(i,k,j-1,n) = c0
# endif
              enddo
            enddo
	  endif
c
c         set del2 = 0 on northern wall
c
          if (jrow-1 .eq. jmt-1) then
            do k=1,km
              do i=is-1,ie+1
# ifdef obc_north
	        del2_vel(i,k,j,n) = del2_vel(i,k,j-1,n)
# else
	        del2_vel(i,k,j,n) = c0
# endif
              enddo
            enddo
          endif
        enddo
c
      enddo  !end of velocity loop
c
# ifdef timing
      call toc ('velocity', 'biharmonic velocity mixing')
# endif
      return
      end
#else
      subroutine delsq_velocity_dummy
      return
      end
#endif
