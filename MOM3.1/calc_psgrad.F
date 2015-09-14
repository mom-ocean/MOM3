      subroutine calc_psgrad (psgrad, uext, vext, js, je, is, ie)
      use arrays_2d
      use memory_window
c
c-----------------------------------------------------------------------
c     compute the surface pressure gradients
c
c     inputs:
c
c     js   = index of starting row
c     je   = index of ending row
c     is   = index of starting longitude
c     ie   = index of ending longitude
c
c     outputs:
c
c     psgrad  = grad(surf press)
c     uext = external mode u (tau+1) for point (ie,je) only
c     vext = external mode v (tau+1) for point (ie,je) only
c-----------------------------------------------------------------------
c
#include "param.h"
#include "emode.h"
#include "grdvar.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"
      dimension psgrad(is:ie,js:je,2)
c
      is1 = max(is,2)
      ie1 = min(ie,imtm1)
c
c     on mixing time steps "ptd" has been multiplied by a factor of
c     two and the time step has to be adjusted also.
c
      if (mod(itt,nmix) .eq. 1) then
        fxa = p5
      else
        fxa = c1
      endif
      r2dtuv = c1/c2dtuv
#ifdef rot_grid
c     f3    = c2dtuv/c2dtsf
      f3    = c2dtuv
      do jrow=js,je
	do i=is1,ie1
          atosp = acor*cori(i,jrow,1)
          f2    = atosp*c2dtuv
#else
      do jrow=js,je
        atosp = acor*cori(jrow,1)
        f2    = atosp*c2dtuv
c       f3    = c2dtuv/c2dtsf
        f3    = c2dtuv
	jrowp1= min(jrow+1,jmt)
	do i=is1,ie1
#endif
          kz = kmu(i,jrow)
          if (kz .ne. 0) then
#if defined rigid_lid_surface_pressure || defined free_surf_on
            uext         = ubar(i,jrow,1)
	    vext         = ubar(i,jrow,2)
            d1           = ps(i+1,jrowp1,1) - ps(i,jrow,1)
            d2           = ps(i+1,jrow,1) - ps(i,jrowp1,1) 
	    psgrad(i,jrow,1)  = (d1 + d2)*dxu2r(i)*csur(jrow)
	    psgrad(i,jrow,2)  = (d1 - d2)*dyu2r(jrow)
#endif
#ifdef stream_function
            diag1        = psi(i+1,jrowp1,1)-psi(i  ,jrow,1)
            diag0        = psi(i  ,jrowp1,1)-psi(i+1,jrow,1)
            uext         = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
            vext         =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)*csur(jrow)
            diag3        = fxa*(ptd(i+1,jrowp1)-ptd(i  ,jrow))
            diag4        = fxa*(ptd(i  ,jrowp1)-ptd(i+1,jrow))
            dubdt        = (diag3+diag4)*dyu2r(jrow)*hr(i,jrow)
            dvbdt        = (diag3-diag4)*dxu2r(i)*hr(i,jrow)*csur(jrow)
            psgrad(i,jrow,1)=r2dtuv*(dubdt + f3*zu(i,jrow,1) + f2*dvbdt)
            psgrad(i,jrow,2)=r2dtuv*(-dvbdt+ f3*zu(i,jrow,2) + f2*dubdt)
#endif
          else
            psgrad(i,jrow,1) = c0 
            psgrad(i,jrow,2) = c0 
          endif
        enddo
      enddo
      return
      end
