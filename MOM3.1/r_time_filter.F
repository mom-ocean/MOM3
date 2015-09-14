      subroutine r_time_filter (joff, js, je, is, ie, nu)
#ifdef robert_time_filter
c
c=======================================================================
c     time filter tracers and internal mode velocity at level tau
c     and write filtered tau time levels to disk
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       nu   = pointer for "tau" disk
c
c     output:
c       t_tf    = time filterd tracers
c       u_tf    = time filterd internal mode velocity
c
c     author:    Tony Rosati  e-mail => ar@gfdl.gov
c=======================================================================
c
      use memory_window
      use sbc_info
# include "param.h"
# include "iounit.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
c
# if !defined max_window
      dimension u_tf(imt,km,jmw,2), t_tf(imt,km,jmw,nt)
      dimension uim(imt,km,jmw,2,-1:0)
# endif
c
c-----------------------------------------------------------------------
c     filter tracers
c-----------------------------------------------------------------------
c
      do n=1,nt
	do j=js,je
	 do k=1,km
	   do i=is,ie
# if defined max_window
	     t(i,k,j,n,tau)  = t(i,k,j,n,tau)
# else
	     t_tf(i,k,j,n)   = t(i,k,j,n,tau)
# endif
     &       + smooth*(0.5*(t(i,k,j,n,taup1) + t(i,k,j,n,taum1))
     &       - t(i,k,j,n,tau))
	   enddo
         enddo
# if defined max_window
	 call setbcx (t(1,1,j,n,tau), imt, km)
# else
	 call setbcx (t_tf(1,1,j,n), imt, km)
# endif
        enddo
      enddo
# if !defined simple_sbc
c
c-----------------------------------------------------------------------
c     if needed, construct the Atmos S.B.C.(surface boundary conditions)
c     averaged over this segment
c     eg: SST and possibly SSS
c-----------------------------------------------------------------------
c
      isst  = index_of_sbc ('sst')
      isss  = index_of_sbc ('sss')
      if (isst .ne. 0 .or. isss .ne. 0) then
#  if defined max_window
        call asbct (joff, js, je, is, ie, isst, isss, t(1,1,1,1,tau))
#  else
        call asbct (joff, js, je, is, ie, isst, isss, t_tf(1,1,1,1))
#  endif
      endif
# endif
c
c-----------------------------------------------------------------------
c     filter internal mode velocity
c-----------------------------------------------------------------------
c
c     remove external mode from velocity at "tau" and "tau-1"
c
# if defined max_window
c
c     note: after this subtraction, "u" will be internal mode in MW
c
      call sub_ext_mode (joff, js, je, is, ie, 'tau'  , u)
      call sub_ext_mode (joff, js, je, is, ie, 'tau-1', u)
# else
      call sub_ext_mode (joff, js, je, is, ie, 'tau'  , uim)
      call sub_ext_mode (joff, js, je, is, ie, 'tau-1', uim)
# endif
      do n=1,2
	do j=js,je
      	  do k=1,km
       	    do i=is,ie
# if defined max_window
	      u(i,k,j,n,tau)   = u(i,k,j,n,tau)
     &        + smooth*(0.5*(u(i,k,j,n,taup1) + u(i,k,j,n,taum1))
     &        - u(i,k,j,n,tau))
# else
	      u_tf(i,k,j,n)   = uim(i,k,j,n,tau)
     &        + smooth*(0.5*(u(i,k,j,n,taup1) + uim(i,k,j,n,taum1))
     &        - uim(i,k,j,n,tau))
# endif
 	    enddo
          enddo
# if defined max_window
	  call setbcx (u(1,1,j,n,tau), imt, km)
# else
	  call setbcx (u_tf(1,1,j,n), imt, km)
# endif
        enddo
      enddo
# if defined max_window
c
c     note: velocities are still internal mode only in the MW
c
# else
c
c     write all newly filtered quantities to disk "tau"
c
      do j=js,je
        jrow   = j + joff
        call putrow (nu, nslab, jrow, u_tf(1,1,j,1), t_tf(1,1,j,1))
      enddo
# endif
c
      return
      end




      subroutine sub_ext_mode (joff, js, je, is, ie, timelev, uim)
c
c=======================================================================
c     subtract external mode for time level "timelev"
c     for robert_time_filter filter only internal mode velocity
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       timelev = "tau" or "tau-1"
c       uim  = internal mode velocity
c     
c=======================================================================
c
      use arrays_2d
      use memory_window
#include "param.h"
#include "emode.h"
#include "grdvar.h"
#include "iounit.h"
#include "mw.h"
#include "switch.h"
#ifdef stream_function
      dimension ext(imt,1:jmw,2)
#endif
      dimension uim(imt,km,jmw,2,-1:1)
      character*(*) timelev
c
      if (timelev .eq. 'tau') then
c
c-----------------------------------------------------------------------
c       subtract external mode "tau" 
c-----------------------------------------------------------------------
c
        do j=js,je
          jrow = j + joff
          if (jrow .lt. jmt) then
#ifdef stream_function
            do i=is,ie
              diag1       = psi(i+1,jrow+1,1) - psi(i  ,jrow,1)
              diag0       = psi(i  ,jrow+1,1) - psi(i+1,jrow,1)
	      ext(i,j,1)  = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
	      ext(i,j,2)  =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)
     &                       *csur(jrow)
            enddo
            do n=1,2
              do k=1,km
                do i=is,ie
                  uim(i,k,j,n,tau) = (u(i,k,j,n,tau) - ext(i,j,n))
     &                             *umask(i,k,j)
                enddo
              enddo
            enddo
#endif
#if defined rigid_lid_surface_pressure || defined free_surf_on
            do n=1,2
	      do k=1,km
                do i=is,ie
                  uim(i,k,j,n,tau) = (u(i,k,j,n,tau) - ubar(i,jrow,n))
     &                             *umask(i,k,j)
                enddo
              enddo
            enddo
#endif
            do n=1,2
	      call setbcx (uim(1,1,j,n,tau), imt, km)
	    enddo
	  endif
	enddo
c
      elseif (timelev .eq. 'tau-1') then
c
c-----------------------------------------------------------------------
c       subtract external mode  "tau-1" 
c-----------------------------------------------------------------------
c
        do j=js,je
          jrow = j + joff
          if (jrow .lt. jmt) then
#ifdef stream_function
            do i=is,ie
              diag1       = psi(i+1,jrow+1,2) - psi(i  ,jrow,2)
              diag0       = psi(i  ,jrow+1,2) - psi(i+1,jrow,2)
	      ext(i,j,1)  = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
	      ext(i,j,2)  =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)
     &                       *csur(jrow)
            enddo
            do n=1,2
              do k=1,km
                do i=is,ie
                  uim(i,k,j,n,taum1) = (u(i,k,j,n,taum1) - ext(i,j,n))
     &                               *umask(i,k,j)
                enddo
              enddo
            enddo
#endif
#if defined rigid_lid_surface_pressure || defined free_surf_on
            do n=1,2
	      do k=1,km
                do i=is,ie
                  uim(i,k,j,n,taum1) = (u(i,k,j,n,taum1)
     &                           - ubarm1(i,jrow,n))*umask(i,k,j)
                enddo
              enddo
            enddo
#endif
            do n=1,2
	      call setbcx (uim(1,1,j,n,taum1), imt, km)
	    enddo
	  endif
	enddo
      else
        write (stdout,'(/a,a/)') 'Error: timelev = ',timelev
	call abort()
      endif
c
#endif
      return
      end
      

