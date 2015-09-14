      subroutine vmixc (joff, js, je, is, ie)
c
c=======================================================================
c     set viscosity coefficient on bottom face of "u" cells
c     set diffusion coefficient on bottom face of "t" cells
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
      use memory_window
#include "param.h"
#include "coord.h"
#include "mw.h"
#include "switch.h"
#include "vmixc.h"
#ifdef redi_diffusion
# include "isopyc.h"
#endif
#if defined constvmix && defined implicitvmix
      dimension temp(imt,km,jmw)
#endif
c
c-----------------------------------------------------------------------
c     bail out if starting row exceeds ending row
c-----------------------------------------------------------------------
c
      if (js .le. je) then
c
c-----------------------------------------------------------------------
c     limit the longitude and latitude indices
c-----------------------------------------------------------------------
c
        istrt = max(2,is)
        iend  = min(imt-1,ie)
        jstrt = max(js-1,1+jbuf)
        jend  = je-1
c
#if defined constvmix 
c
c-----------------------------------------------------------------------
c     constant vertical mixing coefficients
c-----------------------------------------------------------------------
c
        do j=jstrt,jend
          do i=istrt,iend
	    do k=1,kbot-1
	      visc_cbu(i,k,j) = kappa_m
# ifdef bryan_lewis_vertical
	      diff_cbt(i,k,j,1)     = Ahv(k)
	      diff_cbt(i,k,j,nsmix) = Ahv(k)
# else
	      diff_cbt(i,k,j,1)     = kappa_h
	      diff_cbt(i,k,j,nsmix) = kappa_h
# endif
	    enddo
	    do k=kbot,km
	      visc_cbu(i,k,j) = c0
              diff_cbt(i,k,j,1)     = c0
              diff_cbt(i,k,j,nsmix) = c0
	    enddo
	  enddo
        enddo
c
# ifdef implicitvmix
        do ks=1,2
c
c       find density
c
          call statec (t(1,1,1,1,taum1), t(1,1,1,2,taum1), temp(1,1,1)
     &,                jstrt, jend, istrt, iend, ks)
c
c       set large diffusion coefficient between unstable layers
c       (note: viscosity coefficient is not limited but could be here)
c
          do j=jstrt,jend
            do k=ks,kmm1,2
              do i=istrt,iend
                if (temp(i,k,j) .gt. temp(i,k+1,j)) then
                  diff_cbt(i,k,j,1)     = diff_cbt_limit*tmask(i,k+1,j) 
                  diff_cbt(i,k,j,nsmix) = diff_cbt_limit*tmask(i,k+1,j) 
                endif
              enddo
	    enddo
	  enddo
        enddo
# endif
#endif
#ifdef ppvmix
c
c-----------------------------------------------------------------------
c     for momentum and tracers based on the pacanowski & philander
c     richardson mixing scheme (JPO vol 11, #11, 1981).
c-----------------------------------------------------------------------
c
        call ppmix (joff, js, je, istrt, iend)
#endif
#ifdef kppvmix
c
c-----------------------------------------------------------------------
c     compute vertical mixing coefficients based on the k-profile
c     and oceanic planetary boundary layer scheme by 
c     Large, McWilliams, & Doney ( Rev. of Geophys.,32,363_403, 1994).
c-----------------------------------------------------------------------
c
        call kppmix (joff, js, je, istrt, iend)
#endif
#ifdef redi_diffusion
c
c-----------------------------------------------------------------------
c     Add K33 component to vertical diffusion coefficient
c-----------------------------------------------------------------------
c
        do n=1,nsmix
          do j=jstrt,jend
            do i=istrt,iend
              do k=1,kbot
                diff_cbt(i,k,j,n) = diff_cbt(i,k,j,n) + K33(i,k,j)
              enddo
            enddo
          enddo
        enddo  
#endif
c
        do n=1,nsmix
          do j=jstrt,jend
            call setbcx (diff_cbt(1,1,j,n), imt, km)
          enddo
        enddo
c
#ifdef trace_indices
        write (stdout,'(2x,5(a,i4))')
     &   "=> In vmixc: js=",js," je=",je," joff=",joff
     &  ," jstrt=",jstrt," jend=",jend
#endif

      endif

      return
      end
