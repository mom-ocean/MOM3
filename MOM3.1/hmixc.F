      subroutine hmixci (error, cifdef, ifdmax, nifdef
     &,                  htr_mixset, hvel_mixset)
c
c=======================================================================
c     Initialization for horizontal mixing coefficients
c
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
      logical htr_mixset, hvel_mixset, error
      integer io_status
      character*(*) cifdef(ifdmax)
c
#include "param.h"
#include "grdvar.h"
#include "hmixc.h"
#include "iounit.h"
#include "scalar.h"
      namelist /smagnl/ k_smag, diff_c_back, visc_c_back, prandtl
c
c-----------------------------------------------------------------------
c     set constant laplacian horizontal mixing coefficients
c-----------------------------------------------------------------------
c
#if defined velocity_horz_laplacian && defined velocity_horz_mix_const
      visc_cnu  = am
      visc_ceu  = am
#endif
#if defined tracer_horz_laplacian && defined tracer_horz_mix_const
      diff_cnt  = ah
      diff_cet  = ah
#endif
#ifdef full_tensor
      diff_cnt  = ah + kappa_h
      diff_cet  = ah + kappa_h
#endif
c
c-----------------------------------------------------------------------
c     set constant biharmonic horizontal mixing coefficients
c-----------------------------------------------------------------------
c
#if defined velocity_horz_biharmonic && defined velocity_horz_mix_const
      visc_cnu = abs(ambi)
      visc_ceu = abs(ambi)
#endif
#if defined tracer_horz_biharmonic && defined tracer_horz_mix_const
      diff_cet = abs(ahbi)
      diff_cnt = abs(ahbi)
#endif
c
c-----------------------------------------------------------------------
c     items used for computing metric factors in momentum friction 
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
#if defined beta_plane || defined f_plane
c
c       set metric factors = 0 on beta and f plane
c
	am3(jrow)   = c0
	am4(jrow,1) = c0
	am4(jrow,2) = c0
#else
	am3(jrow)   =  (c1-tng(jrow)*tng(jrow))/(radius**2)
	am4(jrow,1) = -c2*sine(jrow)/(radius*csu(jrow)*csu(jrow))
	am4(jrow,2) = -am4(jrow,1)
#endif
      enddo
c
#if defined velocity_horz_mix_var 
c
c-----------------------------------------------------------------------
c     metric factors for variable viscosity coefficients
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        jp1 = min(jrow+1,jmt)
	jm1 = max(jrow-1,1)
        varh1(jrow,1) = -c1/(csu(jrow)*(dyt(jp1)+dyt(jm1)))
        varh1(jrow,2) = -varh1(jrow,1)
        varh2(jrow,1) = -tng(jrow)/(radius*csu(jrow))
        varh2(jrow,2) = -varh2(jrow,1)
        varh3(jrow,1) = c1/(csu(jrow))
        varh3(jrow,2) = -varh3(jrow,1)
        varh4(jrow,1) = sine(jrow)/(radius*csu(jrow))
        varh4(jrow,2) = varh4(jrow,1)
      enddo
c
# if defined velocity_horz_mix_smag
c
c     initial values
c
      diff_c_back = c0
      visc_c_back = c0
      prandtl     = c1
      k_smag      = 0.4
c
c-----------------------------------------------------------------------
c     provide for namelist over-ride of above settings + documentation
c-----------------------------------------------------------------------
c
      call getunit (io, 'namelist.smagnl'
     &,               'formatted sequential rewind')
      read (io,smagnl,IOSTAT=io_status)
      if (io_status .gt. 0) then
        write(stdout,*) 'IOSTAT reading smagnl namelist =',io_status
        call abort()
      endif  
      write (stdout,smagnl)
      call relunit (io)
      call getunit (iodoc, 'document.dta'
     &,               'formatted sequential append')
      write (iodoc, smagnl)
      call relunit (iodoc)
# endif
#endif
      return
      end
      
      
      subroutine hmixc (joff, js, je, is, ie)
c
c=======================================================================
c     set horizontal mixing coeffs on north and east face of T and
c     U cells.
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
#include "grdvar.h"
#include "hmixc.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#ifdef full_tensor
# include "vmixc.h"
#endif
c
c-----------------------------------------------------------------------
c     bail out if starting row exceeds ending row
c-----------------------------------------------------------------------
c
      if (js .le. je) then
c
#if defined velocity_horz_mix_const 
c
c-----------------------------------------------------------------------
c     viscosities are constant in time: set in hmixci
c-----------------------------------------------------------------------
c
#endif
c
#if defined tracer_horz_mix_const
c
c-----------------------------------------------------------------------
c     diffusivities are constant in time: set in hmixci
c-----------------------------------------------------------------------
c
#endif
c
#if defined velocity_horz_mix_var || defined tracer_horz_mix_var
c
c-----------------------------------------------------------------------
c     limit the longitude indices and calculate mixing coefficients
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
      call smagnlc (joff, js, je, istrt, iend)
#endif
c
#ifdef trace_indices
      write (stdout,'(2x,3(a,i4))')
     & "=> In hmixc: js=",js," je=",je," joff=",joff
#endif
c
      endif
      return
      end
