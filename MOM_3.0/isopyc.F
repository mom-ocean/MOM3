      subroutine isopi (error, am, ah)
c
c=======================================================================
c
c      Initialization for isoneutral mixing schemes.  Each scheme
c      can be turned on alone, or in combination with any of 
c      the others.
c
c      -Dredi_diffusion gives either the full or small angle 
c       isoneutral diffusion tensor.  Default is the small tensor. 
c       -Dfull_tensor will use the full Redi tensor.  The full
c       tensor is not recommended for most purposes, and it is 
c       coded only for use with a constant diffusivity and -Dgm_advect.
c
c      -Dgent_mcwilliams gives the Gent_McWilliams mixing scheme, 
c       Default is -Dgm_skew. -Dgm_advect will override the 
c       skew-flux approach to give the older advective
c       formulation. Note: gm_advect without redi_diffusion is not
c       supported. 
c
c      -Dbiharmonic_rm will compute the biharmonic 
c       operator of Roberts and Marshall (1998) as implemented
c       through its skew-flux formulation.  
c
c     input:
c       error  = logical to signal problems
c       slmx   = max slope of isoneutrals
c       ahisop = isoneutral tracer diffusivity(cm**2/sec) 
c       athkdf = GM diffusivity (cm**2/sec)
c       ahbih  = Roberts and Marshall biharmonic diffusivity (cm4/sec).
c       ahsteep= horizontal diffusivity (cm**2/sec) for use
c                in regions of steep isoneutral slopes.
c       del_dm = parameter for the -Ddm_taper tapering scheme.
c       s_dm   = parameter for the -Ddm_taper tapering scheme.
c
c     output:
c       The above input can be reset via namelist
c
c     authors: 
c           R.C.Pacanowski   rcp@gfdl.gov 
c           S.M.Griffies     smg@gfdl.gov 
c=======================================================================
c
      use arrays_2d
      use memory_window
      logical error
c
#include "param.h"
#if defined isoneutralmix
# include "accel.h"
# include "coord.h"
# include "grdvar.h"
# include "iounit.h"
# include "isopyc.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# include "taskrows.h"
# include "vmixc.h"
# ifdef bryan_lewis_horizontal
      common /blisop/ Ahhisop(km)
# endif
c
      namelist /isopyc/ slmx, ahisop, athkdf, del_dm, s_dm, abihrm
     &                 ,ahsteep
c
      write (stdout,*) ' '
      write (stdout,*)
     & 'I S O N E U T R A L M I X   I N I T I A L I Z A T I O N'
c
c-----------------------------------------------------------------------
c     Default parameters 
c-----------------------------------------------------------------------
c
c     Maximum isoneutral slope
c
      slmx  = 1.0/100.0
c
c     Redi isoneutral diffusivity (cm2/sec)
c
      ahisop = 1.e7
c
c     Gent and McWilliams diffusivity (cm2/sec)
c
      athkdf = ahisop
c
c     Roberts and Marshall biharmonic diffusivity (cm4/sec)
c
      abihrm = 1.e19
c
c     Horizontal diffusivity used in regions of steep slopes (cm2/sec)
c
      ahsteep = athkdf
c
c     transition for scaling isoneutral diffusivities with -Ddm_taper
c
      del_dm = 4.0/1000.0
c
c     half width scaling for diffusivity with -Ddm_taper
c
      s_dm = 1.0/1000.0
c
c
c-----------------------------------------------------------------------
c     provide for documentation + namelist over-ride of above settings
c-----------------------------------------------------------------------
c
      call getunit (io, 'namelist.isopyc'
     &,               'formatted sequential rewind')
      read (io,isopyc,end=100)
100   continue
      write (stdout,isopyc)
      call relunit (io)
      call getunit (iodoc, 'document.dta'
     &,               'formatted sequential append')
      write (iodoc, isopyc)
      call relunit (iodoc)
c
c     set reciprocal of maximum isoneutral slope and scaling half width
c
      slmxr = c1/slmx
      s_dmr = c1/s_dm
c
c-----------------------------------------------------------------------
c     check for problems and write some notes to printout
c-----------------------------------------------------------------------
c
# if defined redi_diffusion
#  if defined full_tensor
        write (stdout,*)'==> Note: Full Redi tensor is being used.'
#   if defined partial_cell
        write (stdout,*)
     & '==> Error: option partial_cell has not been extended to work'  
     &,'           with option full_tensor.'
        error = .true.
#   endif
#   if defined nonconst_diffusivity
        write (stdout,*)
     & '==> Error: -Dfull_tensor is not compatible with'
     &  write(stdout,*) '-Dnonconst_diffusivity.'  
        error = .true.
#   endif
#   if defined isotropic_mixed
        write (stdout,*)
     & '==> Error: -Dfull_tensor is not compatible with'
     &  write(stdout,*) '-Disotropic_mixed.'  
        error = .true.
#   endif
#  endif
#  if defined small_tensor
        write (stdout,*)'==> Note: Small angle Redi tensor is used.'
#   if defined gkw_taper
        write (stdout,*)
     & '==> Note: Gerdes, Koberle, and Willebrand taper is used.'
#   endif
#   if defined dm_taper
        write (stdout,*)
     & '==> Note: Danabasoglu and McWilliams taper is used.'
#   endif
#   if defined isotropic_mixed && defined partial_cell
        write (stdout,*)
     & '==> Error: -Disotropic_mixed is not compatible w/ partial_cell.'
        error = .true.
#   endif
#  endif
# endif
c
# if defined gent_mcwilliams
#  if defined gm_skew
        write (stdout,*)'==> Note: GM90 is enabled with skew-fluxes.'
#   ifdef full_tensor
        write (stdout,*)
     & '==> Error: -Dgm_skew is not compatible with -Dfull_tensor.'  
        error = .true.
#   endif
#  endif
#  if defined gm_advect
        write (stdout,*)
     & '==> Warning: GM90 is enabled with advective fluxes.'
        write (stdout,*)
     & 'The recommended formulation is the default -Dgm_skew.'
#  endif
# endif
c
# ifdef biharmonic_rm
      write (stdout,*)
     & '==> Note: Roberts and Marshall biharmonic stirring is on.'  
# endif
c
c-----------------------------------------------------------------------
c     print out tracer diffusion coefficients
c-----------------------------------------------------------------------
c
      call getunit (iodoc, 'document.dta'
     &,               'formatted sequential append')
      write (iodoc,'(a,e14.7)') 'am=', am, 'ah=',ah
     &,     'ahisop=',ahisop, 'athkdf=',athkdf, 'slmxr=',slmxr
      call relunit (iodoc)
      write(stdout,9102) ah
      write(stdout,9103) ahisop
      write(stdout,9104) athkdf
      write(stdout,9105) ahsteep
      write(stdout,9106) abihrm
      write(stdout,9107) dtts
9102  format(' ah(cm2/sec)     = ',e12.6
     &     ,' Constant horizontal background diffusivity ')
9103  format(' ahisop(cm2/sec) = ',e12.6
     &     ,' Constant Redi isoneutral diffusivity')
9104  format(' athkdf(cm2/sec) = ',e12.6
     &     ,' Constant GM diffusivity')
9105  format(' ahsteep(cm2/sec) = ',e12.6
     &     ,' Constant horz diffusivity in steep sloped regions')
9106  format(' abihrm(cm4/sec) = ',e12.6
     &     ,' Constant RM biharmonic diffusivity')
9107  format(' dtts(sec)       = ',e12.6
     &     ,' Tracer time step')
c
# ifdef nonconst_diffusivity
        write (stdout,*)
     & '==> Note: Quasi-geostrophic closure will predict diffusivities.'
      if(ahisop .eq. athkdf) then
        write(stdout,*) 'Redi Isoneutral diffusivity is nonconstant.'
        write(stdout,*) 'GM diffusivity is nonconstant.'
        write(stdout,*) 'Horz diff in steep regions is nonconstant.'
        write(stdout,*) 'Diffusivities are the same nonconstant value.'
      else
        write(stdout,*) 'Redi diffusivity is constant ahisop.'
        write(stdout,*) 'GM diffusivity is nonconstant.'
        write(stdout,*) 'Horz diff in steep regions is same as GM diff.'
      endif     
# endif
c
c-----------------------------------------------------------------------
c     Compute the grid factors which set the maximum slopes available
c     for the mixing schemes. 
c-----------------------------------------------------------------------
c
      ft1 = c1/(4.0*max(ahisop,athkdf)*dtts + epsln)
      ft2 = c1/(4.0*abihrm*dtts + epsln)
      delta_iso1  = dzt(1)*ft1*dxt(1)*cst(jmt/2)
      delta_iso2  = dzt(1)*ft2*(dxt(1)*cst(jmt/2))**3
      i_delta1 = 1
      j_delta1 = 1
      k_delta1 = 1
      i_delta2 = 1
      j_delta2 = 1
      k_delta2 = 1
c
      do jrow=2,jetask-1
        do i=2,imt-1
	  do k=1,kbot
	    delta1a = dxt(i)*cst(jrow)*dzt(k)*ft1
	    delta1b = dyt(jrow)*dzt(k)*ft1
	    if ( delta_iso1 .ge. delta1a
     &      .or. delta_iso1 .ge. delta1b) then
              i_delta1 = i
              j_delta1 = jrow
              k_delta1 = k
	      delta_iso1  = min(delta1a,delta1b)
	    endif 
	    delta2a = (dxt(i)*cst(jrow))**3*dzt(k)*ft2
	    delta2b =  dyt(jrow)**3        *dzt(k)*ft2
	    if ( delta_iso2 .ge. delta2a
     &      .or. delta_iso2 .ge. delta2b) then
              i_delta2 = i
              j_delta2 = jrow
              k_delta2 = k
	      delta_iso2  = min(delta2a,delta2b)
	    endif 
	  enddo
        enddo
      enddo
c
      if (delta_iso1 .lt. p5) then
        s_minus = (c1 - sqrt(c1 - 4.0*delta_iso1**2))/(c2*delta_iso1)
        s_plus  = c1/s_minus
      else
        s_minus = c0
        s_plus  = c0
      endif
c
      write(stdout,'(a)')
      write(stdout,'(a,e14.7)')
     &'The diffusion grid factor delta_iso1 =',delta_iso1
      write(stdout,'(a)')
     &'was determined at the grid point'
      write(stdout,'(a,i4,a,e14.7)')
     &'dxt(',i_delta1,') = ',dxt(i_delta1)
      write(stdout,'(a,i4,a,e14.7)')
     &'dyt(',j_delta1,') = ',dyt(j_delta1)
      write(stdout,'(a,i4,a,e14.7)')
     &'dzt(',k_delta1,') = ',dzt(k_delta1)
      write(stdout,'(a,i4,a,e14.7)')
     &'cst(',j_delta1,') = ',cst(j_delta1)
      write(stdout,*)
     &'Without latitudinal filtering, delta_iso1 is the steepest'
      write(stdout,*)
     &'isoneutral slope available for linear stab of Redi and GM.'
c
# if defined biharmonic_rm
      write(stdout,'(a,e14.7)')
     &'The RM biharmonic grid factor delta_iso2 =',delta_iso2
      write(stdout,'(a)')
     &'was determined at the grid point'
      write(stdout,'(a,i4,a,e14.7)')
     &'dxt(',i_delta2,') = ',dxt(i_delta2)
      write(stdout,'(a,i4,a,e14.7)')
     &'dyt(',j_delta2,') = ',dyt(j_delta2)
      write(stdout,'(a,i4,a,e14.7)')
     &'dzt(',k_delta2,') = ',dzt(k_delta2)
      write(stdout,'(a,i4,a,e14.7)')
     &'cst(',j_delta2,') = ',cst(j_delta2)
      write(stdout,*)
     &'Without latitudinal filtering, delta_iso2 is the steepest'
      write(stdout,*)
     &'isoneutral slope available for linear stab of RM biharmonic.'
# endif
c
# if !defined full_tensor
      write(stdout,'(/a,e14.7/)')
     & 'Maximum allowable isoneutral slope is specified as slmx = ',slmx
# endif
# ifdef full_tensor
      if (delta_iso1 .lt. p5) then
        write(stdout,'(/a,a,/2(a,e14.7)/)')
     & ' The full isoneutral tensor will be rescaled'
     &,' when the slope is in the range:'
     &,' s(-) = ',s_minus,', s_(+) = ',s_plus
      else
        write(stdout,'(/a,a/)') 'Isoneutral slopes in the full tensor'
     &,                  'will not be rescaled since delta_iso1 > 0.5'
      endif
# endif
c
c-----------------------------------------------------------------------
c     store the square root of the tracer timestep acceleration values
c     into variable "dtxsqr" for use in isoneutral mixing
c-----------------------------------------------------------------------
c
      do k=1,km
       dtxsqr(k) = sqrt(dtxcel(k))
      enddo
c
      write (stdout,'(a)') ' Acceleration with depth "dtxsqr(k)"='
      write (stdout,'(5(1x,e12.6))') (dtxsqr(k),k=1,km)
c
c-----------------------------------------------------------------------
c     vertical structure function for the isoneutral diffusion coeff.
c-----------------------------------------------------------------------
c
# ifdef bryan_lewis_horizontal
c
c     the depth dependent diffusion coefficients from blmix.F
c     are used for ahisop.  Note that this overrides ahisop
c     from namelist.isopyc. 
c
      do k=1,kbot
        fzisop(k) = Ahhisop(k)
      enddo
#  ifdef bbl_ag
      fzisop(km) = c0
#  endif
c
# else
      do k=1,kbot
        fzisop(k) = c1
      enddo
#  ifdef bbl_ag
      fzisop(km) = c0
#  endif
# endif
c
      write (stdout,'(a)') '  '
      write (stdout,'(a)') ' Vertical structure function "fzisop(k)"='
      write (stdout,'(5(1x,e12.6))') (fzisop(k),k=1,km)
c
# ifdef bryan_lewis_horizontal 
c
      write(stdout,9110) (ahisop*fzisop(k),k=1,km)
      write(stdout,9111) (athkdf*fzisop(k),k=1,km)
      write(stdout,9112) ah,am
      call getunit (iodoc, 'document.dta'
     &,               'formatted sequential append')
      write (iodoc,'(a,e14.7)') 'am=', am, 'ah=',ah
      write (iodoc,'(a,5(2x,e14.7))') 'ahisop(k)=',ahisop
      write (iodoc,'(a,5(2x,e14.7))') 'athkdf(k)=',athkdf
      write (iodoc,'(a,e14.7)') 'slmxr=',slmxr
      call relunit (iodoc)
c
9110  format(/' depth-varying isoneutral mixing coeff'
     *,       ' ahisop (cm**2/sec) = ',5(2x,e12.6))
9111  format(/' depth-varying GM mixing coeff'
     *,       ' athkdf (cm**2/sec) = ',5(2x,e12.6))
9112  format(/' ah (cm**2/sec) = ',e12.6,' am (cm**2/sec) = ',e12.6)
c
# endif
c
      write(stdout,'(a)')
     &'---------------------------------------------------------------'
c
c-----------------------------------------------------------------------
c     initialize arrays
c-----------------------------------------------------------------------
c
      do j=1,jmw
	do k=1,km
          do i=1,imt
            drho(i,k,j,1)  = c0
            drho(i,k,j,2)  = c0
          enddo
        enddo
      enddo
      do j=2,jmw-1
	do n=1,2
	  do k=1,km
            do i=1,imt
              ddxt(i,k,j,n) = c0
            enddo
          enddo
        enddo
      enddo
      do j=1,jmw-1
	do n=1,2
	  do k=1,km
            do i=1,imt
              ddyt(i,k,j,n) = c0
            enddo
          enddo
        enddo
      enddo
      do j=1,jmw
	do n=1,2
	  do k=0,km
            do i=1,imt
              ddzt(i,k,j,n) = c0
            enddo
          enddo
        enddo
      enddo
      do j=2,jmw-1
	do k=1,km
          do i=1,imt
            do kr=0,1
              do ip=0,1
# if defined redi_diffusion || defined gent_mcwilliams
                Ai_ez(i,k,j,ip,kr)  = c0
# endif
# ifdef biharmonic_rm
                Bi_ez(i,k,j,ip,kr)  = c0
# endif
	      enddo
	    enddo
            do kr=0,1
              do ip=0,1
# if defined redi_diffusion || defined gent_mcwilliams
                Ai_bx(i,k,j,ip,kr) = c0
# endif
# ifdef biharmonic_rm
                Bi_bx(i,k,j,ip,kr) = c0
# endif
	      enddo
	    enddo
            do kr=0,1
              do jq=0,1
# if defined redi_diffusion || defined gent_mcwilliams
                Ai_by(i,k,j,jq,kr) = c0
# endif
# ifdef biharmonic_rm
                Bi_by(i,k,j,jq,kr) = c0
# endif
	      enddo
	    enddo
# ifdef redi_diffusion
            K33(i,k,j) = c0
            K11(i,k,j) = c0
# endif
# ifdef gent_mcwilliams
            adv_vetiso(i,k,j) = c0
# endif
          enddo
        enddo
      enddo
      do j=1,jmw-1
	do k=1,km
          do i=1,imt
            do kr=0,1
              do jq=0,1
# if defined redi_diffusion || defined gent_mcwilliams
                Ai_nz(i,k,j,jq,kr)  = c0
# endif
# ifdef biharmonic_rm
                Bi_nz(i,k,j,jq,kr)  = c0
# endif
	      enddo
	    enddo
# ifdef redi_diffusion
            K22(i,k,j)        = c0
# endif
# ifdef gent_mcwilliams
            adv_vntiso(i,k,j) = c0
# endif
          enddo
        enddo
      enddo
# ifdef gent_mcwilliams
      do j=2,jmw-1
        do k=0,km
          do i=1,imt
            adv_vbtiso(i,k,j) = c0
            adv_fbiso(i,k,j)  = c0
	  enddo
        enddo
      enddo
# endif
c
      do nn=1,6
        redi_gm(nn) = c0
      enddo
# if defined small_tensor && !defined gm_skew
      redi_gm(1) = c1
      redi_gm(2) = c0
# endif
# if !defined small_tensor && defined gm_skew
      redi_gm(1) = c0
      redi_gm(2) = c1
# endif
# if defined small_tensor && defined gm_skew
      redi_gm(1) = c1
      redi_gm(2) = c1
# endif
      if(athkdf .eq. ahisop) then
        redi_gm(3) = c0
        redi_gm(4) = c1
      else
        redi_gm(3) = c1
        redi_gm(4) = c0
      endif
      if(athkdf .eq. ahsteep) then
        redi_gm(5) = c0
        redi_gm(6) = c1
      else
        redi_gm(5) = c1
        redi_gm(6) = c0
      endif
      write(stdout,*)' '
      write(*,*)'The Redi/GM algorithm coefficients are:'
      do nn=1,6
        write(stdout,*)'redi_gm(',nn,') = ',redi_gm(nn)
      enddo
c
# if defined nonconst_diffusivity
      if(init) then 
# endif
        do j=jstask,jetask
          do i=1,imt          
            agm(i,j) = athkdf
          enddo
        enddo
# if defined nonconst_diffusivity
      endif
# endif
c
      return
      end

      subroutine elements (joff, js, je, is, ie)
c
c=======================================================================
c     Estimate drho and normal gradients on faces of T cells
c=======================================================================
c
      use arrays_2d
      use memory_window
# include "param.h"
# include "accel.h"
# include "grdvar.h"
# include "isopyc.h"
# include "mw.h"
# include "state.h"
# include "dens.h"
# include "fdift.h"
c
# ifdef timing
      call tic ('isoneutral ', 'elements')
# endif
c
c-----------------------------------------------------------------------
c     drho_dt and drho_ds at centers of T cells
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=1,kbot
          do i=is,ie
	    tprime = t(i,k,j,1,taum1)-to(k)
	    sprime = t(i,k,j,2,taum1)-so(k)
            drho(i,k,j,1)  = drhodt(tprime, sprime, k)
	    drho(i,k,j,2)  = drhods(tprime, sprime, k)
          enddo
        enddo
        call setbcx (drho(1,1,j,1),  imt, km)
        call setbcx (drho(1,1,j,2),  imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     gradients at bottom face of T cells
c-----------------------------------------------------------------------
c
      do j=js,je
        do n=1,2
          do k=1,kbot
            kp1 = min(k+1,kbot)
            do i=is,ie
	      ddzt(i,k,j,n) = tmask(i,kp1,j)*dz_wtr(i,k,j)*
     &                        (t(i,k,j,n,taum1) - t(i,kp1,j,n,taum1))
            enddo
          enddo
          do i=is,ie
	    ddzt(i,0,j,n) = c0
	  enddo
          call setbcx (ddzt(1,0,j,n), imt, km+1)
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     gradients at eastern face of T cells
c-----------------------------------------------------------------------
c
# if defined full_tensor
      do j=js,je
# else
      do j=max(js-1,1+jbuf),je-1
# endif
        jrow = j + joff
        do n=1,2
          do k=1,kbot
            do i=is,ie
	      ddxt(i,k,j,n) = tmask(i,k,j)*tmask(i+1,k,j)*cstr(jrow)*
     &              dxur(i)*(T_i(i,k,j,n,1) - T_i(i,k,j,n,0))
            enddo
          enddo
          call setbcx (ddxt(1,1,j,n), imt, km)
	enddo
      enddo
c
c-----------------------------------------------------------------------
c     gradients at northern face of T cells
c-----------------------------------------------------------------------
c
      do j=max(js-1,jbuf),je-1
        jrow = j + joff
        do n=1,2
          do k=1,kbot
            kp1 = min(k+1,kbot)
            do i=is,ie
	      ddyt(i,k,j,n) = tmask(i,k,j)*tmask(i,k,j+1)*dyur(jrow)*
     &                 (T_j(i,k,j,n,1) - T_j(i,k,j,n,0))  
            enddo
          enddo
          call setbcx (ddyt(1,1,j,n), imt, km)
	enddo
      enddo
# ifdef timing
      call toc ('isoneutral ', 'elements')
# endif
      return
      end




      subroutine isopyc (joff, js, je, is, ie)
c
c=======================================================================
c
c     Compute isoneutral mixing coefficients.
c
c     input:
c       joff = offset relating row "j" in the MW to latitude "jrow"
c       js   = starting row within the MW for calculations
c       je   = ending row within the MW for calculations
c       is   = starting index longitude within the MW
c       ie   = ending index longitude within the MW
c
c     output:
# if defined redi_diffusion || defined gent_mcwilliams
c       Ai_ez = diffusivity*tapered slope
c               centered on east face of T cells
c       Ai_nz = diffusivity*tapered slope
c               centered on north face of T cells
c       Ai_bx = diffusivity*tapered slope
c               centered on bottom face of T cells
c       Ai_by = diffusivity*tapered slope
c               centered on bottom face of T cells
# endif
c            
# ifdef gent_mcwilliams
c       adv_vetiso = GM advective vel on east face of T cells
c       adv_vntiso = GM advective vel on north face of T cells
c               (Note: this includes the cosine factor as in "adv_vnt")
c       adv_vbtiso = GM advective vel on bottom face of T cells       
# endif
c
# ifdef biharmonic_rm
c     output:
c       Bi_ez = biharmonic diffusivity*tapered slope
c               centered on east face of T cells
c       Bi_nz = biharmonic diffusivity*tapered slope
c               centered on north face of T cells
c       Bi_bx = biharmonic diffusivity*tapered slope
c               centered on bottom face of T cells
c       Bi_by = biharmonic diffusivity*tapered slope
c               centered on bottom face of T cells
# endif
c
c=======================================================================
c
# include "param.h"
# include "switch.h"
c
# ifdef timing
      call tic ('isoneutral ', 'isopyc')
# endif
c
c-----------------------------------------------------------------------
c     set local constants 
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
# ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     & "=> In isopyc: from js=",js," je=",je," joff=",joff
     &," jrows=",js+joff," to ",je+joff
# endif
c
c-----------------------------------------------------------------------
c     estimate drho and gradients on sides of T cells 
c-----------------------------------------------------------------------
c
      call elements (joff, js, je, is, ie)
c
# if defined nonconst_diffusivity
c-----------------------------------------------------------------------
c     compute nonconstant mixing coefficients  
c-----------------------------------------------------------------------
c
      if(diffts) then
#  ifdef hl_diffusivity 
        call hldiffusivity (joff, js, je, is, ie)
#  endif
#  ifdef vmhs_diffusivity 
        call vmhsdiffusivity (joff, js, je, is, ie)
#  endif
        call ncdiffusivity (joff, js, je, is, ie)
      endif
c
# endif
c
c-----------------------------------------------------------------------
c     compute Ai_ez and Bi_ez centered on eastern face of T cells
c-----------------------------------------------------------------------
c
      call ai_east (joff, max(js-1,1+jbuf), je-1, is, ie)
c
c-----------------------------------------------------------------------
c     compute Ai_nz and Bi_nz centered on the northern face of T cells
c-----------------------------------------------------------------------
c
      call ai_north (joff, max(js-1,jbuf), je-1, is, ie)
c
c-----------------------------------------------------------------------
c     evaluate Ai_bx, Ai_by, Bi_bx, and Bi_by
c     centered on bottom face of T cells
c-----------------------------------------------------------------------
c
      call ai_bottom (joff, max(js-1,1+jbuf), je-1, is, ie)
c
# ifdef gent_mcwilliams
c
c-----------------------------------------------------------------------
c     compute GM advective velocities for tracers
c-----------------------------------------------------------------------
c
#  if defined gm_advect
        call isopyc_adv (joff, max(js-1,jbuf), je-1, is, ie)
#  endif
#  if defined gm_skew
      if (snapts .or. vmsfts) then
        call isopyc_adv (joff, max(js-1,jbuf), je-1, is, ie)
      endif
#  endif
c
# endif
c
# ifdef timing
      call toc ('isoneutral ', 'isopyc')
# endif
      return
      end


      subroutine ai_east (joff, js, je, is, ie)
c
c=======================================================================
c     Compute Ai_ez, Bi_ez, and K11 on center of east face of T cell.
c     Note re-scaling factor which reduces diffusivity where abs slope
c     "sxe"  exceeds the critical slope "sc" for the small slope approx.
c     For the full tensor, diffusivity is re-scaled if outside its
c     stable range.
c=======================================================================
c
      use memory_window
# include "param.h"
# include "accel.h"
# include "coord.h"
# include "grdvar.h"
# include "hmixc.h"
# include "isopyc.h"
# include "mw.h"
# include "fdift.h"
# include "isopyc_sf.h"
# if defined nonconst_diffusivity
#  include "switch.h"
# endif
c     
c
      do j=js,je
        jrow = j + joff
        do k=1,kbot
          sc = c1/(slmxr*dtxsqr(k))
	  dzt4r = p5*dzt2r(k)
          kp2 = min(kbot,k+2) 
          do i=2,imtm1
c           
            Ai11 = (redi_gm(3)*ahisop+redi_gm(4)*agm(i,jrow))
c
# if defined small_tensor && defined gm_skew && !defined biharmonic_rm
            if(athkdf.eq.ahisop .and. athkdf.eq.ahsteep) then
    	      K11(i,k,j) = Ai11*tmask(i,k,j)*tmask(i+1,k,j)*fzisop(k)
#  if defined nonconst_diffusivity 
              if(diffts) then
                S2N2 = c0 
                do kr=0,1
                  do ip=0,1
                    S2N2 = S2N2
     &              + drodxe(i,k,j,ip)**2/abs(drodze(i,k,j,ip,kr)-epsln)
	          enddo
	        enddo
                growth2b(i,jrow) = growth2b(i,jrow)
     &          + fznoncst(k)*S2N2*tmask(i,kp2,j)*tmask(i+1,kp2,j)
              endif
#  endif
            else
# endif
c
              Asteep = (redi_gm(5)*ahsteep+redi_gm(6)*agm(i,jrow))
     &                  *fzisop(k)
              Ai13   = (redi_gm(1)*Ai11-redi_gm(2)*agm(i,jrow))
     &                  *fzisop(k)
              Aibih  = abihrm*fzisop(k)
c             
	      sumz = c0
              S2N2 = c0 
              do kr=0,1
                do ip=0,1
                  sxe    = -drodxe(i,k,j,ip)/(drodze(i,k,j,ip,kr)-epsln)
                  abssxe = abs(sxe)
# ifdef dm_taper 
                  fact1 = (abssxe-del_dm)*s_dmr   
                  fact2 = sign(1.0,fact1)*min(19.0,abs(fact1)) 
                  taper = p5*(c1-tanh(fact2))   
# endif          
# ifdef gkw_taper
                  if (abssxe .gt. sc) then
                    taper = (sc/(abssxe + epsln))**2
                  else
                    taper = c1
                  endif
# endif          
# if defined small_tensor
	          sumz = sumz +
#  ifdef partial_cell
     &                 min(delqc(i,k,j,kr),delqc(i+1,k,j,kr))
#  else          
     &                 dzw(k-1+kr)
#  endif         
     &                 *tmask(i,k,j)*tmask(i+1,k,j)
     &                 *max(Asteep,Ai11*taper)
# endif          
# if defined small_tensor || defined gm_skew
	          Ai_ez(i,k,j,ip,kr) = 
     &                 Ai13*taper*sxe*tmask(i,k,j)*tmask(i+1,k,j)
# endif          
# ifdef biharmonic_rm
	          Bi_ez(i,k,j,ip,kr) =
     &                 Aibih*taper*sxe*tmask(i,k,j)*tmask(i+1,k,j)
# endif          
# if defined nonconst_diffusivity
                  S2N2   = S2N2
     &                   + abs(drodze(i,k,j,ip,kr)-epsln)*sxe**2 
# endif          
	        enddo
	      enddo
c
# if defined small_tensor
#  ifdef partial_cell
              K11(i,k,j) = p25*sumz
#  else      
              K11(i,k,j) = dzt4r*sumz
#  endif     
# endif      
# if defined nonconst_diffusivity
              if(diffts) then
                growth2b(i,jrow) = growth2b(i,jrow)
     &          + fznoncst(k)*S2N2*tmask(i,kp2,j)*tmask(i+1,kp2,j)
              endif
# endif      
c             
c             
# if defined small_tensor && defined gm_skew && !defined biharmonic_rm 
            endif
# endif     
c           
c
# ifdef full_tensor
	    sumz = c0
            do kr=0,1
              do ip=0,1
                sxe = abs(drodxe(i,k,j,ip)/(drodze(i,k,j,ip,kr)+epsln))
                if (delta_iso1 .lt. p5 .and.
     &               (s_minus .lt. sxe .and. sxe .lt. s_plus)) then
	          Ai_ez(i,k,j,ip,kr) = Ai11*tmask(i,k,j)*tmask(i+1,k,j)
     &                 *delta_iso1*(sxe + c1/sxe)
	        else
	          Ai_ez(i,k,j,ip,kr) = Ai11*tmask(i,k,j)*tmask(i+1,k,j)
	        endif
	        sumz = sumz + dzw(k-1+kr)*Ai_ez(i,k,j,ip,kr)
     &               *drodze(i,k,j,ip,kr)**2/(drodxe(i,k,j,ip)**2
     &               +  p5*(drodye(i,k,j,ip,0)**2
     &               + drodye(i,k,j,ip,1)**2)
     &               + drodze(i,k,j,ip,kr)**2 + epsln)
              enddo
            enddo
c           
	    sumy = c0
            do jq=0,1
	      facty = csu(jrow-1+jq)*dyu(jrow-1+jq)
	      do ip=0,1
	        sumy = sumy + facty*Ai11*tmask(i,k,j)*tmask(i+1,k,j)
     &               *drodye(i,k,j,ip,jq)**2/(drodxe(i,k,j,ip)**2
     &               + drodye(i,k,j,ip,jq)**2 + epsln
     &               +  p5*(drodze(i,k,j,ip,0)**2
     &               +  drodze(i,k,j,ip,1)**2))
	      enddo
	    enddo
	    K11(i,k,j) = dzt4r*sumz + p25*cstdytr(jrow)*sumy
# endif     
c           
          enddo
        enddo
c       
# if defined redi_diffusion || defined gm_skew
        call setbcx (Ai_ez(1,1,j,0,0), imt, km)
        call setbcx (Ai_ez(1,1,j,1,0), imt, km)
        call setbcx (Ai_ez(1,1,j,0,1), imt, km)
        call setbcx (Ai_ez(1,1,j,1,1), imt, km)
# endif 
# if defined redi_diffusion
        call setbcx (K11(1,1,j), imt, km)
# endif 
# ifdef biharmonic_rm
        call setbcx (Bi_ez(1,1,j,0,0), imt, km)
        call setbcx (Bi_ez(1,1,j,1,0), imt, km)
        call setbcx (Bi_ez(1,1,j,0,1), imt, km)
        call setbcx (Bi_ez(1,1,j,1,1), imt, km)
# endif
      enddo
c     
# ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     &     "=> In ai_east: js=",js," je=",je," joff=",joff
     &     ," jrows=",js+joff," to ",je+joff
# endif
c
      return
      end


      subroutine ai_north (joff, js, je, is, ie)
c     
c=======================================================================
c     Compute "Ai_nz" and "Bi_nz" on center of north face of T cell.
c     Note re-scaling factor which reduces diffusivity where abs slope
c     "syn"  exceeds the critical slope "sc" for the small slope approx.
c     For the full tensor, diffusivity is re-scaled if outside its
c     stable range.
c=======================================================================
c     
      use memory_window
# include "param.h"
# include "accel.h"
# include "coord.h"
# include "grdvar.h"
# include "hmixc.h"
# include "isopyc.h"
# include "mw.h"
# include "fdift.h"
# include "isopyc_sf.h"
# if defined nonconst_diffusivity
#  include "switch.h"
# endif
c     
c     
      do j=js,je
        jrow = j + joff
        do k=1,kbot
          sc = c1/(slmxr*dtxsqr(k))
	  dzt4r = p5*dzt2r(k)
          kp2 = min(kbot,k+2) 
          do i=2,imtm1
c           
            Ai22 = (redi_gm(3)*ahisop+redi_gm(4)*agm(i,jrow))
c
# if defined small_tensor && defined gm_skew && !defined biharmonic_rm
            if(athkdf.eq.ahisop .and. athkdf.eq.ahsteep) then
    	      K22(i,k,j) = Ai22*tmask(i,k,j)*tmask(i,k,j+1)*fzisop(k)
#  if defined nonconst_diffusivity 
              if(diffts) then
                S2N2 = c0 
                do kr=0,1
                  do jq=0,1
                    S2N2 = S2N2
     &              + drodyn(i,k,j,jq)**2/abs(drodzn(i,k,j,jq,kr)-epsln)
	          enddo
	        enddo
                growth2b(i,jrow) = growth2b(i,jrow)
     &          + fznoncst(k)*S2N2*tmask(i,kp2,j)*tmask(i,kp2,j+1)
              endif
#  endif
            else
# endif       
c
              Asteep = (redi_gm(5)*ahsteep+redi_gm(6)*agm(i,jrow))
     &                  *fzisop(k)
              Ai23   = (redi_gm(1)*Ai22-redi_gm(2)*agm(i,jrow))
     &                  *fzisop(k)
              Aibih  = abihrm*fzisop(k)
c
  	      sumz = c0
              S2N2 = c0 
              do kr=0,1
                do jq=0,1
                  syn    = -drodyn(i,k,j,jq)/(drodzn(i,k,j,jq,kr)-epsln)
                  abssyn = abs(syn)
# ifdef dm_taper 
                  fact1 = (abssyn-del_dm)*s_dmr   
                  fact2 = sign(1.0,fact1)*min(19.0,abs(fact1)) 
                  taper = p5*(c1-tanh(fact2))   
# endif          
# ifdef gkw_taper
                  if (abssyn .gt. sc) then
                    taper = (sc/(abssyn + epsln))**2
                  else
                    taper = c1
                  endif
# endif          
# if defined small_tensor
  	          sumz = sumz +
#  ifdef partial_cell
     &                 min(delqc(i,k,j,kr),delqc(i,k,j+1,kr))
#  else           
     &                 dzw(k-1+kr)
#  endif          
     &                 *tmask(i,k,j)*tmask(i,k,j+1)
     &                 *max(Asteep,Ai22*taper)
# endif          
# if defined small_tensor || defined gm_skew
	          Ai_nz(i,k,j,jq,kr) = 
     &                 Ai23*taper*syn*tmask(i,k,j)*tmask(i,k,j+1)
# endif          
# ifdef biharmonic_rm
	          Bi_nz(i,k,j,jq,kr) =
     &                 Aibih*taper*syn*tmask(i,k,j)*tmask(i,k,j+1)
# endif          
# if defined nonconst_diffusivity
                  S2N2   = S2N2
     &                     + abs(drodzn(i,k,j,jq,kr)-epsln)*syn**2 
# endif          
	        enddo
	      enddo
c                
# if defined small_tensor
#  ifdef partial_cell
              K22(i,k,j) = p25*sumz
#  else      
              K22(i,k,j) = dzt4r*sumz
#  endif     
# endif      
# if defined nonconst_diffusivity
            if(diffts) then
               growth2b(i,jrow) = growth2b(i,jrow)
     &         + fznoncst(k)*S2N2*tmask(i,kp2,j)*tmask(i,kp2,j+1)
            endif
# endif  
c             
# if defined small_tensor && defined gm_skew && !defined biharmonic_rm
            endif
# endif     
c           
c           
# ifdef full_tensor
	    sumz = c0
            do kr=0,1
              do jq=0,1
                syn=abs(drodyn(i,k,j,jq)/(drodzn(i,k,j,jq,kr)+epsln))
                if (delta_iso1 .lt. p5 .and.
     &               (s_minus .lt. syn .and. syn .lt. s_plus)) then
	          Ai_nz(i,k,j,jq,kr) = Ai22*tmask(i,k,j)*tmask(i,k,j+1)
     &                 *delta_iso1*(syn + c1/syn)
	        else
	          Ai_nz(i,k,j,jq,kr) = Ai22*tmask(i,k,j)*tmask(i,k,j+1)
	        endif
	        sumz = sumz + dzw(k-1+kr)*Ai_nz(i,k,j,jq,kr)
     &               *drodzn(i,k,j,jq,kr)**2/(
     &               +  p5*(drodxn(i,k,j,0,jq)**2
     &               + drodxn(i,k,j,1,jq)**2)
     &               + drodyn(i,k,j,jq)**2
     &               + drodzn(i,k,j,jq,kr)**2 + epsln)
              enddo
            enddo
c           
	    sumx = c0
	    do ip=0,1
              do jq=0,1
                sumx=sumx+dxu(i-1+ip)*Ai22*tmask(i,k,j)*tmask(i,k,j+1)
     &               *drodxn(i,k,j,ip,jq)**2/(drodxn(i,k,j,ip,jq)**2
     &               + drodyn(i,k,j,jq)**2 + epsln
     &               +    p5*(drodzn(i,k,j,jq,0)**2
     &               +    drodzn(i,k,j,jq,1)**2))
	      enddo
	    enddo
	    K22(i,k,j) = dzt4r*sumz + p25*dxtr(i)*sumx
# endif     
c           
          enddo
        enddo
c       
# if defined redi_diffusion || defined gm_skew
        call setbcx (Ai_nz(1,1,j,0,0), imt, km)
        call setbcx (Ai_nz(1,1,j,1,0), imt, km)
        call setbcx (Ai_nz(1,1,j,0,1), imt, km)
        call setbcx (Ai_nz(1,1,j,1,1), imt, km)
# endif 
# if defined redi_diffusion
        call setbcx (K22(1,1,j), imt, km)
# endif 
# ifdef biharmonic_rm
        call setbcx (Bi_nz(1,1,j,0,0), imt, km)
        call setbcx (Bi_nz(1,1,j,1,0), imt, km)
        call setbcx (Bi_nz(1,1,j,0,1), imt, km)
        call setbcx (Bi_nz(1,1,j,1,1), imt, km)
# endif
      enddo
c     
# ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     &     "=> In ai_north: js=",js," je=",je," joff=",joff
     &     ," jrows=",js+joff," to ",je+joff
# endif
      return
      end



      subroutine ai_bottom (joff, js, je, is, ie)
c
c=======================================================================
c     compute Ai_bx, Ai_by, Bi_bx, Bi_by, and K33 at the center
c     of the bottom face of T cells.
c=======================================================================
c
      use memory_window
# include "param.h"
# include "accel.h"
# include "grdvar.h"
# include "hmixc.h"
# include "isopyc.h"
# include "mw.h"
# include "fdift.h"
# include "isopyc_sf.h"
# if defined nonconst_diffusivity
#  include "switch.h"
# endif
c
c-----------------------------------------------------------------------
c     compute "Ai_bx", "Ai_by", "Bi_bx", "Bi_by, & K33 on bottom face
c     of T cell. Note re-scaling factor to reduce mixing coefficient
c     "Ai" where slopes excede the critical slope "sc" for the small
c     slope approx, gm_skew, or biharmonic_rm.  
c     The full tensor is re-scaled if outside its stable range.
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        do k=1,kbot-1
          sc = c1/(slmxr*dtxsqr(k))
	  kp1   = min(k+1,kbot)
	  kp2   = min(k+2,kbot)
	  fzisopb = p5*(fzisop(k+1)+fzisop(k))
          do i=2,imtm1
c
            Ai31  = fzisopb
     &       *(  redi_gm(1)*(redi_gm(3)*ahisop+redi_gm(4)*agm(i,jrow))
     &         + redi_gm(2)*agm(i,jrow))
            Ai32  = Ai31
            Ai33  = (redi_gm(3)*ahisop+redi_gm(4)*agm(i,jrow))*fzisopb
            Aibih = abihrm*fzisopb
c
# if defined small_tensor || defined gm_skew || defined biharmonic_rm
c
c           eastward slopes at the base of T cells
c
            sumx  = c0
            S2N2x = c0
            do ip=0,1
              do kr=0,1
                sxb = -drodxb(i,k,j,ip,kr)/(drodzb(i,k,j,kr)-epsln)
                abssxb = abs(sxb)
#  ifdef dm_taper
                fact1 = (abssxb-del_dm)*s_dmr 
                fact2 = sign(1.0,fact1)*min(19.0,abs(fact1))
                taper = p5*(c1-tanh(fact2))   
#  endif
#  ifdef gkw_taper 
                if (abssxb .gt. sc) then
                  taper = (sc/(abssxb + epsln))**2
                else
                  taper = c1
                endif
#  endif
#  if defined small_tensor
	        sumx = sumx
     &          + dxu(i-1+ip)*tmask(i,k+1,j)*Ai33*taper*sxb**2
#   ifdef partial_cell
     &          *min(delqc(i-1+ip,k+kr,j,1-kr),delqc(i+ip,k+kr,j,1-kr))
#   endif
#  endif
#  if defined small_tensor || defined gm_skew   
	        Ai_bx(i,k,j,ip,kr) =
     &               Ai31*taper*sxb*tmask(i,k+1,j)
#  endif
#  ifdef biharmonic_rm
	        Bi_bx(i,k,j,ip,kr) =
     &               Aibih*taper*sxb*tmask(i,k+1,j)
#  endif
#  if defined nonconst_diffusivity
                  S2N2x  = S2N2x
     &             + abs(drodzb(i,k,j,kr)-epsln)*sxb**2 
#  endif    
c
              enddo
            enddo
c
c           northward slopes at the base of T cells
c
            sumy  = c0
            S2N2y = c0
            do jq=0,1
	      facty = csu(jrow-1+jq)*dyu(jrow-1+jq)
              do kr=0,1
                syb = -drodyb(i,k,j,jq,kr)/(drodzb(i,k,j,kr)-epsln)
                abssyb = abs(syb)
#  ifdef dm_taper
                fact1 = (abssyb-del_dm)*s_dmr
                fact2 = sign(1.0,fact1)*min(19.0,abs(fact1))
                taper =  p5*(c1-tanh(fact2))  
#  endif
#  ifdef gkw_taper 
                if (abssyb .gt. sc) then
                  taper = (sc/(abssyb + epsln))**2
                else
                  taper = c1
                endif
#  endif
#  if defined small_tensor
	        sumy = sumy
     &          + facty*tmask(i,k+1,j)*Ai33*taper*syb**2
#   ifdef partial_cell
     &          *min(delqc(i,k+kr,j-1+jq,1-kr),delqc(i,k+kr,j+jq,1-kr))
#   endif
#  endif
#  if defined small_tensor || defined gm_skew   
	        Ai_by(i,k,j,jq,kr) =
     &               Ai32*taper*syb*tmask(i,k+1,j)
#  endif
#  ifdef biharmonic_rm
	        Bi_by(i,k,j,jq,kr) =
     &               Aibih*taper*syb*tmask(i,k+1,j)
#  endif
#  if defined nonconst_diffusivity
                  S2N2y  = S2N2y
     &               + abs(drodzb(i,k,j,kr)-epsln)*syb**2 
#  endif 
c
              enddo
            enddo
c
c
#  if defined small_tensor  
            K33(i,k,j) = (dxt4r(i)*sumx + dyt4r(jrow)*cstr(jrow)*sumy)
#   ifdef partial_cell
     &                   /dhwt(i,k,j)
#   endif
#  endif
#  if defined nonconst_diffusivity
            if(diffts) then
               growth2b(i,jrow) = growth2b(i,jrow)
     &           + fznoncst(k)*(S2N2x+S2N2y)*tmask(i,kp2,j)
            endif
#  endif   
c
# endif
c
c
# ifdef full_tensor
c
c           eastward slopes at the base of T cells
c
            sumx = c0
            do ip=0,1
              do kr=0,1
                sxb = abs(drodxb(i,k,j,ip,kr)/(drodzb(i,k,j,kr)+epsln))
                if (delta_iso1 .lt. p5 .and.
     &              (s_minus .lt. sxb .and. sxb .lt. s_plus)) then
	          Ai_bx(i,k,j,ip,kr)  =  Ai33*tmask(i,k+1,j)
     &                                   *delta_iso1*(sxb + c1/sxb)
	        else
	          Ai_bx(i,k,j,ip,kr)  =  Ai33*tmask(i,k+1,j)
	        endif
                sumx = sumx + dxu(i-1+ip)*Ai_bx(i,k,j,ip,kr)
     &            *drodxb(i,k,j,ip,kr)**2/(drodxb(i,k,j,ip,kr)**2
     &       +  p5*(drodyb(i,k,j,0,kr)**2 + drodyb(i,k,j,1,kr)**2)
     &               + drodzb(i,k,j,kr)**2 + epsln)
              enddo
            enddo
c
c           northward slopes at the base of T cells
c
            sumy = c0
            jrow = j + joff
            do jq=0,1
	      facty = csu(jrow-1+jq)*dyu(jrow-1+jq)
              do kr=0,1
                syb = abs(drodyb(i,k,j,jq,kr)/(drodzb(i,k,j,kr)+epsln))
                if (delta_iso1 .lt. p5 .and.
     &              (s_minus .lt. syb .and. syb .lt. s_plus)) then
	          Ai_by(i,k,j,jq,kr)  =  Ai33*tmask(i,k+1,j)
     &                                   *delta_iso1*(syb + c1/syb)
	        else
	          Ai_by(i,k,j,jq,kr)  =  Ai33*tmask(i,k+1,j)
	        endif
	        sumy = sumy + facty*Ai_by(i,k,j,jq,kr)
     &            *drodyb(i,k,j,jq,kr)**2/(
     &           p5*(drodxb(i,k,j,0,kr)**2 + drodxb(i,k,j,1,kr)**2)
     &            + drodyb(i,k,j,jq,kr)**2
     &            + drodzb(i,k,j,kr)**2 + epsln)
              enddo
            enddo
c
            K33(i,k,j) = dxt4r(i)*sumx + dyt4r(jrow)*cstr(jrow)*sumy 
# endif
c
          enddo
        enddo
c
# if defined redi_diffusion || defined gm_skew
        call setbcx (Ai_bx(1,1,j,1,0), imt, km)
        call setbcx (Ai_bx(1,1,j,0,0), imt, km)
        call setbcx (Ai_bx(1,1,j,1,1), imt, km)
        call setbcx (Ai_bx(1,1,j,0,1), imt, km)
        call setbcx (Ai_by(1,1,j,1,0), imt, km)
        call setbcx (Ai_by(1,1,j,0,0), imt, km)
        call setbcx (Ai_by(1,1,j,1,1), imt, km)
        call setbcx (Ai_by(1,1,j,0,1), imt, km)
# endif
# if defined redi_diffusion
        call setbcx (K33(1,1,j), imt, km)
# endif
# ifdef biharmonic_rm
        call setbcx (Bi_bx(1,1,j,1,0), imt, km)
        call setbcx (Bi_bx(1,1,j,0,0), imt, km)
        call setbcx (Bi_bx(1,1,j,1,1), imt, km)
        call setbcx (Bi_bx(1,1,j,0,1), imt, km)
        call setbcx (Bi_by(1,1,j,1,0), imt, km)
        call setbcx (Bi_by(1,1,j,0,0), imt, km)
        call setbcx (Bi_by(1,1,j,1,1), imt, km)
        call setbcx (Bi_by(1,1,j,0,1), imt, km)
# endif
      enddo
c
# if defined nonconst_diffusivity
      if(diffts) then
        if(jrow .eq. jmt-1) then
          do j=1,jmt
            do i=1,imt
              growth2a(i,j) = min(growth2max,growth2b(i,j))
              growth2b(i,j) = c0
            enddo
          enddo
          call fir2d(growth2a,numfltrgrth,imt,jmt)   
          do j=1,jmt
            call setbcx (growth2a(1,j), imt, 1)
          enddo
        endif
      endif
# endif
c
# ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     & "=> In ai_bottom: js=",js," je=",je," joff=",joff
     &," jrows=",js+joff," to ",je+joff
# endif
      return
      end               



      subroutine isoflux (joff, js, je, is, ie, n)
c
c=======================================================================
c     isoneutral diffusive tracer fluxes are computed.
c=======================================================================
c
      use arrays_2d
      use memory_window
# include "param.h"
# include "coord.h"
# include "grdvar.h"
# include "hmixc.h"
# include "isopyc.h"
# include "mw.h"
# include "vmixc.h"
# ifdef debug_tensor
#  include "tmngr.h"
# endif
# include "fdift.h"
# include "isopyc_sf.h"
# ifdef timing
      call tic ('isoneutral ', 'isoflux')
# endif
c
      numtra = min(n,nsmix)
c
c-----------------------------------------------------------------------
c     construct total isoneutral tracer flux at east face of "T" cells 
c-----------------------------------------------------------------------
c
# if defined small_tensor && defined gm_skew && !defined biharmonic_rm
      if (athkdf .ne. ahisop) then     
# endif
        do j=js,je
          jrow = j + joff
          do k=1,kbot
            km1 = max(k-1,1)
  	    dzt4r = p5*dzt2r(k)
            do i=2,imtm1
	      sumz    = c0
	      sumzBih = c0
              do kr=0,1
	        km1kr = max(k-1+kr,1)
	        kpkr  = min(k+kr,kbot)
                do ip=0,1
# if defined small_tensor || defined gm_skew 
	          sumz = sumz + Ai_ez(i,k,j,ip,kr)
#  ifdef partial_cell
     &                *min(delqc(i,k,j,kr),delqc(i+1,k,j,kr))
     &                *dz_wtr(i+ip,km1kr,j)
#  endif
     &             *(t(i+ip,km1kr,j,n,taum1)-t(i+ip,kpkr,j,n,taum1))
# endif
# ifdef full_tensor
	          sumz = sumz - Ai_ez(i,k,j,ip,kr)
     &                *(t(i+ip,km1kr,j,n,taum1)-t(i+ip,kpkr,j,n,taum1))
     &                *drodxe(i,k,j,ip)
     &                *drodze(i,k,j,ip,kr)/(drodxe(i,k,j,ip)**2
     &          +  p5*(drodye(i,k,j,ip,0)**2 + drodye(i,k,j,ip,1)**2)
     &               + drodze(i,k,j,ip,kr)**2 + epsln)
# endif
# ifdef biharmonic_rm
                  rm_laplacian = dxtr(i)*cstr(jrow)**2
     &           *(
     &              (Bi_ez(i+1,k,j,ip,kr)-Bi_ez(i,k,j,ip,kr))
     &                 *dxur(i)*tmask(i+1,k,j)
     &             -(Bi_ez(i,k,j,ip,kr)-Bi_ez(i-1,k,j,ip,kr))
     &                 *dxur(i-1)*tmask(i-1,k,j)
     &            )
     &            + dytr(jrow)*cstr(jrow)
     &           *(
     &             csu(jrow)*(Bi_ez(i,k,j+1,ip,kr)-Bi_ez(i,k,j,ip,kr))
     &                *dyur(jrow)*tmask(i,k,j+1)
     &            -csu(jrow-1)*(Bi_ez(i,k,j,ip,kr)-Bi_ez(i,k,j-1,ip,kr))
     &                *dyur(jrow-1)*tmask(i,k,j-1)
     &            )
	          sumzBih = sumzBih + rm_laplacian
#  ifdef partial_cell
     &                *min(delqc(i,k,j,kr),delqc(i+1,k,j,kr))
     &                *dz_wtr(i+ip,km1kr,j)
#  endif
     &                *(t(i+ip,km1kr,j,n,taum1)-t(i+ip,kpkr,j,n,taum1))
# endif
c
                enddo
	      enddo
# ifdef partial_cell
              flux_x = p25*(sumz+sumzBih)
# else
              flux_x = dzt4r*(sumz+sumzBih)
# endif
c
# ifdef full_tensor
              Ai0 = ahisop*fzisop(k)
	      sumy = c0
              do jq=0,1
	        facty = csu(jrow-1+jq)
	        do ip=0,1
	          sumy = sumy - facty*Ai0*tmask(i,k,j)*tmask(i+1,k,j)
     &       *(t(i+ip,k,j+jq,n,taum1)-t(i+ip,k,j-1+jq,n,taum1))
     &       *drodye(i,k,j,ip,jq)*drodxe(i,k,j,ip)/(drodxe(i,k,j,ip)**2
     &          + drodye(i,k,j,ip,jq)**2 + epsln
     &          +  p5*(drodze(i,k,j,ip,0)**2 + drodze(i,k,j,ip,1)**2))
                enddo
	      enddo
	      flux_x = flux_x + p25*cstdytr(jrow)*sumy
# endif
c
c
              diff_fe(i,k,j) = diff_fe(i,k,j) + flux_x
# ifdef redi_diffusion
     &        + cstdxur(i,j)*(T_i(i,k,j,n,1) - T_i(i,k,j,n,0))
     &             *(K11(i,k,j)
#  if defined isotropic_mixed 
     &        + p25*(  diff_cbt(i,k,j,numtra)     - K33(i,k,j)
     &               + diff_cbt(i,km1,j,numtra)   - K33(i,km1,j)
     &               + diff_cbt(i+1,k,j,numtra)   - K33(i+1,k,j)
     &               + diff_cbt(i+1,km1,j,numtra) - K33(i+1,km1,j)
     &              )
#  endif
     &        )
# endif
c
            enddo
          enddo
          call setbcx (diff_fe(1,1,j), imt, km)
        enddo
c
# if defined small_tensor && defined gm_skew && !defined biharmonic_rm
      else 
        do j=js,je
          do k=1,kbot
            do i=2,imtm1
              diff_fe(i,k,j) = diff_fe(i,k,j)
     &             + K11(i,k,j)*cstdxur(i,j)
     &             *(T_i(i,k,j,n,1) - T_i(i,k,j,n,0))
#  if defined partial_cell
     &             *dhte(i,k,j)
#  endif
            enddo
          enddo
          call setbcx (diff_fe(1,1,j), imt, km)
        enddo
      endif
# endif
c
c
c-----------------------------------------------------------------------
c     construct total isoneutral tracer flux at north face of "T" cells 
c-----------------------------------------------------------------------
c
c    
# if defined small_tensor && defined gm_skew && !defined biharmonic_rm
      if (athkdf .ne. ahisop) then     
# endif
        do j=js-1,je
          jrow = j + joff
          jm1 = max(1,jrow-1)
          do k=1,kbot
	    csu_dzt4r = csu(jrow)*p5*dzt2r(k)
            do i=2,imtm1
	      sumz    = c0
              sumzBih = c0
              do kr=0,1
	        km1kr = max(k-1+kr,1)
	        kpkr = min(k+kr,kbot)
                do jq=0,1
# if defined small_tensor || defined gm_skew 
	         sumz = sumz + Ai_nz(i,k,j,jq,kr)
#  ifdef partial_cell
     &                *min(delqc(i,k,j,kr),delqc(i,k,j+1,kr))
     &                *dz_wtr(i,km1kr,j+jq)
#  endif
     &                 *(t(i,km1kr,j+jq,n,taum1)-t(i,kpkr,j+jq,n,taum1))
# endif
# ifdef full_tensor
	         sumz = sumz - Ai_nz(i,k,j,jq,kr)
     &                 *(t(i,km1kr,j+jq,n,taum1)-t(i,kpkr,j+jq,n,taum1))
     &                 *drodyn(i,k,j,jq)
     &                 *drodzn(i,k,j,jq,kr)/(
     &          +  p5*(drodxn(i,k,j,0,jq)**2 + drodxn(i,k,j,1,jq)**2)
     &          + drodyn(i,k,j,jq)**2 + drodzn(i,k,j,jq,kr)**2 + epsln)
# endif
# ifdef biharmonic_rm
                  rm_laplacian = dxtr(i)*cstr(jrow)**2
     &           *(
     &              (Bi_nz(i+1,k,j,jq,kr)-Bi_nz(i,k,j,jq,kr))
     &                 *dxur(i)*tmask(i+1,k,j)
     &             -(Bi_nz(i,k,j,jq,kr)-Bi_nz(i-1,k,j,jq,kr))
     &                 *dxur(i-1)*tmask(i-1,k,j)
     &            )
     &            + dytr(jrow)*cstr(jrow)
     &           *(
     &             csu(jrow)*(Bi_nz(i,k,j+1,jq,kr)-Bi_nz(i,k,j,jq,kr))
     &                *dyur(jrow)*tmask(i,k,j+1)
     &            -csu(jm1)*(Bi_nz(i,k,j,jq,kr)-Bi_nz(i,k,j-1,jq,kr))
     &                 *dyur(jm1)*tmask(i,k,j-1)
     &            )
	          sumzBih = sumzBih + rm_laplacian
#  ifdef partial_cell
     &                *min(delqc(i,k,j,kr),delqc(i,k,j+1,kr))
     &                *dz_wtr(i,km1kr,j+jq)
#  endif
     &            *(t(i,km1kr,j+jq,n,taum1)-t(i,kpkr,j+jq,n,taum1))

# endif
                enddo
	      enddo
# ifdef partial_cell
              flux_y = csu(jrow)*p25*(sumz+sumzBih)
# else
              flux_y = csu_dzt4r*(sumz+sumzBih)
# endif
c
# ifdef full_tensor
              Ai0 = ahisop*fzisop(k)
	      sumx = c0
	      do ip=0,1
                do jq=0,1
	          sumx = sumx - Ai0*tmask(i,k,j)*tmask(i,k,j+1)
     &               *(t(i+ip,k,j+jq,n,taum1)-t(i-1+ip,k,j+jq,n,taum1))
     &                 *cstr(jrow+jq)
     &               *drodxn(i,k,j,ip,jq)*drodyn(i,k,j,jq)/(
     &             drodxn(i,k,j,ip,jq)**2 + drodyn(i,k,j,jq)**2 + epsln
     &          +  p5*(drodzn(i,k,j,jq,0)**2 + drodzn(i,k,j,jq,1)**2))
	        enddo
	      enddo
	      flux_y = flux_y + p25*csu(jrow)*dxtr(i)*sumx
# endif
c
             diff_fn(i,k,j) = diff_fn(i,k,j) + flux_y 
# ifdef redi_diffusion 
     &        + (T_j(i,k,j,n,1)-T_j(i,k,j,n,0))*csu_dyur(jrow)
     &          * ( K22(i,k,j)
#  if defined isotropic_mixed
     &          + p25*(  diff_cbt(i,k,j,numtra)     - K33(i,k,j)
     &                 + diff_cbt(i,km1,j,numtra)   - K33(i,km1,j)
     &                 + diff_cbt(i,k,j+1,numtra)   - K33(i,k,j+1)
     &                 + diff_cbt(i,km1,j+1,numtra) - K33(i,km1,j+1)
     &                ) 
#  endif
     &            )      
# endif
c
            enddo
          enddo
          call setbcx (diff_fn(1,1,j), imt, km)
        enddo
c
# if defined small_tensor && defined gm_skew && !defined biharmonic_rm
      else 
        do j=js-1,je
          jrow = j + joff
          do k=1,kbot
            do i=2,imtm1
              diff_fn(i,k,j) = diff_fn(i,k,j)
     &        + K22(i,k,j)*csu_dyur(jrow)
     &        *(T_j(i,k,j,n,1)-T_j(i,k,j,n,0))
#  if defined partial_cell
     &        *dhtn(i,k,j)
#  endif
            enddo
          enddo
          call setbcx (diff_fn(1,1,j), imt, km)
        enddo
      endif
# endif
c
c-----------------------------------------------------------------------
c     compute the vertical tracer flux "diff_fbiso" containing the K31
c     and K32 components which are to be solved explicitly. The K33
c     component will be treated implicitly. Note that there are some
c     cancellations of dxu(i-1+ip) and dyu(jrow-1+jq) 
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        do k=1,kbot-1
          do i=2,imtm1
	    sumx    = c0
            sumxBih = c0
            do ip=0,1
              do kr=0,1
# if defined small_tensor || defined gm_skew 
	        sumx = sumx
     &            + Ai_bx(i,k,j,ip,kr)*cstr(jrow)*
#  ifdef partial_cell
     &          min(delqc(i-1+ip,k+kr,j,1-kr),delqc(i+ip,k+kr,j,1-kr))*
#  endif
     &             (T_i(i-1+ip,k+kr,j,n,1) - T_i(i-1+ip,k+kr,j,n,0))
# endif
# ifdef full_tensor
	        sumx = sumx
     &            - Ai_bx(i,k,j,ip,kr)*cstr(jrow)*
     &             (t(i+ip,k+kr,j,n,taum1) - t(i-1+ip,k+kr,j,n,taum1))
     &             *drodxb(i,k,j,ip,kr)
     &             *drodzb(i,k,j,kr)/(drodxb(i,k,j,ip,kr)**2
     &       +  p5*(drodyb(i,k,j,0,kr)**2 + drodyb(i,k,j,1,kr)**2)
     &            + drodzb(i,k,j,kr)**2 + epsln)
# endif
# ifdef biharmonic_rm
                  rm_laplacian = dxtr(i)*cstr(jrow)**2
     &           *(
     &              (Bi_bx(i+1,k,j,ip,kr)-Bi_bx(i,k,j,ip,kr))
     &                 *dxur(i)*tmask(i+1,k,j)
     &             -(Bi_bx(i,k,j,ip,kr)-Bi_bx(i-1,k,j,ip,kr))
     &                 *dxur(i-1)*tmask(i-1,k,j)
     &            )
     &            + dytr(jrow)*cstr(jrow)
     &           *(
     &             csu(jrow)*(Bi_bx(i,k,j+1,ip,kr)-Bi_bx(i,k,j,ip,kr))
     &                *dyur(jrow)*tmask(i,k,j+1)
     &            -csu(jrow-1)*(Bi_bx(i,k,j,ip,kr)-Bi_bx(i,k,j-1,ip,kr))
     &                 *dyur(jrow-1)*tmask(i,k,j-1)
     &            )
	          sumxBih = sumxBih + cstr(jrow)*rm_laplacian
#  ifdef partial_cell
     &          *min(delqc(i-1+ip,k+kr,j,1-kr),delqc(i+ip,k+kr,j,1-kr))
#  endif
     &             *(T_i(i-1+ip,k+kr,j,n,1) - T_i(i-1+ip,k+kr,j,n,0))
# endif
              enddo
	    enddo
c
	    sumy    = c0
            sumyBih = c0
            do jq=0,1
              do kr=0,1
# if defined small_tensor || defined gm_skew 
	        sumy    = sumy 
     &           + Ai_by(i,k,j,jq,kr)*csu(jrow-1+jq)*
#  ifdef partial_cell
     &          min(delqc(i,k+kr,j-1+jq,1-kr),delqc(i,k+kr,j+jq,1-kr))*
#  endif
     &            (T_j(i,k+kr,j-1+jq,n,1)-T_j(i,k+kr,j-1+jq,n,0))
# endif
# ifdef full_tensor
	        sumy    = sumy 
     &           - Ai_by(i,k,j,jq,kr)*csu(jrow-1+jq)*
     &            (t(i,k+kr,j+jq,n,taum1)-t(i,k+kr,j-1+jq,n,taum1))
     &               *drodyb(i,k,j,jq,kr)
     &               *drodzb(i,k,j,kr)/(
     &           p5*(drodxb(i,k,j,0,kr)**2 + drodxb(i,k,j,0,kr)**2)
     &            + drodyb(i,k,j,jq,kr)**2
     &            + drodzb(i,k,j,kr)**2 + epsln)
# endif
# ifdef biharmonic_rm 
                  rm_laplacian = dxtr(i)*cstr(jrow)**2
     &           *(
     &              (Bi_by(i+1,k,j,jq,kr)-Bi_by(i,k,j,jq,kr))
     &                 *dxur(i)*tmask(i+1,k,j)
     &             -(Bi_by(i,k,j,jq,kr)-Bi_by(i-1,k,j,jq,kr))
     &                 *dxur(i-1)*tmask(i-1,k,j)
     &            )
     &            + dytr(jrow)*cstr(jrow)
     &           *(
     &             csu(jrow)*(Bi_by(i,k,j+1,jq,kr)-Bi_by(i,k,j,jq,kr))
     &                *dyur(jrow)*tmask(i,k,j+1)
     &            -csu(jrow-1)*(Bi_by(i,k,j,jq,kr)-Bi_by(i,k,j-1,jq,kr))
     &                 *dyur(jrow-1)*tmask(i,k,j-1)
     &            )
	          sumyBih = sumyBih + csu(jrow-1+jq)*rm_laplacian
#  ifdef partial_cell
     &          *min(delqc(i,k+kr,j-1+jq,1-kr),delqc(i,k+kr,j+jq,1-kr))
#  endif
     &            *(T_j(i,k+kr,j-1+jq,n,1)-T_j(i,k+kr,j-1+jq,n,0))
# endif
              enddo
	    enddo
c
            diff_fbiso(i,k,j) = (dxt4r(i)*(sumx+sumxBih)
     &                        + dyt4r(jrow)*cstr(jrow)*(sumy+sumyBih)) 
# ifdef partial_cell
     &                          /dhwt(i,k,j)
# endif
c
          enddo
        enddo
        do i=2,imtm1
          diff_fbiso(i,0,j)  = c0
          diff_fbiso(i,km,j) = c0
#ifdef bbl_ag
          diff_fbiso(i,km-1,j) = c0
#endif
        enddo
        call setbcx (diff_fbiso(1,0,j), imt, km+1)
      enddo
c
# if defined gm_advect
c
c-----------------------------------------------------------------------
c     compute advective tracer flux at the center of the bottom face of
c     the "T" cells
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=1,kbot-1
          do i=2,imt-1
            adv_fbiso(i,k,j) = adv_vbtiso(i,k,j)*
     &                         (t(i,k,j,n,taum1) + t(i,k+1,j,n,taum1))
          enddo
        enddo
      enddo
c
c     now consider the top and bottom boundaries
c
      do j=js,je
        do i=2,imt-1
          adv_fbiso(i,0,j)  = c0
          adv_fbiso(i,km,j) = c0
        enddo
      enddo
# endif
c
# ifdef trace_indices
      write (stdout,'(2x,9(a,i4))')
     & "=> In isoflux: n=",n," js=",js," je=",je," joff=",joff
     &," diff_fn from row ",js-1," to ",je
     &,", diff_fe from row ",js," to ",je
     &,", diff_fbiso from row ",js," to ",je
# endif
c
# ifdef timing
      call toc ('isoneutral ', 'isoflux')
# endif
      return
      end


# ifdef gent_mcwilliams
      subroutine isopyc_adv (joff, js, je, is, ie)
c
c=======================================================================
c     compute GM eddy-induced transport velocities.
c
c     author:  R.C. Pacanowski e-mail => rcp@gfdl.gov 
c=======================================================================
c
      use arrays_2d
      use memory_window
# include "param.h"
# include "accel.h"
# include "coord.h"
# include "grdvar.h"
# include "isopyc.h"
# include "mw.h"
# include "switch.h"
      dimension top_bc(km), bot_bc(km)
c
      do k=1,km
        top_bc(k) = c1
	bot_bc(k) = c1
      enddo
      top_bc(1)  = c0
      bot_bc(km) = c0
c
c-----------------------------------------------------------------------
c     compute the meridional component of the isoneutral mixing velocity
c     at the center of the northern face of the "T" cells.
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        do k=1,kbot
          sc = c1/(slmxr*dtxsqr(k))
          km1 = max(k-1,1)
          kp1 = min(k+1,kbot)
          do i=1,imt
            Ath0 = agm(i,jrow)*fzisop(k)
c
            at =     (drho(i,k,j,1) + drho(i,k,j+1,1) + drho(i,km1,j,1)
     &              + drho(i,km1,j+1,1))
            bt =     (drho(i,k,j,2) + drho(i,k,j+1,2) + drho(i,km1,j,2)
     &              + drho(i,km1,j+1,2))
            stn = -(at*(ddyt(i,k,j,1) + ddyt(i,km1,j,1))
     &               + bt*(ddyt(i,k,j,2) + ddyt(i,km1,j,2))) /
     &                (at*(ddzt(i,km1,j,1) + ddzt(i,km1,j+1,1))
     &               + bt*(ddzt(i,km1,j,2) + ddzt(i,km1,j+1,2))-epsln)
c
            ab =     (drho(i,k,j,1) + drho(i,k,j+1,1) + drho(i,kp1,j,1)
     &              + drho(i,kp1,j+1,1))
            bb =     (drho(i,k,j,2) + drho(i,k,j+1,2) + drho(i,kp1,j,2)
     &              + drho(i,kp1,j+1,2))
            sbn = -(ab*(ddyt(i,k,j,1) + ddyt(i,kp1,j,1))
     &               + bb*(ddyt(i,k,j,2) + ddyt(i,kp1,j,2))) /
     &                (ab*(ddzt(i,k,j,1) + ddzt(i,k,j+1,1))
     &               + bb*(ddzt(i,k,j,2) + ddzt(i,k,j+1,2))-epsln)
            absstn = abs(stn)
            abssbn = abs(sbn)
# ifdef dm_taper
            fact_t1 =  (absstn-del_dm)*s_dmr
            fact_b1 =  (abssbn-del_dm)*s_dmr
            fact_t =  sign(1.0,fact_t1)*min(19.0,abs(fact_t1))
            fact_b =  sign(1.0,fact_b1)*min(19.0,abs(fact_b1))
	    ath_t = Ath0*tmask(i,k,j)*tmask(i,k,j+1)
     &             *p5*(c1-tanh(fact_t))
	    ath_b = Ath0*tmask(i,kp1,j)*tmask(i,kp1,j+1)
     &             *p5*(c1-tanh(fact_b))
# else
            if (absstn .gt. sc) then
	      ath_t = Ath0*tmask(i,k,j)*tmask(i,k,j+1)
     &              *(sc/(absstn + epsln))**2
	    else
	      ath_t = Ath0*tmask(i,k,j)*tmask(i,k,j+1)
	    endif
            if (abssbn .gt. sc) then
	      ath_b = Ath0*tmask(i,kp1,j)*tmask(i,kp1,j+1)
     &              *(sc/(abssbn + epsln))**2
	    else
	      ath_b = Ath0*tmask(i,kp1,j)*tmask(i,kp1,j+1)
	    endif
# endif
            adv_vntiso(i,k,j) = -(ath_t*stn*top_bc(k) -
     &                            ath_b*sbn*bot_bc(k))*dztr(k)*csu(jrow)
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute the zonal component of the isoneutral mixing velocity
c     at the center of the eastern face of the "T" grid cell.
c-----------------------------------------------------------------------
c
      jstrt = max(js,1+jbuf)
      do j=jstrt,je
        jrow = j + joff 
        do k=1,kbot
          sc = c1/(slmxr*dtxsqr(k))
          km1 = max(k-1,1)
          kp1 = min(k+1,kbot)
          do i=1,imtm1
            Ath0 = agm(i,jrow)*fzisop(k)
c
            at =     (drho(i,k,j,1) + drho(i+1,k,j,1) + drho(i,km1,j,1)
     &              + drho(i+1,km1,j,1))
            bt =     (drho(i,k,j,2) + drho(i+1,k,j,2) + drho(i,km1,j,2)
     &              + drho(i+1,km1,j,2))
            ste =-(at*(ddxt(i,k,j,1) + ddxt(i,km1,j,1))
     &           + bt*(ddxt(i,k,j,2) + ddxt(i,km1,j,2)))
     &          / (at*(ddzt(i,km1,j,1) + ddzt(i+1,km1,j,1))
     &           + bt*(ddzt(i,km1,j,2) + ddzt(i+1,km1,j,2))-epsln)
c
            ab =     (drho(i,k,j,1) + drho(i+1,k,j,1) + drho(i,kp1,j,1)
     &              + drho(i+1,kp1,j,1))
            bb =     (drho(i,k,j,2) + drho(i+1,k,j,2) + drho(i,kp1,j,2)
     &              + drho(i+1,kp1,j,2))
            sbe =-(ab*(ddxt(i,k,j,1) + ddxt(i,kp1,j,1))
     &           + bb*(ddxt(i,k,j,2) + ddxt(i,kp1,j,2))) /
     &            (ab*(ddzt(i,k,j,1) + ddzt(i+1,k,j,1))
     &           + bb*(ddzt(i,k,j,2) + ddzt(i+1,k,j,2))-epsln)

            absste = abs(ste)
            abssbe = abs(sbe)
# ifdef dm_taper
            fact_t1 =  (absste-del_dm)*s_dmr
            fact_b1 =  (abssbe-del_dm)*s_dmr
            fact_t =  sign(1.0,fact_t1)*min(19.0,abs(fact_t1))
            fact_b =  sign(1.0,fact_b1)*min(19.0,abs(fact_b1))
	    ath_t = Ath0*tmask(i,k,j)*tmask(i+1,k,j)
     &             *p5*(c1-tanh(fact_t))
	    ath_b = Ath0*tmask(i,kp1,j)*tmask(i+1,kp1,j)
     &             *p5*(c1-tanh(fact_b))
# else
            if (absste .gt. sc) then
	      ath_t = Ath0*tmask(i,k,j)*tmask(i+1,k,j)
     &              *(sc/(absste + epsln))**2
	    else
	      ath_t = Ath0*tmask(i,k,j)*tmask(i+1,k,j)
	    endif
            if (abssbe .gt. sc) then
	      ath_b = Ath0*tmask(i,kp1,j)*tmask(i+1,kp1,j)
     &              *(sc/(abssbe + epsln))**2
	    else
	      ath_b = Ath0*tmask(i,kp1,j)*tmask(i+1,kp1,j)
	    endif
# endif
            adv_vetiso(i,k,j) = -(ath_t*ste*top_bc(k) - 
     &                            ath_b*sbe*bot_bc(k))*dztr(k)
          enddo
        enddo
      enddo
c
c     set the boundary conditions
c
      do j=jstrt,je
	call setbcx (adv_vetiso(1,1,j), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     compute the vertical component of the isoneutral mixing velocity
c     at the center of the bottom face of the "T" cells, using the
c     continuity equation for the GM velocities
c-----------------------------------------------------------------------
c
      do j=jstrt,je
        do i=1,imt
          adv_vbtiso(i,0,j) = c0
        enddo
      enddo
c
      do j=jstrt,je
        jrow = j + joff
        do k=1,kbot-1
          do i=2,imt
            adv_vbtiso(i,k,j) = dzt(k)*cstr(jrow)*(
     &      (adv_vetiso(i,k,j) - adv_vetiso(i-1,k,j))*dxtr(i) + 
     &      (adv_vntiso(i,k,j) - adv_vntiso(i,k,j-1))*dytr(jrow)) 
          enddo
        enddo
      enddo
c
      do j=jstrt,je
        do k=1,kbot-1
          do i=2,imt
            adv_vbtiso(i,k,j) = adv_vbtiso(i,k,j) + adv_vbtiso(i,k-1,j)
          enddo
        enddo
      enddo
c
      do j=jstrt,je
        jrow = j + joff
        do i=2,imt
          adv_vbtiso(i,kmt(i,jrow),j) = c0
        enddo
      enddo
c
c     set the boundary conditions
c
      do j=jstrt,je
	call setbcx (adv_vbtiso(1,0,j), imt, km+1)
      enddo
c
# ifdef trace_indices
      write (stdout,'(2x,3(a,i4),6(/4x,4(a,i4),a))')
     & "=> In isopyc_adv: js=",js," je=",je," joff=",joff
     &," adv_vntiso calculated on jrows ",js+joff," through ",je+joff
     &," (rows ",js," through ",je," in the MW)"
     &," adv_vetiso calculated on jrows ",jstrt+joff," through ",je+joff
     &," (rows ",jstrt," through ",je," in the MW)"
     &," adv_vbtiso calculated on jrows ",jstrt+joff," through ",je+joff
     &," (rows ",jstrt," through ",je," in the MW)"
# endif
c#  define debug_adv_vel_iso
#  ifdef debug_adv_vel_iso
c
c-----------------------------------------------------------------------
c     to inspect the divergence of "T" cells at point (ipt,jpt)
c     for all levels
c-----------------------------------------------------------------------
c
      sml  = -1.0
      ipt = 61
      jpt = 16
      do j=jstrt,je
        jrow = j + joff
        if (jrow .eq. jpt .and. eorun) then
          write (stdout,*) 'Showing GM advection velocities and'
     &,   'divergence for jrow=',jpt, ' i=',ipt
	  do i=ipt,ipt
            do k=1,kbot
	      divgt =
     &        ((adv_vetiso(i,k,j) - adv_vetiso(i-1,k,j))
     &        *cstr(jrow)*dxtr(i)
     & +      (adv_vntiso(i,k,j) - adv_vntiso(i,k,j-1))
     &        *cstr(jrow)*dytr(jrow)
     & +      (adv_vbtiso(i,k-1,j) - adv_vbtiso(i,k,j))
     &        *dztr(k))*tmask(i,k,j)
              if (abs(divgt) .ge. sml) then
                write (stdout,98) i,k,divgt
     &,                         adv_vetiso(i,k,j), adv_vetiso(i-1,k,j)
     &,                         adv_vntiso(i,k,j), adv_vntiso(i,k,j-1)
     &,                         adv_vbtiso(i,k-1,j), adv_vbtiso(i,k,j)
                write (stdout,*) ' tmask=',tmask(i,k,j)
     &,                        ' kmt=',kmt(i,jrow),' kmu='
     &,                        kmu(i,jrow), ', jrow=',jrow
                write (stdout,*) ' '
	      endif
            enddo
	    write (stdout,*) 'adv_vbtiso(i,0,j) = ',adv_vbtiso(i,0,j)
	    write (stdout,'(///)')
	  enddo
        endif
      enddo
98      format (1x,' i=',i3,', k=',i2,'divtiso=',7e14.7)
#  endif
      return
      end
# endif



#else
      return
      end
#endif

