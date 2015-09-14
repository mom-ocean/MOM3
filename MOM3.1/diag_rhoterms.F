      subroutine density_out (joff, js, je, is, ie, within_window_loop)
#if defined save_density_terms
c=======================================================================
c     Compute terms contributing to the time tendency of 
c     locally referenced potential density.  See manual 
c     for details.
c
c     1. Advection
c     2. vertical diffusion
c     3. horz Laplacian diffusion
c     4. -Dgent_mcwilliams skew-diffusion
c     5. -Dbiharmonic_rm skew-biharmonic diffusion, 
c     6. cabbeling
c     7. thermobaricity
c     8. halobaricity
c
c     Diagnostics written to file rho_terms.dta
c  
c     author:  S. M. Griffies   e-mail=> smg@gfdl.gov
c=======================================================================
c
      use arrays_2d
      use memory_window
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "accel.h"
# include "coord.h"
# include "derived_options.h"
# include "diag.h"
# include "grdvar.h"
# include "hmixc.h"
# include "iounit.h"
# include "isopyc.h"
# include "isopyc_sf.h"
# include "mw.h"
# include "scalar.h"
# include "state.h"
# include "tmngr.h"
# include "dens.h"
# include "switch.h"
# include "taskrows.h"
# include "vmixc.h"
# include "vers.h"
# include "fdift.h"
      character*128 file_name, file_stamp, gvtitle
      logical within_window_loop
      dimension bufzt(km)
      data num_rho /0/
      save num_rho
      save io

c----------------------------------------------------------------------
c     advrho      = contribution to local potential density tendency
c                   due to advection.  Computed in isopyc.F
c                   and commoned inside isopyc.h. 
c     diffvert    = contribution to local potential density tendency 
c                   due to vertical diffusion. 
c     diffhorz    = contribution to local potential density tendency 
c                   due to horizontal diffusion.    
c     diffgmskew  = contribution to local potential density tendency 
c                   due to GM skew diffusion.    
c     diffbihskew = contribution to local potential density tendency 
c                   due to biharmonic skew diffusion.    
c     cabbel      = contribution to local potential density tendency 
c                   due to cabbeling.
c     thermob     = contribution to local potential density tendency 
c                   due to thermobaricity.
c     halob       = contribution to local potential density tendency 
c                   due to halobaricity.
c----------------------------------------------------------------------
# ifdef diagnostics_lomem
      type(domain2D) :: domain
c
c     lomem buffers only within a single mw
c

      real :: buf_advrho(imt,js+joff:je+joff,km)
      real :: diffvert(imt,js+joff:je+joff,km)
      real :: diffhorz(imt,js+joff:je+joff,km)
      real :: diffgmskew(imt,js+joff:je+joff,km)
      real :: diffbihskew(imt,js+joff:je+joff,km)
      real :: cabbel(imt,js+joff:je+joff,km)
      real :: thermob(imt,js+joff:je+joff,km)
      real :: halob(imt,js+joff:je+joff,km)
# endif
# ifdef diagnostics_himem
c
c     domain pointer 
      type(domain2D), pointer :: domain
c
c     himem buffers all the memory windows
c
      real, dimension(:,:,:), save, allocatable :: buf_advrho
      real, dimension(:,:,:), save, allocatable :: diffvert
      real, dimension(:,:,:), save, allocatable :: diffhorz
      real, dimension(:,:,:), save, allocatable :: diffgmskew
      real, dimension(:,:,:), save, allocatable :: diffbihskew
      real, dimension(:,:,:), save, allocatable :: cabbel
      real, dimension(:,:,:), save, allocatable :: thermob
      real, dimension(:,:,:), save, allocatable :: halob
# endif

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_advrho
      type(fieldtype), save :: field_diffvert
      type(fieldtype), save :: field_diffhorz
      type(fieldtype), save :: field_diffgmskew
      type(fieldtype), save :: field_diffbihskew
      type(fieldtype), save :: field_cabbel
      type(fieldtype), save :: field_thermob
      type(fieldtype), save :: field_halob

      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_xt
      type (axistype), save :: axis_yt
      type (axistype), save :: axis_zt
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     Auxiliary fields
c-----------------------------------------------------------------------

      real :: Aiso(0:1,0:1),Athk(0:1,0:1)
      real :: Abih(-1:1,-1:1,0:1,0:1)
      real :: AbihBx(-1:1,-1:1,0:1,0:1)
      real :: AbihBy(-1:1,-1:1,0:1,0:1)
      real :: drhodtdt(imt,km,2:jmw-1)
      real :: drhodtds(imt,km,2:jmw-1)
      real :: drhodsds(imt,km,2:jmw-1)
      real :: ddrhodp(imt,km,2:jmw-1,2)
      real :: horz_fe(imt,km,2:jmw-1,2)
      real :: horz_fn(imt,km,1:jmw-1,2)
      real :: vert_fb(imt,0:km,2:jmw-1,2)
      real :: ison_fe(imt,km,2:jmw-1,2)
      real :: ison_fn(imt,km,1:jmw-1,2)
      real :: ison_fb(imt,0:km,2:jmw-1,2)
      real :: skew_fe(imt,km,2:jmw-1,2)
      real :: skew_fn(imt,km,1:jmw-1,2)
      real :: skew_fb(imt,0:km,2:jmw-1,2)
      real :: bihskew_fe(imt,km,2:jmw-1,2)
      real :: bihskew_fn(imt,km,1:jmw-1,2)
      real :: bihskew_fb(imt,0:km,2:jmw-1,2)
      real :: iso, nonlin
      real :: K11_iso, K22_iso, K11_steep, K22_steep

c-----------------------------------------------------------------------
c     specify file where density terms are written
c-----------------------------------------------------------------------

      file_name = file_stamp ('rho_terms', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, and ocean depth in meters
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Potential density terms' // momver
      bufzt(:) = zt(:)*0.01

c-----------------------------------------------------------------------
c     check whether diagnostic file should be opened
c-----------------------------------------------------------------------

      if (num_rho .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_rho = num_rho + 1

c-----------------------------------------------------------------------
c       axes (on global domain)
c-----------------------------------------------------------------------

        ibeg = mom_domain%x%global%start_index
        iend = mom_domain%x%global%end_index
        jbeg = mom_domain%y%global%start_index
        jend = mom_domain%y%global%end_index

        call mpp_write_meta(
     &      io, axis_xt, 'xt_i', 'degrees_E', 'Longitude of T points',
     &      cartesian='X', domain=mom_domain%x, 
     &      data=xt(ibeg:iend))

        call mpp_write_meta(
     &      io, axis_yt, 'yt_j', 'degrees_N', 'Latitude of T points',
     &      cartesian='Y', domain=mom_domain%y, 
     &      data=yt(jbeg:jend))

        call mpp_write_meta(
     &      io, axis_zt, 'zt_k', 'm', 'Depth of T grid point',
     &      cartesian='Z', sense=-1, data=bufzt)

        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       attributes of data 
c-----------------------------------------------------------------------
c
        call mpp_write_meta(
     &      io, field_advrho, (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'advrho', '(g/cm**3)/s', 'advrho',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

        call mpp_write_meta(
     &      io, field_diffvert, (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'diffvert', '(g/cm**3)/s', 'diffvert',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

        call mpp_write_meta(
     &      io, field_diffhorz, (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'diffhorz', '(g/cm**3)/s', 'diffhorz',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

        call mpp_write_meta(
     &      io, field_diffgmskew, (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'diffgmskew', '(g/cm**3)/s', 'diffgmskew',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

        call mpp_write_meta(
     &      io, field_diffbihskew, (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'diffbihskew', '(g/cm**3)/s', 'diffbihskew',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

        call mpp_write_meta(
     &      io, field_cabbel, (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'cabbel', '(g/cm**3)/s', 'cabbeling',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

        call mpp_write_meta(
     &      io, field_thermob, (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'thermob', '(g/cm**3)/s', 'thermobaricity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

        call mpp_write_meta(
     &      io, field_halob, (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'halob', '(g/cm**3)/s', 'halobaricity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c
c       averaging period 
c
        call mpp_write_meta(
     &      io, field_period, (/axis_t/), 
     &      'period', 'years', 'averaging period',
     &      +0.0E+00, +1.0E+20, -1.0E+34, +1.0E+34)
c
c       figure title 
c
        call mpp_write_meta(io, 'title', cval=gvtitle )

c-----------------------------------------------------------------------
c       write axes now that all meta data has been written
c-----------------------------------------------------------------------
      
        call mpp_write (io,axis_xt)
        call mpp_write (io,axis_yt)
        call mpp_write (io,axis_zt)

c-----------------------------------------------------------------------
c       allocate the buffer variables on data domain if using himem
c-----------------------------------------------------------------------

# ifdef diagnostics_himem
        allocate( buf_advrho(1:imt,jstask:jetask,km) )
        allocate( diffvert(1:imt,jstask:jetask,km) )
        allocate( diffhorz(1:imt,jstask:jetask,km) )
        allocate( diffgmskew(1:imt,jstask:jetask,km) )
        allocate( diffbihskew(1:imt,jstask:jetask,km) )
        allocate( cabbel(1:imt,jstask:jetask,km) )
        allocate( thermob(1:imt,jstask:jetask,km) )
        allocate( halob(1:imt,jstask:jetask,km) )
# endif

c-----------------------------------------------------------------------
c     endif for opening the diagnostic file
c-----------------------------------------------------------------------
      endif

c-----------------------------------------------------------------------
c     point "domain" to the x-y-domain from "mom_domain" 
c-----------------------------------------------------------------------

# ifdef diagnostics_himem
!domain was declared pointer because it needs to inherit linked list
      domain => mom_domain
# else
!domain was not declared pointer
      domain = mom_domain
      domain%y%global = domain%y%compute
      domain%y%compute%start_index = js+joff
      domain%y%compute%end_index   = je+joff
      domain%y%data = domain%y%compute
# endif

c-----------------------------------------------------------------------
c     only do computations when inside the mw. 
c-----------------------------------------------------------------------

      if( within_window_loop )then

c-----------------------------------------------------------------------
c     Initialize to zero 
c-----------------------------------------------------------------------

      do j=js,je
        jrow=j+joff
	do k=1,km
          do i=1,imt
            diffvert(i,jrow,k)     = c0
            diffhorz(i,jrow,k)     = c0
            diffgmskew(i,jrow,k)   = c0
            diffbihskew(i,jrow,k)  = c0
            cabbel(i,jrow,k)       = c0
            thermob(i,jrow,k)      = c0
            halob(i,jrow,k)        = c0
            do n=1,2
              ison_fe(i,k,j,n)    = c0
              skew_fe(i,k,j,n)    = c0
              bihskew_fe(i,k,j,n) = c0
              horz_fe(i,k,j,n)    = c0
            enddo 
          enddo
        enddo
      enddo
      do j=js,je
	do k=0,km
          do i=1,imt
            do n=1,2
              ison_fb(i,k,j,n)    = c0
              skew_fb(i,k,j,n)    = c0
              bihskew_fb(i,k,j,n) = c0
              vert_fb(i,k,j,n)    = c0
            enddo 
          enddo
        enddo
      enddo
      do j=js-1,je
	do k=1,km
          do i=1,imt
            do n=1,2
              ison_fn(i,k,j,n)    = c0
              skew_fn(i,k,j,n)    = c0
              bihskew_fn(i,k,j,n) = c0
              horz_fn(i,k,j,n)    = c0
            enddo 
          enddo
        enddo
      enddo
      do kr=0,1
	do ip=0,1
          Aiso(ip,kr) = c0
          Athk(ip,kr) = c0
          do im=-1,1
            do jn=-1,1
              Abih(im,jn,ip,kr)   = c0
              AbihBx(im,jn,ip,kr) = c0
              AbihBy(im,jn,ip,kr) = c0
            enddo
          enddo
        enddo
      enddo
      if (joff+js .eq. 2) then
c        write(*,*)'joff and js = ',joff,',',js
        do k=1,km
          do i=1,imt
            do n=1,2
              grad_p(i,k,js-1,n) = c0
            enddo 
          enddo 
        enddo 
      endif
c
c-----------------------------------------------------------------------
c     Compute the second derivatives of density with respect to tracers.
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=1,km
          do i=is,ie 
	    tprime = t(i,k,j,1,taum1)-to(k)
	    sprime = t(i,k,j,2,taum1)-so(k)
            drhodtdt(i,k,j) = ddensdtdt(tprime,sprime,k)              
            drhodtds(i,k,j) = ddensdtds(tprime,sprime,k)              
            drhodsds(i,k,j) = ddensdsds(tprime,sprime,k)              
          enddo
        enddo
        call setbcx (drhodtdt(1,1,j),  imt, km)
        call setbcx (drhodtds(1,1,j),  imt, km)
        call setbcx (drhodsds(1,1,j),  imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     Compute the pressure derivative of drhodt and drhods.
c     Compute as vertical average of vertical derivative.
c     Set values on top and bottom equal to those at adjacent
c     interior point.  
c-----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff
        do i=is,ie 

          do k=1,km
            ddrhodp(i,k,j,1) = c0
            ddrhodp(i,k,j,2) = c0
          enddo  

          do k=2,max(kmt(i,jrow)-1,2)
            ddrhodp(i,k,j,1) = -(drho(i,k-1,j,1)-drho(i,k+1,j,1))
     &                /(c1/dz_wtr(i,k-1,j)+c1/dz_wtr(i,k,j))
     &                /(grav*rho0)
            ddrhodp(i,k,j,2) = -(drho(i,k-1,j,2)-drho(i,k+1,j,2))
     &                /(c1/dz_wtr(i,k-1,j)+c1/dz_wtr(i,k,j))
     &                /(grav*rho0)
          enddo

          ddrhodp(i,1,j,1) = ddrhodp(i,2,j,1) 
          ddrhodp(i,1,j,2) = ddrhodp(i,2,j,2)
          if(kmt(i,jrow) .ge. 2) then  
            ddrhodp(i,kmt(i,jrow),j,1) = ddrhodp(i,kmt(i,jrow)-1,j,1) 
            ddrhodp(i,kmt(i,jrow),j,2) = ddrhodp(i,kmt(i,jrow)-1,j,2) 
          endif

        enddo

        call setbcx (ddrhodp(1,1,j,1),imt,km)
        call setbcx (ddrhodp(1,1,j,2),imt,km)
      enddo

c-----------------------------------------------------------------------
c     Vertical diffusion contribution. 
c     Assume vertical diffusion done explicitly.
c     Need to subtract out K33 from diff_cbt(i,k,j,numtra) in order 
c     to remove the isoneutral contribution.  
c-----------------------------------------------------------------------
c
      do j=js,je
        do i=is,ie
          do k=1,km-1
            do n=1,2
              numtra = min(n,nsmix)
              diab_cbt  = diff_cbt(i,k,j,numtra) - K33(i,k,j)
              vert_fb(i,k,j,n) = diab_cbt*dz_wtr(i,k,j)*
     &                         (t(i,k,j,n,taum1) - t(i,k+1,j,n,taum1))
            enddo
          enddo
        enddo
      enddo
      do j=js,je
        do i=is,ie
          do n=1,2
	    kb                = kmt(i,jrow)
            vert_fb(i,0,j,n)  = stf(i,j,n)
            vert_fb(i,kb,j,n) = btf(i,j,n)
          enddo
        enddo
      enddo 
      do j=js,je
        jrow=j+joff
        do k=1,km
          do i=is,ie
            diffvert(i,jrow,k) = c0
            do n=1,2
              diffvert(i,jrow,k) = diffvert(i,jrow,k)  
     &          + drho(i,k,j,n)
     &          *(diff_fb(i,k-1,j) - diff_fb(i,k,j))*dz_tr(i,k,j)
            enddo 
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     Compute the fluxes from 
c     1. -Dredi_diffusion assuming -Dsmall_tensor and full cells.
c     2. -Dgent_mcwilliams assuming -Dgm_skew.
c     3. -Dbiharmonic_rm skew-diffusion.
c
c     NOTE: Cannot use any fluxes from isopyc.F since the isopyc.F
c     calculation combined a number of processes into an  
c     overall flux. Of interest here is the contribution from each 
c     process separately.  Also, cannot use K11 or K22 since they 
c     may have contributions from ahsteep, which is strictly not part
c     of the isoneutral diffusion tensor and a represents diabatic 
c     process. The contribution from ahsteep will be placed inside 
c     horz_fe and horz_fn, as will the contribution from 
c     -Disotropic_mixed.  K33 from isopyc.F can be used
c     since it represents just contributions from isoneutral diffusion.
c-----------------------------------------------------------------------
c     
c     Eastern face of a T-box
c     
      do j = js,je
        jrow = j+joff  
        do k = 1,km
          sc = c1/(slmxr*dtxsqr(k))
          dzt4r = p5*dzt2r(k)
          km1 = max(k-1,1)
          do i=2,imtm1
c           
c           Diffusion coefficients 
c           
            Ai0  = (redi_gm(3)*ahisop+redi_gm(4)*agm(i,jrow))*fzisop(k)
            Ath0 = agm(i,jrow)*fzisop(k)
            Abih0= abihrm*fzisop(k)
            Astp = (redi_gm(5)*ahsteep+redi_gm(6)*agm(i,jrow))*fzisop(k)
c
            sumziso   = c0
            sumzsteep = c0
            do kr=0,1
              km1kr = max(k-1+kr,1)
              kpkr  = min(k+kr,km)
              do ip=0,1
                sxe    = -drodxe(i,k,j,ip)
     &                    /(drodze(i,k,j,ip,kr)+epsln)
                abssxe = abs(sxe)
#  ifdef dm_taper
                fact1 = (abssxe-del_dm)*s_dmr   
                fact2 = sign(1.0,fact1)*min(19.0,abs(fact1)) 
                taper = p5*(c1-tanh(fact2))   
#  endif        
#  ifdef gkw_taper
                if (abssxe .gt. sc) then
                  taper = (sc/(abssxe + epsln))**2
                else
                  taper = c1
                endif
#  endif        
	        sumziso = sumziso +
     &                 Ai0*taper*tmask(i,k,j)*tmask(i+1,k,j)*
     &                 dzw(k-1+kr)
                if(Astp .gt. Ai0*taper) then
                  sumzsteep = sumzsteep + (Astp - Ai0*taper)
     &                 *tmask(i,k,j)*tmask(i+1,k,j)*
     &                 dzw(k-1+kr)
                endif
c
                Aiso(ip,kr)=Ai0*taper*sxe*tmask(i,k,j)*tmask(i+1,k,j)
                Athk(ip,kr)=Ath0*taper*sxe*tmask(i,k,j)*tmask(i+1,k,j)
c               
c               The biharmonic scheme needs
c               to reach out to i=+-1 and j=+-1.
c               
                do im = -1,1
                  do jn = -1,1
                    sxe   = -drodxe(i+im,k,j+jn,ip)
     &                   /(drodze(i+im,k,j+jn,ip,kr)+epsln)
                    abssxe = abs(sxe)
#  ifdef dm_taper   
                    fact1 = (abssxe-del_dm)*s_dmr   
                    fact2 = sign(1.0,fact1)*min(19.0,abs(fact1)) 
                    taper = p5*(c1-tanh(fact2))   
#  endif            
#  ifdef gkw_taper  
                    if (abssxe .gt. sc) then
                      taper = (sc/(abssxe + epsln))**2
                    else
                      taper = c1
                    endif
#  endif            
                    Abih(im,jn,ip,kr)=Abih0*taper*sxe
     &                   *tmask(i+im,k,j+jn)*tmask(i+im+1,k,j+jn)
                  enddo
                enddo
c               
              enddo
            enddo
c
	    K11_iso   = dzt4r*sumziso
#   ifdef partial_cell
	    K11_steep = p25*sumzsteep
#   else
    	    K11_steep = dzt4r*sumzsteep
#   endif
c
c           X-component of flux for the two active tracers
c           
            do n=1,2
              numtra = min(n,nsmix)
              sumzdiff = c0
              sumzskew = c0
              sumzbih  = c0
              do kr=0,1
                km1kr = max(k-1+kr,1)
                kpkr  = min(k+kr,km)
                do ip=0,1
                  sumzdiff = sumzdiff + Aiso(ip,kr)*
     &                 (t(i+ip,km1kr,j,n,taum1)-t(i+ip,kpkr,j,n,taum1))
                  sumzskew = sumzskew + Athk(ip,kr)*
#  ifdef partial_cell
     &                 min(delqc(i,k,j,kr),delqc(i+1,k,j,kr))*
#  endif          
     &                 (t(i+ip,km1kr,j,n,taum1)-t(i+ip,kpkr,j,n,taum1))
                  rm_laplacian = dxtr(i)*cstr(jrow)**2
     &                 *(
     &                 (Abih(1,0,ip,kr)-Abih(0,0,ip,kr))
     &                 *dxur(i)*tmask(i+1,k,j)
     &                 -(Abih(0,0,ip,kr)-Abih(-1,0,ip,kr))
     &                 *dxur(i-1)*tmask(i-1,k,j)
     &                 )
     &                 + dytr(jrow)*cstr(jrow)
     &                 *(
     &                 csu(jrow)*(Abih(0,1,ip,kr)-Abih(0,0,ip,kr))
     &                 *dyur(jrow)*tmask(i,k,j+1)
     &                 -csu(jrow-1)*(Abih(0,0,ip,kr)-Abih(0,-1,ip,kr))
     &                 *dyur(jrow-1)*tmask(i,k,j-1)
     &                 )
	          sumzBih = sumzBih + rm_laplacian*
#  ifdef partial_cell
     &                 min(delqc(i,k,j,kr),delqc(i+1,k,j,kr))*
#  endif          
     &                 (t(i+ip,km1kr,j,n,taum1)-t(i+ip,kpkr,j,n,taum1))
                enddo
              enddo
c             
              flux_xdiff = sumzdiff*dzt4r
# ifdef partial_cell
              flux_xskew = sumzskew*p25
              flux_xbih  = sumzbih*p25
# else        
              flux_xskew = sumzskew*dzt4r
              flux_xbih  = sumzbih*dzt4r
# endif       
c             
              ison_fe(i,k,j,n) =  flux_xdiff*tmask(i+1,k,j)
     &        + cstdxur(i,j)*(T_i(i,k,j,n,1) - T_i(i,k,j,n,0))
     &          *K11_iso
              horz_fe(i,k,j,n) = 
     &        + cstdxur(i,j)*(T_i(i,k,j,n,1) - T_i(i,k,j,n,0))
     &          *(K11_steep
#  ifdef isotropic_mixed
     &        + p25*(  diff_cbt(i,k,j,numtra)     - K33(i,k,j)
     &               + diff_cbt(i,km1,j,numtra)   - K33(i,km1,j)
     &               + diff_cbt(i+1,k,j,numtra)   - K33(i+1,k,j)
     &               + diff_cbt(i+1,km1,j,numtra) - K33(i+1,km1,j)
     &              )
#  endif
     &        )
              skew_fe(i,k,j,n)    = flux_xskew*tmask(i+1,k,j)
              bihskew_fe(i,k,j,n) = flux_xbih*tmask(i+1,k,j)
c             
            enddo
          enddo
        enddo
        call setbcx (ison_fe(1,1,j,1), imt, km)
        call setbcx (ison_fe(1,1,j,2), imt, km)
        call setbcx (skew_fe(1,1,j,1), imt, km)
        call setbcx (skew_fe(1,1,j,2), imt, km)
        call setbcx (bihskew_fe(1,1,j,1), imt, km)
        call setbcx (bihskew_fe(1,1,j,2), imt, km)
        call setbcx (horz_fe(1,1,j,1), imt, km)
        call setbcx (horz_fe(1,1,j,2), imt, km)
      enddo
c     
c     Northern face of a T-box
c     
      do j=js-1,je
        jrow = j + joff
        do k=1,km
          sc = c1/(slmxr*dtxsqr(k))
          dzt4r = p5*dzt2r(k)
          csu_dzt4r = csu(jrow)*p5*dzt2r(k)
          do i=2,imtm1
c           
c           Diffusion coefficients 
c           
            Ai0  = (redi_gm(3)*ahisop+redi_gm(4)*agm(i,jrow))*fzisop(k)
            Ath0 = agm(i,jrow)*fzisop(k)
            Abih0= abihrm*fzisop(k)
            Astp = (redi_gm(5)*ahsteep+redi_gm(6)*agm(i,jrow))*fzisop(k)
c
            sumziso   = c0
            sumzsteep = c0
            do kr=0,1
              km1kr = max(k-1+kr,1)
              kpkr = min(k+kr,km)
              do jq=0,1
                syn= -drodyn(i,k,j,jq)
     &               /(drodzn(i,k,j,jq,kr)+epsln)
                abssyn  = abs(syn)
# ifdef dm_taper
                fact1 = (abssyn-del_dm)*s_dmr
                fact2 = sign(1.0,fact1)*min(19.0,abs(fact1))
                taper = p5*(c1-tanh(fact2)) 
#  endif        
#  ifdef gkw_taper
                if (abssyn .gt. sc) then
                  taper = (sc/(abssyn + epsln))**2
                else
                  taper = c1
                endif
#  endif        
	        sumziso = sumziso +
     &                 Ai0*taper*tmask(i,k,j)*tmask(i,k,j+1)*
     &                 dzw(k-1+kr)
                if(Astp .gt. Ai0*taper) then
                  sumzsteep = sumzsteep + (Astp - Ai0*taper)
     &                 *tmask(i,k,j)*tmask(i,k,j+1)*
     &                 dzw(k-1+kr)
                endif

                Aiso(jq,kr)=Ai0*taper*syn*tmask(i,k,j)*tmask(i,k,j+1)
                Athk(jq,kr)=Ath0*taper*syn*tmask(i,k,j)*tmask(i,k,j+1)
c               
c               The biharmonic scheme needs to 
c               reach out to i+-1 and j+-1.
c               
                do im = -1,1
                  do jn = -1,1
                    syn   = -drodyn(i+im,k,j+jn,jq)
     &                       /(drodzn(i+im,k,j+jn,jq,kr)+epsln)
                    abssyn = abs(syn)
#  ifdef dm_taper   
                    fact1 = (abssyn-del_dm)*s_dmr   
                    fact2 = sign(1.0,fact1)*min(19.0,abs(fact1)) 
                    taper = p5*(c1-tanh(fact2))   
#  endif            
#  ifdef gkw_taper  
                    if (abssyn .gt. sc) then
                      taper = (sc/(abssyn + epsln))**2
                    else
                      taper = c1
                    endif
#  endif            
                    Abih(im,jn,ip,kr)=Abih0*taper*syn
     &                   *tmask(i+im,k,j+jn)*tmask(i+im,k,j+1+jn)
                  enddo
                enddo
c
              enddo
            enddo
c
	    K22_iso   = dzt4r*sumziso
#   ifdef partial_cell
	    K22_steep = p25*sumzsteep
#   else
	    K22_steep = dzt4r*sumzsteep
#   endif
c
c           
c           Y-component of flux for the two active tracers
c           
            do n=1,2
              numtra = min(n,nsmix)
              sumzdiff =c0
              sumzskew =c0 
              sumzbih  =c0
              do kr=0,1
                km1kr = max(k-1+kr,1)
                kpkr = min(k+kr,km)
                do jq=0,1
                  sumzdiff = sumzdiff + Aiso(jq,kr)*
     &                 (t(i,km1kr,j+jq,n,taum1)-t(i,kpkr,j+jq,n,taum1))
                  sumzskew = sumzskew + Athk(jq,kr)*
#  ifdef partial_cell
     &                 min(delqc(i,k,j,kr),delqc(i,k,j+1,kr))*
#  endif          
     &                 (t(i,km1kr,j+jq,n,taum1)-t(i,kpkr,j+jq,n,taum1))
                  rm_laplacian = dxtr(i)*cstr(jrow)**2
     &                 *(
     &                 (Abih(1,0,jq,kr)-Abih(0,0,jq,kr))
     &                 *dxur(i)*tmask(i+1,k,j)
     &                 -(Abih(0,0,jq,kr)-Abih(-1,0,jq,kr))
     &                 *dxur(i-1)*tmask(i-1,k,j)
     &                 )
     &                 + dytr(jrow)*cstr(jrow)
     &                 *(
     &                 csu(jrow)*(Abih(0,1,jq,kr)-Abih(0,0,jq,kr))
     &                 *dyur(jrow)*tmask(i,k,j+1)
     &                 -csu(jrow-1)*(Abih(0,0,jq,kr)-Abih(0,-1,jq,kr))
     &                 *dyur(jrow-1)*tmask(i,k,j-1)
     &                 )
	          sumzBih = sumzBih + rm_laplacian*
#  ifdef partial_cell
     &                min(delqc(i,k,j,kr),delqc(i,k,j+1,kr))*
#  endif
     &                (t(i,km1kr,j+jq,n,taum1)-t(i,kpkr,j+jq,n,taum1))
                enddo
              enddo
c             
              flux_ydiff = csu_dzt4r*sumzdiff
# ifdef partial_cell
              flux_yskew = csu(jrow)*p25*sumzskew
              flux_ybih  = csu(jrow)*p25*sumzbih
# else        
              flux_yskew = csu_dzt4r*sumzskew
              flux_ybih  = csu_dzt4r*sumzbih
# endif       
c             
              ison_fn(i,k,j,n) =  flux_ydiff*tmask(i,k,j+1)  
     &             + csu_dyur(jrow)*(T_j(i,k,j,n,1)-T_j(i,k,j,n,0))
     &               *K22_iso
              horz_fn(i,k,j,n) = 
     &          csu_dyur(jrow)*(T_j(i,k,j,n,1)-T_j(i,k,j,n,0))
     &          * (K22_steep
#  ifdef isotropic_mixed
     &          + p25*(  diff_cbt(i,k,j,numtra)     - K33(i,k,j)
     &                 + diff_cbt(i,km1,j,numtra)   - K33(i,km1,j)
     &                 + diff_cbt(i,k,j+1,numtra)   - K33(i,k,j+1)
     &                 + diff_cbt(i,km1,j+1,numtra) - K33(i,km1,j+1)
     &                ) 
#  endif
     &            )
              skew_fn(i,k,j,n)    = flux_yskew*tmask(i,k,j+1) 
              bihskew_fn(i,k,j,n) = flux_ybih*tmask(i,k,j+1)
c             
            enddo
          enddo
        enddo
        call setbcx (ison_fn(1,1,j,1), imt, km)
        call setbcx (ison_fn(1,1,j,2), imt, km)
        call setbcx (skew_fn(1,1,j,1), imt, km)
        call setbcx (skew_fn(1,1,j,2), imt, km)
        call setbcx (bihskew_fn(1,1,j,1), imt, km)
        call setbcx (bihskew_fn(1,1,j,2), imt, km)
        call setbcx (horz_fn(1,1,j,1), imt, km)
        call setbcx (horz_fn(1,1,j,2), imt, km)
      enddo
c
c     Bottom face of T-box.
c
      do j=js,je
        jrow = j + joff
        do k=1,kmm1
          sc = c1/(slmxr*dtxsqr(k))
          kp1   = min(k+1,km)
          fzisopb = p5*(fzisop(k+1)+fzisop(k))
          do i=2,imtm1
            Ai0  = (redi_gm(3)*ahisop+redi_gm(4)*agm(i,jrow))*fzisop(k)
            Ath0 = agm(i,jrow)*fzisopb
            Abih0= abihrm*fzisopb
c              
c           Eastward slopes at the base of T cells
c
            do ip=0,1
              do kr=0,1
                sxb = -drodxb(i,k,j,ip,kr)
     &                /(drodzb(i,k,j,kr)+epsln)
                abssxb = abs(sxb)
# ifdef dm_taper
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
                Aiso(ip,kr)  =  Ai0*taper*sxb*tmask(i,k+1,j)
                Athk(ip,kr)  =  Ath0*taper*sxb*tmask(i,k+1,j)
c
c               The biharmonic scheme needs to 
c               reach out to i+-1 and j+-1.
c               
                do im = -1,1
                  do jn = -1,1
                    sxb   = -drodxb(i+im,k,j+jn,ip,kr)
     &                /(drodzb(i+im,k,j+jn,kr)+epsln)
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
                    AbihBx(im,jn,ip,kr)=Abih0*taper*sxb
     &                   *tmask(i+im,k,j+jn)*tmask(i+1+im,k,j+jn)
                  enddo
                enddo
c
              enddo
            enddo
c
c           Northward slopes at the base of T cells
c
            do jq=0,1
              facty = csu(jrow-1+jq)*dyu(jrow-1+jq)
              do kr=0,1
                syb   =  -drodyb(i,k,j,jq,kr)
     &                   /(drodzb(i,k,j,kr)+epsln)
                abssyb = abs(syb)
# ifdef dm_taper  
                fact1  = (abssyb-del_dm)*s_dmr
                fact2  = sign(1.0,fact1)*min(19.0,abs(fact1))
                taper  = p5*(c1-tanh(fact2))
# endif
# ifdef gkw_taper
                if (abssyb .gt. sc) then
                  taper = (sc/(abssyb + epsln))**2
                else
                  taper = c1
                endif
# endif
  
                Aiso(jq,kr) = Ai0*taper*syb*tmask(i,k+1,j)
                Athk(jq,kr) = Ath0*taper*syb*tmask(i,k+1,j)
c
c               The biharmonic scheme needs to 
c               reach out to i+-1 and j+-1.
c               
                do im = -1,1
                  do jn = -1,1
                    syb = -drodyb(i+im,k,j+jn,jq,kr)
     &                   /(drodzb(i+im,k,j+jn,kr)+epsln)
                    abssyb = abs(syb)
# ifdef dm_taper   
                    fact1 = (abssyb-del_dm)*s_dmr   
                    fact2 = sign(1.0,fact1)*min(19.0,abs(fact1)) 
                    taper = p5*(c1-tanh(fact2))   
# endif            
# ifdef gkw_taper  
                    if (abssyb .gt. sc) then
                      taper = (sc/(abssyb + epsln))**2
                    else
                      taper = c1
                    endif
# endif            
                    AbihBy(im,jn,jq,kr)=Abih0*taper*syb
     &                   *tmask(i+im,k,j+jn)*tmask(i+1+im,k,j+jn)
                  enddo
                enddo
c
              enddo
            enddo
c
c           Z-component of flux for the two active tracers
c
            do n=1,2
              sumxdiff = c0
              sumxskew = c0
              sumxbih  = c0
              sumydiff = c0
              sumyskew = c0
              sumybih  = c0
c
              do ip=0,1
                do kr=0,1
                  sumxdiff = sumxdiff + Aiso(ip,kr)*cstr(jrow)*
     &           (T_i(i,k+kr,j,n,ip) - T_i(i-1,k+kr,j,n,ip))
                 sumxskew = sumxskew + Athk(ip,kr)*cstr(jrow)*
# ifdef partial_cell
     &           min(delqc(i-1+ip,k+kr,j,1-kr),delqc(i+ip,k+kr,j,1-kr))*
# endif
     &           (T_i(i,k+kr,j,n,ip) - T_i(i-1,k+kr,j,n,ip))
                  rm_laplacian = dxtr(i)*cstr(jrow)**2
     &               *(
     &               (AbihBx(1,0,ip,kr)-AbihBx(0,0,ip,kr))
     &               *dxur(i)*tmask(i+1,k,j)
     &               -(AbihBx(0,0,ip,kr)-AbihBx(-1,0,ip,kr))
     &               *dxur(i-1)*tmask(i-1,k,j)
     &               )
     &               + dytr(jrow)*cstr(jrow)
     &               *(
     &               csu(jrow)*(AbihBx(0,1,ip,kr)-AbihBx(0,0,ip,kr))
     &               *dyur(jrow)*tmask(i,k,j+1)
     &               -csu(jrow-1)*(AbihBx(0,0,ip,kr)-AbihBx(0,-1,ip,kr))
     &               *dyur(jrow-1)*tmask(i,k,j-1)
     &               )            
	          sumxBih = sumxBih + cstr(jrow)*rm_laplacian*
# ifdef partial_cell
     &          min(delqc(i-1+ip,k+kr,j,1-kr),delqc(i+ip,k+kr,j,1-kr))*
# endif
     &          (T_i(i,k+kr,j,n,ip) - T_i(i-1,k+kr,j,n,ip))
                enddo
              enddo
c
              do jq=0,1
                do kr=0,1
                  sumydiff = sumydiff + Aiso(jq,kr)*csu(jrow-1+jq)*
     &            (T_j(i,k+kr,j,n,jq)-T_j(i,k+kr,j-1,n,jq))
                  sumyskew = sumyskew + Athk(jq,kr)*csu(jrow-1+jq)*
# ifdef partial_cell
     &          min(delqc(i,k+kr,j-1+jq,1-kr),delqc(i,k+kr,j+jq,1-kr))*
# endif
     &            (T_j(i,k+kr,j,n,jq)-T_j(i,k+kr,j-1,n,jq))
                  rm_laplacian = dxtr(i)*cstr(jrow)**2
     &               *(
     &               (AbihBy(1,0,jq,kr)-AbihBy(0,0,jq,kr))
     &               *dxur(i)*tmask(i+1,k,j)
     &               -(AbihBy(0,0,jq,kr)-AbihBy(-1,0,jq,kr))
     &               *dxur(i-1)*tmask(i-1,k,j)
     &               )
     &               + dytr(jrow)*cstr(jrow)
     &               *(
     &               csu(jrow)*(AbihBy(0,1,jq,kr)-AbihBy(0,0,jq,kr))
     &               *dyur(jrow)*tmask(i,k,j+1)
     &               -csu(jrow-1)*(AbihBy(0,0,jq,kr)-AbihBy(0,-1,jq,kr))
     &               *dyur(jrow-1)*tmask(i,k,j-1)
     &               )            
	          sumyBih = sumyBih + cstr(jrow)*rm_laplacian*
# ifdef partial_cell
     &          min(delqc(i,k+kr,j-1+jq,1-kr),delqc(i,k+kr,j+jq,1-kr))*
# endif
     &            (T_j(i,k+kr,j,n,jq)-T_j(i,k+kr,j-1,n,jq))
                enddo
              enddo
c
c             z-flux component
c             Note that the isoneutral diffusive flux through the bottom
c             includes the diagonal term for purposes of diagnostics. 
c             For the actual calculation in the tracer equation, 
c             K33 is done implicitly in time.
c
              ison_fb(i,k,j,n) =
     &                tmask(i,k+1,j)*(
     &                dxt4r(i)*sumxdiff
     &                + dyt4r(jrow)*cstr(jrow)*sumydiff
     &                + K33(i,k,j)*dzwr(k)
     &                * (t(i,k,j,n,taum1) - t(i,k+1,j,n,taum1))
     &                )
              skew_fb(i,k,j,n) =
     &                tmask(i,k+1,j)*(
     &                dxt4r(i)*sumxskew
     &                + dyt4r(jrow)*cstr(jrow)*sumyskew
     &                )
# ifdef partial_cell
     &                /dht(i,k,j)
# endif
              bihskew_fb(i,k,j,n) =
     &                tmask(i,k+1,j)*(
     &                dxt4r(i)*sumxBih
     &                + dyt4r(jrow)*cstr(jrow)*sumyBih
     &                )
# ifdef partial_cell
     &                /dht(i,k,j)
# endif
c
            enddo
          enddo
        enddo
c
        do i=2,imtm1
          do n=1,2
            ison_fb(i,0,j,n)     = c0
            ison_fb(i,km,j,n)    = c0
            skew_fb(i,0,j,n)     = c0
            skew_fb(i,km,j,n)    = c0
            bihskew_fb(i,0,j,n)  = c0
            bihskew_fb(i,km,j,n) = c0
          enddo
        enddo
c
        call setbcx (ison_fb(1,0,j,1), imt, km+1)
        call setbcx (ison_fb(1,0,j,2), imt, km+1)
        call setbcx (skew_fb(1,0,j,1), imt, km+1)
        call setbcx (skew_fb(1,0,j,2), imt, km+1)
        call setbcx (bihskew_fb(1,0,j,1), imt, km+1)
        call setbcx (bihskew_fb(1,0,j,2), imt, km+1)
      enddo
c
c-----------------------------------------------------------------------
c     Horizontal diffusion contribution. This term is present only  
c     when there is nonzero horizontal background diffusion
c     or when there is an added contribution from isotropic_mixed
c     or ahsteep. 
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=1,km
            do i=is-1,ie              
              do n=1,2
                horz_fe(i,k,j,n) = horz_fe(i,k,j,n)
     &               + diff_cet*cstdxur(i,j)*
     &                 (T_i(i,k,j,n,1) - T_i(i,k,j,n,0))
              enddo
            enddo
          enddo
        enddo
        do j=js-1,je
          jrow = j + joff
          do k=1,km
            do i=is,ie
              do n=1,2
                horz_fn(i,k,j,n) = horz_fn(i,k,j,n)
     &               + diff_cnt*csu_dyur(jrow)*
     &                 (T_j(i,k,j,n,1) - T_j(i,k,j,n,0))
              enddo
            enddo
          enddo
        enddo
c
        do j=js,je
          jrow = j + joff
          do k=1,kbot
            do i=is,ie 
              diffhorz(i,jrow,k) = c0
              do n=1,2
                diffhorz(i,jrow,k) = diffhorz(i,jrow,k) +
     &               drho(i,k,j,n)
     &               *(
     &                  (horz_fe(i,k,j,n)*tmask(i+1,k,j)
     &                  -horz_fe(i-1,k,j,n)*tmask(i-1,k,j))*dx_tr(i,k,j)
     &                 +(horz_fn(i,k,j,n)*tmask(i,k,j+1)
     &                  -horz_fn(i,k,j-1,n)*tmask(i,k,j-1))*dy_tr(i,k,j)
     &               )
              enddo
            enddo
          enddo
        enddo 
c
c
c-----------------------------------------------------------------------
c     Construct contribution from GM.
c     Do so assuming GM implemented with skew-diffusion,
c     even if the isopyc.F code actually uses the advective flux. 
c-----------------------------------------------------------------------
c       
# if defined gent_mcwilliams
          do j=js,je
            jrow = j + joff
            do k=1,km
              do i=is,ie
                diffgmskew(i,jrow,k) = c0
                do n=1,2
                  diffgmskew(i,jrow,k) = diffgmskew(i,jrow,k) 
     &             + drho(i,k,j,n)
     &             *(
     &               (skew_fe(i,k,j,n)-skew_fe(i-1,k,j,n))*dx_tr(i,k,j)
     &             + (skew_fn(i,k,j,n)-skew_fn(i,k,j-1,n))*dy_tr(i,k,j)
     &             + (skew_fb(i,k-1,j,n)-skew_fb(i,k,j,n))*dz_tr(i,k,j)
     &             )
                enddo
              enddo
            enddo
          enddo
# else
          do j=js,je
            jrow = j + joff
            do k=1,km
              do i=is,ie
                diffgmskew(i,jrow,k) = c0
              enddo
            enddo
          enddo
# endif
c
c-----------------------------------------------------------------------
c     Construct contribution from -Dbiharmonic_rm.
c-----------------------------------------------------------------------
c       
# if defined biharmonic_rm
          do j=js,je
            jrow = j + joff
            do k=1,km
              do i=is,ie
                diffbihskew(i,jrow,k) = c0
                do n=1,2
                  diffbihskew(i,jrow,k) = diffbihskew(i,jrow,k) 
     &             + drho(i,k,j,n)
     &             *(
     &               (bihskew_fe(i,k,j,n)-bihskew_fe(i-1,k,j,n))
     &                 *dx_tr(i,k,j)
     &             + (bihskew_fn(i,k,j,n)-bihskew_fn(i,k,j-1,n))
     &                 *dy_tr(i,k,j)
     &             + (bihskew_fb(i,k-1,j,n)-bihskew_fb(i,k,j,n))
     &                 *dz_tr(i,k,j)
     &             )
                enddo
              enddo
            enddo
          enddo
# else
          do j=js,je
            jrow = j + joff
            do k=1,km
              do i=is,ie
                diffbihskew(i,jrow,k) = c0
              enddo
            enddo
          enddo
# endif
c
c-----------------------------------------------------------------------
c     Construct the contribution from cabbeling 
c     Note the extra minus sign is to absorb the 
c     minus convention used to compute the diffusive flux components.
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        do k=1,km
          do i=is,ie
c
            iso =
     &          p5*(
     &          ison_fe(i,k,j,1)
     &          *(t(i+1,k,j,1,taum1)-t(i,k,j,1,taum1))*cstdxur(i,j)
     &          +
     &          ison_fe(i-1,k,j,1)
     &          *(t(i,k,j,1,taum1)-t(i-1,k,j,1,taum1))*cstdxur(i-1,j)
     &             )
     &          +
     &          p5*(
     &          ison_fn(i,k,j,1)
     &          *(t(i,k,j+1,1,taum1)-t(i,k,j,1,taum1))*csu_dyur(jrow)
     &          +
     &          ison_fn(i,k,j-1,1)
     &          *(t(i,k,j,1,taum1)-t(i,k,j-1,1,taum1))*csu_dyur(jrow-1)
     &             )
     &          +
     &          p5*(
     &          ison_fb(i,k,j,1)
     &          *(t(i,k,j,1,taum1)-t(i,k+1,j,1,taum1))*dzwr(k)
     &          +
     &          ison_fb(i,k-1,j,1)
     &          *(t(i,k-1,j,1,taum1)-t(i,k,j,1,taum1))*dzwr(k-1)
     &             )
c
            nonlin = 
     &         (
     &          drhodtdt(i,k,j)
     &        - 2.0*drhodtds(i,k,j)*drho(i,k,j,1)/(epsln+drho(i,k,j,2)) 
     &        + drhodsds(i,k,j)*(drho(i,k,j,1)/(epsln+drho(i,k,j,2)))**2
     &         )
c
            cabbel(i,jrow,k) = -iso * nonlin
c 
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     Construct the contributions from thermobaricity and halobaricity.
c     Note the extra minus sign is to absorb the 
c     minus convention used to compute the diffusive flux components.
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        do k=1,km
          do i=is,ie
c
c           thermobaricity
c
            termX = 
     &      (grad_p(i,k,j,1)   + grad_p(i,k,j-1,1)  )*ison_fe(i,k,j,1)
     &     +(grad_p(i-1,k,j,1) + grad_p(i-1,k,j-1,1))*ison_fe(i-1,k,j,1)
c
            termY = 
     &      (grad_p(i,k,j,2)   + grad_p(i-1,k,j,2)  )*ison_fn(i,k,j,1)
     &     +(grad_p(i,k,j-1,2) + grad_p(i-1,k,j-1,2))*ison_fn(i,k,j-1,1)
               termZ = 
     &         -rho0*grav
     &         *(ison_fb(i,k-1,j,1)+ison_fb(i,k,j,1))
c
            thermob(i,jrow,k) = -ddrhodp(i,k,j,1)
     &              *(p25*termX + p25*termY + p5*termZ)
c
c           halobaricity
c
            halob(i,jrow,k) = -(drho(i,k,j,1)*ddrhodp(i,k,j,2))/
     &                        (drho(i,k,j,2)*ddrhodp(i,k,j,1)+epsln)
     &                        *thermob(i,jrow,k)
c 
          enddo
        enddo
      enddo
c
# ifdef save_density_terms_debug
      do j=js,je
        do k=1,km
          do i=is,ie
            jrow = j + joff
            write(*,*)'----------------------------------------------- '
            write(*,*)
     &      'drhodtdt(',i,',',k,',',jrow,')  = ',drhodtdt(i,k,j)
            write(*,*)
     &      'drhodtds(',i,',',k,',',jrow,')  = ',drhodtds(i,k,j)
            write(*,*)
     &      'drhodsds(',i,',',k,',',jrow,')    = ',drhodsds(i,k,j)
            write(*,*)
     &      'ddrhodp(',i,',',k,',',jrow,',1)   = ',ddrhodp(i,k,j,1)
            write(*,*)
     &      'ddrhodp(',i,',',k,',',jrow,',2)   = ',ddrhodp(i,k,j,2)
            write(*,*)
     &      'grad_p(',i,',',k,',',jrow,',1)    = ',grad_p(i,k,j,1)
            write(*,*)
     &      'grad_p(',i,',',k,',',jrow,',2)    = ',grad_p(i,k,j,2)
            write(*,*)
     &      'grad_p(',i-1,',',k,',',jrow,',1)  = ',grad_p(i-1,k,j,1)
            write(*,*)
     &      'grad_p(',i-1,',',k,',',jrow,',2)  = ',grad_p(i-1,k,j,2)
            write(*,*)
     &      'grad_p(',i,',',k,',',jrow-1,',1)  = ',grad_p(i,k,j-1,1)
            write(*,*)
     &      'grad_p(',i,',',k,',',jrow-1,',2)  = ',grad_p(i,k,j-1,2)
            write(*,*)
     &      'ison_fe(',i,',',k,',',jrow,',1)   = ',ison_fe(i,k,j,1)
            write(*,*)
     &      'ison_fe(',i,',',k,',',jrow,',2)   = ',ison_fe(i,k,j,2)
            write(*,*)
     &      'ison_fe(',i-1,',',k,',',jrow,',1) = ',ison_fe(i-1,k,j,1)
            write(*,*)
     &      'ison_fe(',i-1,',',k,',',jrow,',2) = ',ison_fe(i-1,k,j,2)
            write(*,*)
     &      'ison_fn(',i,',',k,',',jrow,',1)   = ',ison_fn(i,k,j,1)
            write(*,*)
     &      'ison_fn(',i,',',k,',',jrow,',2)   = ',ison_fn(i,k,j,2)
            write(*,*)
     &      'ison_fn(',i,',',k,',',jrow-1,',1) = ',ison_fn(i,k,j-1,1)
            write(*,*)
     &      'ison_fn(',i,',',k,',',jrow-1,',2) = ',ison_fn(i,k,j-1,2)
            write(*,*)
     &      'ison_fb(',i,',',k,',',jrow,',1)   = ',ison_fb(i,k,j,1)
            write(*,*)
     &      'ison_fb(',i,',',k,',',jrow,',2)   = ',ison_fb(i,k,j,2)
            write(*,*)
     &      'ison_fb(',i,',',k-1,',',jrow,',1) = ',ison_fb(i,k-1,j,1)
            write(*,*)
     &      'ison_fb(',i,',',k-1,',',jrow,',2) = ',ison_fb(i,k-1,j,2)
            write(*,*)
     &      'cabbel(',i,',',jrow,',',k,')      = ',cabbel(i,jrow,k) 
            write(*,*)
     &      'thermob(',i,',',jrow,',',k,')     = ',thermob(i,jrow,k) 
            write(*,*)
     &      'halob(',i,',',jrow,',',k,')       = ',halob(i,jrow,k) 
          enddo
        enddo
      enddo
# endif

c
c     compute the advrho buffer, and add mask values for Ferret
c
      do k=1,km
        do j=js,je
          jrow = j+joff
          do i=1,imt
            buf_advrho(i,jrow,k) = advrho(i,k,j)
     &          + (1-tmask(i,k,j))*field_advrho%missing
            diffvert(i,jrow,k) = diffvert(i,jrow,k)
     &          + (1-tmask(i,k,j))*field_diffvert%missing
            diffhorz(i,jrow,k) = diffhorz(i,jrow,k)
     &          + (1-tmask(i,k,j))*field_diffhorz%missing
            diffgmskew(i,jrow,k) = diffgmskew(i,jrow,k)
     &          + (1-tmask(i,k,j))*field_diffgmskew%missing
            diffbihskew(i,jrow,k) = diffbihskew(i,jrow,k)
     &          + (1-tmask(i,k,j))*field_diffbihskew%missing
            cabbel(i,jrow,k) = cabbel(i,jrow,k)
     &          + (1-tmask(i,k,j))*field_cabbel%missing
            thermob(i,jrow,k) = thermob(i,jrow,k)
     &          + (1-tmask(i,k,j))*field_thermob%missing
            halob(i,jrow,k) = halob(i,jrow,k)
     &          + (1-tmask(i,k,j))*field_halob%missing
          enddo
	enddo
      enddo

c----------------------------------------------------------------------
c     Finsished with computation within the mw 
c----------------------------------------------------------------------
      endif

c----------------------------------------------------------------------
c     Data is written only outside window loop for himem
c     It is written within window loop for lomem
c----------------------------------------------------------------------

# if defined diagnostics_himem
      if (.not. within_window_loop) then
# endif
        call mpp_write(io, field_advrho, domain, buf_advrho, timrec)
        call mpp_write(io, field_diffvert, domain, diffvert, timrec)
        call mpp_write(io, field_diffhorz, domain, diffhorz, timrec)
        call mpp_write(io, field_diffgmskew, domain, diffgmskew, timrec)
        call mpp_write(io, field_diffbihskew, domain,diffbihskew,timrec)
        call mpp_write(io, field_cabbel, domain, cabbel, timrec)
        call mpp_write(io, field_thermob, domain, thermob, timrec)
        call mpp_write(io, field_halob, domain, halob, timrec)

c-----------------------------------------------------------------------
c       write out averaging period 
c-----------------------------------------------------------------------
        call mpp_write(io, field_period, period, timrec)

# if defined diagnostics_himem
      endif
# endif

      if(js+joff .eq. jscomp) then
          write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8,a,i8)')
     &'   ==>  Local potential density terms written to file '
     &,   file_name, 'on ts = ', itt,' ', stamp
     &,'   between latitudes:',jscomp,' and ', jecomp
     &,'   and from the surface down to depth level ',km
      endif

      return
      end
#else
      return
      end
#endif



