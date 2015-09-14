      subroutine baroclinic (joff, js, je, is, ie)
c
c=======================================================================
c     compute internal mode velocity components for rows js through je
c     in the MW.
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c=======================================================================
c
      use arrays_2d
      use memory_window       
      use sbc_info
      use sbc_ocean
#include "param.h"
      parameter (istrt=2, iend=imt-1)
#if defined bbl_ag
      dimension bdrag(imt,jmw)
# include "bbl.h"
#endif
#ifdef neptune
# include "cnep.h"
#endif
#include "coord.h"
#include "grdvar.h"
#include "hmixc.h"
#include "index.h"
#include "emode.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "vmixc.h"
      dimension baru(imt,2:jmw-1,2)
#include "fdifm.h"
c
c-----------------------------------------------------------------------
c     bail out if starting row exceeds ending row
c-----------------------------------------------------------------------
c
      if (js .le. je) then
c
# ifdef timing
      call tic ('ocean', 'baroclinic')
      call tic ('baroclinic', 'misc')
# endif
c
c-----------------------------------------------------------------------
c     limit the longitude indices based on those from the argument list
c     Note: these are currently bypassed. istrt and iend are set as
c           parameters to optimize performance
c-----------------------------------------------------------------------
c
c      istrt = max(2,is)
c      iend  = min(imt-1,ie)
c
c-----------------------------------------------------------------------
c     solve for one component of velocity at a time
c     n = 1 => zonal component
c     n = 2 => meridional component
c-----------------------------------------------------------------------
c
      do n=1,2
c
#if !defined linearized_advection
c
c-----------------------------------------------------------------------
c       calculate 2*advective flux (for speed) across east and north 
c       face of U cells.
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=1,km
            do i=istrt-1,iend
	      adv_fe(i,k,j) = adv_veu(i,k,j)*(u(i,  k,j,n,tau) + 
     &                                        u(i+1,k,j,n,tau))
            enddo
          enddo
        enddo
        do j=js-1,je
          do k=1,km
            do i=istrt,iend
	      adv_fn(i,k,j) = adv_vnu(i,k,j)*(u(i,k,j  ,n,tau) + 
     &                                        u(i,k,j+1,n,tau))
            enddo
          enddo
        enddo
#endif
c
c-----------------------------------------------------------------------
c       Diffusive momentum flux across east and north faces of U cells
c-----------------------------------------------------------------------
c
c       east face
c
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt-1,iend              
#if defined velocity_horz_laplacian
              diff_fe(i,k,j) = (u(i+1,k,j,n,taum1) - u(i,k,j,n,taum1)                        
# if defined neptune
     &                        - unep(i+1,jrow,n)*umask(i+1,k,j)
     &                        + unep(i,jrow,n)*umask(i,k,j)
# endif
#endif
#if defined velocity_horz_biharmonic
              diff_fe(i,k,j) = (del2_vel(i+1,k,j,n) - del2_vel(i,k,j,n)
#endif
     &                         )
#if defined velocity_horz_mix_const
     &                          *am_csudxtr(i+1,j)
#endif
#if defined velocity_horz_mix_var
     &                          *visc_ceu(i,k,j)*csudxtr(i+1,j)
#endif
#if defined partial_cell
     &                          *dhue(i,k,j)          
#endif
            enddo
          enddo
        enddo
c
c       north face
c
        do j=js-1,je
          jrow = j + joff
	  jp1 = min(jrow+1,jmt-1) - joff
#if defined velocity_horz_mix_const
          am_cst_dytr = visc_cnu*cst_dytr(jrow+1)
#endif
          do k=1,km
            do i=istrt,iend
#if defined velocity_horz_laplacian
              diff_fn(i,k,j) = (u(i,k,j+1,n,taum1) - u(i,k,j,n,taum1)
# if defined neptune
     &                        - unep(i,jrow+1,n)*umask(i,k,j+1) 
     &                        + unep(i,jrow,n)*umask(i,k,j)
# endif
#endif
#if defined velocity_horz_biharmonic
              diff_fn(i,k,j) = (del2_vel(i,k,jp1,n) - del2_vel(i,k,j,n)
#endif
     &                         )
#if defined velocity_horz_mix_const
     &                         *am_cst_dytr
#endif
#if defined velocity_horz_mix_var
     &                         *visc_cnu(i,k,j)*cst_dytr(jrow+1) 
#endif
#if defined partial_cell
     &                         *dhun(i,k,j)
#endif
            enddo
          enddo
        enddo
c
#if defined velocity_horz_mix_var && defined velocity_horz_laplacian
c
c-----------------------------------------------------------------------
c       Compute metric term due to variable viscosity. 
c-----------------------------------------------------------------------
c
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
# if defined am_cosine || defined am_taper_highlats
              a_lambda = c0
# else 
              a_lambda = (am_varh(i,k,j-1) - am_varh(i-1,k,j-1)
     &                  + am_varh(i,k,j)   - am_varh(i-1,k,j))
     &                  *dxu2r(i)
# endif
              a_phi    = (am_varh(i,k,j)   - am_varh(i,k,j-1)
     &                  + am_varh(i-1,k,j) - am_varh(i-1,k,j-1))
     &                  *dyu2r(jrow)
c
              varh_met(i,k,j) = a_phi*(
     &              varh3(jrow,n)*dxmetr(i)
     &             *(u(i+1,k,j,3-n,taum1)-u(i-1,k,j,3-n,taum1))
     &             + varh4(jrow,n)*u(i,k,j,n,taum1)
     &                                )
# if !defined am_cosine && !defined am_taper_highlats
     &                         +a_lambda*(
     &              varh1(jrow,n)
     &             *(u(i,k,j+1,3-n,taum1)-u(i,k,j-1,3-n,taum1))
     &             +varh2(jrow,n)*u(i,k,j,3-n,taum1)
     &                                   )
# endif
            enddo
          enddo
        enddo
#endif
c
#if defined velocity_horz_mix_var && defined velocity_horz_biharmonic
c
c-----------------------------------------------------------------------
c       Compute metric term due to variable viscosity. 
c-----------------------------------------------------------------------
c
        do j=js,je
          jrow = j + joff
          jp1 = min(jrow+1,jmt-1) - joff
          do k=1,km
            do i=istrt,iend
# if defined am_cosine || defined am_taper_highlats
              a_lambda = c0
# else 
              a_lambda = (am_varh(i,k,j-1) - am_varh(i-1,k,j-1)
     &                  + am_varh(i,k,j)   - am_varh(i-1,k,j))
     &                  *dxu2r(i)
# endif
              a_phi    = (am_varh(i,k,j)   - am_varh(i,k,j-1)
     &                  + am_varh(i-1,k,j) - am_varh(i-1,k,j-1))
     &                  *dyu2r(jrow)
              varh_met(i,k,j) = a_phi*(
     &              varh3(jrow,n)*dxmetr(i)
     &             *(del2_vel(i+1,k,j,3-n)-del2_vel(i-1,k,j,3-n))
     &             + varh4(jrow,n)*del2_vel(i,k,j,n)
     &                                )
# if !defined am_cosine && !defined am_taper_highlats
     &                         +a_lambda*(
     &              varh1(jrow,n)
     &             *(del2_vel(i,k,jp1,3-n)-del2_vel(i,k,j-1,3-n))
     &             + varh2(jrow,n)*del2_vel(i,k,j,3-n)
     &                                   )
# endif
            enddo
          enddo
        enddo
#endif
c
c-----------------------------------------------------------------------
c       calculate 2*advective flux (2 absorbed for speed) on
c       bottom face of U cells. 
c       Also compute the diffusive flux on bottom face of U cells
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=1,kmm1
            do i=istrt,iend
#if !defined linearized_advection
	      adv_fb(i,k,j)  = adv_vbu(i,k,j)*(u(i,k,  j,n,tau) +
     &                                         u(i,k+1,j,n,tau))
#endif
              diff_fb(i,k,j) = visc_cbu(i,k,j)*
     &                         (u(i,k,j,n,taum1) - u(i,k+1,j,n,taum1))
#ifdef partial_cell
     &                      /min(dhwt(i,k,j),   dhwt(i+1,k,j)
     &,                           dhwt(i,k,j+1), dhwt(i+1,k,j+1))
#else
     &                        *dzwr(k)
#endif
            enddo
          enddo
        enddo
#ifdef partial_cell
c
c       compute momentum sink by partial-cells
c
        if (n .eq. 1) then
          do j=js,je
            jrow = j + joff
            do k=1,km
              do i=istrt,iend              
                pc_sink(i,k,j) =
# if defined velocity_horz_mix_var
     &           - visc_ceu(i,k,j)*         
# else
     &           - visc_ceu*         
# endif
     &          (csudxur(i,j)/dhu(i,k,j))*
     &          (csudxtr(i+1,j)*(dhu(i,k,j)-dhue(i,k,j))
     &         + csudxtr(i,j)*(dhu(i,k,j)-dhue(i-1,k,j)))
# if defined velocity_horz_mix_var
     &         - visc_cnu(i,k,j)*
# else
     &         - visc_cnu*
# endif
     &          (csudyur(jrow)/dhu(i,k,j))*
     &          (cst_dytr(jrow+1)*(dhu(i,k,j)-dhun(i,k,j))
     &         + cst_dytr(jrow)*(dhu(i,k,j)-dhun(i,k,j-1)))
              enddo
            enddo
            call setbcx (pc_sink(1,1,j), imt, km)
          enddo
        endif
#endif
c
c-----------------------------------------------------------------------
c       set surface and bottom vert b.c. on U cells for mixing
c       and advection.
c       note: the b.c. at adv_fb(i,k=bottom,j) is set by the above code.
c             However, it is not set when k=km so it is set below.
c             adv_fb(i,km,j) is always zero (to within roundoff).
c-----------------------------------------------------------------------
c
        do j=js,je
          jrow = j + joff
          do i=istrt,iend
	    kb              = kmu(i,jrow)
# ifdef bbl_ag
	    kz              = max(kmu(i,jrow),1)
	    bmf(i,j,n) = (Kvbbltop*(u(i,kz,j,n,taum1)-u(i,km,j,n,taum1))
     &                    *dztr(km))*umask(i,km,j)
# endif
            diff_fb(i,0,j)  = smf(i,j,n)
            diff_fb(i,kb,j) = bmf(i,j,n)
	    adv_fb(i,0,j)   = adv_vbu(i,0,j)*(u(i,1,j,n,tau) +
     &                                        u(i,1,j,n,tau)) 
	    adv_fb(i,km,j)  = adv_vbu(i,km,j)*u(i,km,j,n,tau) 
          enddo
        enddo
c
#ifdef source_term
c
c-----------------------------------------------------------------------
c       set source term for U cells
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=1,km
            do i=istrt,iend
	      source(i,k,j) = c0
            enddo
          enddo
        enddo
#endif
c
c-----------------------------------------------------------------------
c       solve for the internal mode part of du/dt at center of 
c       U cells by neglecting the surface pressure gradients. use
c       statement functions to represent each component of the 
c       calculation. 
c-----------------------------------------------------------------------
c
        do j=js,je
          jrow = j + joff
          do i=istrt,iend
            zu(i,jrow,n) = c0
#if defined explicit_free_surface && !defined linearized_advection
            zu(i,jrow,n) =  (p5*adv_fb(i,0,j)
# if defined explicit_fresh_water_flux  
     &                       +sffu(i,jrow)*u(i,1,j,n,tau)
# endif  
     &                      )*umask(i,1,j)
#endif
          enddo
          do k=1,kbot
            do i=istrt,iend
              coritemp = CORIOLIS(i,k,j,n)
# if defined explicit_free_surface
	      zu(i,jrow,n) = zu(i,jrow,n)-coritemp*umask(i,k,j)
#  if defined partial_cell
     &                       *dhu(i,k,j)
#  else
     &                       *dzt(k)
#  endif
# endif
	      u(i,k,j,n,taup1) = 
     &            (DIFF_Ux(i,k,j) + DIFF_Uy(i,k,j) + DIFF_Uz(i,k,j) 
     &             + DIFF_metric(i,k,j,n)
#  if !defined linearized_advection
     &             - ADV_Ux(i,k,j) - ADV_Uy(i,k,j) - ADV_Uz(i,k,j)
     &             + ADV_metric(i,k,j,n) 
#  endif
     &             - grad_p(i,k,j,n) + coritemp
#  ifdef source_term
     &             + source(i,k,j)
#  endif
     &              )*umask(i,k,j)
            enddo
          enddo
#ifdef obc
# include "baroclinic_obc.inc"
#endif
        enddo
c
#if defined bbl_ag
c
c-----------------------------------------------------------------------
c       solve for the internal mode part of du/dt at center of U-cells
c       in the bbl (level km) 
c-----------------------------------------------------------------------
c
        if (n .eq. 1) then
          do j=js,je
            do i=istrt,iend
              speed = sqrt(u(i,km,j,1,taum1)**2 + u(i,km,j,2,taum1)**2)
	      bdrag(i,j) = cdbotbbl*speed*dztr(km)
            enddo
	  enddo
	endif
        do j=js,je
          jrow = j + joff
          do k=km,km
            do i=istrt,iend
              coritemp = CORIOLIS(i,k,j,n)
# if defined explicit_free_surface
	      zu(i,jrow,n) = zu(i,jrow,n)-coritemp*umask(i,k,j)
#  if defined partial_cell
     &                       *dhu(i,k,j)
#  else
     &                       *dzt(k)
#  endif
# endif
	      u(i,k,j,n,taup1) = (
     &             - raydrag*u(i,k,j,n,taum1)
     &             - bdrag(i,j)*u(i,k,j,n,taum1)
     &             - grad_p(i,k,j,n) + coritemp
     &             + bmf(i,j,n)*dztr(k)
     &              )*umask(i,k,j)
            enddo
          enddo
        enddo
#endif
c       
#ifdef implicitvmix
c
c-----------------------------------------------------------------------
c       add in du/dt component due to implicit vertical diffusion
c-----------------------------------------------------------------------
c
        call ivdifu (joff, js, je, istrt, iend, n)
#endif
c
c-----------------------------------------------------------------------
c       construct diagnostics associated with velocity component "n" 
c-----------------------------------------------------------------------
c
        call diagc1 (joff, js, je, istrt, iend, n)
c
c-----------------------------------------------------------------------
c       construct the vertical average of du/dt for forcing
c       the barotropic equation
c-----------------------------------------------------------------------
c
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              zu(i,jrow,n) = zu(i,jrow,n) + u(i,k,j,n,taup1)
#ifdef partial_cell
     &                       *dhu(i,k,j)
#else
     &                       *dzt(k)
#endif
            enddo
          enddo
        enddo
c
        do j=js,je
          jrow = j + joff
          do i=istrt,iend
#if !defined explicit_free_surface 
            zu(i,jrow,n) = zu(i,jrow,n)*hr(i,jrow)
#endif
          enddo
        enddo
#ifdef symmetry
        do j=js,je
          jrow = j + joff
          if (jrow .eq. jmtm1 .and. n .eq. 2) then
            do i=istrt,iend
              zu(i,jrow,2)   =  c0
	      zu(i,jrow+1,2) = -zu(i,jrow-1,2)
	      zu(i,jrow+1,1) =  zu(i,jrow-1,1)
            enddo
          endif
        enddo
#endif
c
c-----------------------------------------------------------------------
c       end of velocity component "n" loop
c-----------------------------------------------------------------------
c
      enddo
c
c-----------------------------------------------------------------------
c     compute "tau+1" velocities accounting for implicit part of the
c     coriolis term if treated implicitly. velocities are in error by an
c     arbitrary constant related to neglecting the unknown surface
c     pressure gradients
c-----------------------------------------------------------------------
c
#if defined damp_inertial_oscillation
      do j=js,je
        jrow  = j + joff
# ifdef rot_grid
        do k=1,km
          do i=istrt,iend
          fx    = c2dtuv*acor*cori(i,jrow,1)
          detmr = c1/(c1 + fx*fx)
# else
        fx    = c2dtuv*acor*cori(jrow,1)
        detmr = c1/(c1 + fx*fx)
        do k=1,km
          do i=istrt,iend
# endif            
            t1 = (u(i,k,j,1,taup1) + fx*u(i,k,j,2,taup1))*detmr
            t2 = (u(i,k,j,2,taup1) - fx*u(i,k,j,1,taup1))*detmr
            u(i,k,j,1,taup1) = u(i,k,j,1,taum1) + c2dtuv*t1
            u(i,k,j,2,taup1) = u(i,k,j,2,taum1) + c2dtuv*t2
          enddo
	enddo
      enddo
#else
      do n=1,2
	do j=js,je
          do k=1,km
            do i=istrt,iend
              u(i,k,j,n,taup1) = u(i,k,j,n,taum1)
     &                            + c2dtuv*u(i,k,j,n,taup1)
            enddo
	  enddo
	enddo
      enddo
#endif
#if defined firfil
c
c-----------------------------------------------------------------------
c     filter velocity components at high latitudes
c-----------------------------------------------------------------------
c
      if (istrt .eq. 2 .and. iend .eq. imt-1) then
        call filuv (joff, js, je)
      else
        write (stdout,'(a)') 
     &  'Error: filtering requires is=2 and ie=imt-1 in baroclinic'
        call abort()
      endif
#endif
      do j=js,je
        call setbcx (u(1,1,j,1,taup1), imt, km)
        call setbcx (u(1,1,j,2,taup1), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     subtract vertical mean to get internal mode.
c-----------------------------------------------------------------------
c
      do n=1,2
        do j=js,je
          do i=istrt,iend
            baru(i,j,n) = c0
          enddo
        enddo
        do j=js,je
          do k=1,km
            do i=istrt,iend
              baru(i,j,n) = baru(i,j,n) + u(i,k,j,n,taup1)
#ifdef partial_cell
     &                                    *dhu(i,k,j)
#else
     &                                    *dzt(k)
#endif
            enddo
	  enddo
        enddo
        do j=js,je
          jrow  = j + joff
          do i=istrt,iend
#if defined explicit_free_surface && defined time_averaged_surface
            baru(i,j,n) = baru(i,j,n)*depthur(i,jrow)
#else
            baru(i,j,n) = baru(i,j,n)*hr(i,jrow)
#endif
	  enddo
        enddo
        do j=js,je
          do k=1,km
            do i=istrt,iend
              u(i,k,j,n,taup1) = u(i,k,j,n,taup1)
     &                          - umask(i,k,j)*baru(i,j,n)
            enddo
	  enddo
          call setbcx (u(1,1,j,n,taup1), imt, km)
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     construct diagnostics involving internal mode velocity at "tau+1"
c-----------------------------------------------------------------------
c
      call diagc2 (joff, js, je, is, ie)
c
#if defined fourfil
c
c-----------------------------------------------------------------------
c     filter velocity components at high latitudes
c-----------------------------------------------------------------------
c
      if (istrt .eq. 2 .and. iend .eq. imt-1) then
        call filuv (joff, js, je)
      else
        write (stdout,'(a)') 
     &  'Error: filtering requires is=2 and ie=imt-1 in baroclinic'
        call abort()
      endif
#endif
      do j=js,je
        call setbcx (u(1,1,j,1,taup1), imt, km)
        call setbcx (u(1,1,j,2,taup1), imt, km)
      enddo
c
#if defined coupled
c
c-----------------------------------------------------------------------
c     if needed, construct the Atmos S.B.C.(surface boundary conditions)
c     averaged over this segment      
c     reference: Pacanowski, R.C., Effect of Equatorial Currents
c                on Surface Stress (JPO, Vol 17, No. 6, June 1987)
c-----------------------------------------------------------------------
c
# if defined robert_time_filter
c
c  note: the following call to "asbcu" will not preserve answers across
c        restarts when option robert_time_filter is enabled. The reason
c        is that on non-restart time steps, "ascbu" is working with
c        unfiltered values of u(tau). However, filtered values of u(tau)
c        are written to a restart file. Therefore, when starting from
c        a restart file, filtered values of u(tau) are loaded from disk
c        and used by "asbcu" on the first time step of a restart.
c
# endif
      iu = index_of_sbc ('surface u')
      iv = index_of_sbc ('surface v')
      if (iu .ne. 0 .and. iv .ne. 0) then
        call asbcu (joff, js, je, istrt, iend, iu, iv)
      endif
#endif
c
#if defined debug_baroclinic
      do j=js,je
	jrow = j+joff
        write(stdout,*) 'In baroclinic: j=',j,' jrow=',jrow
        call print_checksum (u(1,1,j,1,taup1), imt, km,'ck u(taup1) = ')
        call print_checksum (u(1,1,j,2,taup1), imt, km,'ck v(taup1) = ')
      enddo
#endif
#ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     & "=> In baroclinic: js=",js," je=",je," joff=",joff
     &,", calculating latitude rows: ",js+joff," to ",je+joff
#endif
#ifdef timing
      call toc ('baroclinic', 'misc')
      call toc ('ocean', 'baroclinic')
#endif

      endif

      return
      end


      subroutine pressure_gradients (joff, js, je)
c
c-----------------------------------------------------------------------
c     construct the hydrostatic pressure gradients: 1 = dp/dx; 2 = dp/dy
c-----------------------------------------------------------------------
c
      use arrays_2d
      use memory_window       
      use sbc_info
      use sbc_ocean
#include "param.h"
      parameter (istrt=2, iend=imt-1)
#include "coord.h"
#include "grdvar.h"
#include "mw.h"
#include "scalar.h"
#if defined bbl_ag
# include "bbl.h"
      dimension px_hor(imt,0:jmw), py_hor(imt,0:jmw)
#endif
c
      js_press = max(js,1+jbuf)-1
      je_press = je-1
c
      do j=js_press,je_press
        jrow = j + joff
        fxa  = rho0r*csur(jrow)
        fxb  = rho0r*dyu2r(jrow)
	do k=1,kbot
          do i=istrt-1,iend
            t1              = pressure(i+1,k,j+1) - pressure(i  ,k,j)
            t2              = pressure(i  ,k,j+1) - pressure(i+1,k,j)
            grad_p(i,k,j,1) = (t1-t2)*fxa*dxu2r(i)*umask(i,k,j)
            grad_p(i,k,j,2) = (t1+t2)*fxb*umask(i,k,j)
          enddo
        enddo
      enddo
c
#if defined partial_cell
c
c     add in correction due to variable thickness partial cell
c
      grav_rho0r = grav*rho0r
      do j=js_press,je_press
        jrow = j + joff
        fxa = grav_rho0r*csur(jrow)
        fxb = grav_rho0r*dyu4r(jrow)
        do i=istrt-1,iend
          k = kmu(i,jrow)
	  if (k .gt. 1) then
	    rhoij     = rho(i,k,j)
	    rhoip1    = rho(i+1,k,j)
	    rhojp1    = rho(i,k,j+1)
	    rhoip1jp1 = rho(i+1,k,j+1)
            grad_p(i,k,j,1) = grad_p(i,k,j,1) - umask(i,k,j)*fxa
     &         *dxu4r(i)*((rhoip1 + rhoij)*(ztp(i+1,k,j) - ztp(i,k,j))
     &           +(rhoip1jp1 + rhojp1)*(ztp(i+1,k,j+1) - ztp(i,k,j+1)))
            grad_p(i,k,j,2) = grad_p(i,k,j,2) - umask(i,k,j)*fxb*(  
     &                    (rhojp1 + rhoij)*(ztp(i,k,j+1) - ztp(i,k,j))
     &           +(rhoip1jp1 + rhoip1)*(ztp(i+1,k,j+1) - ztp(i+1,k,j)))
     &                        
          endif
	enddo
      enddo
#endif
c
      do j=js_press,je_press
        call setbcx (grad_p(1,1,j,1), imt, km)
        call setbcx (grad_p(1,1,j,2), imt, km)
      enddo
c
#if defined bbl_ag
c
c-----------------------------------------------------------------------
c     Horizontal pressure gradients on U points in BBL (k=km)
c-----------------------------------------------------------------------
c
      jstbe = max(js,1+jbuf)
      jsu   = jstbe - 1
c
c     zonal pressure difference
c
      p5grav = p5*grav
      do j=js,je
        jrow = j+joff
        do i=istrt,iend
          kb     = max(kmt(i,jrow),1)
          kb1    = max(kmt(i+1,jrow),1)
c
c         zonal pressure difference "px_hor", to be used only in case of
c         detrainement, between
c         BBL cell and interior cell (when kb <> kb1)
c         BBL cell and adjacent BBL cell (when kb = kb1)
c         
          kbhor  = min(kb,kb1)
          px_hor(i,j) = 0.0
          if (kbhor .gt. 1) then
            if (kb .lt. kb1) then
              pcorr =-p5grav*((rho(i+1,kbhor+1,j) + rho(i,kbhor+1,j))*
# if defined partial_cell
     &                        (ztp(i+1,kbhor+1,j) - zt_bbl(i,j)))
# else
     &                        (zt(kbhor+1) - zt_bbl(i,j)))
# endif
	    elseif (kb .gt. kb1) then
              pcorr =-p5grav*((rho(i+1,kbhor+1,j) + rho(i,kbhor+1,j))*
# if defined partial_cell
     &                        (zt_bbl(i+1,j) - ztp(i,kbhor+1,j)))
# else
     &                        (zt_bbl(i+1,j) - zt(kbhor+1)))
# endif
	    else
              pcorr =-p5grav*((rho(i+1,kbhor+1,j) + rho(i,kbhor+1,j))*
     &                        (zt_bbl(i+1,j) - zt_bbl(i,j)))
	    endif
c
            px_hor(i,j) = pressure(i+1,kbhor+1,j) - 
     &                    pressure(i  ,kbhor+1,j) + pcorr
          endif
c
c         determine detrainment
c         oppsignx = 1.0 => no detrainment
c         oppsignx = 0.0 => detrainment
c
          if (kb1 .gt. kb) then
            if (rho(i+1,km,j) .gt. rho(i,km,j)) then
              oppsignx(i,j) = 0.0
            else
              oppsignx(i,j) = 1.0
            endif
          elseif (kb1 .lt. kb) then
            if (rho(i+1,km,j) .lt. rho(i,km,j)) then
              oppsignx(i,j) = 0.0
            else
              oppsignx(i,j) = 1.0
            endif
          else
            oppsignx(i,j) = 1.0
          endif
#ifdef bbl_always_detrain
          if (kb1 .ne. kb) oppsignx(i,j) = 0.0
#endif
#ifdef bbl_never_detrain
          oppsignx(i,j) = 1.0
#endif
c
# if defined beckmann_doescher
          px_bbl(i,j) = px_hor(i,j)
# else
c
c         The pressure difference is within the BBL or between the BBL
c         and the interior if detraining.
c
          rhom = (rho(i,km,j)+rho(i+1,km,j))*p5
          px_bbl(i,j) = (pressure(i+1,km,j)-pressure(i,km,j)
     &                   -grav*rhom*(zt_bbl(i+1,j)-zt_bbl(i,j))
     &                  )*tmask(i+1,km,j)*oppsignx(i,j)
     &                   + px_hor(i,j)*(1.0-oppsignx(i,j))
# endif
        enddo
        call setbcx (px_bbl(1,j), imt, 1)
        call setbcx (oppsignx(1,j), imt, 1)
      enddo
c
c     meridional pressure difference
c
      do j=jsu,je-1
        jrow = max(j + joff,1)
        do i=istrt,iend
          kb  = max(kmt(i,jrow),1)
          kb1 = max(kmt(i,jrow+1),1)
c
c         meridional pressure difference "py_hor", to be used only in case of
c         detrainement,  between
c         BBL cell and interior cell (when kb <> kb1)
c         BBL cell and adjacent BBL cell (when kb = kb1)
c
          kbhor = min(kb,kb1)
          py_hor(i,j) = 0.0
          if (kbhor .gt. 1) then
c
            if (kb .lt. kb1) then
              pcorr =-p5grav*((rho(i,kbhor+1,j+1)+rho(i,kbhor+1,j))*
# if defined partial_cell
     &                        (ztp(i,kbhor+1,j+1) - zt_bbl(i,j)))
# else
     &                        (zt(kbhor+1) - zt_bbl(i,j)))
# endif
	    elseif (kb .gt. kb1) then
              pcorr =-p5grav*((rho(i,kbhor+1,j+1)+rho(i,kbhor+1,j))*
# if defined partial_cell
     &                        (zt_bbl(i,j+1) - ztp(i,kbhor+1,j)))
# else
     &                        (zt_bbl(i,j+1) - zt(kbhor+1)))
# endif
	    else
              pcorr =-p5grav*((rho(i,kbhor+1,j+1)+rho(i,kbhor+1,j))*
     &                        (zt_bbl(i,j+1) - zt_bbl(i,j)))
	    endif
c
            py_hor(i,j) = pressure(i,kbhor+1,j+1)
     &                  - pressure(i,kbhor+1,j) + pcorr   
          endif
c
c         determine detrainment
c         oppsigny = 1.0 => no detrainment
c         oppsigny = 0.0 => detrainment
c
          if (kb1 .gt. kb) then
            if (rho(i,km,j+1) .gt. rho(i,km,j)) then
              oppsigny(i,j) = 0.0
            else
              oppsigny(i,j) = 1.0
            endif
          elseif (kb1 .lt. kb) then
            if (rho(i,km,j+1) .lt. rho(i,km,j)) then
              oppsigny(i,j) = 0.0
            else
              oppsigny(i,j) = 1.0
            endif
          else
            oppsigny(i,j) = 1.0
          endif
#ifdef bbl_always_detrain
          if (kb1 .ne. kb) oppsigny(i,j) = 0.0
#endif
#ifdef bbl_never_detrain
          oppsigny(i,j) = 1.0
#endif
c
# if defined beckmann_doescher
          py_bbl(i,j) = py_hor(i,j)
# else
c
c         The pressure difference is within the BBL or between the BBL
c         and the interior if detraining.
c
          rhom = (rho(i,km,j+1)+rho(i,km,j))*p5
          py_bbl(i,j) = (pressure(i,km,j+1)-pressure(i,km,j)
     &                   -grav*rhom*(zt_bbl(i,j+1)-zt_bbl(i,j))
     &                  )*tmask(i,km,j+1)*oppsigny(i,j)
     &                  + py_hor(i,j)*(1.0-oppsigny(i,j))
# endif
        enddo
        call setbcx (py_bbl(1,j), imt, 1)
        call setbcx (oppsigny(1,j), imt, 1)
      enddo
c
c-----------------------------------------------------------------------
c     construct the horizontal pressure gradient on U points in the BBL
c-----------------------------------------------------------------------
c
      p5_rho0r = p5*rho0r
      do j=jsu,je-1
        jrow = j + joff
        do i=istrt-1,iend
          grad_p(i,km,j,1) = (px_bbl(i,j)+px_bbl(i,j+1))
     &         *p5_rho0r*csudxur(i,j)
        enddo
        call setbcx (grad_p(1,km,j,1), imt, 1)
      enddo
      do j=jsu,je-1
        jrow = max(j + joff,1)
        do i=istrt-1,iend
          grad_p(i,km,j,2) = (py_bbl(i,j)+py_bbl(i+1,j))
     &         *p5_rho0r*dyur(jrow)
        enddo
        call setbcx (grad_p(1,km,j,2), imt, 1)
      enddo
#endif
      return
      end



      subroutine diagc1 (joff, js, je, is, ie, n)
c
c-----------------------------------------------------------------------
c     construct diagnostics which don`t require internal mode velocity
c     at "tau+1" for each velocity component "n"
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       n    = (1,2) = (u,v) velocity component
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      use mpp_io_mod
      use memory_window
#include "param.h"
#include "coord.h"
#include "diag.h"
#include "grdvar.h"
#include "hmixc.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "vmixc.h"
      dimension temp(imt,km)
c
#ifdef timing
        if (n .eq. 1) then
	  call tic ('baroclinic', 'diagnostics')
	else
	  call ticr ('baroclinic', 'diagnostics')
	endif
#endif
c
#if defined save_mixing_coeff
c
c-----------------------------------------------------------------------
c     diagnostic: estimate mixing coefficients on east, north, and
c                 bottom face of U cells from the flux
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (cmixts .and. n .eq. 1 .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'mixing coeff')
# endif
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=2,imt-1
# if defined velocity_horz_mix_var
              ce(i,k,j,1) = visc_ceu(i,k,j)*umask(i,k,j)*umask(i+1,k,j)
              cn(i,k,j,1) = visc_cnu(i,k,j)*umask(i,k,j)*umask(i,k,j+1)
# endif
# if defined velocity_horz_laplacian && defined velocity_horz_mix_const 
              ce(i,k,j,1) = am*umask(i,k,j)*umask(i+1,k,j)
              cn(i,k,j,1) = am*umask(i,k,j)*umask(i,k,j+1)
#  endif
# if defined velocity_horz_biharmonic && defined velocity_horz_mix_const 
              ce(i,k,j,1) = ambi*umask(i,k,j)*umask(i+1,k,j)
              cn(i,k,j,1) = ambi*umask(i,k,j)*umask(i,k,j+1)
#  endif
            enddo
          enddo
        enddo
        do j=js,je
          jrow = j + joff
          do k=1,km-1
            do i=2,imt-1
              cb(i,k,j,1) = visc_cbu(i,k,j)*umask(i,k,j)*umask(i,k+1,j)
            enddo
          enddo
          do i=2,imt-1
            cb(i,km,j,1) = 0.0
          enddo
        enddo
c
        do j=js,je
          call setbcx (ce(1,1,j,1), imt, km)
          call setbcx (cn(1,1,j,1), imt, km)
          call setbcx (cb(1,1,j,1), imt, km)
        enddo
# ifdef timing
        call toc ('diagnostic', 'mixing coeff')
# endif
      endif
#endif
c
#ifdef time_step_monitor
c
c-----------------------------------------------------------------------
c     diagnostic: accumulate global kinetic energy on "tau" velocity
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c               (based on diagnostic by M. Cox)
c-----------------------------------------------------------------------
c
      if (tsits .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'time step monitor')
# endif
        do j=js,je
          jrow = j + joff
          fx = rho0*p5*csu(jrow)*dyu(jrow)
# ifdef symmetry
          if (jrow .eq. jmtm1) fx = fx*p5
# endif
          do k=1,km
            do i=is,ie
	      weight    = fx*dxu(i)
# ifdef partial_cell
     &                             *dhu(i,k,j)
# else
     &                             *dzt(k)
# endif
              temp(i,k) = u(i,k,j,n,tau)**2*weight
            enddo
            do i=is,ie
              ektot(k,jrow) = ektot(k,jrow) + temp(i,k)
            enddo
          enddo
        enddo
# ifdef timing
        call toc ('diagnostic', 'time step monitor')
# endif
      endif
#endif
c
#ifdef energy_analysis
c
c-----------------------------------------------------------------------
c     diagnostic: integrate work done by the r.h.s. terms in the
c                  momentum equations.
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c               (this is not done the same way as in MOM 1)
c-----------------------------------------------------------------------
c
      if (glents .and. eots) call ge1 (joff, js, je, is, ie, n)
#endif
c
#ifdef term_balances
c
c-----------------------------------------------------------------------
c     diagnostic: integrate r.h.s. terms in the momentum equations
c                 over specified regional volumes
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (trmbts .and. eots) call utb1 (joff, js, je, is, ie, n)
#endif
#ifdef save_xbts
c
c-----------------------------------------------------------------------
c     diagnostic: accumulate r.h.s terms in the momentum equation
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (xbtperts .and. eots) call uxbt1 (joff, js, je, n)
#endif
#ifdef timing
        call toc ('baroclinic', 'diagnostics')
#endif
      return
      end






      subroutine diagc2 (joff, js, je, is, ie)
c
c-----------------------------------------------------------------------
c     construct diagnostics requiring internal mode velocity at "tau+1"
c     and those not dependent on velocity component fluxes.
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      use arrays_2d
#include "param.h"
#include "scalar.h"
#include "switch.h"
c
#ifdef timing
      call ticr ('baroclinic', 'diagnostics')
#endif
#ifdef energy_analysis
c
c-----------------------------------------------------------------------
c     diagnostic: integrate work done by du/dt in the momentum equations
c                 the external mode part at "tau+1" will be 
c                 accounted for after the external mode is solved.
c                 also, integrate the work done by buoyancy.
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c               (this is not done the same way as in MOM 1)
c-----------------------------------------------------------------------
c
      if (glents .and. eots) then
        call ge2 (joff, js, je, is, ie, c2dtuv, grav, rho0r)
      endif
#endif
c
#ifdef term_balances
c
c-----------------------------------------------------------------------
c     diagnostic: add du/dt and implicit coriolis terms to the integrals
c                 over specified volumes. the external mode parts will
c                 be accounted for after the external mode is solved.
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (trmbts .and. eots) then
        call utb2 (joff, js, je, is, ie, c2dtuv, acor)
      endif
#endif
c
#ifdef save_xbts
c
c-----------------------------------------------------------------------
c     diagnostic: accumulate du/dt and implicit coriolis terms from the 
c                 momentum equations
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (xbtperts .and. eots) call uxbt2 (joff, js, je, c2dtuv, acor)
#endif
#ifdef timing
      call toc ('baroclinic', 'diagnostics')
#endif
      return
      end



#if !defined simple_sbc
      subroutine asbcu (joff, js, je, is, ie, iu, iv)
c
c-----------------------------------------------------------------------
c     construct the Atmos S.B.C.(surface boundary conditions)
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       iu   = index for u component 
c       iv   = index for v component 
c
c     reference: Pacanowski, R.C., Effect of Equatorial Currents
c                on Surface Stress (JPO, Vol 17, No. 6, June 1987)
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      use arrays_2d
      use memory_window
      use sbc_info
      use sbc_ocean
      integer tlev
# include "param.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# ifdef timing
      call tic ('baroclinic', 'constructing SBC')
# endif
c
c     initialize S.B.C. at the beginning of each ocean segment
c     (do not alter values in land)
c
      if (eots .and. osegs .and. iu .ne. 0 .and. iv .ne. 0) then
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) then
              sbcocn(i,jrow,iu) = c0
              sbcocn(i,jrow,iv) = c0
            endif
          enddo
        enddo
      endif
c
c     accumulate surface currents for the Atmos S.B.C. every time step
c
      if (eots .and. iu .ne. 0 .and. iv .ne. 0) then
        tlev = tau
        if (euler2) tlev = taum1
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            sbcocn(i,jrow,iu) = sbcocn(i,jrow,iu) + p25*(
     &                          u(i,1,j,1,tlev) + u(i-1,1,j,1,tlev) 
     &                        + u(i,1,j-1,1,tlev) + u(i-1,1,j-1,1,tlev))
            sbcocn(i,jrow,iv) = sbcocn(i,jrow,iv) + p25*(
     &                          u(i,1,j,2,tlev) + u(i-1,1,j,2,tlev) 
     &                        + u(i,1,j-1,2,tlev) + u(i-1,1,j-1,2,tlev))
          enddo
        enddo
      endif
c
c     average the surface currents for the Atmos S.B.C. at the end
c     of each ocean segment. (do not alter values in land)
c
      if (eots .and. osege .and. iu .ne. 0 .and. iv .ne. 0) then
        rts = c1/ntspos
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) then
              sbcocn(i,jrow,iu) = rts*sbcocn(i,jrow,iu)
	      sbcocn(i,jrow,iv) = rts*sbcocn(i,jrow,iv)
            endif
          enddo
        enddo
      endif
# ifdef timing
      call toc ('baroclinic', 'constructing SBC')
# endif
      return
      end
#endif





#ifdef implicitvmix
      subroutine ivdifu (joff, js, je, is, ie, n)
c
c-----------------------------------------------------------------------
c     solve vertical diffusion of velocity implicitly
c
c     input:
c       joff  = offset relating "j" in the MW to latitude "jrow"
c       js    = starting row in the MW
c       je    = ending row in the MW
c       is    = starting longitude index in the MW
c       ie    = ending longitude index in the MW
c       n     = velocity component
c       twodt = (2*dtuv, dtuv) on (leapfrog, mixing) time steps
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      use arrays_2d
      use memory_window
# include "param.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# include "vmixc.h"
      dimension twodt(km)
c
# ifdef timing
        call tic ('baroclinic', 'implicit vert diff')
# endif
c
c     construct the "tau+1" velocity without implicit vertical diffusion
c
      do j=js,je
        do k=1,km
          do i=is,ie
            u(i,k,j,n,taup1) = u(i,k,j,n,taum1)+c2dtuv*u(i,k,j,n,taup1)
          enddo
        enddo
      enddo
# if defined save_xbts || defined energy_analysis || defined term_balances
c
c     store terms to compute implicit vertical diffusion on
c     diagnostic time steps
c
      if ((xbtperts .or. glents .or. trmbts) .and. eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = u(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif
# endif
c
c     add in the implicit vertical diffusion
c
      do k=1,km
        twodt(k) = c2dtuv
      enddo
# ifdef tcvmix
      call invtri (up1(1,1,1,n), smf(1,1,n), bmf(1,1,n), vvca(???)
     &,            twodt, kmu, umask(1,1,1), is, ie, joff, js, je)
      if (imt-1 .gt. 1) then
        write (stdout,*) 'Error: ivdif.F is not converted for tcvmix'
	call abort()
      endif
# else
      call invtri (u(1,1,1,n,taup1), smf(1,1,n), bmf(1,1,n)
     &, visc_cbu(1,1,2), twodt, kmu, umask(1,1,1), is, ie
     &, joff, js, je, 'baroclinic')
# endif
c
      r2dtuv = c1/c2dtuv
# if defined save_xbts || defined energy_analysis || defined term_balances
c
c     compute residual implicit vertical diffusion for diagnostics
c
      if ((xbtperts .or. glents .or. trmbts) .and. eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = r2dtuv*(u(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
          enddo
        enddo
      endif
# endif
c
c     convert back to time change of velocity
c
      do j=js,je
        do k=1,km
          do i=is,ie
            u(i,k,j,n,taup1) =r2dtuv*(u(i,k,j,n,taup1)-u(i,k,j,n,taum1))
          enddo
        enddo
      enddo
# ifdef timing
        call toc ('baroclinic', 'implicit vert diff')
# endif
      return
      end
#endif

