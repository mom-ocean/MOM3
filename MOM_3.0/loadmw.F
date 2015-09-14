      subroutine loadmw (joff, js, je, is, ie, num1, nu)
c
c=======================================================================
c     load variables into the MW for rows "js" through "je"
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row for loading jrows into the MW
c       je   = ending row for loading jrows into the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       num1 = "tau-1" latitude disk unit
c       nu   = "tau" latitude disk unit
c     
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
      use arrays_2d
      use memory_window
#include "param.h"
#ifdef bbl_ag
# include "bbl.h"
#endif
#include "coord.h"
#include "emode.h"
#include "grdvar.h"
#include "hmixc.h"
#include "iounit.h"
#ifdef isoneutralmix
# include "isopyc.h"
#endif
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "taskrows.h"
#include "tmngr.h"
      dimension rho_full(imt,km,jmw)
      save taum1old
c
#ifdef timing
      call tic ('loadmw', 'loading the MW')
#endif
c
c-----------------------------------------------------------------------
c     limit the longitude indices
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
#ifdef trace_indices
      if (js .le. je) then
        write (stdout,'(1x,6(a,i4))')
     &   "=> In loadmw: itt=",itt,", loading latitude rows "
     &,  js+joff," to ",je+joff," into MW rows js=",js," to je=",je
     &,  ", joff=",joff
      endif
#endif
c
c-----------------------------------------------------------------------
c     for all MW`s after the first, move the MW northward by copying
c     data from the last few rows into the first few rows.
c     (Refer to the manual for specifics)
c-----------------------------------------------------------------------
c
      if (js+joff .ne. jstask) then
        call move_memory_window (istrt-1, iend+1)
#if defined bbl_ag
      else
        do i=1,imt
	  py_bbl(i,1)      = c0
	  oppsigny(i,1)    = c0
	enddo
# ifdef detrain
        do n=1,nt
          do i=1,imt
            adv_fdetn(i,1,n) = c0
	  enddo
	enddo
# endif
	do k=1,km
	  do i=1,imt
	    adv_vet(i,k,2) = c0
	    adv_vnt(i,k,1) = c0
	    adv_vnu(i,k,1) = c0
	  enddo
	enddo
	do k=0,km
	  do i=1,imt
	    adv_vbt(i,k,jmw) = c0
	    adv_vbu(i,k,jmw-1) = c0
	  enddo
	enddo
# if defined max_window && !defined second_order_window
        if (je+joff .eq. jmt) then
          do i=1,imt
	    px_bbl(i,jmw-1) = c0
	  enddo
	  do k=0,km
	    do i=1,imt
	      adv_vbu(i,k,jmw-2) = c0
	    enddo
	  enddo
	endif
# endif
#endif
      endif
#if defined max_window
      if (first) then
#endif
c
c-----------------------------------------------------------------------
c     construct grid factors to minimize computations 
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        do i=istrt-1,iend+1
	  csudxur(i,j)    = csur(jrow)*dxur(i)
	  csudxu2r(i,j)   = csur(jrow)*dxur(i)*p5
	  csudxtr(i,j)    = csur(jrow)*dxtr(i)
	  cstdxtr(i,j)    = cstr(jrow)*dxtr(i)
	  cstdxt2r(i,j)   = cstr(jrow)*dxtr(i)*p5
	  cstdxur(i,j)    = cstr(jrow)*dxur(i)
#if defined velocity_horz_mix_const
	  am_csudxtr(i,j) = visc_ceu*csur(jrow)*dxtr(i)
#endif
#if defined tracer_horz_mix_const
	  ah_cstdxur(i,j) = diff_cet*cstr(jrow)*dxur(i)
#endif
	enddo
      enddo
c
c-----------------------------------------------------------------------
c     construct T cell and U cell land/sea masks
c-----------------------------------------------------------------------
c
#ifdef timing
      call tic ('loadmw', 'getvar+mask+ext mode')
#endif
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt-1,iend+1
              if (kmt(i,jrow) .ge. k) then
                tmask(i,k,j) = c1
              else
                tmask(i,k,j) = c0
              endif
              if (kmu(i,jrow) .ge. k) then
                umask(i,k,j) = c1
              else
                umask(i,k,j) = c0
              endif
            enddo
          enddo
#ifdef bbl_ag
          do i=istrt-1,iend+1
            if (kmt(i,jrow) .gt. 0
# ifdef bbl_mask
     &           .and.mask_t_bbl(i,jrow).eq.1
# endif
     &           ) then
              tmask(i,km,j) = c1
            endif
            if (kmu(i,jrow) .gt. 0
# ifdef bbl_mask
     &           .and.mask_u_bbl(i,jrow).eq.1
# endif
     &           ) then
              umask(i,km,j) = c1
            endif
	  enddo
#endif
        enddo
c
#ifdef partial_cell
c
c-----------------------------------------------------------------------
c     construct thickness of T cells and U cells
c     (an aberration of the canonical indexing is used here)
c-----------------------------------------------------------------------
c
      do j=js,je+1
        jrow = j + joff
	minj = min(jrow,jmt)
        do k=1,km
          do i=1,imt
	    dht(i,k,j) = dzt(k)
	  enddo
	enddo
	do i=2,imt-1
	  kb = kmt(i,minj)
	  if (kb .gt. 1) then
	    dht(i,kb,j) = htp(i,minj) - zw(kb-1)
	  endif
	enddo
# ifdef cyclic
	call setbcx (dht(1,1,j), imt, km)
# else
c       do not set zero boundaries on dht
# endif
      enddo
c
      do j=js,je
        do k=0,km
          do i=1,imt
            dhwt(i,k,j) = dzw(k)
	  enddo
	enddo
        do k=1,km
	  do i=1,imt
            ztp(i,k,j) = zt(k)
	  enddo
	enddo
# ifdef isoneutralmix
        do k=1,km
	  do i=1,imt
            delqc(i,k,j,0) = c2*fracdz(k,0)*dht(i,k,j)
            delqc(i,k,j,1) = c2*fracdz(k,1)*dht(i,k,j)
	  enddo
	enddo
# endif
        jrow = j + joff
	do i=2,imt-1
	  kb = kmt(i,jrow)
	  if (kb .gt. 1) then
	    ztp(i,kb,j) = zw(kb-1) + fracdz(kb,0)*dht(i,kb,j)
            dhwt(i,kb-1,j) = ztp(i,kb,j) - zt(kb-1)
#ifdef bbl_ag
	    ztp_bbl  = htp(i,jrow) + p5*dzt(km)
            do k=kb+1,km
	      ztp(i,k,j)  = ztp_bbl
	    enddo
	    dhwt_bbl = ztp_bbl - ztp(i,kb,j)
            do k=kb,km
	      dhwt(i,k,j) = dhwt_bbl
	    enddo
# endif
	  endif
	enddo
        do k=1,km
          do i=2,imt-1
	    dhu(i,k,j) = min(dht(i,k,j), dht(i+1,k,j), dht(i,k,j+1)
     &,                      dht(i+1,k,j+1))
	  enddo
          dhu(1,k,j)   = dht(1,k,j)
          dhu(imt,k,j) = dht(imt,k,j)
	enddo
# ifdef cyclic
	call setbcx (ztp(1,1,j), imt, km)
	call setbcx (dhu(1,1,j), imt, km)
	call setbcx (dhwt(1,0,j), imt, km+1)
# else
c       do not set zero boundaries on dhwt or dhu
# endif
      enddo
#endif
c
#ifdef bbl_ag
      do j=js,je
        jrow = j + joff
        do i=1,imt
          kb = max(kmt(i,jrow),1)
# ifdef partial_cell
          zt_bbl(i,j)  = ztp(i,kb+1,j)
          dzw_bbl(i,j) = ztp(i,kb+1,j) - ztp(i,kb,j)
# else
          zt_bbl(i,j)  = zw(kb) + p5*dzt(km)
          dzw_bbl(i,j) = zw(kb) - zt(kb) + p5*dzt(km)
# endif
        enddo
      enddo
#endif
#if defined max_window
      endif
#endif
c
#if defined explicit_free_surface_nonlinear && defined partial_cell 
c
c-----------------------------------------------------------------------
c     construct time dependent thickness of surface T cells and U cells
c-----------------------------------------------------------------------
c
      call calc_thickness (joff, js, je)
c
#endif
c
c-----------------------------------------------------------------------
c     read data from "tau-1" and "tau" disk units into the MW
c     or if MW is wide open, copy data to proper time level and add
c     external mode to "tau" only since "tau-1" was done last timestep
c-----------------------------------------------------------------------
c
#if defined max_window
        if (leapfrog) then
          call add_ext_mode (joff, js, je, istrt, iend, 'tau')
	  if (first) then
	    call add_ext_mode (joff, js, je, istrt, iend, 'tau-1')
	  else
# if defined robert_time_filter
	    call add_ext_mode (joff, js, je, istrt, iend, 'tau-1')
# endif
	  endif
	elseif (euler2) then
          tau   = taup1
	  taup1 = taum1old
          call add_ext_mode (joff, js, je, istrt, iend, 'tau')
	elseif (euler1 .or. forward) then
          taum1old = taum1
	  taum1    = tau
          call add_ext_mode (joff, js, je, istrt, iend, 'tau')
	endif
#else
        if (leapfrog .or. euler2) then
c
c         read "tau","tau-1" disk data into "tau","tau-1" MW positions
c
          call getvar (joff, js, je, istrt-1, iend+1, num1, nu)
	elseif (forward .or. euler1) then
c
c         read "tau" disk data into "tau" and "tau-1" MW positions
c
          call getvar (joff, js, je, istrt-1, iend+1, nu, nu)
	endif
c
c       add external mode to both since only internal modes are on disk
c
        call add_ext_mode (joff, js, je, istrt, iend, 'tau-1')
        call add_ext_mode (joff, js, je, istrt, iend, 'tau')
#endif
#ifdef timing
      call toc ('loadmw', 'getvar+mask+ext mode')
#endif
c
#ifdef symmetry
c
c-----------------------------------------------------------------------
c     set symmetry condition at latitude of U cell row jmtm1 
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        if (jrow .eq. jmt) then
          do k=1,km
	    do i=istrt-1,iend+1
	      u(i,k,j,1,tau)     = u(i,k,j-2,1,tau)
	      u(i,k,j,1,taum1)   = u(i,k,j-2,1,taum1)
	      u(i,k,j-1,2,tau)   = c0
	      u(i,k,j-1,2,taum1) = c0
	      u(i,k,j,2,tau)     = -u(i,k,j-2,2,tau)
	      u(i,k,j,2,taum1)   = -u(i,k,j-2,2,taum1)
	    enddo
	  enddo
	  do n=1,nt
            do k=1,km
	      do i=istrt-1,iend+1
	        t(i,k,j,n,tau)   = t(i,k,j-1,n,tau)
	        t(i,k,j,n,taum1) = t(i,k,j-1,n,taum1)
	      enddo
	    enddo
	  enddo
        endif
      enddo
#endif
#if !defined second_order_window
c
c-----------------------------------------------------------------------
c     when jrow=1, initalize buffer rows below jrow=1
c-----------------------------------------------------------------------
c
      if (js+joff .eq. 1) then
        do jj=1,js-1
          do k=1,km
	    do i=istrt-1,iend+1
	      tmask(i,k,jj) = c0
	      umask(i,k,jj) = c0
	      do n=1,nt
	        t(i,k,jj,n,tau)   = c0
	        t(i,k,jj,n,taum1) = c0
	      enddo
	      do n=1,2
	        u(i,k,jj,n,tau)   = c0
	        u(i,k,jj,n,taum1) = c0
	      enddo
# if defined partial_cell
              dht(i,k,jj) = dzt(k)
              dhu(i,k,jj) = dzt(k)
# endif
	    enddo
	  enddo
# if defined partial_cell
          do k=1,km
            do i=istrt-1,iend+1
	      ztp(i,k,jj) = zt(k)
	    enddo
	  enddo
          do k=0,km
	    do i=istrt-1,iend+1
              dhwt(i,k,jj) = dzw(k)
	    enddo
	  enddo
# endif
	enddo
      endif
      if (je+joff .eq. jmt) then
c
c-----------------------------------------------------------------------
c       when jrow=jmt, initialize buffer rows above jrow=jmt
c-----------------------------------------------------------------------
c
        do jj=je+1,jmw
          do k=1,km
	    do i=istrt-1,iend+1
	      tmask(i,k,jj) = c0
	      umask(i,k,jj) = c0
	      do n=1,nt
	        t(i,k,jj,n,tau)   = c0
	        t(i,k,jj,n,taum1) = c0
	      enddo
	      do n=1,2
	        u(i,k,jj,n,tau)   = c0
	        u(i,k,jj,n,taum1) = c0
	      enddo
# if defined partial_cell
              dht(i,k,jj) = dzt(k)
              dhu(i,k,jj) = dzt(k)
# endif
	    enddo
	  enddo
# if defined partial_cell
          do k=1,km
            do i=istrt-1,iend+1
	      ztp(i,k,jj) = zt(k)
	    enddo
	  enddo
          do k=0,km
	    do i=istrt-1,iend+1
              dhwt(i,k,jj) = dzw(k)
	    enddo
	  enddo
# endif
	enddo
      endif
#endif
#if defined pressure_depth
c
c     construct full density at "tau-1" for use in calculating pressure
c     at depth of T cell grid points. This pressure will be used in the  
c     UNESCO equation of state to construct the density anomoly.     
c     
      do j=js,je
	do k=1,km
	  do i=1,imt
	    rho_full(i,k,j) = rho(i,k,j) + rho0
	  enddo
	enddo
      enddo
#endif
#ifdef timing
      call tic ('loadmw', 'state')
#endif
c
c-----------------------------------------------------------------------
c     compute density anomoly (rho=density-rho0) at T cell grid points
c-----------------------------------------------------------------------
c
      call state (t(1,1,1,1,tau), t(1,1,1,2,tau), rho, rho_full
     &,           js, je, istrt-1, iend+1, joff)
c
c-----------------------------------------------------------------------
c     pressure = hydrostatic pressure anomoly = (rho - rho0)*grav*z 
c     for use in horizontal pressure gradient calculations
c-----------------------------------------------------------------------
c
      call hydrostatic_pressure (pressure, rho, js, je, joff)
c
#if defined explicit_free_surface
      do j=js,je
        jrow=j+joff
        do i=1,imt
          rhosurf(i,jrow) = (rho(i,1,j) + rho0)*tmask(i,1,j)      
        enddo
      enddo
#endif    
#ifdef timing
      call toc ('loadmw', 'state')
#endif
c
c-----------------------------------------------------------------------
c       open boundaries: extrapolate tracer and velocity values onto  
c       boundary to prevent diffusion
c-----------------------------------------------------------------------
c
#if defined obc
# include "loadmw_obc.inc"
#endif
c
#ifdef timing
      call toc ('loadmw', 'loading the MW')
#endif
c
#ifdef debug_loadmw
      do j=js,je
        jrow = j + joff
        write(stdout,*) 'In loadmw:  jrow=',jrow
        call print_checksum (u(1,1,j,1,tau), imt, km, 'ck u(tau) = ') 
        call print_checksum (u(1,1,j,2,tau), imt, km, 'ck v(tau) = ') 
        call print_checksum (t(1,1,j,1,tau), imt, km, 'ck T(tau) = ') 
        call print_checksum (t(1,1,j,2,tau), imt, km, 'ck S(tau) = ') 
c
        call print_checksum (u(1,1,j,1,taum1), imt, km,'ck u(tau-1) = ')
        call print_checksum (u(1,1,j,2,taum1), imt, km,'ck v(tau-1) = ')
        call print_checksum (t(1,1,j,1,taum1), imt, km,'ck T(tau-1) = ')
        call print_checksum (t(1,1,j,2,taum1), imt, km,'ck S(tau-1) = ')
      enddo        
#endif
      return
      end      




      subroutine move_memory_window (iss, iee)
c
c=======================================================================
c     move the MW up (northward) by copying data from the last two rows
c     into the first two rows.
c     (last 3 rows if using fourth order schemes)
c
c     input:
c      is = starting longitude index in the MW
c      ie = ending longitude index in the MW
c      note: iss,iee bypassed to optimize performance
c            is,ie set by parameter statement below
c     
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
      use memory_window
#include "param.h"
#ifdef bbl_ag
# include "bbl.h"
#endif
#include "hmixc.h"
#ifdef isoneutralmix
# include "isopyc.h"
#endif
#include "mw.h"
#include "rinum.h"
#include "vmixc.h"
      parameter (is=1, ie=imt)
c
#ifdef timing
      call tic ('loadmw', 'moving MW northward')
#endif
c
      do j=1,num_rows_to_copy
        jfrom = jmw - num_rows_to_copy + j
	jto   = j
c
c-----------------------------------------------------------------------
c       copy quantities with rows dimensioned (1:jmw)
c-----------------------------------------------------------------------
c
	do i=is,ie
	  csudxur(i,jto)    = csudxur(i,jfrom)
	  csudxu2r(i,jto)   = csudxu2r(i,jfrom)
	  csudxtr(i,jto)    = csudxtr(i,jfrom)
	  cstdxtr(i,jto)    = cstdxtr(i,jfrom)
	  cstdxt2r(i,jto)   = cstdxt2r(i,jfrom)
	  cstdxur(i,jto)    = cstdxur(i,jfrom)
#if defined velocity_horz_mix_const
	  am_csudxtr(i,jto) = am_csudxtr(i,jfrom)
#endif
#if defined tracer_horz_mix_const
	  ah_cstdxur(i,jto) = ah_cstdxur(i,jfrom)
#endif
        enddo
c
#ifdef bbl_ag
        do i=is,ie
          px_bbl(i,jto)   = px_bbl(i,jfrom)
          oppsignx(i,jto) = oppsignx(i,jfrom)
        enddo
#endif
        do k=1,km
	  do i=is,ie
	    do n=1,nt
	      t(i,k,jto,n,taum1) = t(i,k,jfrom,n,taum1)
	      t(i,k,jto,n,tau)   = t(i,k,jfrom,n,tau)
	    enddo
	    do n=1,2
	      u(i,k,jto,n,taum1) = u(i,k,jfrom,n,taum1)
	      u(i,k,jto,n,tau)   = u(i,k,jfrom,n,tau)
	    enddo
	    tmask(i,k,jto)    = tmask(i,k,jfrom)
	    umask(i,k,jto)    = umask(i,k,jfrom)
	    adv_vnt(i,k,jto)  = adv_vnt(i,k,jfrom)
	    rho(i,k,jto)      = rho(i,k,jfrom)
	    pressure(i,k,jto) = pressure(i,k,jfrom)
          enddo
        enddo
c
        do i=is,ie
	  do n=1,nt
	    stf(i,jto,n) = stf(i,jfrom,n)
	    btf(i,jto,n) = btf(i,jfrom,n)
	  enddo
	  smf(i,jto,1) = smf(i,jfrom,1)
	  smf(i,jto,2) = smf(i,jfrom,2)
	  bmf(i,jto,1) = bmf(i,jfrom,1)
	  bmf(i,jto,2) = bmf(i,jfrom,2)
#ifdef bbl_ag
          pbot(i,jto)    = pbot(i,jfrom)
          rhoz(i,jto)    = rhoz(i,jfrom)
          uz2(i,jto)     = uz2(i,jfrom)
          zt_bbl(i,jto)  = zt_bbl(i,jfrom)
          dzw_bbl(i,jto) = dzw_bbl(i,jfrom)
#endif
	enddo
#ifdef partial_cell
        do k=1,km
	  do i=is,ie
	    dht(i,k,jto)  = dht(i,k,jfrom)
	    dhu(i,k,jto)  = dhu(i,k,jfrom)
	    ztp(i,k,jto)  = ztp(i,k,jfrom)
# ifdef isoneutralmix
	    delqc(i,k,jto,0)  = delqc(i,k,jfrom,0)
	    delqc(i,k,jto,1)  = delqc(i,k,jfrom,1)
# endif
          enddo
        enddo
        do k=0,km
	  do i=is,ie
	    dhwt(i,k,jto)  = dhwt(i,k,jfrom)
          enddo
        enddo
#endif
c
c-----------------------------------------------------------------------
c       copy quantities with rows dimensioned (1:jmw-1)
c-----------------------------------------------------------------------
c
        if (jfrom .le. jmw-1) then
          do k=1,km
	    do i=is,ie
	      adv_vnu(i,k,jto) = adv_vnu(i,k,jfrom)
            enddo
          enddo
          do n=1,2
            do k=1,km
	      do i=is,ie
  	        grad_p(i,k,jto,n) = grad_p(i,k,jfrom,n)
              enddo
            enddo
          enddo
#ifdef tracer_advection_fct
          do n=1,nt
            do k=1,km
	      do i=is,ie
                anti_fn(i,k,jto,n)  = anti_fn(i,k,jfrom,n)
                R_plusY(i,k,jto,n)  = R_plusY(i,k,jfrom,n)
                R_minusY(i,k,jto,n) = R_minusY(i,k,jfrom,n)
# ifdef fct_3d
                R_plus3(i,k,jto,n)  = R_plus3(i,k,jfrom,n)
                R_minus3(i,k,jto,n) = R_minus3(i,k,jfrom,n)
# endif
              enddo
            enddo
          enddo
#endif
#ifdef trajectories
          do k=0,km
	    do i=is,ie
	      adv_vbu(i,k,jto) = adv_vbu(i,k,jfrom)
            enddo
          enddo
#endif
#if defined velocity_horz_mix_var
          do k=1,km
	    do i=is,ie
	      am_varh(i,k,jto)  = am_varh(i,k,jfrom)
	      visc_cnu(i,k,jto) = visc_cnu(i,k,jfrom)
            enddo
          enddo
#endif
#if defined tracer_horz_mix_var
          do k=1,km
	    do i=is,ie
	      diff_cnt(i,k,jto) = diff_cnt(i,k,jfrom)
            enddo
          enddo
#endif
#ifdef bbl_ag
          do i=is,ie
            py_bbl(i,jto)      = py_bbl(i,jfrom)
            oppsigny(i,jto)    = oppsigny(i,jfrom)
          enddo
# if defined detrain
	  do n=1,nt
            do i=is,ie
	      adv_fdetn(i,jto,n) = adv_fdetn(i,jfrom,n)
	    enddo
          enddo
# endif
#endif
#if defined fourth_order_tracer_advection || defined quicker
          do k=1,km
	    do i=is,ie
	      do n=1,nt
	        adv_f4n(i,k,jto,n) = adv_f4n(i,k,jfrom,n)
	      enddo
            enddo
          enddo
#endif
c
        endif
c
c-----------------------------------------------------------------------
c       copy quantities with rows dimensioned (2:jmw)
c-----------------------------------------------------------------------
c
        if (jto .ge. 2) then
          do k=1,km
	    do i=is,ie
	      adv_vet(i,k,jto) = adv_vet(i,k,jfrom)
            enddo
          enddo
          do k=0,km
	    do i=is,ie
	      adv_vbt(i,k,jto) = adv_vbt(i,k,jfrom)
            enddo
          enddo
#ifdef bbl_ag
          do i=is,ie
            Kbbl(i,jto) = Kbbl(i,jfrom)
          enddo
# if defined detrain
	  do n=1,nt
            do i=is,ie
	      adv_fdets(i,jto,n) = adv_fdets(i,jfrom,n)
	      adv_fdete(i,jto,n) = adv_fdete(i,jfrom,n)
	      adv_fdetw(i,jto,n) = adv_fdetw(i,jfrom,n)
	    enddo
          enddo
# endif
#endif
	endif
c
c
c-----------------------------------------------------------------------
c       copy quantities with rows dimensioned (2:jmw-1)
c-----------------------------------------------------------------------
c
        if (jto .ge. 2 .and. jfrom .le. jmw-1) then
          do k=1,km
	    do i=is,ie
	      visc_cbu(i,k,jto) = visc_cbu(i,k,jfrom)
	      diff_cbt(i,k,jto,1) = diff_cbt(i,k,jfrom,1)
	      diff_cbt(i,k,jto,nsmix) = diff_cbt(i,k,jfrom,nsmix)
	      adv_veu(i,k,jto)  = adv_veu(i,k,jfrom)
# if defined velocity_horz_mix_var
              varh_cu(i,k,jto)  = varh_cu(i,k,jfrom)
              visc_ceu(i,k,jto) = visc_ceu(i,k,jfrom)
# endif
# if defined tracer_horz_mix_var
              diff_cet(i,k,jto) = diff_cet(i,k,jfrom)
# endif
            enddo
          enddo
c
# ifdef tracer_advection_fct
c         row "jmw" is not used when "jmw" < "jmt" so it is not copied
c         treat these arrays as if dimensioned (2:jmw-1)
          do n=1,nt
            do k=1,km
	      do i=is,ie
                anti_fe(i,k,jto,n)  = anti_fe(i,k,jfrom,n)
#  ifdef fct_dlm2
                t_lo(i,k,jto,n) = t_lo(i,k,jfrom,n)
#  endif
              enddo
            enddo
            do k=0,km
	      do i=is,ie
                anti_fb(i,k,jto,n) = anti_fb(i,k,jfrom,n)
              enddo
            enddo
          enddo
# endif
c
# if !defined trajectories
          do k=0,km
	    do i=is,ie
	      adv_vbu(i,k,jto) = adv_vbu(i,k,jfrom)
            enddo
          enddo
# endif
c
# if defined velocity_horz_biharmonic
          do n=1,2
            do k=1,km
	      do i=is,ie
	        del2_vel(i,k,jto,n) = del2_vel(i,k,jfrom,n)
              enddo
            enddo
          enddo
# endif
# if defined tracer_horz_biharmonic
          do n=1,nt
            do k=1,km
	      do i=is,ie
	        del2_tracer(i,k,jto,n) = del2_tracer(i,k,jfrom,n)
              enddo
            enddo
          enddo
# endif
c
# ifdef bbl_ag
          do i=is,ie
            entrain(i,jto)  = entrain(i,jfrom)
          enddo
# endif
	endif
c
#ifdef isoneutralmix
c
c-----------------------------------------------------------------------
c       copy quantities with rows dimensioned (1:jmw)
c-----------------------------------------------------------------------
c
        do k=1,km
	  do i=is,ie
            drho(i,k,jto,1)  = drho(i,k,jfrom,1)
            drho(i,k,jto,2)  = drho(i,k,jfrom,2)
# if defined full_tensor
	    ddxt(i,k,jto,1) = ddxt(i,k,jfrom,1)
	    ddxt(i,k,jto,2) = ddxt(i,k,jfrom,2)
# endif
          enddo
        enddo
        do k=0,km
	  do i=is,ie
	    ddzt(i,k,jto,1) = ddzt(i,k,jfrom,1)
	    ddzt(i,k,jto,2) = ddzt(i,k,jfrom,2)
          enddo
        enddo
c
	  do i=is,ie
# if defined hl_diffusivity
	    hllength(i,jto) = hllength(i,jfrom)
	    hldiff(i,jto)   = hldiff(i,jfrom)
# endif
# if defined vmhs_diffusivity
	    vmhslength(i,jto) = vmhslength(i,jfrom)
	    vmhsdiff(i,jto)   = vmhsdiff(i,jfrom)
# endif
          enddo
c
c-----------------------------------------------------------------------
c       copy quantities with rows dimensioned (1:jmw-1)
c-----------------------------------------------------------------------
c
        if (jfrom .le. jmw-1) then
          do k=1,km
	    do i=is,ie
# if defined redi_diffusion
              K22(i,k,jto)    = K22(i,k,jfrom)
# endif
	      ddyt(i,k,jto,1) = ddyt(i,k,jfrom,1)
	      ddyt(i,k,jto,2) = ddyt(i,k,jfrom,2)
              do kr=0,1
                do jq=0,1
# if defined redi_diffusion || defined gm_skew
                  Ai_nz(i,k,jto,jq,kr) = Ai_nz(i,k,jfrom,jq,kr)
# endif
# if defined biharmonic_rm
                  Bi_nz(i,k,jto,jq,kr) = Bi_nz(i,k,jfrom,jq,kr)
# endif
	        enddo
	      enddo
# ifdef gent_mcwilliams
	      adv_vntiso(i,k,jto)  = adv_vntiso(i,k,jfrom)
# endif
            enddo
          enddo
	endif
c
c-----------------------------------------------------------------------
c       copy quantities with rows dimensioned (2:jmw-1)
c-----------------------------------------------------------------------
c
        if (jto .ge. 2 .and. jfrom .le. jmw-1) then
          do k=1,km
	    do i=is,ie
              do kr=0,1
                do ip=0,1
#  if defined redi_diffusion || defined gm_skew
                  Ai_ez(i,k,jto,ip,kr)  = Ai_ez(i,k,jfrom,ip,kr)
#  endif
#  if defined biharmonic_rm
                  Bi_ez(i,k,jto,ip,kr)  = Bi_ez(i,k,jfrom,ip,kr)
#  endif
	        enddo
	      enddo
              do kr=0,1
                do ip=0,1
#  if defined redi_diffusion || defined gm_skew
                  Ai_bx(i,k,jto,ip,kr) = Ai_bx(i,k,jfrom,ip,kr)
#  endif
#  if defined biharmonic_rm
                  Bi_bx(i,k,jto,ip,kr) = Bi_bx(i,k,jfrom,ip,kr)
#  endif
	        enddo
	      enddo
              do kr=0,1
                do jq=0,1
#  if defined redi_diffusion || defined gm_skew
                  Ai_by(i,k,jto,jq,kr) = Ai_by(i,k,jfrom,jq,kr)
#  endif
#  if defined biharmonic_rm
                  Bi_by(i,k,jto,jq,kr) = Bi_by(i,k,jfrom,jq,kr)
#  endif
	        enddo
	      enddo
#  if defined redi_diffusion
	      K11(i,k,jto)   = K11(i,k,jfrom)
	      K33(i,k,jto)   = K33(i,k,jfrom)
#   if !defined full_tensor
	      ddxt(i,k,jto,1) = ddxt(i,k,jfrom,1)
	      ddxt(i,k,jto,2) = ddxt(i,k,jfrom,2)
#   endif
#  endif
#  ifdef gent_mcwilliams
	      adv_vetiso(i,k,jto) = adv_vetiso(i,k,jfrom)
#  endif
            enddo
          enddo
#  ifdef gent_mcwilliams
          do k=0,km
	    do i=is,ie
	      adv_vbtiso(i,k,jto) = adv_vbtiso(i,k,jfrom)
	      adv_fbiso(i,k,jto)  = adv_fbiso(i,k,jfrom)
            enddo
          enddo
#  endif
        endif
# endif
c
#ifdef ppvmix
c
c-----------------------------------------------------------------------
c       copy quantities with rows dimensioned (1:jmw)
c-----------------------------------------------------------------------
c
        do k=1,km
	  do i=is,ie
	    rhom1z(i,k,jto)  = rhom1z(i,k,jfrom)
          enddo
        enddo
c
c-----------------------------------------------------------------------
c       copy quantities with rows dimensioned (1:jmw-1)
c-----------------------------------------------------------------------
c
        if (jfrom .le. jmw-1) then
          do k=1,km
	    do i=is,ie
	      riu(i,k,jto)  = riu(i,k,jfrom)
            enddo
          enddo
        endif
c
c-----------------------------------------------------------------------
c       copy quantities with rows dimensioned (2:jmw-1)
c-----------------------------------------------------------------------
c
        if (jto .ge. 2 .and. jfrom .le. jmw-1) then
          do k=1,km
	    do i=is,ie
	      rit(i,k,jto)  = rit(i,k,jfrom)
            enddo
          enddo
        endif
c
#endif
#ifdef kppvmix
c
c-----------------------------------------------------------------------
c       copy quantities with rows dimensioned (1:jmw)
c-----------------------------------------------------------------------
c
        do k=1,km
          do i=is,ie
            rhom1z(i,k,jto) = rhom1z(i,k,jfrom)
            dVsq(i,k,jto)   = dVsq(i,k,jfrom)
          enddo
        enddo
c
c-----------------------------------------------------------------------
c       copy quantities with rows dimensioned (1:jmw-1)
c-----------------------------------------------------------------------
c
        if (jfrom .le. jmw-1) then
          do k=1,km
	    do i=is,ie
	      riu(i,k,jto)  = riu(i,k,jfrom)
            enddo
          enddo
        endif
c
c-----------------------------------------------------------------------
c       copy quantities with rows dimensioned (2:jmw-1)
c-----------------------------------------------------------------------
c
        if (jto .ge. 2 .and. jfrom .le. jmw-1) then
          do k=1,km
            do i=is,ie
	      rit(i,k,jto)     = rit(i,k,jfrom)
              dbloc(i,k,jto)   = dbloc(i,k,jfrom)
              dbsfc(i,k,jto)   = dbsfc(i,k,jfrom)
              talpha(i,k,jto)  = talpha(i,k,jfrom)
              sbeta(i,k,jto)   = sbeta(i,k,jfrom)
              alphaDT(i,k,jto) = alphaDT(i,k,jfrom)
              betaDS(i,k,jto)  = betaDS(i,k,jfrom)
              blmc(i,k,jto,1)  = blmc(i,k,jfrom,1)
              blmc(i,k,jto,2)  = blmc(i,k,jfrom,2)
              blmc(i,k,jto,3)  = blmc(i,k,jfrom,3)
              ghats(i,k,jto)   = ghats(i,k,jfrom)
            enddo
          enddo
	  do i=is,ie
            Rib(i,1,jto)   = Rib(i,1,jfrom)
            Rib(i,2,jto)   = Rib(i,2,jfrom)
	    Bo(i,jto)      = Bo(i,jfrom)
	    Bosol(i,jto)   = Bosol(i,jfrom)
	    jwtype(i,jto)  = jwtype(i,jfrom)
	    bfsfc(i,jto)   = bfsfc(i,jfrom)
	    ustar(i,jto)   = ustar(i,jfrom)
	    caseA(i,jto)   = caseA(i,jfrom)
	    sigma(i,jto)   = sigma(i,jfrom)
	    rhosfc(i,jto)  = rhosfc(i,jfrom)
	    kbl(i,jto)     = kbl(i,jfrom)
	    hbl(i,jto)     = hbl(i,jfrom)
	  enddo
	  do n=1,3
	    do i=is,ie
              dkm1(i,jto,n)  = dkm1(i,jfrom,n)
              gat1(i,jto,n)  = gat1(i,jfrom,n)
              dat1(i,jto,n)  = dat1(i,jfrom,n)
	    enddo
	  enddo
        endif
#endif
c
#ifdef trace_indices
        write (stdout,'(4x,2(a,i4))')
     & "=> In movemw: copying variables on row ",jfrom," to row ",jto
#endif
c
      enddo
#ifdef timing
      call toc ('loadmw', 'moving MW northward')
#endif
      return
      end



      subroutine getvar (joff, js, je, is, ie, num1, nu)
c
c=======================================================================
c     read prognostic quantities from disk units "num1" (tau-1) and
c     "nu" (tau) into the MW for rows "js" through "je"
c
c     input:
c       joff = offset between "j" in the MW and latitude "jrow"
c       js = starting row in the MW
c       je = ending row in the MW
c       is = starting longitude index in the MW
c       ie = ending longitude index in the MW
c       num1 = "tau-1" latitude disk unit
c       nu   = "tau" latitude disk unit
c     
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
      use memory_window
#include "param.h"
#include "mw.h"
c
#ifdef timing
      call tic ('getvar', 'from disk')
#endif
c
      do j=js,je
        jrow   = j + joff
c
c       read "tau-1" data into the "tau-1" portion of the MW
c
        call getrow (num1, nslab, jrow, u(1,1,j,1,taum1)
     &,                                 t(1,1,j,1,taum1))
        if (nu .ne. num1) then
c
c         read "tau" data into the "tau" portion of the MW 
c
          call getrow (nu, nslab, jrow, u(1,1,j,1,tau)
     &,                                 t(1,1,j,1,tau))
        else
c
c         copy "tau-1" data into "tau" portion of the MW
c
	  do k=1,km
	    do i=is,ie
	      u(i,k,j,1,tau) = u(i,k,j,1,taum1)
	      u(i,k,j,2,tau) = u(i,k,j,2,taum1)
	    enddo
	  enddo
	  do n=1,nt
	    do k=1,km
	      do i=is,ie
	        t(i,k,j,n,tau) = t(i,k,j,n,taum1)
	      enddo
	    enddo
	  enddo
        endif
      enddo
c
#ifdef timing
      call toc ('getvar', 'from disk')
#endif
c
      return
      end


      subroutine add_ext_mode (joff, js, je, is, ie, timelev)
c
c=======================================================================
c     add external mode to velocity for time level "timelev"
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       timelev = "tau" or "tau-1"
c     
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
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
      character*(*) timelev
c
      if (timelev .eq. 'tau') then
c
c-----------------------------------------------------------------------
c       add external mode to "tau" velocity
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
                  u(i,k,j,n,tau) = (u(i,k,j,n,tau) + ext(i,j,n))
     &                             *umask(i,k,j)
                enddo
              enddo
            enddo
#endif
#if defined rigid_lid_surface_pressure || defined free_surf_on
            do n=1,2
	      do k=1,km
                do i=is,ie
                  u(i,k,j,n,tau) = (u(i,k,j,n,tau) + ubar(i,jrow,n))
     &                             *umask(i,k,j)
                enddo
              enddo
            enddo
# if defined explicit_free_surface_nonlinear && defined partial_cell
c
c           diagnose ubar from full velocity
c
            do n=1,2
              do i=is,ie
                ubar(i,jrow,n) = c0
	        do k=1,km
                  ubar(i,jrow,n) = ubar(i,jrow,n)
     &                           + u(i,k,j,n,tau)*dhu(i,k,j)
                enddo
                ubar(i,jrow,n)   = ubar(i,jrow,n)*depthur(i,jrow)
              enddo
	      call setbcx (ubar(1,jrow,n), imt, 1)
            enddo
# endif
#endif
            do n=1,2
	      call setbcx (u(1,1,j,n,tau), imt, km)
	    enddo
	  endif
	enddo
c
      elseif (timelev .eq. 'tau-1') then
c
c-----------------------------------------------------------------------
c       add external mode to "tau-1" velocity
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
                  u(i,k,j,n,taum1) = (u(i,k,j,n,taum1) + ext(i,j,n))
     &                               *umask(i,k,j)
                enddo
              enddo
            enddo
#endif
#if defined rigid_lid_surface_pressure || defined free_surf_on
            do n=1,2
	      do k=1,km
                do i=is,ie
                  u(i,k,j,n,taum1) = (u(i,k,j,n,taum1)+ubarm1(i,jrow,n))
     &                           *umask(i,k,j)
                enddo
              enddo
            enddo
#endif
            do n=1,2
	      call setbcx (u(1,1,j,n,taum1), imt, km)
	    enddo
	  endif
	enddo
      else
        write (stdout,'(/a,a/)') 'Error: timelev = ',timelev
	call abort()
      endif
c
      return
      end
      


      subroutine putmw (joff, js, je, nup1)
c
c=======================================================================
c     write prognostic quantities from MW to disk unit "nup1" (tau+1)
c     for rows "js" to "je".
c     
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
      use memory_window
#include "param.h"
#include "mw.h"
c
#ifdef timing
      call tic ('putmw', 'to disk')
#endif
c
c     write all newly computed quantities to disk "tau+1"
c
      do j=js,je
        jrow   = j + joff
        call putrow (nup1, nslab, jrow, u(1,1,j,1,taup1)
     &,                                 t(1,1,j,1,taup1))
      enddo
c
#ifdef timing
      call toc ('putmw', 'to disk')
#endif
c
      return
      end


      
