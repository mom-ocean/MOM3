      subroutine expl_freesurf (js, je)

      use mom3_mpp_mod
      use arrays_2d
# include "param.h"
# include "coord.h"
# if defined explicit_polar_filter && defined firfil
#  include "cpolar.h"
# endif
# include "cprnts.h"
# include "emode.h"
# include "grdvar.h"
# include "hmixc.h"
# include "index.h"
# include "scalar.h"
# include "switch.h"
# include "taskrows.h"
# include "tmngr.h"
# ifdef obc
#  include "cobc.h"
#  ifndef orlanski
#   include "obc_data.h"
#  endif
# endif
# if defined xlandmix
#  include "xland.h"
# endif
#if defined explicit_free_surface && defined time_averaged_surface
c
c=======================================================================
c
c     Compute the barotropic fields:
c     vertically averaged velocity ubar and free surface height etat  
c
c     Two main options
c      
c     A. -Dexplicit_free_surface_linear
c     Keeps the ocean volume constant.  Does not conserve tracers. 
c     Similar to the equations solved with Killworth etal
c     1991 and Dukowicz and Smith 1993.
c
c     B. -Dexplicit_free_surface_nonlinear
c     Top ocean cell undulates in time using partial cell methods. 
c     Total tracers constant if no boundary fluxes and dtts=dtuv.
c     Requires -Dpartial_cell
c
c     Basics of the algorithm:
c
c     1. The baroclinic equations have been previously solved for
c     time taup1 = T0+dtuv = itt*dtuv. The vertically integrated
c     forcing terms zu, centered at baroclinic time tau=T0, are known.  
c
c     2. Using zu as forcing terms, held fixed at time tau=T0,
c     the barotropic equations are integrated for N=2*dtuv/dtsf 
c     barotropic timesteps t, each of length dtsf. 
c     The barotropic integration goes from baroclinic 
c     time T0 to T0+2*dtuv (tau to tau+2).
c
c     3. To update the ocean surface height onto the baroclinic
c     time step taup1=T0+dtuv, a leap-frog is used.  However, the 
c     eta(taum1) value is taken from the time averaged eta(taum1),
c     instead of using a Robert time filter.  This approach
c     leads to a smoother and more stable solution.  If the
c     nonlinear free surface is used, then eta(taup1) is used 
c     to update the top cell thicknesses dht(i,1,j) & dhu(i,1,j).
c
c     4. For the nonlinear free surface, the updated full velocity is
c     computed in loadmw.  The corrected updated barotropic velocity
c     is diagnosed from the updated full velocity and eta(taup1).
c     This corrected ubar is used to initialize the next barotropic
c     time step.  
c
c     For spin-ups to equlibrium in coarse models, long tracer 
c     time steps are available, so far as the scheme remains
c     numerically stable.  However, conservation is manifest
c     only with the nonlinear free surface and with dtts=dtuv.  
c
c     authors:
c             R. C. Pacanowski e-mail=> rcp@gfdl.gov    
c             S. M. Griffies   e-mail=> smg@gfdl.gov 
c     (mpp items)  V. Balaji   e-mail=> vb@gfdl.gov
c                  C. Kerr     e-mail=> ck@gfdl.gov
c      
c=======================================================================
c     Some fields local to this subroutine
c
c     eta    = sea surface elevation in cm at barotropic time t
c     etap   = sea surface elevation in cm at barotropic time t+1
c     hubar  = depth over u-cells * barotropic velocity
c-----------------------------------------------------------------------
c
# ifdef obc
      dimension etarel(imt,js:je)
# endif
      dimension eta(imt,js:je), etap(imt,js:je)
      dimension hubar(imt,js:je,2)
      dimension utempm(2)
# if defined explicit_eta_laplacian
      dimension tempeta(imt,js:je)
      dimension weighte(js:je),weightn(js:je)
# endif
# if defined explicit_polar_filter
      dimension hu0(imt,js:je,2)
      dimension dhudt(imt,js:je,2), dedt(imt,js:je)
# endif
# if defined explicit_polar_filter && defined firfil
      dimension tempik(imt,2)
# endif
# ifdef timing
      call tic ('barotropic', 'expl_freesurf')
# endif
c
c-----------------------------------------------------------------------
c     some scalars 
c-----------------------------------------------------------------------
c
      mscan      = 0
      gravrho0   = grav*rho0  
      gravrho0r  = c1/gravrho0
c
c-----------------------------------------------------------------------
c     initialization for barotropic equations
c-----------------------------------------------------------------------
c
      do jrow=js,je
	do i=1,imt
c
c         shift (tau) into (tau-1) for next baroclinic time step
c
          ubarm1(i,jrow,1) = ubar(i,jrow,1)
          ubarm1(i,jrow,2) = ubar(i,jrow,2)
          ps(i,jrow,2)     = ps(i,jrow,1) 
c
c         save leap-frog etat(tau) in thicktau
c         save time mean etat(tau) temporarily in etat(..,1)
c
          thicktau(i,jrow)   = max(epsln,dzt(1)+etat(i,jrow,1))
          etat(i,jrow,1)     = etat(i,jrow,3)
c
c         initialize fields used in the barotropic time stepping 
c
          eta(i,jrow)      = etat(i,jrow,3)
	  etap(i,jrow)     = c0
          hubar(i,jrow,1)  = depthu(i,jrow)*ubar(i,jrow,1)
          hubar(i,jrow,2)  = depthu(i,jrow)*ubar(i,jrow,2)
# if defined explicit_polar_filter
          hu0(i,jrow,1)    = hubar(i,jrow,1)
          hu0(i,jrow,2)    = hubar(i,jrow,2)
# endif
c
c         initialize ubar for temporary accumulation of hubar
c
          ubar(i,jrow,1) = hubar(i,jrow,1) 
          ubar(i,jrow,2) = hubar(i,jrow,2) 
c
# ifdef obc
          etarel(i,jrow) = c0
#  ifdef obc_west
          umsk(1,jrow)   = c0
          umsk(2,jrow)   = c0
#  endif
#  ifdef obc_east
          umsk(imt,jrow)   = c0
          umsk(imtm1,jrow) = c0
#  endif
# endif
        enddo
      enddo
c
# ifdef symmetry
c
c     set meridional symmetry on vertically integrated forcing
c
      if (je .eq. jmt) then   
        do i=1,imt
	  zu(i,je-1,2) = c0
        enddo
      endif
# endif
# ifdef obc_south
      if (js .eq. 1) then   
        do i=1,imt
	  umsk(i,js)   = c0
  	  umsk(i,js+1) = c0
        enddo
      endif
# endif
# ifdef obc_north
#  ifdef obc_north_anywhere
      if (je .eq. jobcn+1) then   
#  else
      if (je .eq. jmt) then   
#  endif
        do i=1,imt
	  umsk(i,je)     = c0
 	  umsk(i,je-1)   = c0
        enddo
      endif
# endif
c
c=======================================================================
c     integrate barotropic equations for "nb" barotropic time steps
c     using an Euler forward-backward scheme from baroclinic time T0
c     to baroclinic time T0+2*dtuv.  
c=======================================================================
c
      nb = 2*nint(dtuv/dtsf)
      do itbt=1,nb
c
c-----------------------------------------------------------------------
c       time step free surface height using a forward barotropic step
c-----------------------------------------------------------------------
c
        do jrow=jscomp,jecomp
	  do i=2,imtm1
	    conv_hu       = -dxt2r(i)*dytr(jrow)*cstr(jrow)
     &                    *( hubar(i  ,jrow,1  )*dyu(jrow)
     &                      +hubar(i,jrow-1,1  )*dyu(jrow-1)
     &                    -( hubar(i-1,jrow,1  )*dyu(jrow)
     &                      +hubar(i-1,jrow-1,1)*dyu(jrow-1))
     &                    +( hubar(i,jrow,2)*dxu(i) + 
     &                       hubar(i-1,jrow,2)*dxu(i-1))*csu(jrow)
     &                    -( hubar(i,jrow-1,2)*dxu(i) + 
     &                       hubar(i-1,jrow-1,2)*dxu(i-1))*csu(jrow-1))
	    etap(i,jrow) = eta(i,jrow) + dtsf*(
     &                                         conv_hu
# ifdef explicit_fresh_water_flux
     &                                         +sfft(i,jrow)
# endif
# ifdef xlandmix_eta
     &                                         +source_eta(i,jrow)
# endif
     &                                         )
	  enddo
# ifdef cyclic
	  etap(  1,jrow) = etap(imtm1,jrow)
          etap(imt,jrow) = etap(2,jrow)
# else
	  etap(  1,jrow) = c0
          etap(imt,jrow) = c0
# endif
        enddo
c
# ifdef obc_south
      jrow = 2
      var = -dyu(jrow+1)/dtsf
      do i=2,imtm1
c
         var1(i,1,1)= eta(i,jrow+2) - eta(i,jrow+1)
         cmax = -c2*sqrt(grav*zw(kmt(i,jrow)))
         if (var1(i,1,1) .eq. c0) then
           c1ps(i) = cmax
         else
           c1ps(i) = var*(etap(i,jrow+1)-eta(i,jrow+1))/var1(i,1,1)
           if (c1ps(i) .ge. c0)   c1ps(i) = c0
           if (c1ps(i) .lt. cmax) c1ps(i) = cmax
         endif
         etap(i,jrow)  = eta(i,jrow) 
     &           - c1ps(i)*dtsf*dyur(jrow)*(eta(i,jrow+1)-eta(i,jrow))
      enddo
# endif
# ifdef obc_north
#  ifdef obc_north_anywhere
      jrow = jobcn
#  else
      jrow = jmtm1
#  endif
      var = dyu(jrow-1)/dtsf
      do i=2,imtm1
         var1(i,1,1)= eta(i,jrow-1) - eta(i,jrow-2)
         cmax = c2*sqrt(grav*ABS(zw(kmt(i,jrow))))
         if (var1(i,1,1) .eq. c0) then
           c1pn(i) = cmax
         else
           c1pn(i) = -var*(etap(i,jrow-1)-eta(i,jrow-1))/var1(i,1,1)
           if (c1pn(i) .le. c0)   c1pn(i) = c0
           if (c1pn(i) .gt. cmax) c1pn(i) = cmax
         endif
         etap(i,jrow)  = eta(i,jrow) 
     &             - c1pn(i)*dtsf*dyur(jrow)*(eta(i,jrow)-eta(i,jrow-1))
      enddo
# endif
# ifdef obc_west
      i = 2
      do jrow=js,je
         var = -dxu(i+1)*csu(jrow)/dtsf
         var2(1,jrow,1)= etap(i+2,jrow) - etap(i+1,jrow)
         cmax = -c2*sqrt(grav*zw(kmt(i,jrow)))
         if (var2(1,jrow,1) .eq. c0) then
           c1pw(jrow) = cmax
         else
          c1pw(jrow) = var*(etap(i+1,jrow)-eta(i+1,jrow))/var2(1,jrow,1)
           if (c1pw(jrow) .ge. c0)   c1pw(jrow) = c0
           if (c1pw(jrow) .lt. cmax) c1pw(jrow) = cmax
         endif
         etap(i,jrow)  = eta(i,jrow) - c1pw(jrow)*dtsf*
     &                   dxur(i)*csur(jrow)*(eta(i+1,jrow)-eta(i,jrow))
      enddo
# endif
# ifdef obc_east
      i = imtm1
      do jrow=js,je
         var = dxu(i-1)*csu(jrow)/dtsf
         var2(1,jrow,1)= etap(i-1,jrow) - etap(i-2,jrow)
         cmax = c2*sqrt(grav*zw(kmt(i,jrow)))
         if (var2(1,jrow,1) .eq. c0) then
           c1pe(jrow) = cmax
         else
           c1pe(jrow) = -var*(etap(i-1,jrow)-eta(i-1,jrow))
     &                       /var2(1,jrow,1)
           if (c1pe(jrow) .le. c0)   c1pe(jrow) = c0
           if (c1pe(jrow) .gt. cmax) c1pe(jrow) = cmax
         endif
         etap(i,jrow)  = eta(i,jrow) - c1pe(jrow)*dtsf*
     &                   dxur(i)*csur(jrow)*(eta(i,jrow)-eta(i-1,jrow))
      enddo
# endif
c
# if defined explicit_polar_filter 
c
c-----------------------------------------------------------------------
c       spatially filter time tendency of eta at high latitudes
c-----------------------------------------------------------------------
c
        rdtsf = c1/dtsf
        do jrow=jscomp,jecomp
          if ((jrow.le.jft1 .or. jrow.ge.jft2) .and. jrow.ge.jfrst) then
	    do i=1,imt
	      dedt(i,jrow)    = (etap(i,jrow)-eta(i,jrow))*rdtsf   
	    enddo
#  if defined firfil
            call fast_fir (dedt(1,jrow), imt, 1, 1, emsk(1,jrow)
     &,                    numflt(jrow))
#  else
            call fileta (dedt(1,jrow), jrow)
#  endif
	    do i=1,imt
	      etap(i,jrow)    = eta(i,jrow) + dtsf*dedt(i,jrow)  
	    enddo
          endif
	enddo
# endif
c
# ifdef symmetry
c
c-----------------------------------------------------------------------
c       set meridional symmetry
c-----------------------------------------------------------------------
c
        if (je .eq. jmt) then
          do i=1,imt
	    etap(i,je) = etap(i,je-2)
	  enddo
        endif
# endif
c
c-----------------------------------------------------------------------
c       diagnose surface pressure and absorb rho0
c-----------------------------------------------------------------------
c
        ps(:,:,1) = etap(:,:)*grav*rhosurf(:,:)/rho0
# ifdef parallel_1d
        call mpp_update_domains( ps, mom_domain )
# endif
c
c-----------------------------------------------------------------------
c       time step hubar using updated ps.
c       this is the "backward" step in the forward-backward scheme.
c-----------------------------------------------------------------------
c
#if defined obc        
# if defined orlanski
#  ifdef obc_south
        jrow  = 2
        etamin =  1.e99
        etamax = -1.e99
        do i=2,imtm1
           etamin = MIN(etamin,eta(i,jrow+1))
           etamax = MAX(etamax,eta(i,jrow+1))
        enddo
        if (ABS(etamin) .gt. ABS(etamax)) then
           etashift = etamax
        else
           etashift = etamin
        endif
        do i=2,imtm1
           etarel(i,jrow) = eta(i,jrow+1)-(etashift-obc_etanull)
        enddo
#  endif
#  ifdef obc_north
        jrow  = jmtm1
        etamin =  1.e99
        etamax = -1.e99
        do i=2,imtm1
           etamin = MIN(etamin,eta(i,jrow-1))
           etamax = MAX(etamax,eta(i,jrow-1))
        enddo
        if (ABS(etamin) .gt. ABS(etamax)) then
           etashift = etamax
        else
           etashift = etamin
        endif
        do i=2,imtm1
           etarel(i,jrow) = eta(i,jrow-1)-(etashift-obc_etanull)
        enddo
#  endif
#  ifdef obc_west
        i  = 2
        etamin =  1.e99
        etamax = -1.e99
        do jrow=2,jmtm1
           etamin = MIN(etamin,eta(i+1,jrow))
           etamax = MAX(etamax,eta(i+1,jrow))
        enddo
        if (ABS(etamin) .gt. ABS(etamax)) then
           etashift = etamax
        else
           etashift = etamin
        endif
        do jrow=2,jmtm1
           etarel(i,jrow) = eta(i+1,jrow)-(etashift-obc_etanull)
        enddo
#  endif
#  ifdef obc_east
        i  = imtm1
        etamin =  1.e99
        etamax = -1.e99
        do j=2,jmtm1
           etamin = MIN(etamin,eta(i-1,jrow))
           etamax = MAX(etamax,eta(i-1,jrow))
        enddo
        if (ABS(etamin) .gt. ABS(etamax)) then
           etashift = etamax
        else
           etashift = etamin
        endif
        do j=2,jmtm1
           etarel(i,jrow) = eta(i-1,jrow)-(etashift-obc_etanull)
        enddo
#  endif
# else
        call addobceta(js,je,etarel)
# endif
#endif
c
        do jrow=jscomp,jecomp
# ifndef rot_grid
	  fc = cori(jrow,1)*p5
	  fx = fc*dtsf
	  fy = c1/(c1 + fx**2)
# endif
	  do i=2,imtm1
# ifdef rot_grid
	    fc = cori(i,jrow,1)*p5
	    fx = fc*dtsf
	    fy = c1/(c1+fx**2)
# endif
c
c           surface pressure gradients + forcing + explicit coriolis
c
	    utempm(1) = hubar(i,jrow,1) + dtsf*(
     &                  - depthu(i,jrow)*csur(jrow)*dxu2r(i)
     &                  *( ps(i+1,jrow+1,1) + ps(i+1,jrow,1)
     &                   - ps(i  ,jrow+1,1) - ps(i  ,jrow,1) )
     &                  + zu(i,jrow,1)+fc*hubar(i,jrow,2)
     &                  )
	    utempm(2) = hubar(i,jrow,2) + dtsf*(
     &                  - depthu(i,jrow)*dyu2r(jrow)
     &                  *( ps(i,jrow+1,1) + ps(i+1,jrow+1,1)
     &                   - ps(i,jrow,1) - ps(i+1,jrow,1) )
     &                  + zu(i,jrow,2)-fc*hubar(i,jrow,1)
     &                  )
c       
#ifdef obc
c
c           add reference water level
c     
#  ifdef obc_south
            if (jrow .eq. 2 ) then
               utempm(2) = utempm(2) - depthu(i,jrow)*grav*dyur(jrow)
     &              *dtsf*(eta(i,jrow+1) - etarel(i,jrow)) * obcs
           endif
#  endif
#  ifdef obc_north
#   ifdef obc_north_anywhere
            if (jrow .eq. jobcn-1) then
#   else
            if (jrow .eq. jmt-1) then
#   endif
               utempm(2) = utempm(2) - depthu(i,jrow)*grav*dyur(jrow)
     &              *dtsf*(etarel(i,jrow) - eta(i,jrow-1)) * obcn
            endif
#  endif
#  ifdef obc_west
            if (i .eq. 2) then
               utempm(1) = utempm(1) - depthu(i,jrow)*grav*dxur(i)
     &              *csur(jrow)*dtsf*(eta(i+1,jrow) - etarel(i,jrow)) 
     &              * obcw
            endif
#  endif
#  ifdef obc_east
            if (i .eq. imt-1) then
               utempm(1) = utempm(1)  - depthu(i,jrow)*grav*dxur(i)
     &              *csur(jrow)*dtsf*(etarel(i,jrow) - eta(i-1,jrow))
     &              * obce
            endif
# endif
#endif
c
c           solve for hubar using implicit coriolis piece
c
	    hubar(i,jrow,1) = umsk(i,jrow)*(utempm(1) + fx*utempm(2))*fy
	    hubar(i,jrow,2) = umsk(i,jrow)*(utempm(2) - fx*utempm(1))*fy
	  enddo
# ifdef cyclic
	  hubar(  1,jrow,1) = hubar(imtm1,jrow,1)
	  hubar(  1,jrow,2) = hubar(imtm1,jrow,2)
	  hubar(imt,jrow,1) = hubar(2,jrow,1)
	  hubar(imt,jrow,2) = hubar(2,jrow,2)
# endif
        enddo
# ifdef symmetry
        if (je .eq. jmt) then 
	  do i=1,imt
	    hubar(i,je  ,1) = hubar(i,je-2,1)
	    hubar(i,je-1,2) = c0
	  enddo
        endif
# endif
c
# if defined explicit_polar_filter 
c
c-----------------------------------------------------------------------
c       spatially filter time tendency of hubar at high latitudes
c-----------------------------------------------------------------------
c
        rdtsf = c1/dtsf
        do jrow=jscomp,jecomp
          if ((jrow.le.jfu1 .or. jrow.ge.jfu2) .and. jrow.ge.jfrst) then
	    do i=1,imt
	      dhudt(i,jrow,1) = (hubar(i,jrow,1)-hu0(i,jrow,1))*rdtsf
     &                        *umsk(i,jrow)   
	      dhudt(i,jrow,2) = (hubar(i,jrow,2)-hu0(i,jrow,2))*rdtsf   
     &                        *umsk(i,jrow)   
	    enddo
#  if defined firfil
            fx = sign(1.0,phi(jrow))
            do i=1,imt
              tempik(i,1) = -fx*dhudt(i,jrow,1)*spsin(i)
     &                             - dhudt(i,jrow,2)*spcos(i)
              tempik(i,2) =  fx*dhudt(i,jrow,1)*spcos(i)
     &                             - dhudt(i,jrow,2)*spsin(i)
	    enddo
c
	    call fast_fir (tempik, imt, 1, 2, umsk(1,jrow),numflu(jrow))
c
            do i=1,imt
              dhudt(i,jrow,1) = fx*(-tempik(i,1)*spsin(i)
     &                             + tempik(i,2)*spcos(i))
              dhudt(i,jrow,2) = -tempik(i,1)*spcos(i)
     &                         - tempik(i,2)*spsin(i)
            enddo
#  else
            call filhat (dhudt(1,jrow,1), dhudt(1,jrow,2), jrow)
#  endif
	    do i=1,imt
	      hubar(i,jrow,1) = hu0(i,jrow,1) + dtsf*dhudt(i,jrow,1)   
	      hubar(i,jrow,2) = hu0(i,jrow,2) + dtsf*dhudt(i,jrow,2)   
	    enddo
          endif
	enddo
# endif
c
# ifdef parallel_1d
        call mpp_update_domains (hubar, mom_domain)
# endif
c
c-----------------------------------------------------------------------
c       accumulate hubar using ubar as a temporary field 
c       accumulate etat in etat(,,3)
c       shift etap = eta(t+1) into eta(t) for next barotropic step
c-----------------------------------------------------------------------
c
        do jrow=js,je
	  do i=1,imt
            ubar(i,jrow,1) = ubar(i,jrow,1) + hubar(i,jrow,1)
            ubar(i,jrow,2) = ubar(i,jrow,2) + hubar(i,jrow,2)
            etat(i,jrow,3) = etat(i,jrow,3) + etap(i,jrow)
            eta(i,jrow)    = etap(i,jrow)           
# if defined explicit_polar_filter 
            hu0(i,jrow,1)  = hubar(i,jrow,1)
            hu0(i,jrow,2)  = hubar(i,jrow,2)
# endif
          enddo
        enddo
c
      enddo  !end of barotropic time loop
c
c-----------------------------------------------------------------------
c     Compute time averaged ubar and eta centered at baroclinic 
c     step taup1 = T0+dtuv. 
c-----------------------------------------------------------------------
c    
      rnb = 1.0/(float(nb)+1.0)      
      do jrow=js,je
	do i=1,imt
          ubar(i,jrow,1) = ubar(i,jrow,1)*depthur(i,jrow)*rnb
          ubar(i,jrow,2) = ubar(i,jrow,2)*depthur(i,jrow)*rnb
          etat(i,jrow,3) = etat(i,jrow,3)*rnb
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     Update etat(,,1) to baroclinic time T0+dtuv using a leap-frog.  
c     convU, sff, and source_eta are at baroclinic time tau=T0.  
c     Also update etat(,,2) to be the time mean etat(tau).
c     No need to Robert time filter since the etat(taum1) value is
c     a time average.
c-----------------------------------------------------------------------
c
      do jrow=jscomp,jecomp
        eta(:,jrow) = etat(:,jrow,1)
	do i=2,imtm1
	  etat(i,jrow,1) = etat(i,jrow,2) + c2dtuv*(
     &                                        convU(i,jrow) 
#  ifdef explicit_fresh_water_flux
     &                                       +sfft(i,jrow)
#  endif
#  ifdef xlandmix_eta
     &                                       +source_eta(i,jrow)
#  endif
     &                                       )*emsk(i,jrow)
	enddo
#  ifdef cyclic
	etat(  1,jrow,1) = etat(imtm1,jrow,1)
        etat(imt,jrow,1) = etat(2,jrow,1)
#  else
	etat(  1,jrow,1) = c0
        etat(imt,jrow,1) = c0
#  endif
        etat(:,jrow,2)   = eta(:,jrow)
      enddo
c
# if defined explicit_eta_laplacian
c     
c-----------------------------------------------------------------------
c     Smooth free surface height with a laplacian. 
c     Do smoothing only once per baroclinic time step.
c     The laplacian conserves total volume.  
c
c     eta_mix is a mixing namelist parameter. If it is too large, 
c     then smoothing can adversely affect other fields, especially 
c     the vertical velocity, and it will overly dissipate the free 
c     surface height. Tests with a global 3degree model indicate 
c     eta_mix can be order 1e7. Cosine tapering is used to avoid 
c     numerical stability problems in the high latitudes. 
c
c     If neta=2, do smoothing on etat(tau).  This will 
c     keep etat(taup1) as a pure leap-frog, which helps with
c     conservation properties of total tracer quantity. 
c     The alternative is neta=1, which does eta-laplacian
c     on etat(taup1).  
c-----------------------------------------------------------------------
c
# ifdef parallel_1d
        call mpp_update_domains (etat, mom_etat_domain)
# endif
      neta=2

      do jrow=js,je
	do i=1,imt
	  tempeta(i,jrow) = etat(i,jrow,neta)
        enddo
      enddo 

      wfactor = c2*dtuv*eta_mix
      do jrow=js,je
        weighte(jrow) = wfactor*cst(jrow)
        weightn(jrow) = wfactor*csu(jrow)
      enddo

      do jrow=jscomp,jecomp
        do i=2,imtm1
          fei   = weighte(jrow)*(tempeta(i+1,jrow)-tempeta(i,jrow))
     &          *cstr(jrow)*dxur(i)*emsk(i+1,jrow)
          feim1 = weighte(jrow)*(tempeta(i,jrow)-tempeta(i-1,jrow))
     &          *cstr(jrow)*dxur(i-1)*emsk(i-1,jrow)
          fnj   = weightn(jrow)*(tempeta(i,jrow+1)-tempeta(i,jrow))
     &          *csu_dyur(jrow)*emsk(i,jrow+1)
          fnjm1 = weightn(jrow-1)*(tempeta(i,jrow)-tempeta(i,jrow-1))
     &          *csu_dyur(jrow-1)*emsk(i,jrow-1)
          smooth_lap = (fei-feim1)*cstr(jrow)*dxtr(i)
     &                +(fnj-fnjm1)*cstr(jrow)*dytr(jrow) 
          etat(i,jrow,neta) = etat(i,jrow,neta)+smooth_lap*emsk(i,jrow)  
         enddo
      enddo

# ifdef parallel_1d
        call mpp_update_domains (etat, mom_etat_domain)
# endif
c
      do jrow=js,je
#  ifdef cyclic
        etat(  1,jrow,neta) = etat(imtm1,jrow,neta)
        etat(imt,jrow,neta) = etat(    2,jrow,neta)
#  else 
        etat(  1,jrow,neta) = c0
        etat(imt,jrow,neta) = c0
#  endif
      enddo
c 
# endif
c
c-----------------------------------------------------------------------
c     update surface pressure to new baroclinic time step taup1
c-----------------------------------------------------------------------
c
      do jrow=js,je
	do i=1,imt
          ps(i,jrow,1) = etat(i,jrow,1)*grav*rhosurf(i,jrow)
        enddo
      enddo 
c
c-----------------------------------------------------------------------
c     update boundaries 
c-----------------------------------------------------------------------
c
# ifdef parallel_1d
      call mpp_update_domains (etat, mom_etat_domain)
      call mpp_update_domains (ps, mom_domain)
# endif

# ifdef xlandmix_eta
c     zero source_eta so it can be set in xland.F for next time step:
      do jrow=js,je
	do i=1,imt
          source_eta(i,jrow)=c0
        enddo
      enddo
# endif

# if defined time_step_monitor
c
c-----------------------------------------------------------------------
c     construct scalar diagnostics
c-----------------------------------------------------------------------
c     
      if(tsits) then
        call diag_fs
      endif
# endif
      
# ifdef timing
      call toc ('barotropic', 'expl_freesurf')
# endif
      return
      end
#endif



#if defined explicit_free_surface && !defined time_averaged_surface
c
c=======================================================================
c     calculate the barotropic velocities and the sea surface elevation
c     by sub-cycling the barotropic equations over one baroclinic time
c     step
c
c     js = global starting latitude index for the task
c     je = global ending latitude index for the task
c     
c     adopted from
c
c     P.D. Killworth, D.Stainforth, D.J. Webb and S.M. Paterson
c     IOS Deacon Laboratory Report No. 270 (1989)
c     and JPO vol 21, 1333-1348 (1991)
c
c     author: M. Schmidt      e-mail=> Martin.Schmidt@io-warnemuende.de
c=======================================================================
c
c
c     additional variables:
c     etad   - temporary field, contains eta
c     etagrd - time change of the sea surface elevation
c     uhat   - vertically integrated barotropic velocity
c
# include "switch.h"
      real, dimension(:,:), allocatable :: eta, etad, etagrd
      real, dimension(:,:,:), allocatable :: uhat
      dimension utempm(2)
# ifdef timing
      call tic ('barotropic', 'expl_freesurf')
# endif
c
      mscan      = 0
      gravrrho0r = rho0r/grav
c
c-----------------------------------------------------------------------
c     initialize fields
c-----------------------------------------------------------------------
c
      allocate (eta   (imt,js:je))
      allocate (etad  (imt,js:je))
      allocate (etagrd(imt,js:je))
      allocate (uhat  (imt,js:je,2))
c
c-----------------------------------------------------------------------
c     set eta, etad, etagrd 
c-----------------------------------------------------------------------
c
      do jrow=js,je
	do i=1,imt
	  eta(i,jrow)    = ps(i,jrow,1)*gravrrho0r
# ifndef explicit_dpdc
          etad(i,jrow)   = eta(i,jrow)
# endif
          etagrd(i,jrow) = c0
	enddo
      enddo
c
c-----------------------------------------------------------------------
c     and compute transports for the horizontal divergence of u
c-----------------------------------------------------------------------
c
      do jrow=js,je
        do i=1,imt
	  uhat(i,jrow,1) = ubar(i,jrow,1)*h(i,jrow)
	  uhat(i,jrow,2) = ubar(i,jrow,2)*h(i,jrow)
	enddo
      enddo
c
c-----------------------------------------------------------------------
c     shift ubar and ps
c-----------------------------------------------------------------------
c
      do jrow=js,je
	do i=1,imt
          ubarm1(i,jrow,1) = ubar(i,jrow,1)
          ubarm1(i,jrow,2) = ubar(i,jrow,2)
          ps(i,jrow,2)     = ps(i,jrow,1)
        enddo
      enddo
# ifdef trace_indices
      write (stdout,*) 'expl: ubar  -> ubarm1'
      write (stdout,*) 'expl: ps(1) -> ps(2)'
# endif
# ifdef explhmix
c
c-----------------------------------------------------------------------
c     start loop to subtract friction term from the forcing terms zu
c-----------------------------------------------------------------------
c
      do m=1,2
	msig = 2 * m - 3
	do jrow=jscomp,jecomp
	  do i=2,imtm1
	    fx = am*(
     &       bhd(jrow)*dxur(i)*
     &      (  dxtr(i+1)*( ubar(i+1,jrow,m) - ubar(i  ,jrow,m))
     &       - dxtr(i  )*( ubar(i  ,jrow,m) - ubar(i-1,jrow,m)) )
     &       + chd(jrow)*dytr(jrow+1)
     &                  *( ubar(i,jrow+1,m) - ubar(i,jrow  ,m))
     &       - dhd(jrow)*dytr(jrow  )
     &                  *( ubar(i,jrow  ,m) - ubar(i,jrow-1,m))
     &       + ghd(jrow)*ubar(i,jrow,m)
     &       + msig*hhd(jrow)*dxur(i)
     &                  *( ubar(i+1,jrow,3-m)-ubar(i-1,jrow,3-m)))
c
	    zu(i,jrow,m) = zu(i,jrow,m) - fx 
	  enddo
	enddo
      enddo
#  ifdef parallel_1d
        call mpp_update_domains (zu, mom_domain)
#  endif
#  if defined debug_expl_freesurf
      write (stdout,*) 'expl_freesurf ck2d: js=',js,' je=',je,' pn=',pn
      write (stdout,*) 'ck2d zu1=',ck2d(zu(1,js,1), js, je)
      write (stdout,*) 'ck2d zu2=',ck2d(zu(1,js,2), js, je)
#  endif
# endif
# ifdef symmetry
c
c-----------------------------------------------------------------------
c  set meridional forcing on symmetry row to zero
c-----------------------------------------------------------------------
c
      if (je .eq. jmt) then   
        do i=1,imt
	  zu(i,je-1,2) = c0
        enddo
      endif
# endif
c
c=======================================================================
c     sub-cycle barotropic equation for "nb" barotropic time steps
c=======================================================================
c
      nb = nint(dtuv/dtsf)
c
      do itbt=1,nb
# if defined debug_expl_freesurf
c      write (stdout,*) 'itbt=',itbt
# endif
# ifdef explicit_dpdc
c
c-----------------------------------------------------------------------
c       remove local checkerboards with del-plus minus del-cross filter
c       set pseudo-eta values on land for dlpl
c       note: nippts(1) is used to store the number of perimeter land
c             points
c-----------------------------------------------------------------------
c
        do n=1,nippts(1)
	  i = iperm(n)
	  jrow = jperm(n)
	  eta(i,jrow) = ( emsk(i+1,jrow) * eta(i+1,jrow)
     &                  + emsk(i-1,jrow) * eta(i-1,jrow)
     &                  + emsk(i,jrow+1) * eta(i,jrow+1)
     &                  + emsk(i,jrow-1) * eta(i,jrow-1))*splr(n)
        enddo
c
        if (je .eq. jmt) then
	  do i=2,imtm1
            eta(i,je) = eta(i,je-1)
	  enddo
        endif
c      
        if (js .eq. 1) then
	  do i=2,imtm1
	    eta(i,1) = eta(i,2)
	  enddo
        endif
#  ifdef parallel_1d
        call mpp_update_domains (eta, mom_domain)
#  endif
        do jrow=js,je
#  ifdef cyclic
	  eta(  1,jrow) = eta(imtm1,jrow)
	  eta(imt,jrow) = eta(    2,jrow)
#  else
	  eta(  1,jrow) = eta(    2,jrow)
	  eta(imt,jrow) = eta(imtm1,jrow)
#  endif
        enddo
        do jrow=js,je
	  do i=1,imt
	    etad(i,jrow) = eta(i,jrow)
	  enddo
        enddo
# endif
# ifdef explicit_eb
c
c-----------------------------------------------------------------------
c       start of loop for euler backward step
c-----------------------------------------------------------------------
c
        do id=1,2
c
# else
          id = 1
# endif
c-----------------------------------------------------------------------
c         compute etagrd in a forward step
c-----------------------------------------------------------------------
c
# ifdef explicit_dpdc
c         set wght to an appropriate value. Do some experiments.
          wght = c1
          fx = dtsf*wght*grav*zt(km)
# endif
          do jrow=jscomp,jecomp
	    do i=2,imtm1
# ifdef explicit_dpdc
	      dlpl = eta(i,jrow+1) + eta(i+1,jrow)
     &             + eta(i-1,jrow) + eta(i,jrow-1) - c4*eta(i,jrow)
c
	      dlcr = p5* ( emsk(i+1,jrow+1)*eta(i+1,jrow+1)
     &                +emsk(i-1,jrow+1)*eta(i-1,jrow+1)
     &                +emsk(i+1,jrow-1)*eta(i+1,jrow-1)
     &                +emsk(i-1,jrow-1)*eta(i-1,jrow-1)
     &                - ( emsk(i+1,jrow+1) + emsk(i-1,jrow+1)
     &                  + emsk(i+1,jrow-1) + emsk(i-1,jrow-1))
     &                       *eta(i,jrow))
# endif
	      etagrd(i,jrow) = dxt2r(i)*dytr(jrow)
     &            *( uhat(i  ,jrow,1  )*dyu(jrow)
     &              +uhat(i,jrow-1,1  )*dyu(jrow-1)
     &            -( uhat(i-1,jrow,1  )*dyu(jrow)
     &              +uhat(i-1,jrow-1,1)*dyu(jrow-1))
     &            +( uhat(i,jrow,2  )*dxu(i) + uhat(i-1,jrow,2  )
     &              *dxu(i-1))*csu(jrow)
     &            -( uhat(i,jrow-1,2)*dxu(i) + uhat(i-1,jrow-1,2)
     &              *dxu(i-1))*csu(jrow-1))
# ifdef explicit_dpdc
     &             -fx*dxtr(i)*dytr(jrow)*( dlpl - dlcr )
# endif
	    enddo
          enddo
# ifdef parallel_1d
          call mpp_update_domains (etagrd, mom_domain)
# endif
# if defined debug_expl_freesurf
c      write (stdout,*) 'ck2d etagrd=',ck2d(etagrd(1,js), js, je)
# endif
# ifdef explicit_efb
c
c-----------------------------------------------------------------------
c         compute etad from etagrd and eta.
c-----------------------------------------------------------------------
c
          do jrow=js,je
	    do i=1,imt
	      etad(i,jrow) = emsk(i,jrow) * (eta(i,jrow)
     &                      - cstr(jrow)*dtsf*etagrd(i,jrow))
	    enddo
#  ifdef cyclic
	    etad(  1,jrow) = etad(imtm1,jrow)
            etad(imt,jrow) = etad(2,jrow)
#  endif
          enddo
#  ifdef symmetry
          if (je .eq. jmt) then
            do i=1,imt
	      etad(i,je) = etad(i,je-2)
	    enddo
          endif
#  endif
# endif
c
c-----------------------------------------------------------------------
c         start loop to get uhat from etad and ut
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c         calculate  uhat
c-----------------------------------------------------------------------
c
          do jrow=jscomp,jecomp
# ifndef rot_grid
	    fc = cori(jrow,1)*p5
	    fx = fc*dtsf
	    fy = c1/(c1 + fx**2)
# endif
	    do i=2,imtm1
# ifdef rot_grid
	      fc = cori(i,jrow,1)*p5
	      fx = fc*dtsf
	      fy = c1/(c1+fx**2)
# endif
c
c-----------------------------------------------------------------------
c             calculate the surface pressure gradients
c-----------------------------------------------------------------------
c
	      utempm(1) = - h(i,jrow)* grav*csur(jrow)*dxu2r(i)
     &                    *( etad(i+1,jrow+1) + etad(i+1,jrow)
     &                     - etad(i  ,jrow+1) - etad(i  ,jrow) )
	      utempm(2) = - h(i,jrow)*grav*dyu2r(jrow)
     &                    *( etad(i,jrow+1) + etad(i+1,jrow+1)
     &                     - etad(i,jrow  ) - etad(i+1,jrow  ) )
	      do m=1,2
	        msig = 2*m-3
	        mc   = 3-m
# ifdef explhmix
c
c-----------------------------------------------------------------------
c               calculate friction term
c-----------------------------------------------------------------------
c
	        fric = am*h(i,jrow)*(bhd(jrow)*dxur(i)*
     &    ( dxtr(i+1)*( ubar(i+1,jrow,m) - ubar(i  ,jrow,m))
     &    - dxtr(i  )*( ubar(i  ,jrow,m) - ubar(i-1,jrow,m)) )
     &       + chd(jrow)*dytr(jrow+1)
     &                     *( ubar(i,jrow+1,m)- ubar(i,jrow  ,m))
     &       - dhd(jrow)*dytr(jrow  )
     &                     *( ubar(i,jrow  ,m)- ubar(i,jrow-1,m))
     &       + ghd(jrow)*ubar(i,jrow,m)
     &       + msig*hhd(jrow)*dxur(i)
     &              *(ubar(i+1,jrow,mc) - ubar(i-1,jrow,mc)))
# endif
c
c-----------------------------------------------------------------------
c               add in forcing, explicit coriolis term and friction
c-----------------------------------------------------------------------
c
	        utemp = utempm(m) + ( zu(i,jrow,m)
     &                             - h(i,jrow)*msig*fc*ubar(i,jrow,mc) )
# ifdef explhmix
     &                  + fric
# endif
c
	        utempm(m) = ubar(i,jrow,m)*h(i,jrow) + dtsf*utemp
	      enddo
c
c-----------------------------------------------------------------------
c             add in implicit coriolis term
c-----------------------------------------------------------------------
c
	      uhat(i,jrow,1) = (utempm(1) + fx*utempm(2))*fy
	      uhat(i,jrow,2) = (utempm(2) - fx*utempm(1))*fy
	    enddo
          enddo
# ifdef symmetry
c
c-----------------------------------------------------------------------
c         set symmetry conditions at uhat and vhat
c-----------------------------------------------------------------------
c
          if (je .eq. jmt) then 
	    do i=1,imt
	      uhat(i,je  ,1) = uhat(i,je-2,1)
	      uhat(i,je-1,2) = c0
	    enddo
          endif
# endif
c
# ifdef cyclic
          do jrow=jscomp,jecomp
	    uhat(  1,jrow,1) = uhat(imtm1,jrow,1)
	    uhat(  1,jrow,2) = uhat(imtm1,jrow,2)
	    uhat(imt,jrow,1) = uhat(2,jrow,1)
	    uhat(imt,jrow,2) = uhat(2,jrow,2)
          enddo
# endif
# ifdef explicit_eb
c
c-----------------------------------------------------------------------
c         compute etad from etagrd and eta.
c-----------------------------------------------------------------------
c
          do jrow=jscomp,jecomp
	    do i=1,imt
	      etad(i,jrow) = emsk(i,jrow) * (eta(i,jrow)
     &                       - cstr(jrow)*dtsf*etagrd(i,jrow))
	    enddo
#  ifdef cyclic
	    etad(  1,jrow) = etad(imtm1,jrow)
	    etad(imt,jrow) = etad(2,jrow)
#  endif
          enddo
#  ifdef symmetry
          if (je .eq. jmt) then
	    do i=1,imt
	      etad(i,je) = etad(i,je-2)
	    enddo
          endif
#  endif
# endif
c
# ifdef explicit_fresh_water_flux
c
c-----------------------------------------------------------------------
c         add the fresh water flux to eta
c         and accumulate the total fresh water flux
c-----------------------------------------------------------------------
c
#  ifdef explicit_eb
          if (id.eq.2) then
#  endif
          do jrow=jscomp,jecomp
	    do i=1,imt
              etad(i,jrow  ) = etad(i,jrow  ) + sfft(i,jrow)*dtsf
            enddo
          enddo
#  ifdef explicit_eb
          endif 
#  endif
# endif
# ifdef parallel_1d
#  ifdef explicit_eb
          call mpp_update_domains (etad, mom_domain)
#  endif
          call mpp_update_domains (uhat, mom_domain)
# endif
# if defined debug_expl_freesurf
c          write (stdout,*) 'ck2d etad=',ck2d(etad(1,js), js, je)
c          write (stdout,*) 'ck2d uhat1=',ck2d(uhat(1,js,1), js, je)
c          write (stdout,*) 'ck2d uhat2=',ck2d(uhat(1,js,2), js, je)
# endif
c
# ifdef explicit_eb
c
c-----------------------------------------------------------------------
c       end of loop for euler backward step
c-----------------------------------------------------------------------
c
        enddo
# endif
# if defined explicit_polar_filter 
c
c-----------------------------------------------------------------------
c       filter uhat, vhat, and eta at high latitudes
c-----------------------------------------------------------------------
c
        do jrow=jscomp,jecomp
          if ((jrow.le.jfu1 .or. jrow.ge.jfu2) .and. jrow.ge.jfrst) then
            call filhat (uhat(1,jrow,1), uhat(1,jrow,2), jrow)
          endif
          if ((jrow.le.jft1 .or. jrow.ge.jft2) .and. jrow.ge.jfrst) then
            call fileta (etad(1,jrow), jrow)
          endif
        enddo
# endif
c
c-----------------------------------------------------------------------
c       compute ubar and set eta equal to etad and add fresh water flux
c-----------------------------------------------------------------------
c
        do jrow=js,je
	  do i=1,imt
            ubar(i,jrow,1) = uhat(i,jrow,1)*hr(i,jrow)
            ubar(i,jrow,2) = uhat(i,jrow,2)*hr(i,jrow)
            eta (i,jrow  ) = etad(i,jrow  )           
          enddo
        enddo
      enddo
c
c=======================================================================
c     end of sub-cycle barotropic equation loop
c=======================================================================
c
c-----------------------------------------------------------------------
c     analyze the mass balance
c-----------------------------------------------------------------------
c     
      if (tsits) then     
c       
        do m=1,4
          etaana(m) = c0
        enddo
        do jrow=jscomp,jecomp
          fu = c1 
# ifdef symmetry
          if (jrow.eq.jmtm1) fu = p5
# endif   
# if defined cyclic
          do i=2,imtm1
# else      
            do i=1,imt
#endif        
              fx = dxt(i)*cst(jrow)*dyt(jrow)
# ifdef explicit_fresh_water_flux
              etaana(1) = etaana(1) + fx*sfft(i,jrow)*dtsf*nb
# endif       
              etaana(2) = etaana(2) + fx*eta(i,jrow)
              etaana(4) = etaana(4) + fx*eta(i,jrow)**2
              etaana(3) = etaana(3) + fu*dxu(i)*csu(jrow)*dyu(jrow)
     &             *h(i,jrow)*(ubar(i,jrow,1)**2 + ubar(i,jrow,2)**2)
            enddo
          enddo
          etaana(4) = etaana(4)*p5*grav*rho0
          etaana(3) = etaana(3)*p5*rho0
c         
        endif
c-----------------------------------------------------------------------
c     put results in ps
c-----------------------------------------------------------------------
c
      do jrow=js,je
	do i=1,imt
          ps(i,jrow,1) = eta(i,jrow)*grav*rho0
        enddo
      enddo
c
#  if defined debug_expl_freesurf
      write (stdout,*) 'ck2d ps(,,1)=',ck2d(ps(1,js,1), js, je)
      write (stdout,*) 'ck2d uhat1=',ck2d(uhat(1,js,1), js, je)
      write (stdout,*) 'ck2d uhat2=',ck2d(uhat(1,js,2), js, je)
#  endif
# ifdef trace_indices
      write (stdout,*) 'expl: new ubar and ps(1) ready'
# endif
      deallocate(eta, etad, etagrd, uhat)
# ifdef timing
      call toc ('barotropic', 'expl_freesurf')
# endif
      return
      end

#endif
#if !defined explicit_free_surface
      return
      end
#endif


#if defined explicit_polar_filter 
      subroutine fileta (eta, jrow)
c
c=======================================================================
c     subroutine fileta sets up input needed for fourier filtering
c     (when the "fourfil" option is defined) -or- symmetric finite
c     impulse response filtering (when the "firfil" option is defined)
c     of eta at the specified high latitude row "jrow".
c=======================================================================
c
      use arrays_2d
      use memory_window
# include "param.h"
# include "grdvar.h"
# include "index.h"
# include "mw.h"
      dimension eta(imt)
      dimension tempik(imt)
c
# ifdef timing
      call tic ('eta', 'filtering')
# endif
c
c---------------------------------------------------------------------
c     fourier filter eta at high latitudes
c---------------------------------------------------------------------
c
      jj = jrow-jfrst+1
c
      if (phi(jrow) .gt. c0) then
        jj = jj - jskpt + 1
        refcos = refcosn
      else
        refcos = refcoss
      endif
c
c    if previous strips were of same length, do not recompute
c    fourier coeffs
c
      isave = 0
      ieave = 0
      do l=1,lsegf
        k=1
        if (istf(jj,l,k) .ne. 0) then
          is    = istf(jj,l,k)
          ie    = ietf(jj,l,k)
          iredo = 0
          if (is.ne.isave .or. ie.ne.ieave) then
            iredo = -1
            isave = is
            ieave = ie
            im = ie-is+1
# ifdef cyclic
            if (im.ne.imtm2 .or. kmt(1,jrow).lt.k) then
              m = 1
              n = nint(im*cst(jrow)/refcos)
            else
              m = 3
              n = nint(im*cst(jrow)*0.5/refcos)
            endif
# else
            m = 1
            n = nint(im*cst(jrow)/refcos)
# endif
          endif
c
          idx  = iredo+1
          ism1 = is-1
          iea  = ie
          if (ie .ge. imt) iea = imtm1
          do i=is,iea
            tempik(i-ism1) = eta(i)
          enddo
          if (ie .ge. imt) then
            ieb = ie-imtm2
            ii  = imtm1-is
            do i=2,ieb
              tempik(i+ii) = eta(i)
            enddo
          endif
c
# ifdef fourfil
          call filtr (tempik(1), im, m, n, idx)
# endif
# ifdef firfil
          call filfir (tempik(1), im, m, numflt(jrow))
# endif
c
          do i=is,iea
            eta(i) = tempik(i-ism1)
          enddo
          if (ie .ge. imt) then
            do i=2,ieb
              eta(i) = tempik(i+ii)
            enddo
          endif
        endif
      enddo
# ifdef cyclic
      eta(1)   = eta(imt-1)
      eta(imt) = eta(2)
# else
      eta(1)   = c0
      eta(imt) = c0
# endif
c
# ifdef timing
      call toc ('eta', 'filtering')
# endif
      return
      end
#endif


      subroutine init_surface
#if defined explicit_free_surface
c
c=======================================================================
c     Initialize some terms for free surface.
c=======================================================================
c
      use arrays_2d
      use mpp_io_mod      
#include "param.h"
#include "coord.h"      
#include "diag.h"
#include "emode.h"
#include "scalar.h"
#include "switch.h"
#include "taskrows.h"
#ifdef xlandmix_eta
# include "xland.h"
#endif
c
      do m=1,4 
        etaana(m) = c0
      enddo
      volume_tot = c0
      salt_tot   = c0
      do jrow=jstask,jetask
        do i=1,imt
          umsk(i,jrow)       = min(kmu(i,jrow),1)
          emsk(i,jrow)       = min(kmt(i,jrow),1)
          depthu(i,jrow)     = h(i,jrow)
          depthur(i,jrow)    = hr(i,jrow)
# if defined snapshots || defined time_averages          
          psiu(i,jrow)       = c0
          psiv(i,jrow)       = c0
# endif
# if defined time_step_monitor && defined explicit_free_surface_nonlinear
          surftracer(i,jrow,:,:) = c0
# endif
          convU(i,jrow)      = c0
          rhosurf(i,jrow)    = rho0*emsk(i,jrow)
          etau(i,jrow)       = c0
          thicktau(i,jrow)   = dzt(1)
# if defined xlandmix_eta
          source_eta(i,jrow) = c0
# endif   
        enddo
      enddo
      if(init) then
        do jrow=jstask,jetask
          do i=1,imt       
            etat(i,jrow,:) = c0
            ps(i,jrow,:)   = c0
          enddo
        enddo
      endif
c
#endif  
      return
      end

      subroutine calc_convU (joff, js, je)
#if defined explicit_free_surface
c
c=======================================================================
c     Compute surface vertical velocity centered on bottom of T-cell
c     Also diagnose "streamfunctions" for vertically integrated flow. 
c=======================================================================
c
      use arrays_2d
#include "param.h"
#include "emode.h"
#include "coord.h"
#include "grdvar.h"
#include "switch.h"
#include "taskrows.h"

c
      do j=js,je
        jrow = j + joff
        jrowm1=max(1,jrow-1)
        do i=2,imtm1
          factor = cstr(jrow)*dxt2r(i)*dytr(jrow)
          convU(i,jrow) = -factor
     &         *(
     &         dyu(jrow)*(depthu(i,jrow)*ubar(i,jrow,1)
     &         -depthu(i-1,jrow)*ubar(i-1,jrow,1))
     &         + dyu(jrowm1)*(depthu(i,jrowm1)*ubar(i,jrowm1,1)
     &         -depthu(i-1,jrowm1)*ubar(i-1,jrowm1,1))
     &         + dxu(i)*(csu(jrow)*depthu(i,jrow)*ubar(i,jrow,2)
     &         -csu(jrowm1)*depthu(i,jrowm1)*ubar(i,jrowm1,2))
     &         + dxu(i-1)*(csu(jrow)*depthu(i-1,jrow)*ubar(i-1,jrow,2)
     &         -csu(jrowm1)*depthu(i-1,jrowm1)*ubar(i-1,jrowm1,2))
     &         )       
        enddo
#  ifdef cyclic
        convU(  1,jrow) = convU(imtm1,jrow)
        convU(imt,jrow) = convU(    2,jrow)
#  else 
        convU(  1,jrow) = c0
        convU(imt,jrow) = c0
#  endif
      enddo
c
#  if defined snapshots || defined time_averages       
      if (snapts .or. first .or. timavgperts) then
c         
c       Two "streamfunctions" (Sv) for depthu*ubar 
c
        do i=1,imt
          do j=js,je
            jrow = j + joff
            if (jrow .EQ. 1) then 
              psiu(i,jrow) = c0
            else
              psiu(i,jrow) = psiu(i,jrow-1)
     &                     -dyu(jrow)*ubar(i,jrow,1)*depthu(i,jrow)
     &                     *umsk(i,jrow)*1.e-12 
           endif 
          enddo
        enddo

        do j=js,je
          jrow = j + joff
          psiv(1,jrow)  = psiu(1,jrow)  
          do i = 2,imt          
            psiv(i,jrow)= psiv(i-1,jrow)
     &                   +dxu(i)*csu(jrow)*ubar(i,jrow,2)*depthu(i,jrow)
     &                   *umsk(i,jrow)*1.e-12
          enddo
        enddo
c
#   if defined explicit_psi_normalize
c
c       Normalize by subtracting value defined over America
c       This code will not work on parallel machines. 
c
        if(j+joff .eq. jmt-1) then
          normi = indp (300.0, xu, imt) 
          normj = indp (-20.0, yu, jmt) 
          psiunorm = psiu(normi,normj)
          psivnorm = psiv(normi,normj)
          do jrow=1,jmt
            do i=1,imt
              psiu(i,jrow) = (psiu(i,jrow)-psiunorm)*umsk(i,jrow)
              psiv(i,jrow) = (psiv(i,jrow)-psivnorm)*umsk(i,jrow)
            enddo
          enddo
        endif 
#   endif

      endif 
#  endif

#endif  
      return
      end


      subroutine calc_thickness (joff, js, je)
# if defined explicit_free_surface_nonlinear && defined partial_cell 
c
c=======================================================================
c     Update top grid cell thickness as well as ocean depth and volume
c     k=1 throughout this routine.  Uses partial cells.
c=======================================================================
c
      use mpp_io_mod            
      use arrays_2d
#include "param.h"
#include "coord.h"
#include "mw.h"
#include "diag.h"
#ifdef isoneutralmix
# include "isopyc.h"
#endif
#include "grdvar.h"
#include "scalar.h"
#include "switch.h"
#include "taskrows.h"
#ifdef bbl_ag
# include "bbl.h"
#endif
c
      k=1
c
c     thickness of upper t-cell at baroclinic time tau
c
      do j=js,je+1
        jrow = j + joff
        minj = min(jrow,jmt) 
        do i=1,imt
          dht(i,k,j) = max(epsln,dzt(k)+etat(i,minj,1))
        enddo
# ifdef cyclic
        dht(1,k,j)   = dht(imt-1,k,j)
        dht(imt,k,j) = dht(2,k,j)
# endif   
      enddo
c       
c     thickness of upper u-cell, isoneutral geometric factor,
c     and depth of u-cell ocean.
c
      do j=js,je
        jrow = max(1,j + joff)
# ifdef isoneutralmix
        do i=1,imt
          delqc(i,k,j,0) = c2*fracdz(k,0)*dht(i,k,j)
          delqc(i,k,j,1) = c2*fracdz(k,1)*dht(i,k,j)
        enddo
# endif   
        do i=2,imt-1
          dhu(i,k,j)     = min(dht(i,k,j), dht(i+1,k,j), dht(i,k,j+1)
     &,                        dht(i+1,k,j+1))
          etau(i,jrow)   = (dhu(i,k,j)-dzt(k))*umsk(i,jrow)  
          depthu(i,jrow) = h(i,jrow) + etau(i,jrow)
          if(umsk(i,jrow) .eq. c1) then  
            depthur(i,jrow) = c1/depthu(i,jrow)
          else
            depthur(i,jrow) = c0
          endif
        enddo
# ifdef cyclic
        dhu(1,k,j)        = dhu(imt-1,k,j)
        dhu(imt,k,j)      = dhu(2,k,j)
        etau(1,jrow)      = etau(imt-1,jrow)
        etau(imt,jrow)    = etau(2,jrow)
        depthu(1,jrow)    = depthu(imt-1,jrow)
        depthu(imt,jrow)  = depthu(2,jrow)
        depthur(1,jrow)   = depthur(imt-1,jrow)
        depthur(imt,jrow) = depthur(2,jrow)
# else    
        dhu(1,k,j)        = dht(1,k,j)
        dhu(imt,k,j)      = dht(imt,k,j)
# endif   
      enddo
c
c     time dependent volume of the ocean for T and U cells
c
      do j=js,je      
        if(j+joff .eq. jecomp) then

          do jrow = 1,jmt
            tcellv_j(jrow) = c0
            ucellv_j(jrow) = c0
          enddo
          do jrow=jscomp,jecomp
            tcellv_tmp = c0 
            ucellv_tmp = c0
            do i=2,imtm1
              if (kmt(i,jrow) .gt. 0) then
                areat_tmp  = cst(jrow)*dxt(i)*dyt(jrow) 
                tcellv_tmp = tcellv_tmp
     &             + areat_tmp*(etat(i,jrow,1)+htp(i,jrow))
# ifdef bbl_ag
                tcellv_tmp = tcellv_tmp + areat_tmp*dzt(km)
#  ifdef bbl_mask
     &           *mask_t_bbl(i,jrow)
#  endif
# endif       
              endif
              if (kmu(i,jrow) .gt. 0) then
              areau_tmp  = csu(jrow)*dxu(i)*dyu(jrow)
              ucellv_tmp = ucellv_tmp + areau_tmp*depthu(i,jrow)
# ifdef bbl_ag 
              ucellv_tmp = ucellv_tmp + areau_tmp*dzt(km)
#  ifdef bbl_mask
     &           *mask_u_bbl(i,jrow)
#  endif
# endif       
              endif
            enddo
            tcellv_j(jrow) = tcellv_tmp
            ucellv_j(jrow) = ucellv_tmp
         enddo

       endif
      enddo
#endif
      return
      end


#if defined explicit_free_surface && defined time_step_monitor
      subroutine diag_fs
c
c=======================================================================
c     Compute free surface scalar diagnostics
c=======================================================================
c
      use mpp_io_mod            
      use arrays_2d
#include "param.h"
#include "coord.h"
#include "diag.h"
#include "emode.h"
#include "grdvar.h"
#include "scalar.h"
#include "taskrows.h"
c
      negpoints = 0
      do m=1,4
        etaana(m) = c0
      enddo
c
      do jrow=jscomp,jecomp
        fu = c1 
# ifdef symmetry
        if (jrow.eq.jmtm1) fu = p5
# endif 
# if defined cyclic
        do i=2,imtm1
# else    
        do i=1,imt
# endif     
          t_area = dxt(i)*cst(jrow)*dyt(jrow)*emsk(i,jrow)
c       
c         energetics as well as fresh water volume fluxes
c
# ifdef explicit_fresh_water_flux
          etaana(1) = etaana(1) + t_area*sfft(i,jrow)*1.e-12
# endif     
          etaana(2) = etaana(2) + t_area*etat(i,jrow,1)
          etaana(4) = etaana(4) + t_area*etat(i,jrow,1)**2
          etaana(3) = etaana(3) + fu*dxu(i)*csu(jrow)*dyu(jrow)
     &         *(
     &         (depthu(i,jrow)*ubar(i,jrow,1))**2
     &         +(depthu(i,jrow)*ubar(i,jrow,2))**2
     &         )
c           
# if defined explicit_free_surface_nonlinear
c           
c         flag negative top cell thicknesses
c           
          thickness = dzt(1)+etat(i,jrow,1) 
          if(thickness .le. c0 .and. negpoints .le. 50) then
            negpoints = negpoints+1
            if(negpoints .eq. 1) then
              write(stdout,*)'Error: have negative cell thickness!!!' 
            endif
            write(stdout,*)'dzt(1)+etat(',i,',',jrow,') = ',thickness
          endif
          if(negpoints .eq. 50) then
            write(stdout,*)'More than 50 negative thickness cells!!!'
            write(stdout,*)'Model now aborting from diag_fs'
            call abort()
          endif 
# endif
          
        enddo
      enddo

# if defined explicit_free_surface_nonlinear     
c           
c     now that have the taup1 value for eta,
c     recompute surface level contribution to tbar
c
      do n=1,nt 
        do jrow=jscomp,jecomp
          cosdyt  = cst(jrow)*dyt(jrow)         
          tbar(1,n,jrow) = c0
          do i=2,imt-1
            darea_t = dxt(i)*cosdyt*emsk(i,jrow)
            thick1  = dzt(1) + etat(i,jrow,1)
            thick2  = thicktau(i,jrow)
            thk_tracer= thick1*surftracer(i,jrow,n,2)
     &                 +thick2*surftracer(i,jrow,n,1)
            tbar(1,n,jrow) = tbar(1,n,jrow) + p5*darea_t*thk_tracer
          enddo
        enddo
      enddo
# endif
      
      etaana(4) = etaana(4)*p5*grav*rho0
      etaana(3) = etaana(3)*p5*rho0

      return
      end
#endif

