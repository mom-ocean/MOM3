#ifdef kppvmix
      subroutine kppmix (joff, js, je, istrt, iend)
c
c-----------------------------------------------------------------------
c
c     Main driver subroutine for kpp vertical mixing scheme and 
c     interface to greater ocean model
c     eqn. numbers are referenced to :
c     Large, McWilliams, & Doney ( Rev. of Geophys.,32,363_403, 1994).
c
c
c     written  by: bill large,    june  6, 1994
c     modified by: jan morzel,    june 30, 1994
c                  bill large,  august 11, 1994
c
c-----------------------------------------------------------------------
c     
      use arrays_2d
      use memory_window
      use sbc_info
      use sbc_ocean
# include "param.h"
# include "coord.h"
# include "cprnts.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# include "tmngr.h"
# include "vmixc.h"
c
c
c
c model
c      integer km,kmp1        ! number of vertical levels
c      integer imt            ! number of horizontal grid points
c      real zt(1:km)          ! vertical grid                   (cm)
c      real dzt(1:km)         ! layer thicknesses               (cm)
c input
c     cori                   ! Coriolis parameter              (1/s)
c      integer jwtype(imt,jmw)    ! Jerlov water type        (1 to 5)
c                              (only used ifdef "sdflxsw")
c
c output
c      real ghats(imt,km,jmw)     ! nonlocal transport         (s/cm^2)
c      real hbl(imt,jmw)          ! boundary layer depth (cm)
c
c local
c      real bfsfc(imt,jmw)      ! surface buoyancy forcing    (cm^2/s^3)
c      real ws(imt,jmw)         ! scalar velocity scale
c      real wm(imt,jmw)         ! momentum velocity scale 
c      real caseA(imt,jmw)      ! = 1 in case A; =0 in case B
c      real stable(imt,jmw)     ! = 1 in stable forcing; =0 in unstable
c      real dkm1(imt,jmw,3)     ! boundary layer diff_cbt at kbl-1 level
c      real blmc(imt,km,jmw,3)  ! boundary layer mixing coefficients
c      real sigma(imt,jmw)      ! normalized depth (d / hbl)
c      real rhosfc(imt,jmw)     ! potential density of sfc layer(g/cm^3)
c      real talpha(imt,km,jmw)  ! d(rho)/ d(pot.temperature)  (g/cm^3/C)
c      real sbeta(imt,km,jmw)   ! d(rho)/ d(salinity)       (g/cm^3/PSU)
c      real alphaDT(imt,km,jmw) ! alpha * DT  across interfaces (g/cm^3)
c      real betaDS(imt,km,jmw)  ! beta  * DS  across interfaces (g/cm^3)
c      real ustar(imt,jmw)      ! surface friction velocity       (cm/s)
c      real Bo(imt,jmw)         ! surface turb buoy. forcing  (cm^2/s^3)
c      real Bosol(imt,jmw)      ! radiative buoy forcing      (cm^2/s^3)
c      real dbloc(imt,km,jmw)   ! local delta buoy at interfaces(cm/s^2)
c      real dbsfc(imt,km,jmw)   ! delta buoy w/ respect to sfc  (cm/s^2)
c      integer kbl(imt,jmw)     ! index of first grid level below hbl
c
c     dbsfchmxl         = critical dbsfc value that is used to
c                         compute mixed layer depth = 0.03    (cm/s2)
      data dbsfchmxl / 0.03 /
c
# ifdef timing
      call tic ('vmixc', 'kppmix')
# endif
c
c-----------------------------------------------------------------------
c     compute gradient Ri within memory window
c-----------------------------------------------------------------------
c
      call rinum (joff, js, je, istrt, iend)
c
c-----------------------------------------------------------------------
c     compute vertical difference of velocity squared
c-----------------------------------------------------------------------
c
      do j=js,je
	do k=1,kbot
	  do i=istrt-1,iend+1
	    dVsq(i,k,j) = (u(i,1,j,1,tau) - u(i,k,j,1,tau))**2 
     &                  + (u(i,1,j,2,tau) - u(i,k,j,2,tau))**2
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     set starting and ending indices for computation of mixing
c     coefficients diff_cbt and visc_cbu
c-----------------------------------------------------------------------
c
      jstrt = max(js-1,1+jbuf)
      jend  = je-1
c
c-----------------------------------------------------------------------
c
c     density related quantities
c     --------------------------
c     based on mom equation of state, which computes normalized
c     density  drho{k,kr} for layer k with respect to reference layer
c     kr, i.e. drho{T(k),S(k), tr(kr), sr(kr), zt(kr)}, as the
c     difference of drho=rho{k,kr}-rhor(kr):
c
c     density of surface layer                                   (g/cm3)
c           rho   = drho{T(1),S(1),zt(1)} + 1. + rhor(zt(1))
c     local buoyancy gradient at km interfaces:                  (cm/s2)
c           dbloc = g/rho{k+1,k+1} * [ drho{k,k+1}-drho{k+1,k+1} ]
c     buoyancy difference with respect to "zref",i.e. the surface(cm/s2)
c           dbsfc = g * [ drho{1,k}/rho{1,k} - drho{k,k}/rho{k,k} ]
c     thermal expansion coefficient without 1/rho factor       (g/cm3/C)
c           talpha= d(rho{k,k})/d(T(k))
c     salt expansion coefficient without 1/rho factor        (g/cm3/PSU)
c           sbeta = d(rho{k,k})/d(S(k))
c-----------------------------------------------------------------------
c
      call statekpp (joff, jstrt, jend, istrt, iend)
c
c-----------------------------------------------------------------------
c     kinematic surface fluxes on the "t-grid" :
c     --------------------------------------------------------
c
c     wusfc = kinematic zonal velocity sfc flux on jp           (cm2/s2)
c           = -smf(i,j,1)
c     wvsfc = kinematic merid.velocity ..                       (cm2/s2)
c           = -smf(i,j,2)
c     wtsfc = kinematic temperature    ..                       (C*cm/s)
c           = -stf(i,j,1)
c     wssfc = kinematic salinity       ..                     (PSU*cm/s)
c           = -stf(i,j,2)
c     wtsol = kinematic solar temp.    ..                       (C*cm/s)
c           = -sbcocn(i,jrow,ishortw) 
c
c     friction velocity, turbulent and radiative sfc buoyancy forcing
c     ---------------------------------------------------------------
c
c
c     ustar(i) = sqrt(sqrt(wusfc(i)**2 + wvsfc(i)**2))           (cm/s)
c     Bo(i)    = -g*(talpha(i,j)*wtsfc(i)-sbeta(i,j)*wssfc(i))/rho(i)
c     Bosol(i) = -g*(talpha(i,j)*wtsol(i)                    )/rho(i)
c                                                              (cm2/s3)
c-----------------------------------------------------------------------
c
# ifdef shortwave
      ishortw = index_of_sbc ('short wave')
      if (ishortw .eq. 0) then
        write (stdout,*)'=>Error: short wave S.B.C. is not setup'
	call abort()
      endif
# endif
      do j=jstrt,jend
	jrow = j + joff
	do i=istrt,iend
c
c         the ustar field is needed on the "T-grid" 
c
          active_cells = umask(i,1,j) + umask(i-1,1,j) +
     &         umask(i,1,j-1) + umask(i-1,1,j-1)+epsln
          smftu = (smf(i,j,1)   + smf(i-1,j,1)
     &         + smf(i,j-1,1) + smf(i-1,j-1,1))
     &         /active_cells
          smftv = (smf(i,j,2)   + smf(i-1,j,2)
     &         + smf(i,j-1,2) + smf(i-1,j-1,2))
     &         /active_cells
          ustar(i,j) = sqrt(sqrt( smftu**2 + smftv**2 ))
          Bo(i,j)    = grav * (
     &         talpha(i,1,j) * stf(i,j,1)
     &         -sbeta (i,1,j) * stf(i,j,2)
     &         )
     &                   / (epsln + rhosfc(i,j))
# ifdef shortwave
          Bosol(i,j) = grav * talpha(i,1,j) * sbcocn(i,jrow,ishortw)
     &                 / (epsln + rhosfc(i,j))
# endif
	enddo
        call setbcx (ustar(1,j), imt, 1)
        call setbcx (Bo(1,j)   , imt, 1)
# ifdef shortwave
        call setbcx (Bosol(1,j), imt, 1)
# endif
      enddo
c
c-----------------------------------------------------------------------
c     compute interior mixing coefficients everywhere, due to constant 
c     internal wave activity, static instability, and local shear 
c     instability.
c-----------------------------------------------------------------------
c
      call ri_iwmix (joff, jstrt, jend, istrt, iend)  
c
c-----------------------------------------------------------------------
c     add double diffusion if desired
c-----------------------------------------------------------------------
c
      if (ldd) then
        call ddmix (joff, jstrt, jend, istrt, iend)  
      endif
c 
c-----------------------------------------------------------------------
c     boundary layer mixing coefficients: diagnose new b.l. depth
c-----------------------------------------------------------------------
c
      call bldepth (joff, jstrt, jend, istrt, iend)
c 
c-----------------------------------------------------------------------
c     boundary layer diffusivities
c-----------------------------------------------------------------------
c
      call blmix_kpp (joff, jstrt, jend, istrt, iend) 
c 
c-----------------------------------------------------------------------
c     enhance diffusivity at interface kbl - 1
c-----------------------------------------------------------------------
c
      call enhance  (joff, jstrt, jend, istrt, iend)
c 
c-----------------------------------------------------------------------
c     combine interior and b.l. coefficients and nonlocal term
c-----------------------------------------------------------------------
c
      do j=jstrt,jend
        jrow = j + joff
        do ki=1,kbot-1
          do i=istrt,iend
            if (ki .lt. kbl(i,j)) then
              visc_cbu(i,ki,j)   = blmc(i,ki,j,1)
              diff_cbt(i,ki,j,1) = blmc(i,ki,j,2)
              diff_cbt(i,ki,j,2) = blmc(i,ki,j,3)
            else
              ghats(i,ki,j)=0.0
            endif
	    visc_cbu(i,ki,j) = visc_cbu(i,ki,j)
     &           *umask(i,min(ki+1,kbot),j)
	    diff_cbt(i,ki,j,1) = diff_cbt(i,ki,j,1)
     &           *tmask(i,min(ki+1,kbot),j)
	    diff_cbt(i,ki,j,2) = diff_cbt(i,ki,j,2)
     &           *tmask(i,min(ki+1,kbot),j)
          enddo
        enddo
        do i=istrt,iend
	  visc_cbu(i,kbot,j) = c0
	  diff_cbt(i,kbot,j,1) = c0
	  diff_cbt(i,kbot,j,2) = c0
#ifdef bbl_ag
	  visc_cbu(i,km,j) = c0
	  diff_cbt(i,km,j,1) = c0
	  diff_cbt(i,km,j,2) = c0
#endif
        enddo
	call setbcx (visc_cbu(1,1,j), imt, km)
	call setbcx (diff_cbt(1,1,j,1), imt, km)
	call setbcx (diff_cbt(1,1,j,2), imt, km)
      enddo
c
      if (snapts.or.prxzts) then
c-----------------------------------------------------------------------
c       compute mixed layer depth (cm), which is defined as the depth (>0)
c       where the buoyancy difference with respect to the surface level is
c       equal to 0.03 (cm/s2).
c       use dbsfc, which is defined on zt levels
c       hxml is used for diagnostics only.            
c-----------------------------------------------------------------------
c       

        do j=jstrt,jend
          jrow = j + joff
          do i=istrt,iend
            if (kmt(i,jrow).eq.0) then
              hmxl(i,jrow) = 0.0
            else
              hmxl(i,jrow) = zt(kmt(i,jrow))
            end if
          enddo
        enddo
        do j=jstrt,jend
          jrow = j + joff
          do k=2,kbot
            km1 = k-1  
            do i=istrt,iend
              if (kmt(i,jrow).eq.0) then
                hmxl(i,jrow) = 0.0
              else
                if ( dbsfc(i,k,j) .ge. dbsfchmxl 
     &               .and. hmxl(i,jrow) .eq. zt(kmt(i,jrow)) ) then
                  hmxl(i,jrow) = zt(km1) - (zt(km1) - zt(k))
     &                 * (dbsfchmxl-dbsfc(i,km1,j))
     &                 / (dbsfc(i,k,j) - dbsfc(i,km1,j) + epsln)
                endif
                hmxl(i,jrow) = hmxl(i,jrow) * tmask(i,1,j)
              endif
            enddo
          enddo
        enddo
        
      end if
c      
c-----------------------------------------------------------------------
c      boundary layer depth
c-----------------------------------------------------------------------
c
       do j=jstrt,jend
         jrow = j + joff
         do i=istrt,iend
	   hblt(i,jrow) = hbl(i,j) * tmask(i,1,j)
         enddo
       enddo
c
c-----------------------------------------------------------------------
c      when using kmix vertical mixing scheme, save global array
c      of boundary layer depth (cm) and mixed-layer depth (cm).
c-----------------------------------------------------------------------
c

# ifdef matrix_sections
      if (prxzts .and. eots) then
c
        call diagvmix (joff, jstrt, jend)
c
        do j=jstrt,jend
          jrow = j + joff
	  if (jrow .eq. jmt-1) then
	    isml = indp (slonxy, xt, imt)
	    ieml = indp (elonxy, xt, imt)
	    jsml = indp (slatxy, yt, jmt)
	    jeml = indp (elatxy, yt, jmt)
c
            scl = 100.
            write (stdout,7100) 'Mixed Layer Depth (m)'
     &,     itt, xt(isml), xt(ieml), yt(jsml), yt(jeml), scl
c
	    call matrix (hmxl(1,1), imt, isml, ieml, -jsml, -jeml, scl)
c
7100       format(1x,a30,1x,'ts=',i7
     &,    ', lon:',f6.2,' ==> ',f6.2,',   lat:',f6.2,' ==> ',f6.2
     &,    ', scaling=',1pg10.3)
	  endif
        enddo
      endif
# endif
# ifdef trace_indices
      write (*,'(2x,7(a,i4))')
     & "=> In kppmix: js=",js," je=",je," joff=",joff
     &," jstrt=",jstrt," jend=",jend," jrow=",jstrt+joff
     &," to ",jend+joff
# endif
# ifdef timing
      call toc ('vmixc', 'kppmix')
# endif

      return
      end


      subroutine  bldepth (joff, js, je, is, ie)
c
c=======================================================================
c
c
c     the oceanic planetray boundary layer depth, hbl, is determined as
c     the shallowest depth where the bulk richardson number is
c     equal to the critical value, Ricr.
c
c     bulk richardson numbers are evaluated by computing velocity and
c     buoyancy differences between values at zt(kl) and surface
c     reference values.
c     in this configuration, the reference values are equal to the
c     values in the surface layer.  
c     when using a very fine vertical grid, these values should be 
c     computed as the vertical average of velocity and buoyancy from 
c     the surface down to epsilon*zt(kl).
c
c     when the bulk richardson number at k exceeds Ricr, hbl is
c     linearly interpolated between grid levels zt(k) and zt(k-1).
c
c     The water column and the surface forcing are diagnosed for 
c     stable/ustable forcing conditions, and where hbl is relative 
c     to grid points (caseA), so that conditional branches can be 
c     avoided in later subroutines.
c
      use arrays_2d
      use memory_window
# include "param.h"
# include "coord.h"
# ifdef shortwave
#  include "cshort.h"
# endif
# include "mw.h"
# include "vmixc.h"

c  model  
c      integer km,kmp1        ! number of vertical levels
c      integer imt            ! number of horizontal grid points
c      real zt(1:km)          ! vertical grid                   (cm)
c      real dzt(1:km)         ! layer thicknesses               (cm)
c
c  input
c      real dbloc(imt,km,jmw)   ! local delta buoyancy         (cm/s^2)
c      real dbsfc(imt,km,jmw) ! delta buoyancy w/ respect to sfc(cm/s)^2
c      real ustar(imt,jmw)      ! surface friction velocity     (cm/s)
c      real Bo(imt,jmw)   ! surface turbulent buoyancy forcing(cm^2/s^3)
c      real Bosol(imt,jmw)      ! radiative buoyancy forcing (cm^2/s^3)
c      real cori                ! Coriolis parameter            (1/s)
c      integer jwtype(imt,jmw)  ! Jerlov water type           (1 to 5)
c
c  output
c      real hbl(imt,jmw)        ! boundary layer depth              (cm)
c      real bfsfc(imt,jmw) !Bo+radiation absorbed to d=hbf*hbl(cm^2/s^3)
c      real stable(imt,jmw)     ! =1 in stable forcing; =0 unstable
c      real caseA(imt,jmw)      ! =1 in case A, =0 in case B 
c      integer kbl(imt,jmw)     ! index of first grid level below hbl 
c 
c  local
c      real dVsq(imt,km,jmw)    ! (velocity shear re sfc)^2   (cm/s)^2
      real Ritop               ! numerator of bulk Richardson Number
c      real Rib(imt,2,jmw)      ! Bulk Richardson number
c      real sigma(imt,jmw)      ! normalized depth (d/hbl)
c      real wm(imt,jmw),ws(imt,jmw) ! turbulent velocity scales  (cm/s)
c
c=======================================================================
c
c
c find bulk Richardson number at every grid level until > Ricr
c
c note: the reference depth is -epsilon/2.*zt(k), but the reference
c       u,v,t,s values are simply the surface layer values,
c       and not the averaged values from 0 to 2*ref.depth,
c       which is necessary for very fine grids(top layer < 2m thickness)
c note: max values when Ricr never satisfied are
c       kbl(i)=kmt(i,jrow) and hbl(i)=zt(i,kmt(i,jrow))
c
c
c-----------------------------------------------------------------------
c     set local constants
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
c     initialize hbl and kbl to bottomed out values
c
        do j=js,je
	  jrow   = j + joff
	  do i=istrt,iend
            Rib(i,1,j) = 0.0
            kbl(i,j)    = MAX(kmt(i,jrow),2)
            hbl(i,j)    = zt(kbl(i,j))
          enddo
        enddo
c
      do j=js,je
        jrow   = j + joff
c
c       indices for array Rib(i,k,j), the bulk Richardson number.
c       
        ka = 1
        ku = 2
c
	do kl=2,kbot
c
c         compute bfsfc = sw fraction at hbf * zt
c
# ifdef shortwave
c         This is commented out because water type is a constant in MOM:
c         call swfrac(imt,1.0,zt(kl),jwtype,bfsfc)
# endif

          do i=istrt,iend
c
c           compute bfsfc= Bo + radiative contribution down to hbf * hbl
c
            bfsfc(i,j)  = Bo(i,j) 
# ifdef shortwave
     &                  + Bosol(i,j) * (pen(0) - pen(kl-1))
# endif
            stable(i,j) = 0.5 + SIGN( 0.5, bfsfc(i,j) )
            sigma(i,j)  = stable(i,j) * 1. + (1.-stable(i,j)) * epsilon
          enddo
c
c-----------------------------------------------------------------------
c         compute velocity scales at sigma, for hbl = zt(kl):
c-----------------------------------------------------------------------
c
          iwscale_use_hbl_eq_zt=1
          call wscale (joff, j, j, istrt, iend
     &,                iwscale_use_hbl_eq_zt, zt(kl))


	  do i=istrt,iend
c
c-----------------------------------------------------------------------
c           compute the turbulent shear contribution to Rib
c           eqn. (23)
c-----------------------------------------------------------------------
c
            bvsq =0.5*
     &            ( dbloc(i,kl-1,j) / (dzw(kl-1))+
     &              dbloc(i,kl  ,j) / (dzw(kl)) )
            Vtsq =  zt(kl) * ws(i,j) * sqrt(abs(bvsq)) * Vtc
c 
c-----------------------------------------------------------------------
c           compute bulk Richardson number at new level
c           note: Ritop needs to be zero on land and ocean bottom
c           points so that the following if statement gets triggered
c           correctly. otherwise, hbl might get set to (big) negative
c           values, that might exceed the limit for the "exp" function
c           in "swfrac"
c
c           eqn. (21)
c
c
c           numerator of bulk richardson number on grid levels
c           --------------------------------------------------
c           note: land and ocean bottom values need to be set to zero
c           so that the subroutine "bldepth" works correctly
c
c-----------------------------------------------------------------------
c
            Ritop = (zt(kl)-zt(1)) * dbsfc(i,kl,j) * tmask(i,kl,j)
	    active_cells = umask(i,kl,j) + umask(i-1,kl,j) +
     &                     umask(i,kl,j-1) + umask(i-1,kl,j-1)+epsln
            dVsq_avg_tcell=
     &           (dVsq(i,kl,j) + dVsq(i-1,kl,j)
     &           + dVsq(i,kl,j-1) + dVsq(i-1,kl,j-1))
     &           /active_cells
            Rib(i,ku,j) = Ritop / ( dVsq_avg_tcell + Vtsq + epsln )

            if((kbl(i,j).eq.kmt(i,jrow)).and.(Rib(i,ku,j).gt.Ricr)) then
c
c              linearly interpolate to find hbl where Rib = Ricr
c
               hbl(i,j) = zt(kl-1) + (dzw(kl-1)) *
     &                  (Ricr - Rib(i,ka,j))
     &                 /(Rib(i,ku,j)-Rib(i,ka,j) + epsln)
               kbl(i,j) = kl
            endif
         enddo
c
         ksave = ka
         ka    = ku
         ku    = ksave
c
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     find stability and buoyancy forcing for boundary layer
c-----------------------------------------------------------------------
c
# ifdef shortwave
c      call swfrac(imt,1.0,hbl,jwtype,bfsfc)
# endif
c
      do j=js,je
        do i = istrt,iend
          bfsfc(i,j)  = Bo(i,j) 
# ifdef shortwave
     &                + Bosol(i,j) * (pen(0) - pen(kbl(i,j)))
# endif
          stable(i,j) = 0.5 + SIGN( 0.5, bfsfc(i,j) )
c         ensures bfsfc never=0
          bfsfc(i,j)  = bfsfc(i,j) + stable(i,j) * epsln
        enddo
      enddo
c
# ifndef kpp_dont_check_ekman_moninobukov
c
c-----------------------------------------------------------------------
c        check hbl limits for hekman or hmonob
c        eqn. (24)
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow   = j + joff
        do i = istrt,iend
          if (bfsfc(i,j) .gt. 0.0) then
             hekman = cekman * ustar(i,j) / (abs(cori(jrow,1))+epsln)
             hmonob = cmonob * ustar(i,j)*ustar(i,j)*ustar(i,j) 
     &                /vonk / (bfsfc(i,j)+epsln) 
             hlimit = stable(i,j)     * AMIN1(hekman,hmonob) +
     &                (stable(i,j)-1.) * (zt(kbot))
             hbl(i,j) = AMIN1(hbl(i,j),hlimit)
             hbl(i,j) = AMAX1(hbl(i,j),zt(1))
          endif
          kbl(i,j) = kmt(i,jrow)
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     find new kbl
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow   = j + joff
	do kl=2,kbot
          do i = istrt,iend
            if((kbl(i,j).eq.kmt(i,jrow)).and.(zt(kl).gt.hbl(i,j))) then
              kbl(i,j) = kl
            endif
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     find stability and buoyancy forcing for final hbl values
c-----------------------------------------------------------------------
c
#  ifdef shortwave
c     call swfrac(imt,1.0,hbl,jwtype,bfsfc)
#  endif
      do j=js,je
        do i = istrt,iend
          bfsfc(i,j)  = Bo(i,j) 
#  ifdef shortwave
     &                + Bosol(i,j) * (pen(0) - pen(kbl(i,j)))
#  endif
          stable(i,j) = 0.5 + SIGN( 0.5, bfsfc(i,j) )
          bfsfc(i,j)  = bfsfc(i,j) + stable(i,j) * epsln 
        enddo
      enddo
# endif
c
c-----------------------------------------------------------------------
c     determine caseA and caseB
c     (if hbl is below (deeper than) the mid point of level kbl
c     then caseA=0  else caseA=1)
c-----------------------------------------------------------------------
c
      do j=js,je
        do i = istrt,iend
          caseA(i,j)  = 0.5 + 
     &        SIGN( 0.5,zt(kbl(i,j)) - 0.5*dzt(kbl(i,j)) - hbl(i,j) )
        enddo
      enddo
      return
      end

c *********************************************************************

      subroutine wscale(joff, js, je, istrt, iend
     &,                 iwscale_use_hbl_eq_zt, zt_kl)
c
c     compute turbulent velocity scales.
c     use a 2D-lookup table for wm and ws as functions of ustar and
c     zetahat (=vonk*sigma*hbl*bfsfc).
c
c     note: the lookup table is only used for unstable conditions
c     (zehat.le.0), in the stable domain wm (=ws) gets computed
c     directly.
c
      use memory_window
# include "param.h"
# include "mw.h"
# include "vmixc.h"
# include "switch.h"


c  model

c  input
c      real sigma(imt,jmw)      ! normalized depth (d/hbl)
c      real hbl(imt,jmw)        ! boundary layer depth (cm)
c      real ustar(imt,jmw)      ! surface friction velocity    (cm/s)
c      real bfsfc(imt,jmw)      ! total surface buoyancy flux (cm^2/s^3)

c  output
c      real wm(imt,jmw),ws(imt,jmw) ! turbulent velocity scales at sigma

c local
      real zehat           ! = zeta *  ustar**3

c
c-----------------------------------------------------------------------
c     use lookup table for zehat < zmax only; otherwise use
c     stable formulae
c-----------------------------------------------------------------------
c
      do j=js,je
        do i=istrt,iend
          if (iwscale_use_hbl_eq_zt .eq. 1) then
            zehat = vonk * sigma(i,j) * zt_kl * bfsfc(i,j)
          else
            zehat = vonk * sigma(i,j) * hbl(i,j) * bfsfc(i,j)
          end if

         if (zehat.le.zmax) then

            zdiff  = zehat-zmin
            iz = int( zdiff/deltaz )
            iz = min( iz , nni )
            iz = max( iz , 0  )
            izp1=iz+1
            
            udiff  = ustar(i,j)-umin
            ju = int( udiff/deltau)
            ju = min( ju , nnj )
            ju = max( ju , 0  )
            jup1=ju+1

            zfrac = zdiff/deltaz - float(iz)
            ufrac = udiff/deltau - float(ju)

            fzfrac= 1.-zfrac
            wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
            wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
            wm(i,j) = (1.-ufrac)* wbm          + ufrac*wam

            was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
            wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
            ws(i,j) = (1.-ufrac)* wbs          + ufrac*was
          else
            u3    = ustar(i,j)*ustar(i,j)*ustar(i,j)
            wm(i,j) = vonk * ustar(i,j) * u3
     &           / ( u3 + conc1*zehat + epsln )
            ws(i,j) = wm(i,j)
          endif
        enddo
      enddo
      return
      end

 
      subroutine ri_iwmix (joff, js, je, is, ie)   
c
c=======================================================================
c
c     compute interior viscosity diffusivity coefficients due 
c     to shear instability (dependent on a local richardson number),
c     to background internal wave activity, and 
c     to static instability (local richardson number < 0).
c
c
c
c     inputs:
c
c      joff   = offset between rows in the MW and latitude rows
c      js     = starting row for loading variables to calculate
c               coefficients.
c      je     = ending row for loading variables to calculate
c               coefficients. calculations end at je-1
c      is     = starting index for calculating coefficients in the
c               longitude direction
c      ie     = ending index for calculating coefficients in the
c               longitude direction
c      km     = number of vertical levels
c      grav   = gravity (cm/sec**2)
c      umask  = land/sea mask on "u" grid (land=0.0, sea=1.0)
c      tmask  = land/sea mask on "t" grid (land=0.0, sea=1.0)
c      fricmx = max viscosity (cm**2/sec)
c      wndmix = min viscosity at bottom of 1st level to simulate
c               missing high frequency windstress components (cm**2/sec)
c      visc_cbu_back = background "visc_cbu" (cm**2/sec)
c      diff_cbt_back = background "diff_cbt" (cm**2/sec)
c      visc_cbu_limit = largest "visc_cbu" in regions of gravitational
c                      instability (cm**2/sec)
c      diff_cbt_limit = largest "diff_cbt" in regions of gravitational
c                       instability (cm**2/sec)
c
c
c
c     outputs:
c
c      visc_cbu = viscosity coefficient at bottom of "u" cells (cm**2/s)
c      diff_cbt = diffusion coefficient at bottom of "t" cells (cm**2/s)
c
c=======================================================================
c
      use memory_window
# include "param.h"
# include "mw.h"
# include "rinum.h"
# include "vmixc.h"
c      real dbloc(imt,km,jmw)
c
c
c
c-----------------------------------------------------------------------
c     set local constants
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
c-----------------------------------------------------------------------
c     diffusion and viscosity coefficients on bottom of "T" cells
c-----------------------------------------------------------------------
c
      do j=js,je
	do k=1,kbot-1
	  do i=istrt,iend
c  
c-----------------------------------------------------------------------
c           evaluate function of Ri for shear instability eqn. (28b&c)
c-----------------------------------------------------------------------
c
            Rigg  = AMAX1( rit(i,k,j) , 0.0)
            ratio = AMIN1( Rigg/Riinfty , 1.0 )
            frit  = (1. - ratio*ratio)
            frit  = frit*frit*frit
            Rigg  = AMAX1( riu(i,k,j) , 0.0)
            ratio = AMIN1( Rigg/Riinfty , 1.0 )
            friu  = (1. - ratio*ratio)
            friu  = friu*friu*friu
c
c-----------------------------------------------------------------------
c           evaluate function of Ri for convection eqn. (28a)
c-----------------------------------------------------------------------
c
            fcont  = 0.5 * ( abs(rit(i,k,j)) - rit(i,k,j) )
            fcont  = fcont / (AMAX1(fcont,epsln))
            fconu  = 0.5 * ( abs(riu(i,k,j)) - riu(i,k,j) )
            fconu  = fconu / (AMAX1(fconu,epsln))
c     
c-----------------------------------------------------------------------
c           mixing due to internal wave activity and static instability 
c           eqn. (29)
c-----------------------------------------------------------------------
c
            visc_cbu(i,k,j)   = visc_cbu_iw + fconu * visc_con_limit
            diff_cbt(i,k,j,1) = diff_cbt_iw + fcont * diff_con_limit
            diff_cbt(i,k,j,2) = diff_cbt_iw + fcont * diff_con_limit
c
c-----------------------------------------------------------------------
c           add contrubution due to shear instability
c-----------------------------------------------------------------------
c
            if (lri) then
              visc_cbu(i,k,j)   = visc_cbu(i,k,j)  +visc_cbu_limit*friu
              diff_cbt(i,k,j,1) = diff_cbt(i,k,j,1)+diff_cbt_limit*frit
              diff_cbt(i,k,j,2) = diff_cbt(i,k,j,2)+diff_cbt_limit*frit
            endif
#ifdef bryan_lewis_vertical
c
c-----------------------------------------------------------------------
c           add Bryan-Lewis mixing if wanted
c-----------------------------------------------------------------------
c
            diff_cbt(i,k,j,1) = diff_cbt(i,k,j,1) + Ahv(k)
            diff_cbt(i,k,j,2) = diff_cbt(i,k,j,2) + Ahv(k)
#endif
	   enddo
        enddo
	call setbcx (visc_cbu(1,1,j), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c       set lateral bc
c-----------------------------------------------------------------------
c
      do j=js,je
	call setbcx (diff_cbt(1,1,j,1), imt, km)
        call setbcx (diff_cbt(1,1,j,2), imt, km)
      enddo
c
      return
      end

c *********************************************************************

      subroutine ddmix (joff, js, je, istrt, iend)
c
c=======================================================================
c     Rrho dependent interior flux parameterization.
c     Add double-diffusion diffusivities to Ri-mix values at blending
c     interface and below.
c=======================================================================
c
      use memory_window
# include "param.h"
# include "mw.h"
# include "vmixc.h"

c
c=======================================================================
c
c     inputs:
c
c      joff   = offset between rows in the MW and latitude rows
c      js     = starting row for loading variables to calculate
c               coefficients.
c      je     = ending row for loading variables to calculate
c               coefficients. calculations end at je-1
c      is     = starting index for calculating coefficients in the
c               longitude direction
c      ie     = ending index for calculating coefficients in the
c               longitude direction
c      km     = number of vertical levels
c      real talpha(imt,km,jmw)   ! d(rho)/ d(pot.temperature) (g/cm^3/C)
c      real sbeta(imt,km,jmw)    ! d(rho)/ d(salinity)     (g/cm^3/PSU)
c
c     outputs:
c
c      visc_cbu = viscosity coefficient at bottom of "u" cells (cm**2/s)
c      diff_cbt = diffusion coefficient at bottom of "t" cells (cm**2/s)
c
c=======================================================================
c
c
c local
c      real alphaDT(imt,km,jmw)   ! alpha * DT  across interfaces
c      real betaDS(imt,km,jmw)    ! beta  * DS  across interfaces
      real diffdd                ! double diffusion diffusivity scale
      real prandtl               ! prandtl number

c
c-----------------------------------------------------------------------
c
c     for double diffusion
c     --------------------
c     temperature and salt contributions of buoyancy gradients
c     on interfaces
c
      do j=js,je
	do k=1,kbot-1
	kp1 = k+1
	  do i=istrt,iend
	    alphaDT(i,k,j) = 0.5 * (talpha(i,k,j) + talpha(i,kp1,j) )
     &                        * (t(i,k,j,1,tau) - t(i,kp1,j,1,tau))
            betaDS (i,k,j) = 0.5 * (sbeta (i,k,j) + sbeta(i,kp1,j)  )
     &                        *(t(i,k,j,2,tau) - t(i,kp1,j,2,tau))
	  enddo
        enddo
      enddo
      do j=js,je
	do i=istrt,iend
	  alphaDT(i,kbot,j) = 0.0
	  betaDS (i,kbot,j) = 0.0
#ifdef bbl_ag
	  alphaDT(i,km,j) = 0.0
	  betaDS (i,km,j) = 0.0
#endif
	enddo
      enddo
c-----------------------------------------------------------------------
      do j=js,je
	do ki=1,kbot
          do i = istrt,iend
c
c-----------------------------------------------------------------------
c           salt fingering case
c           eqn. (31)
c-----------------------------------------------------------------------
c
            if ((alphaDT(i,ki,j) .gt. betaDS(i,ki,j)) .and.
     &          (betaDS (i,ki,j) .gt. 0.         )) then

              Rrho       = MIN(alphaDT(i,ki,j) / betaDS(i,ki,j), Rrho0)
c             diffdd     = dsfmax*(1.0-((Rrho-1)/(Rrho0-1))**2)**pexp2
              diffdd     =         1.0-((Rrho-1)/(Rrho0-1))**2
              diffdd     = dsfmax*diffdd*diffdd*diffdd
	      diff_cbt(i,ki,j,1) = diff_cbt(i,ki,j,1) + 0.7*diffdd
              diff_cbt(i,ki,j,2) = diff_cbt(i,ki,j,2) + diffdd
c
c-----------------------------------------------------------------------
c             diffusive convection eqn. (32)
c-----------------------------------------------------------------------
c
            else if ((alphaDT(i,ki,j) .lt. 0.0).and.
     &                (betaDS(i,ki,j) .lt. 0.0)
     &               .and. (alphaDT(i,ki,j) .gt. betaDS(i,ki,j)) ) then
              viscosity_molecular=1.5e-2
              Rrho    = alphaDT(i,ki,j) / betaDS(i,ki,j) 
              diffdd  = viscosity_molecular
     &             *9.0*0.101*exp(4.6*exp(-0.54*(1/Rrho-1)))
c
              prandtl = 0.15*Rrho
              if (Rrho.gt.0.5) prandtl = (1.85-0.85/Rrho)*Rrho
c
	      diff_cbt(i,ki,j,1) = diff_cbt(i,ki,j,1) + diffdd
              diff_cbt(i,ki,j,2) = diff_cbt(i,ki,j,2) + prandtl*diffdd
c 
            endif
          enddo
        enddo
      enddo
      return
      end

c *********************************************************************

      subroutine blmix_kpp (joff, js, je, istrt, iend)
c
c     mixing coefficients within boundary layer depend on surface
c     forcing and the magnitude and gradient of interior mixing below
c     the boundary layer ("matching").
c
c     caution: if mixing bottoms out at hbl = zgrid(km) then
c     fictitious layer at kmp1 is needed with small but finite width 
c     hwide(kmp1) (eg. epsln = 1.e-20).
c
      use arrays_2d
      use memory_window
# include "param.h"
# include "coord.h"
# include "mw.h"
# include "vmixc.h"
# include "switch.h"

c
c=======================================================================
c
c     inputs:
c
c      joff   = offset between rows in the MW and latitude rows
c      js     = starting row for loading variables to calculate
c               coefficients.
c      je     = ending row for loading variables to calculate
c               coefficients. calculations end at je-1
c      is     = starting index for calculating coefficients in the
c               longitude direction
c      ie     = ending index for calculating coefficients in the
c               longitude direction
c      km     = number of vertical levels
c
c      real ustar(imt,jmw)    ! surface friction velocity         (cm/s)
c      real bfsfc(imt,jmw)    ! surface buoyancy forcing     (cm^2/s^3)
c      real hbl(imt,jmw)      ! boundary layer depth              (cm)
c      real stable(imt,jmw)   ! = 1 in stable forcing
c      real caseA(imt,jmw)    ! = 1 in case A
c      integer kbl(imt,jmw)   ! index of first grid level below hbl
c
c output
c
c     outputs:
c
c      visc_cbu = viscosity coefficient at bottom of "u" cells (cm**2/s)
c      diff_cbt = diffusion coefficient at bottom of "t" cells (cm**2/s)
c
c=======================================================================

c      real dkm1(imt,jmw,3)    ! boundary layer diff_cbt at kbl-1 level
c      real blmc(imt,km,jmw,3) ! boundary layer mixing coeff.(cm**2/s)
c      real ghats(imt,km,jmw)  ! nonlocal scalar transport
c
c  local
c      real gat1(imt,jmw,3)
c      real dat1(imt,jmw,3)
c      real sigma(imt,jmw)            ! normalized depth (d / hbl)
c      real ws(imt,jmw), wm(imt,jmw)  ! turbulent velocity scales (cm/s)
c
      do j=js,je
c
c-----------------------------------------------------------------------
c       compute velocity scales at hbl
c-----------------------------------------------------------------------
c
        do i = istrt,iend
          sigma(i,j) = stable(i,j) * 1.0 + (1.-stable(i,j)) * epsilon
        enddo
c
        iwscale_use_hbl_eq_zt = 0
        zt_kl_dummy=0.0
        call wscale (joff, j, j, istrt, iend
     &,             iwscale_use_hbl_eq_zt, zt_kl_dummy)
c 
        do i = istrt,iend 
          kn    = ifix(caseA(i,j)+epsln) *(kbl(i,j) -1) +
     &            (1-ifix(caseA(i,j)+epsln)) * kbl(i,j)
          knm1 = max(kn-1,1)
          knp1 = min(kn+1,kbot)
c
c-----------------------------------------------------------------------
c         find the interior viscosities and derivatives at hbl(i) 
c         eqn. (18)
c-----------------------------------------------------------------------
c
          delhat = 0.5*dzt(kn) + zt(kn) - hbl(i,j)
          R      = 1.0 - delhat / dzt(kn)
	  dvdzup = (visc_cbu(i,knm1,j) - visc_cbu(i,kn,j))/dzt(kn)
	  dvdzdn = (visc_cbu(i,kn,j) - visc_cbu(i,knp1,j))/dzt(knp1)
          viscp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     &                        R  * (dvdzdn + abs(dvdzdn)) )
c
	  dvdzup = (diff_cbt(i,knm1,j,2) - diff_cbt(i,kn,j,2))/dzt(kn)
	  dvdzdn = (diff_cbt(i,kn,j,2) - diff_cbt(i,knp1,j,2))/dzt(knp1)
          difsp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     &                        R  * (dvdzdn + abs(dvdzdn)) )
c
	  dvdzup = (diff_cbt(i,knm1,j,1) - diff_cbt(i,kn,j,1))/dzt(kn)
	  dvdzdn = (diff_cbt(i,kn,j,1) - diff_cbt(i,knp1,j,1))/dzt(knp1)
          diftp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     &                        R  * (dvdzdn + abs(dvdzdn)) )
c
          visch  = visc_cbu(i,kn,j)   + viscp * delhat
          difsh  = diff_cbt(i,kn,j,2) + difsp * delhat
          difth  = diff_cbt(i,kn,j,1) + diftp * delhat
c
          f1 = stable(i,j) * conc1 * bfsfc(i,j) / (ustar(i,j)**4+epsln) 
c
          gat1(i,j,1) = visch / (hbl(i,j)+epsln) / (wm(i,j)+epsln)
          dat1(i,j,1) = -viscp / (wm(i,j)+epsln) + f1 * visch
          dat1(i,j,1) = min(dat1(i,j,1),0.) 
c
          gat1(i,j,2) = difsh  / (hbl(i,j)+epsln) / (ws(i,j)+epsln)
          dat1(i,j,2) = -difsp / (ws(i,j)+epsln) + f1 * difsh 
          dat1(i,j,2) = min(dat1(i,j,2),0.) 
c  
          gat1(i,j,3) = difth /  (hbl(i,j)+epsln) / (ws(i,j)+epsln)
          dat1(i,j,3) = -diftp / (ws(i,j)+epsln) + f1 * difth 
          dat1(i,j,3) = min(dat1(i,j,3),0.) 
        enddo
      enddo
c
      blmc = 0.0
      do j=js,je
	do ki=1,kbot
c
c-----------------------------------------------------------------------
c         compute turbulent velocity scales on the interfaces
c-----------------------------------------------------------------------
c
          do i  = istrt,iend
            sig     = (zt(ki) + 0.5 * dzt(ki)) / (hbl(i,j)+epsln)
            sigma(i,j) = stable(i,j)*sig 
     &                 + (1.-stable(i,j))*AMIN1(sig,epsilon)
          enddo

          iwscale_use_hbl_eq_zt=0
          zt_kl_dummy=0.0
          call wscale(joff, j, j, istrt, iend
     &,               iwscale_use_hbl_eq_zt, zt_kl_dummy)
c
c-----------------------------------------------------------------------
c         compute the dimensionless shape functions at the interfaces
c         eqn. (11)
c-----------------------------------------------------------------------
c
          do i = istrt,iend
	    if (ki .lt. kbl(i,j)) then
              sig = (zt(ki) + 0.5 * dzt(ki)) / (hbl(i,j)+epsln)
              a1 = sig - 2.
              a2 = 3.-2.*sig
              a3 = sig - 1.
c
              Gm = a1
#if !defined kpp_no_matching
     &             + a2 * gat1(i,j,1) + a3 * dat1(i,j,1) 
#endif
              Gs = a1
#if !defined kpp_no_matching
     &             + a2 * gat1(i,j,2) + a3 * dat1(i,j,2)
#endif
              Gt = a1
#if !defined kpp_no_matching
     &             + a2 * gat1(i,j,3) + a3 * dat1(i,j,3)
#endif
c
c-----------------------------------------------------------------------
c             compute boundary layer diffusivities at the interfaces
c             eqn. (10)
c-----------------------------------------------------------------------
c
              blmc(i,ki,j,1) = hbl(i,j) * wm(i,j) * sig * (1.+sig*Gm)
              blmc(i,ki,j,2) = hbl(i,j) * ws(i,j) * sig * (1.+sig*Gt)
              blmc(i,ki,j,3) = hbl(i,j) * ws(i,j) * sig * (1.+sig*Gs)
c
c-----------------------------------------------------------------------
c             nonlocal transport term = ghats * <ws>o (eqn. 20)
c-----------------------------------------------------------------------
c
              ghats(i,ki,j) = (1.-stable(i,j))
#if !defined kpp_nonlocal_trans_no_match
     &             * cg 
     &             / (ws(i,j) * hbl(i,j) + epsln)
#else
     &             * cg 
c             Nondimensional nonlocal transport term under this option 
c             should not be multiplied by diffusivity when using it to form the
c             source term for the tracer equations.  The tracer source term 
c             should be: \partial_z(ghats*stf), where ghats < 1.
c             Non dimensional G(sigma) function;  needs to satisfy
c             G(sig=0)=G(sig=1)=0; G'(sig=0)=1; G'(sig=1)=0:
     &             * ( sig*(1.0+sig*(sig-2.0)) )
#endif
            endif
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     find diffusivities at kbl-1 grid level 
c-----------------------------------------------------------------------
c
      do j=js,je
        do i= istrt,iend
          sig      =  zt(kbl(i,j)-1)  / (hbl(i,j)+epsln)
          sigma(i,j) =stable(i,j) * sig 
     &              + (1.-stable(i,j)) * MIN(sig,epsilon)
        enddo
c
        iwscale_use_hbl_eq_zt=0
        zt_kl_dummy=0.0
        call wscale(joff, j, j, istrt, iend
     &,             iwscale_use_hbl_eq_zt, zt_kl_dummy)
c
        do i = istrt,iend
          sig = zt(kbl(i,j)-1) / (hbl(i,j)+epsln)
          a1= sig - 2.
          a2 = 3.-2.*sig
          a3 = sig - 1.
          Gm = a1 + a2 * gat1(i,j,1) + a3 * dat1(i,j,1)
          Gs = a1 + a2 * gat1(i,j,2) + a3 * dat1(i,j,2)
          Gt = a1 + a2 * gat1(i,j,3) + a3 * dat1(i,j,3)
          dkm1(i,j,1) = hbl(i,j) * wm(i,j) * sig * (1. + sig * Gm)
          dkm1(i,j,2) = hbl(i,j) * ws(i,j) * sig * (1. + sig * Gs)
          dkm1(i,j,3) = hbl(i,j) * ws(i,j) * sig * (1. + sig * Gt)
        enddo
      enddo
      return
      end

c ******************************************************************

      subroutine enhance (joff, js, je, istrt, iend) 
c
c enhance the diffusivity at the kbl-.5 interface
c
      use memory_window
# include "param.h"
# include "coord.h"
# include "mw.h"
# include "vmixc.h"
c
c input
c      integer kbl(imt,jmw)       ! grid above hbl
c      real hbl(imt,jmw)          ! boundary layer depth          (cm)
c      real dkm1(imt,jmw,3)       ! bl diffusivity at kbl-1 grid level
c      real caseA(imt,jmw)        ! = 1 in caseA, = 0 in case B
c
c input/output
c      real ghats(imt,km,jmw)     ! nonlocal transport     (s/cm**2)
c                              modified ghats at kbl(i)-1 interface
c output
c      real blmc(imt,km,jmw,3)  ! enhanced bound. layer mixing coeff.
c
c local
      real delta             ! fraction hbl lies beteen zt neighbors
c
      do j=js,je
	do ki=1,kbot-1
          do i = istrt,iend

            if(ki .eq. (kbl(i,j) - 1) ) then

              delta = (hbl(i,j)-zt(ki)) / (zt(ki+1)-zt(ki))

              dkmp5 = caseA(i,j) * visc_cbu(i,ki,j) 
     &               + (1.-caseA(i,j)) * blmc(i,ki,j,1)
              dstar = (1.-delta)**2 * dkm1(i,j,1) + delta**2 * dkmp5
              blmc(i,ki,j,1) = (1.-delta) * visc_cbu(i,ki,j) 
     &                        + delta * dstar

c             temperature:
              dkmp5 = caseA(i,j) * diff_cbt(i,ki,j,1) 
     &               + (1.-caseA(i,j)) * blmc(i,ki,j,2)
              dstar = (1.-delta)**2 * dkm1(i,j,3) + delta**2 * dkmp5    
              blmc(i,ki,j,2) = (1.-delta) * diff_cbt(i,ki,j,1) 
     &                      + delta * dstar

c             salinity:   
              dkmp5 = caseA(i,j) * diff_cbt(i,ki,j,2) 
     &               + (1.-caseA(i,j)) * blmc(i,ki,j,3)
              dstar = (1.-delta)**2 * dkm1(i,j,2) + delta**2 * dkmp5
              blmc(i,ki,j,3) = (1.-delta) * diff_cbt(i,ki,j,2) 
     &                      + delta * dstar
            
              ghats(i,ki,j) = (1.-caseA(i,j)) * ghats(i,ki,j)

            endif

          enddo
        enddo
      enddo
      return
      end

c***********************************************************************

      subroutine kppmixi (error, cifdef, ifdmax, nifdef, vmixset
     &,                   jstask,jetask)
c
      use arrays_2d

      logical vmixset, error
      character*(*) cifdef(ifdmax)
c

# include "param.h"
# include "vmixc.h"
# include "rinum.h"
c
c=======================================================================
c     Initialization for the kpp vertical mixing scheme
c     Large, McWilliams, & Doney ( Rev. of Geophys.,32,363_403, 1994).
c
c     input:
c       dzt    = thickness of vertical levels (cm)
c       km     = number of vertical levels
c       yt     = latitude of grid points (deg)
c       jmt    = number of latitudes
c       dtxcel = time step accelerator as a function of level
c       dtts   = density time step (sec)
c       dtuv   = internal mode time step (sec)
c       error  = logical to signal problems
c       cifdef = array of character strings for listing enabled "ifdefs"
c       ifdmax = size of "cifdef"
c       nifdef = current number of enabled "ifdefs"
c       vmixset= logical to determine if a vertical mixing scheme was
c                chosen
c
c     output:
c       lri = logical switch for shear instability mixing
c       ldd = logical switch for double-diffusive mixing
c       visc_cbu_limit = visc max due to shear instability  (cm**2/sec)
c       diff_cbt_limit = diffusivity ..                     (cm**2/sec)
c       visc_cbu_iw  = visc background due to internal waves(cm**2/sec)
c       diff_cbt_iw  = diffusivity ..                       (cm**2/sec)
c       visc_con_limit = visc due to convective instability (cm**2/sec)
c       diff_con_limit = diffusivity ..                     (cm**2/sec)
c       Vtc = non-dimensional constant used in calc. bulk Ri
c       cg  = constant used in calc.nonlocal transport term
c       wmt = turbulent velocity scale for momentum
c       wst = turbulent velocity scale for scaler
c       cifdef = array of character strings for listing enabled "ifdefs"
c       nifdef = incremented by 1 if this routine is called
c       error  = true if some inconsistancy was found
c       vmixset= true
c
c=======================================================================
c
c
c local
      real zehat                        ! = zeta *  ustar**3
      real zeta                         ! = stability parameter d/L
       
      namelist /kppmix/ lri,ldd 
     &,                diff_cbt_iw, visc_cbu_iw
     &,                visc_cbu_limit, diff_cbt_limit
     &,                visc_con_limit,diff_con_limit
     &,                concv,Ricr
c
      write (stdout,'(/,20x,a,/)')
     & 'K P P V M I X    I N I T I A L I Z A T I O N'
c
c-----------------------------------------------------------------------
c     initialize variables (all mixing units are cm**2/sec.)
c-----------------------------------------------------------------------
c
      lri = .true.
      ldd = .false.
      visc_cbu_limit = 50.0
      diff_cbt_limit = 50.0
      visc_cbu_iw    = 1.0
      diff_cbt_iw    = 0.1
      visc_con_limit = 1000.0
      diff_con_limit = 1000.0
c
c-----------------------------------------------------------------------
c     provide for namelist over-ride of above settings + documentation
c-----------------------------------------------------------------------
c
      call getunit (io, 'namelist.kppmix'
     &,             'formatted sequential rewind')
      read (io,kppmix,end=100)
100   continue
c
c-----------------------------------------------------------------------
c     set no-flux condition on density difference across bottom level
c-----------------------------------------------------------------------
c
      do j=1,jmw
        do i=1,imt
	  rhom1z(i,kbot,j) = c0
#ifdef bbl_ag
	  rhom1z(i,km,j) = c0
#endif
        enddo
      enddo
c
c--------------------------------------------------------------------
c     Initialize the boundary layer and mixed layer depths
c--------------------------------------------------------------------
c
      do j=jstask,jetask
        do i=1,imt
	   hblt(i,j) = c0
	   hmxl(i,j) = c0
        enddo
      enddo
c
c--------------------------------------------------------------------
c     Initialize viscosity and diffusivity
c--------------------------------------------------------------------
c
      do j=2,jmw-1
        do i=1,imt
#ifdef bbl_ag
c         values are initialized all the way into the bbl at k=km
#endif
          do k=1,km
            visc_cbu(i,k,j)     = c0
	    ghats(i,k,j)        = c0
            do n=1,nsmix
              diff_cbt(i,k,j,n) = c0
            enddo
          enddo
        enddo 
      enddo
c
c-----------------------------------------------------------------------
c     initialize some constants for kmix subroutines, and initialize
c     for kmix subroutine "wscale" the 2D-lookup table for wm and ws
c     as functions of ustar and zetahat (=vonk*sigma*hbl*bfsfc).
c-----------------------------------------------------------------------
c
c  
c-----------------------------------------------------------------------
c define some non-dimensional constants  
c-----------------------------------------------------------------------
c
c     Vtc used in eqn. 23
c
      Vtc     = concv * sqrt(0.2/concs/epsilon) / vonk**2 / Ricr
c
c     cg = cs in eqn. 20
c
      cg      = cstar * vonk * (concs * vonk * epsilon)**(1./3.)
c
c-----------------------------------------------------------------------
c construct the wm and ws lookup tables (eqn. 13 & B1)
c-----------------------------------------------------------------------
c
      deltaz = (zmax-zmin)/(nni+1) 
      deltau = (umax-umin)/(nnj+1)
      
      do i=0,nni+1
         zehat = deltaz*(i) + zmin
         do j=0,nnj+1
            usta = deltau*(j) + umin
            zeta = zehat/(usta**3+epsln)

            if(zehat.ge.0.) then
               wmt(i,j) = vonk*usta/(1.+conc1*zeta)
               wst(i,j) = wmt(i,j)
            else
               if(zeta.gt.zetam) then
                  wmt(i,j) = vonk* usta * (1.-conc2*zeta)**(1./4.)
               else
                  wmt(i,j) = vonk* (conam*usta**3-concm*zehat)**(1./3.)
               endif
               if(zeta.gt.zetas) then
                  wst(i,j) = vonk* usta * (1.-conc3*zeta)**(1./2.)
               else
                  wst(i,j) = vonk* (conas*usta**3-concs*zehat)**(1./3.)
               endif
            endif   
         enddo   
      enddo
c
c-----------------------------------------------------------------------
c     add character string to "ifdef option list" indicating which 
c     kpp options are enabled
c-----------------------------------------------------------------------
c
      nifdef = nifdef + 1
      cifdef(nifdef) = 'kppvmix     '
# ifdef kpp_nonlocal_trans_no_match
      nifdef = nifdef + 1
      cifdef(nifdef) = 'kpp_nonlocal_trans_no_match  '
# endif
# ifdef kpp_cstar_5     
      nifdef = nifdef + 1
      cifdef(nifdef) = 'kpp_cstar_5  '
# endif
c-----------------------------------------------------------------------
c     check for problems
c-----------------------------------------------------------------------
c
# if defined kppvmix && !defined implicitvmix
        write (stdout,'(/,(1x,a))')
     & '==> Error:  "kppvmix" must use "implicitvmix"   '
        error = .true.
# endif
# ifdef bryan_lewis_vertical
        write (stdout,'(/,(1x,a/1x,a/1x,a/1x,a))')
     & '==> Warning: "bryan_lewis_vertical" tracer diffus coefficients'
     &,'              will  be added to "kppvmix" diffus coefficients  '
     &,'              Note that diff_cbt_iw is being reset to zero  '
     &,'              while visc_cbu_iw is unchanged.               '
      diff_cbt_iw =  0.0
# endif
      if (vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "kppvmix" cannot be enabled because another   '
     &,'            vertical mixing scheme has been enabled           '
        error = .true.
      else
        vmixset = .true.
      endif
c
c     write out namelist values
c
      write (stdout,kppmix)
      call relunit (io)
      call getunit (iodoc, 'document.dta'
     &,             'formatted sequential append')
      write (iodoc, kppmix)
      call relunit (iodoc)
      return
      end



      block data kmixbdta
# include "param.h"
# include "vmixc.h"
c
c     data statements needed for "kmix" vertical mixing subroutines
c
c-----------------------------------------------------------------------
c     parameters for several subroutines
c-----------------------------------------------------------------------
c
      data epsilon                 /   0.1    /
      data vonk                    /   0.4    /
      data conc1                   /   5.0    /
      data conam,concm,conc2,zetam /   1.257  ,  8.380, 16.0, - 0.2 /
      data conas,concs,conc3,zetas / -28.86   , 98.96 , 16.0, - 1.0 /
c
c-----------------------------------------------------------------------
c     parameters for subroutine "bldepth"
c-----------------------------------------------------------------------
c
c     Ricr    = critical bulk Richardson Number            = 0.3
c     hbf     = fraction of bounadry layer depth to 
c               which absorbed solar radiation 
c               contributes to surface buoyancy forcing    = 1.0
c
      data Ricr      /  0.3  /
      data cekman    /  0.7  /
      data cmonob    /  1.0  /
      data concv     /  1.8  /
      data hbf       /  1.0  /
c
c-----------------------------------------------------------------------
c     limits for lookup table of wm and ws in subroutine "kmixinit"
c     (zmin and zmax are in cm3/s3, umin and umax are in cm/s)
c-----------------------------------------------------------------------
c
      data zmin,zmax / -4.e-1 ,  0.0  /
      data umin,umax /  0.    ,   4.0 /
c
c-----------------------------------------------------------------------
c     parameters for subroutine "ri_iwmix"
c     to compute vertical mixing coefficients below boundary layer:
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c     Riinfty = local Richardson Number limit 
c               for shear instability     
c-----------------------------------------------------------------------
c
      data Riinfty   /  0.8   /
c
c-----------------------------------------------------------------------
c     visc_cbu_limit = visc max due to shear instability (cm**2/sec)
c     diff_cbt_limit = diffusivity ..                    (cm**2/sec)
c     visc_cbu_iw    = visc background due to internal waves (cm**2/sec)
c     diff_cbt_iw    = diffusivity ..                        (cm**2/sec)
c     visc_con_limit = visc due to convective instability (cm**2/sec)
c     diff_con_limit = diffusivity ..                     (cm**2/sec)
c-----------------------------------------------------------------------
c
c      data visc_cbu_limit   /   50.0  /
c      data diff_cbt_limit   /   50.0  /
c      data visc_cbu_iw      /   10.0  / 
c      data diff_cbt_iw      /    0.3  / 
c      data visc_con_limit   / 1000.0  / 
c      data diff_con_limit   / 1000.0  / 
c
c-----------------------------------------------------------------------
c     parameters for subroutine "ddmix"
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c     to compute additional diffusivity due to double diffusion:
c     Rrho0   = limit for double diffusive density ratio
c     dsfmax  = maximum diffusivity in case of salt fingering (cm2/s)
c-----------------------------------------------------------------------
c
      data Rrho0     /    1.9   /
      data dsfmax    /   10.0   /
c
c-----------------------------------------------------------------------
c     parameters for subroutine "blmix_kpp"
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c     to compute mixing within boundary layer:
c
c     cstar   = proportionality coefficient
c               for nonlocal transport
c-----------------------------------------------------------------------
c
#ifdef kpp_cstar_5
      data cstar     /   5.0   /
#else
      data cstar     /   10.0   /
#endif
      end

#else
      subroutine kppmix
      return
      end
#endif




