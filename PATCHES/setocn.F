      subroutine setocn
c
c=======================================================================
c
c     set up everything which must be done only once per run
c
c    (mpp items)  v. balaji     e-mail=> vb@gfdl.gov
c=======================================================================
c

      use mom3_mpp_mod
      use mpp_io_mod
      use datamod
      use arrays_2d
      use memory_window
# ifdef posdef
      use posdef_mod
# endif
#ifdef restorst
      use sbc_info
#endif
#if defined sponges
      use xsponge_mod
#else
      logical annlev
#endif
c
      logical error, vmixset, htr_mixset, hvel_mixset
      parameter (ifdmax=200)
      character*30 cifdef(ifdmax)
c
#include "param.h"
#include "accel.h"
#if defined bbl_ag
# include "bbl.h"
#endif
#ifdef neptune
# include "cnep.h"
#endif
#include "coord.h"
#if defined fourfil || defined firfil
# include "filtr.h"
# include "cpolar.h"
#endif
#include "cprnts.h"
#ifdef shortwave
# include "cshort.h"
#endif
#include "diag.h"
#ifdef trajectories
# include "diag_ptraj.h"
#else
      logical initpt
#endif
#include "diag_region.h"
#include "diag_stab.h"
#include "docnam.h"
#include "emode.h"
#include "grdvar.h"
#include "hmixc.h"
#include "index.h"
#include "iounit.h"
#ifdef isoneutralmix
# include "isopyc.h"
#endif
#include "mw.h"
#include "scalar.h"
#include "state.h"
#include "switch.h"
#include "taskrows.h"
#include "tmngr.h"
#include "vers.h"
#include "vmixc.h"
#ifdef xlandmix
# include "xland.h"
#endif
#if defined max_window && defined parallel_1d
      real,dimension(imt*km,jmw,nt) :: tmp_var_t
      real,dimension(imt*km,jmw,2) :: tmp_var_u
      pointer (ptr_t,tmp_var_t)
      pointer (ptr_u,tmp_var_u)
#endif
      integer, dimension(:,:), allocatable :: kmz
      real, dimension(:,:), allocatable :: a2d
#ifdef rot_grid
      real, dimension(:,:), allocatable :: corbuf
#endif
#if !defined partial_cell
      real, dimension(:,:), allocatable :: htp, hup
#endif
      real, dimension(nt) :: dampts
c
      namelist /contrl/ init, runlen, rununits, segtim, restrt, initpt
     &,                 taux0,tauy0
      namelist /tsteps/ dtts, dtuv, dtsf
      namelist /riglid/ mxscan, tolrsf, tolrsp, tolrfs
     &,                 land_mass_a, land_mass_b
      namelist /mixing/ am, ah, ambi, ahbi
     &,                 kappa_m, kappa_h, cdbot
     &,                 spnep, senep
     &,                 aidif
     &,                 ncon, nmix, eb, acor, dampts, annlev
     &,                 smooth
     &,                 rjfrst, filter_reflat_n, filter_reflat_s
     &,                 eta_mix, neta
      namelist /diagn/  tsiint, tavgint, itavg, tmbint, tmbper, itmb
     &,                 stabint, zmbcint, glenint, trmbint, itrmb
     &,                 vmsfint, gyreint,igyre, extint, prxzint, trajint
     &,                 exconvint, dspint, dspper, snapint
     &,                 timavgint, timavgper, cmixint 
     &,                 prlat, prslon, prelon, prsdpt, predpt
     &,                 slatxy, elatxy, slonxy, elonxy
     &,                 cflons, cflone, cflats, cflate, cfldps, cfldpe
     &,                 maxcfl, xbtint, xbtper, crossint, pressint
     &,                 fctint, tyzint, rhoint
#if defined snapshots_regions 
     &,                 snapint_reg
#endif
      namelist /io/     expnam, iotavg 
     &,                 iotmb, iotrmb
     &,                 ioglen, iovmsf, iogyre
     &,                 ioprxz, ioext, iodsp, iotsi, iozmbc, iotraj
     &,                 ioxbt
      namelist /ictime/ year0, month0, day0, hour0, min0, sec0
     &,                 ryear, rmonth, rday, rhour, rmin, rsec
     &,                 refrun, refinit, refuser, eqyear, eqmon
#if defined bbl_ag
      namelist /bbl/ raydrag, Kvbbltop, cdbotbbl
     &,              entrain_kbbl, entrain_kbbl_constant
#endif
c
c-----------------------------------------------------------------------
c                                INITIALIZATIONS
c
c     grep xxx *.h to find the description (where xxx is what you want)
c
c-----------------------------------------------------------------------
#ifdef timing
      call tic ('driver', 'setocn')
#endif
c
c-----------------------------------------------------------------------
c     time at initial conditions.
c-----------------------------------------------------------------------
c
      year0  =1
      month0 =1
      day0   =1
      hour0  =0
      min0   =0
      sec0   =0
c
c-----------------------------------------------------------------------
c     set whether calculations for logical switches are referenced to
c
c     refrun  = T ==>  the start of each run
c     refinit = T ==>  initial condition time given by:
c                     year0, month0, day0, hour0, min0, sec0
c     refuser = T ==>  user specified reference time
c                     ryear, rmonth, rday, rhour, rmin, rsec
c
c     choose by setting one of the above to TRUE
c-----------------------------------------------------------------------
c
      refrun  = .false.
      refinit = .true.
      refuser = .false.
c
c     set the date and time for referencing switches
c     (only needed if rfuser = .T.)
c
      ryear  = 1900
      rmonth = 1
      rday   = 1
      rhour  = 0
      rmin   = 0
      rsec   = 0
c
c-----------------------------------------------------------------------
c     set type of calendar
c
c     if (eqyear and eqmon) then thirty_day_months
c     if (eqyear and !eqmon) then no_leap (feb always has 28 days)
c     if (!eqyear and !eqmon) then julian
c-----------------------------------------------------------------------
c
      eqyear   = .false.
      eqmon    = .false.
      call set_calendar (eqyear, eqmon)
c
c     stability diagnostic
c
      call stabi
c
c     misc variables
c
      error       = .false.
      vmixset     = .false.
      htr_mixset  = .false.
      hvel_mixset = .false.
      momver      = ' {MOM 3.0}'
      expnam      = ' MOM test case  '
      init        =.true.
      initpt      =.true.
c
c     initialize ifdef list
c
      do n=1,ifdmax
        cifdef(n) = '                              '
      enddo
c
      nifdef = 0
      cifdef(1) = 'no ifdefs enabled (turned on)'
c
c     integration time related variables 
c
      runlen   = 5.0
      rununits = 'days'
      segtim   = 1.0
      dtts     = -3600.0
      dtuv     = -3600.0
      dtsf     = -3600.0
      do k=1,km
        dtxcel(k) = 1.0
      enddo
c
c     control for leapfrog/mixing time steps
c
      nmix     = 17
      eb       = .true.
c
c     smoothing coefficient for robert time filter
c
      smooth = 0.01
c
c     tolerance settings for poisson solvers
c
      tolrsf   = 1.0e8
      tolrsp   = 5.0e4
      tolrfs   = 5.0e4
c
c     specified flow between two land masses
c
      land_mass_a = 0
      land_mass_b = 0
c
c     vertical and horizontal mixing + Newtonian damping variables
c
      kappa_h        = 1.0
      kappa_m        = 20.0
      am             = 1.0e9
      ah             = 2.0e7
      ambi           = 1.0e23
      ahbi           = 5.0e22
      cdbot          = 0.0
      visc_cbu_limit = 1.0e6
      diff_cbt_limit = 1.0e6

      neta           = 2
#ifdef explicit_eta_laplacian
      eta_mix        = 1.0e7
#endif
#ifdef explicit_eta_dpdc
      eta_mix        = 1.0
#endif
c
c     damping time scales for surface tracers in test case
c
      do n=1,nt
        dampts(n)    = 50.0
      enddo
      annlev         = .false.
c
c     number of convective passes for standard explicit convection
c
      ncon     = 1
c
c     physical constants
c
      rho0     = 1.035
      rho0r    = c1/rho0
      grav     = 980.6
      radius   = 6371.0e5
c
c     time centering for equations
c
#if defined implicitvmix || defined redi_diffusion
      aidif    = 1.0
#else
      aidif    = 0.0
#endif
c
#ifdef rigid_lid_surface_pressure
      alph     = 1.0
      gam      = 0.0
      theta    = 1.0
#endif
c
#ifdef implicit_free_surface
      alph     = 0.3333333
      gam      = 0.3333333
      theta    = 0.5
#endif
c
c     initialize neptune
c
      spnep = 0.0
      senep = 0.0
c
c     initialize polar filtering latitudes
c
      rjfrst = -81.0
      filter_reflat_s = -70.0
      filter_reflat_n = 70.0
c
c     initialize simple windstess for idealized equatorial studies
c
      taux0 = 0.0
      tauy0 = 0.0
c
#if defined bbl_ag
c
c     Coefficients for bottom boundary layer
c
      raydrag      = 2.0e-5
      cdbotbbl     = 3.0e-3
      Kvbbltop     = 200.0
c
c
c     Parameter to determine the vertical diffusion between the bbl
c     and interior:
c
      entrain_kbbl = 0.0
c
c     A Constant background vertical diffusion:
c
      entrain_kbbl_constant = 0.0
#endif
c
c-----------------------------------------------------------------------
c     write MOM version number and experiment name
c-----------------------------------------------------------------------
c
      write (stdout,'(/1x,a,a/)') 'GFDL MOM version = ',momver
      write (stdout,'(/1x,a,a/)') 'Experiment name = ',expnam
c
c-----------------------------------------------------------------------
c     provide for change in above presets using "namelist"
c-----------------------------------------------------------------------
c
c
      call getunit (ioun, 'namelist.contrl'
     &,             'formatted sequential rewind')
      read  (ioun, contrl, IOSTAT=io_status)
      write (stdout,contrl)
      call relunit (ioun)
      if (io_status .gt. 0) then
        write(stdout,*) '==> ERROR reading namelist.contrl'
     &,  ' IOSTAT =',io_status
        call abort()
      endif  
c
      call getunit (ioun, 'namelist.tsteps'
     &,             'formatted sequential rewind')
      read  (ioun, tsteps)
      write (stdout,tsteps)
      call relunit (ioun)
c
      call getunit (ioun, 'namelist.riglid'
     &,             'formatted sequential rewind')
      read  (ioun, riglid)
      write (stdout,riglid)
      call relunit (ioun)
c
      call getunit (ioun, 'namelist.mixing'
     &,             'formatted sequential rewind')
      read  (ioun, mixing)
      write (stdout,mixing)
      call relunit (ioun)
c
      call getunit (ioun, 'namelist.diagn'
     &,             'formatted sequential rewind')
      read  (ioun, diagn)
      write (stdout,diagn)
      call relunit (ioun)
c
      call getunit (ioun, 'namelist.io'
     &,             'formatted sequential rewind')
      read  (ioun, io)
      write (stdout,io)
      call relunit (ioun)
c
      call getunit (ioun, 'namelist.ictime'
     &,             'formatted sequential rewind')
      read  (ioun, ictime)
      write (stdout,ictime)
      call relunit (ioun)
c
#if defined bbl_ag
      call getunit (ioun, 'namelist.bbl', 'formatted sequential rewind')
      read  (ioun, bbl)
      write (stdout,bbl)
      call relunit (ioun)
#endif
c
      call getunit (iodoc, 'document.dta'
     &,             'formatted sequential rewind')
      write (iodoc, contrl)
      write (iodoc, tsteps)
      write (iodoc, riglid)
      write (iodoc, mixing)
      write (iodoc, diagn)
      write (iodoc, io)
      write (iodoc, ictime)
#if defined bbl_ag
      write (iodoc, bbl)
#endif
      write (iodoc,'(a,i6)') 'imt=', imt, 'jmt=', jmt, 'km=', km
     &,                       'nt=', nt, 'nvar=', nvar
c
c     user specified tracer names are placed into "trname" here.
c
      do m=1,nt
        trname(m) = '**unknown***'
      enddo
      trname(1) = 'temp        ' ! 'potentl_temp'
      trname(2) = 'salinity    '
      write (iodoc,'(a,a)') 'trname(1)=', trname(1)
     &,                     'trname(2)=', trname(2)
      call relunit (iodoc)
c
c-----------------------------------------------------------------------
c     open i/o units needed for prognostic variables
c     nbuf = number of i/o buffers (currently not used) 
c-----------------------------------------------------------------------
c
#if defined ramdrive 
      call getunitnumber (kflds)
      call getunitnumber (latdisk(1))
      call getunitnumber (latdisk(2))
#endif
c
      nbuf = 2
      call ostart (kflds, 'kflds', nkflds*nwds, nwds, nbuf)
      call ostart (latdisk(1), 'latdisk1', jmt*nslab, nslab, nbuf)
      call ostart (latdisk(2), 'latdisk2', jmt*nslab, nslab, nbuf)
c
c-----------------------------------------------------------------------
c     set up the grids in x (longitude), y (latitude), and z (depth)
c     corresponding to Arakawa "b" grid system
c-----------------------------------------------------------------------
c
      call grids
c
c-----------------------------------------------------------------------
c     compute density coefficients based on depth of grid points
c-----------------------------------------------------------------------
c
      call eqstate (zt, km, ro0, to, so, c, tmink, tmaxk, smink, smaxk)
c
      cksum = checksum(c, km, 9)
      write (stdout,*)
     & '      Checksum for density coefficients = ', cksum
c
#if defined linearized_advection || defined equatorial_thermocline
c
c     set the initial idealized stratification as a function of
c     temperature only. Salinity is fixed at 35 psu (tbar(k,2)=0) but 
c     passive tracers (if nt>2) are 1.0 at k=1 and 0.0 for k>1 
c
      tzero = 7.5
      tone  = 10.0
      z0    = 30.0e2
      bigl  = 80.0e2
      write (stdout,'(//a/)')
     & 'Initial Tracer profile as a function of depth'
      do k=1,km
        tbarz(k,1) = tzero*(1.0-tanh((zt(k)-bigl)/z0)) + 
     &                tone*(1.0-zt(k)/zt(km))
# if defined salinity_psu
        tbarz(k,2) = 35.0
# else
        tbarz(k,2) = c0        
# endif
        if (nt .gt. 2) then
	  do n=3,nt
	    if (k .eq. 1) then
	      tbarz(k,n) = c1
	    else
	      tbarz(k,n) = c0
	    endif
	  enddo
	endif
        write (stdout,'(1x,"k=",i3,", zt(k)=",f8.1
     &,        "cm,  tbarz(k,n),n=1,nt =",10e14.7)')
     &      k, zt(k),(tbarz(k,n),n=1,nt)
      enddo
#endif
c
c-----------------------------------------------------------------------
c     set the earth`s rotation rate to computer precision
c-----------------------------------------------------------------------
c
      omega  = pi/43082.0
      call getunit (iodoc, 'document.dta'
     &,             'formatted sequential append')
      write (iodoc,'(a,e14.7)') 'omega=', omega
      call relunit (iodoc)
c
#if !defined max_window
c
c     time level indices are fixed in time
c
      taum1 = -1
      tau   =  0
      taup1 = +1
#endif
c
c-----------------------------------------------------------------------
c     allocate 2-D arrays for entire domain on each processor
c-----------------------------------------------------------------------
c
      allocate (kmt(imt,jmt))
      allocate (kmu(imt,jmt))
      allocate (map(imt,jmt))
      allocate (kmz(imt,jmt))
      allocate (h(imt,jmt))
      allocate (hr(imt,jmt))
      allocate (htp(imt,jmt))
      allocate (hup(imt,jmt))
c
c-----------------------------------------------------------------------
c     set up topography and derived fields
c-----------------------------------------------------------------------
c
      call topog (kmt, map, xt, yt, zt, xu, yu, zw, dzt, imt, jmt, km
     &,           htp)
      call calc_kmu (kmt, htp, imt, kmu, hup, 1, jmt, jmt)
      call depth_u (kmu, hup, imt, zw, km, h, hr, 1, jmt)
c
#ifdef stability_tests
c
c-----------------------------------------------------------------------
c     convert starting and ending longitudes for the stability tests
c     to nearest model grid points.
c-----------------------------------------------------------------------
c
      if (stabint .ge. c0) then
        iscfl  = max(indp (cflons, xt, imt), 2)
        cflons = xt(iscfl)
        iecfl  = min(indp (cflone, xt, imt), imt-1)
        cflone = xt(iecfl)
        jscfl  = max(indp (cflats, yt, jmt), 2)
        cflats = yt(jscfl)
        jecfl  = min(indp (cflate, yt, jmt), jmt-1)
        cflate = yt(jecfl)
        kscfl  = indp (cfldps, zt, km)
        cfldps = zt(kscfl)
        kecfl  = indp (cfldpe, zt, km)
        cfldpe = zt(kecfl)
      endif
#endif
c
#ifdef restorst
c
c-----------------------------------------------------------------------
c      damp surface tracers back to data. schematically, the restoring
c      term will be = (dzt(1)/(dampts(n)*86400)) * (data - t)
c      where dampts(n) is the time scale (days) for the nth tracer
c-----------------------------------------------------------------------
c
      write (stdout,'(/a/)') 'Restoring coefficient for surface tracers'
      allocate (damp_factor(nt))
      do n=1,nt
        write (stdout,'(/,1x,a,i2,a,a,/,1x,1pe14.7,a,1pe14.7,a)')
     &  ' Surface tracer #',n,' will be damped back to data using a'
     &, ' Newtonian restoring time scale of '
     &,  dampts(n),' days. and a level thickness =', dzt(1),' cm.'
        damp_factor(n) = dzt(1)*secday/dampts(n)
        if (n .eq. 1) then
          cwatts = 4.186e4
          dampc = cwatts*damp_factor(n)
          write (stdout,'(1x,a,1pe14.7,a)')
     &  ' resulting in a damping coefficient of ', dampc
     &, ' Watts/M**2/tracer unit'
        elseif (n .eq. 2) then
          ckg = 10.0
          dampc = ckg*damp_factor(n)
          write (stdout,'(1x,a,1pe14.7,a)')
     &  ' resulting in a damping coefficient of ', dampc
# if defined salinity_psu
     &, ' kg/m**2/s/psu'
# else
     &, ' kg/m**2/s/model salinity units'
# endif
        else
          ckg = 10.0
          dampc = ckg*damp_factor(n)
          write (stdout,'(1x,a,1pe14.7,a)')
     &  ' resulting in a damping coefficient of ', dampc
     &, ' kg/M**2/sec/tracer unit'
        endif
      enddo
      call getunit (iodoc, 'document.dta'
     &,             'formatted sequential append')
      write (iodoc,'(a,(1x,e14.7))') 'dampts=', dampts
      write (iodoc,'(a,(1x,e14.7))') 'dampdz=',dzt(1)
      call relunit (iodoc)
#endif
#ifdef shortwave
c
c-----------------------------------------------------------------------
c     Solar Shortwave energy penetrates below the ocean surface. Clear
c     water assumes energy partitions between two exponentials as 
c     follows:
c     
c     58% of the energy decays with a 35 cm e-folding scale
c     42% of the energy decays with a 23 m e-folding scale
c
c     if the thickness of the first ocean level "dzt(1)" is 50 meters,
c     then shortwave penetration wouldn't matter. however, for 
c     dzt(1) = 10 meters, the effect can be significant. this can be
c     particularly noticable in the summer hemisphere.
c
c     Paulson and Simpson (1977 Irradiance measurements in the upper
c                               ocean JPO 7, 952-956)
c     Also see ... Jerlov (1968 Optical oceanography. Elsevier)
c                  A General Circulation Model for Upper Ocean
c                  Simulation (Rosati and Miyakoda JPO vol 18,Nov 1988)
c-----------------------------------------------------------------------
c
      write (stdout,'(/a/)') 
     &' => Shortwave penetration is a double exponential as follows:'
      rpart  = 0.58
      efold1 = 35.0e0
      efold2 = 23.0e2
      rscl1  = 1.0/efold1
      rscl2  = 1.0/efold2
      call getunit (iodoc, 'document.dta'
     &,             'formatted sequential append')
      write (iodoc,'(a,e14.7)') 'rpart=', rpart, 'efold1=',efold1
     &, 'efold2=',efold2
      call relunit (iodoc)
c
c     note that pen(0) is set 0.0 instead of 1.0 to compensate for the
c     shortwave part of the total surface flux in "stf(i,1)"
c
      pen(0) = c0
c
      do k=1,km
        swarg1 = -min(zw(k)*rscl1,70.0)
        swarg2 = -min(zw(k)*rscl2,70.0)
        pen(k) = rpart*exp(swarg1) + (c1-rpart)*exp(swarg2)
	divpen(k) = (pen(k-1) - pen(k))/dzt(k)
	write (stdout,9200) k, zw(k)*1.e-2, pen(k), zt(k)*1.e-2
     &,                     divpen(k)
      enddo
      write (stdout,*) ' '
9200  format (1x,'k=',i3,' zw(k)=',f8.2,'(m)  pen(k)=',e14.7
     &,       ' zt(k)=',f8.2,'(m)  divpen(k)=',e14.7)
#endif
c
c-----------------------------------------------------------------------
c     compute surface area and volume of ocean (T cells and U cells)
c     (note that areas are defined at each level)
c-----------------------------------------------------------------------
c
      do k=1,km
        tcella(k) = c0
        ucella(k) = c0
      enddo
      ocnp   = c0
      tcellv = c0
      ucellv = c0
      do jrow=1,jmt
        tcellv_j(jrow) = c0
        ucellv_j(jrow) = c0
      enddo
c
c     this comment directive turns off autotasking on the YMP
c     for the following loop
c
cfpp$ noconcur l
      do jrow=2,jmt-1
        do i=2,imtm1
          if (kmt(i,jrow) .gt. 0) then
	    do k=1,kmt(i,jrow)
              tcella(k) = tcella(k) + cst(jrow)*dxt(i)*dyt(jrow)
	    enddo
            tcellv_j(jrow) = tcellv_j(jrow) + cst(jrow)*dxt(i)*dyt(jrow)
#if defined partial_cell
     &                                                  *htp(i,jrow)
#else
     &                                                  *zw(kmt(i,jrow))
#endif
#ifdef bbl_ag
	    tcella(km) = tcella(km) + cst(jrow)*dxt(i)*dyt(jrow)
            tcellv_j(jrow) = tcellv_j(jrow)
     &                       + cst(jrow)*dxt(i)*dyt(jrow)*dzt(km)
#endif
            ocnp   = ocnp + float(kmt(i,jrow))
          endif
          if (kmu(i,jrow) .gt. 0) then
	    do k=1,kmu(i,jrow)
              ucella(k) = ucella(k) + csu(jrow)*dxu(i)*dyu(jrow)
	    enddo
            ucellv_j(jrow) = ucellv_j(jrow) + csu(jrow)*dxu(i)*dyu(jrow)
#if defined partial_cell
     &                                                  *hup(i,jrow)
#else
     &                                                  *zw(kmu(i,jrow))
#endif
#ifdef bbl_ag
	    ucella(km) = ucella(km) + csu(jrow)*dxu(i)*dyu(jrow)
            ucellv_j(jrow) = ucellv_j(jrow)
     &                       + csu(jrow)*dxu(i)*dyu(jrow)*dzt(km)
#endif
          endif
        enddo
        tcellv = tcellv + tcellv_j(jrow)
        ucellv = ucellv + ucellv_j(jrow)
      enddo
c
      write (stdout,9341) tcella(1), tcellv, ucella(1), ucellv
c
#if defined tracer_averages || defined term_balances
# define averages_balances
#endif
c
#ifdef meridional_overturning_basin
c-----------------------------------------------------------------------
c     read the basin mask to divide the ocean into basins for 
c     meridional transport calculation
c---------------------------------------------------------------------
c
      write(stdout,*)
     &    " Reading basin mask for meridional overturning"
     &,   " from file northmasks."
c
      call merid1st
c
#endif
c
#if defined averages_balances || defined gyre_components
c
c---------------------------------------------------------------------
c     set the horizontal regional masks and names to be used when
c     computing averages on the T grid in subroutine "sethr".
c     also set the vertical regional masks and names for use (along
c     with the horizontal ones) in term balance calculations for
c     tracers & momentum. For term balance calculations the number
c     of masks is the product of horizontal & vertical regions.
c---------------------------------------------------------------------
c
      allocate (mskhr(imt,jmt))
c
# ifdef readrmsk
c
c    read in horizontal & vertical regional masks ("mskhr" & "mskvr")
c    and names ("hregnm" & "vregnm") on unit iormsk
c
      call getunit (iormsk, 'region_masks'
     &,             'formatted sequential rewind')
      call reg1st (iormsk, .true., .false., .false., .false., .true.)
      call relunit (iormsk)
# else
c
c     set up the horizontal regions:
c     specify "mskhr" and "hregnm" values rather than reading them in
c     from a file (arbitrarily defaulted below for 5 zonal bands).
c     note: there must be "nhreg" calls ... one for each horizontal
c     region. The form is:
c     call sethr (region number, starting lon, ending lon, starting lat,
c                 ending lat)
c     where the lon and lat limits are for the edges of the T cells
c     "sethr" will fit the region using the nearest model grid points
c     and the region number will be increased by one for each call
c
      do jrow=1,jmt
        do i=1,imt
          mskhr(i,jrow) = 0
        enddo
      enddo
c
      n = 1
      start_lat = -90.0
      end_lat   = -60.1
      call sethr (n, xu(1), xu(imtm1), start_lat, end_lat)
c
      start_lat = -60.1
      end_lat   = -24.9
      call sethr (n, xu(1), xu(imtm1), start_lat, end_lat)
c
      start_lat = -24.9
      end_lat   = 24.9
      call sethr (n, xu(1), xu(imtm1), start_lat, end_lat)
c
      start_lat = 24.9
      end_lat   = 59.9
      call sethr (n, xu(1), xu(imtm1), start_lat, end_lat)
c
      start_lat = 59.9
      end_lat   = 90.0
      call sethr (n, xu(1), xu(imtm1), start_lat, end_lat)
c
      n = n - 1
      if (n .ne. nhreg) then
        write (stdout,*)
     &   '=>Error: parameter "nhreg" was specified as ',nhreg,' but '
     &,  '         ',n,' horizontal regions were set up in setocn.'
     &,  '         "nhreg" must equal the number of regions set up.' 
	error = .true.
      endif
c
      do jrow=1,jmt
        mskhr(1,jrow)   = 0
        mskhr(imt,jrow) = 0
      enddo
c
c     set up the vertical regions:
c     specify "mskvr" and "vregnm" values. also the starting & ending
c     levels for the vertical regions (arbitrarily defaulted for two
c     vertical regions). Note: there must be "nvreg" calls... one for
c     each vertical region. 
c     The form is:
c     call setvr (region number, starting depth, ending depth)
c     where the depth limits are for the edges of the T cells
c     "setvr" will fit the region using the nearest model grid points
c     and the region number will be increased by one for each call
c
      n = 1
      call setvr (n, 0.0e2, 245.0e2)
      call setvr (n, 245.0e2, 1570.e2)
c
      n = n - 1
      if (n .gt. nvreg) then
        write (stdout,*) "=>Error: increase parameter 'nvreg' "
	write (stdout,*) "         stop in setocn"
	call abort()
      endif
      if (n .ne. nvreg) then
        write (stdout,*)
     &   '=>Error: parameter "nvreg" was specified as ',nvreg,' but '
     &,  '         ',n,' vertical regions were set up in setocn.'
     &,  '         "nvreg" must equal the number of regions set up.'
        error = .true. 
      endif
c
      do k=1,km
        mskvr(k) = 0
      enddo
      do n=1,nvreg
        ks = llvreg(n,1)
        ke = llvreg(n,2)
        do k=ks,ke
          mskvr(k) = n
        enddo
      enddo
# endif
# if defined tracer_averages
c
c     print out regional mask info for tracer average diagnostic
c
      if ((iotavg .eq. stdout .or. iotavg .lt. 0) .and. itavg) then
        if (pn .eq. 1) then
          call reg1st (stdout, .true., .false., .false., .true.,.false.)
! t.sft
#  ifdef writrmsk
          write(stdout,*) 'print region_masks'
          call getunit (iormsk, 'region_masks'
     &,             'formatted sequential rewind')
          call reg1st(iormsk, .true., .false., .false., .true., .false.)
          call relunit (iormsk)
#  endif
        endif
      endif
# endif
c
c-----------------------------------------------------------------------
c     compute the surface area and volume of the ocean regions. index 0
c     refers to the sum of all regions.
c     (values used in subroutine region are done in terms of meters,
c     rather than centimeters)
c-----------------------------------------------------------------------
c
      areag = c0
      volgt  = c0
c
      do k=1,km
        volgk(k) = c0
      enddo
c
      do n=0,numreg
        areat(n) = c0
        areau(n) = c0
        volt(n)  = c0
        volu(n)  = c0
      enddo
      do mask=1,nhreg
        areab(mask) = c0
        volbt(mask) = c0
        do k=1,km
          volbk(mask,k) = c0
        enddo
      enddo
c
      do jrow=2,jmt-1
        do i=2,imtm1
          mask = mskhr(i,jrow)
          kz   = kmt(i,jrow)
          if (kz .gt. 0 .and. mask .gt. 0) then
            boxat = cst(jrow) * dxt(i) * dyt(jrow)
            if (kmu(i,jrow) .ne. 0) then
              boxau = csu(jrow) * dxu(i) * dyu(jrow)
            else
              boxau = c0
            endif
            areat(0) = areat(0) +  boxat
            areau(0) = areau(0) +  boxau
	    areab(mask)  = areab(mask) + boxat * 1.e-4
	    do k=1,kz
              thick_t = dzt(k)
              thick_u = dzt(k)
# ifdef partial_cell
              if (k .eq. kz .and. k .gt. 1) then
                thick_t = htp(i,jrow) - zw(k-1)
              endif
              if (k .eq. kmu(i,jrow) .and. k .gt. 1) then
                thick_u = hup(i,jrow) - zw(k-1)
              endif
# endif
              volbk(mask,k) = volbk(mask,k) + boxat * thick_t * 1.e-6
              dvolt   = boxat * thick_t
              if (kmu(i,jrow) .ge. k) then
                dvolu   = boxau * thick_u
              else
                dvolu   = c0
              endif
              n = nhreg*(mskvr(k)-1) + mask
	      if (n .gt. 0) then
                volt(0) = volt(0) +  dvolt
                volu(0) = volu(0) +  dvolu
                volt(n) = volt(n) +  dvolt
                volu(n) = volu(n) +  dvolu
                do nv=1,nvreg
                  ks = llvreg(nv,1)
                  if (k .eq. ks) then
                    areat(n) = areat(n) +  boxat
                    areau(n) = areau(n) +  boxau
                  endif
                enddo
	      endif
            enddo
          endif
        enddo
      enddo
c
      do mask=0,numreg
        if (areat(mask) .eq. c0) then
          rareat(mask) = c0
        else
          rareat(mask) = c1/areat(mask)
        endif
c
        if (areau(mask) .eq. c0) then
          rareau(mask) = c0
        else
          rareau(mask) = c1/areau(mask)
        endif
c
        if (volt(mask) .eq. c0) then
          rvolt(mask) = c0
        else
          rvolt(mask) = c1/volt(mask)
        endif
c
        if (volu(mask) .eq. c0) then
          rvolu(mask) = c0
        else
          rvolu(mask) = c1/volu(mask)
        endif
      enddo
c
      do mask=1,nhreg
        do k=1,km
          volbt(mask) = volbt(mask) + volbk(mask,k)
          volgk(k) = volgk(k) + volbk(mask,k)
        enddo
        areag = areag + areab(mask)
        volgt = volgt + volbt(mask)
      enddo
c
      if (iotavg .eq. stdout .or. iotavg .lt. 0) then
        call reg1st (stdout, .true., .true., .true., .false., .false.)
      endif
c
c     compute and print statistics for regions
c
      sum = c0
      do n=1,numreg
        sum = sum + volt(n)
      enddo
      sum    = 100.0*sum/tcellv
      pctocn = 100.0*ocnp/float((imt-2)*(jmt-2)*km)
      diffa  = 100.0 * (c1 - (tcella(1) - 10000.0*areag)/tcella(1))
c
      write (stdout,9342) diffa, numreg, sum, pctocn
9342  format ('  the horizontal regional masks cover',f8.3
     &, '% of the total ocean surface area.'/
     &, '  there are ', i6, ' regions over which tracer & '
     &, 'momentum balances will be computed',/,'  accounting for '
     &, f6.2,'% of the total ocean volume.'/
     &, 1x,f6.2,'% of the grid points lie within the ocean.'/)
#endif
c
#if defined averages_balances || defined gyre_components
c
c     distribute mskhr across processors
c
      do jrow=1,jmt
        do i=1,imt
	  kmz(i,jrow) = mskhr(i,jrow)
	enddo
      enddo
      deallocate (mskhr)
      allocate (mskhr(imt,jstask:jetask))
      do jrow=jstask,jetask
        do i=1,imt
          mskhr(i,jrow) = kmz(i,jrow)
        enddo
      enddo
#endif
c
c-----------------------------------------------------------------------
c     find all land mass perimeters for poisson solvers
c-----------------------------------------------------------------------
c
      call isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle
     &,                    imt, jmt, km, mnisle, maxipp)
#ifdef zero_island_flow
c
c     set no-net-flow between land_mass_a and land_mass_b by combining
c     their perimeters and re-ordering related quantities. logically,
c     the number of distinct land masses will be reduced by one.
c     physically, the original number of distinct land masses remains.
c
c     no-flow between other land masses may be set by duplicating this
c     code and working with the new land mass numbers.
c
      write (stdout,'(/,1x,a,i2,a,i2,/)')
     & '=>Warning: Enforcing zero net transport between land mass'
     &, land_mass_a,' and land mass ',land_mass_b 
c
      land_1 = min(land_mass_a, land_mass_b)
      land_2 = max(land_mass_a, land_mass_b)
      if (land_2 .gt. nisle .or. land_1 .lt. 1) then
        write (stdout,'(/a/a,i2,a,i2/a)')
     & '=>Error: specify valid land masses from "map(i,j)" numbers'
     &,'         land_mass_a = ',land_mass_a, ' land_mass_b = '
     &,  land_mass_b
     &,'         set "land_mass_a" and "land_mass_b" in namelist'
        call abort()
      endif
      if (land_2 .eq. land_1) then
        write (stdout,'(/a/a)')
     & '=>Error: two distinct land mass numbers must be specified'
     &,'         set "land_mass_a" and "land_mass_b" in namelist'
        call abort()
      endif
c
c     equate land_1 and land_2 and re-number the land masses
c
      do j=1,jmt
        do i=1,imt
          if (map(i,j) .eq.  land_2) map(i,j) = land_1
	  if (map(i,j) .gt.  land_2) map(i,j) = map(i,j)-1
	  if (map(i,j) .lt. -land_2) map(i,j) = map(i,j)+1
        enddo
      enddo
c
c     combine perimeter points of land_1 and land_2 and re-order all
c
      num = 0
      do m=1,nisle
        if (m .le. land_1 .or. m .gt. land_2) then
	  isle = m
	elseif (m .eq. land_1+1) then
	  isle = land_2
	elseif (m .gt. land_1+1 .and. m .le. land_2) then
	  isle = m-1
        endif
        do n=1,nippts(isle)
	  num = num + 1
          ipermt(num,1) = iperm(iofs(isle)+n)
          ipermt(num,2) = jperm(iofs(isle)+n)
        enddo
      enddo
      do n=1,maxipp
        iperm(n) = ipermt(n,1)
        jperm(n) = ipermt(n,2)
      enddo
c
c     re-calculate length of perimeters and re-order 
c
      nippts(land_1) = nippts(land_1) + nippts(land_2)
      do isle=land_2,nisle-1
        nippts(isle) = nippts(isle+1)
      enddo
      nisle = nisle - 1
      do isle=2,nisle
        iofs(isle) = nippts(isle-1) + iofs(isle-1)
      enddo
c
# if !defined skip_island_map
      write (stdout,'(/a/)')
     &   'Use -Dskip_island_map to eliminate the following map'
      write (stdout,'(/a,a,i3,a,i3/)')
     &   'Re-mapped land masses after specifying zero net flow between '
     &,  'disconnected land masses with numbers =',land_1,' & ',land_2
      call showmap (map, imt, jmt, linewidth)
# else
      write (stdout,'(/a/)')
     &   'Remove -Dskip_island_map to show island re-mapping'
# endif
#endif
c
c     set mask for land mass perimeters on which to perform calculations
c     imask(-n) = .false.  [no equations ever on dry land mass n]
c     imask(0)  = .true.   [equations at all mid ocean points]
c     imask(n)  = .true./.false [controls whether there will be
c                                equations on the ocean perimeter of
c                                land mass n]
c     note: land mass 1 is the northwest-most land mass
c           usually includes the "north pole", and at low resolutions,
c           the "main continent"
c     for the numbering of the other landmasses, see generated map(i,j)
c           by selecting option -Dshow_details
c
      do isle=-mnisle,mnisle
        if (isle .ge. 0 .and. isle .le. nisle) then
          imask(isle) = .true.
        else
          imask(isle) = .false.
        end if
      end do
#ifdef symmetry
c
c     do not perform island integrals for land masses whose perimeters
c     touch the equator in models symmetric about the equator
c
      do isle=1,nisle
        do n=1,nippts(isle)
          j = jperm(iofs(isle)+n)
          if (j .eq. jmt-1) then 
	    imask(isle) = .false.
	  endif
        end do
      end do
#endif
c
c     user-specified changes to island mask
c       imask(1) = .true.
c       imask(2) = .true.
c
c     there are problems if imask is set .true. for a nonexistent
c     island.
c
c     print diagnostic information
c
      do isle=-mnisle,mnisle
        if (imask(isle)) then
	  if (isle .eq. 0) then
            write (stdout,'(a)') 
     &      '=> calculations enabled for mid ocean points'
	  else
            write (stdout,'(2a,i3)') 
     &       '=> calculations enabled for ocean '
     &,                      'perimeter of land mass',isle
	  end if
	end if
      end do
      do isle=0,nisle
        if (.not. imask(isle)) then
            write (stdout,'(2a,i3)') 
     &      '=> calculations disabled for ocean '
     &,                      'perimeter of land mass',isle
	end if
      end do
c
c     imain is the land mass on which dpsi is normalized to 0
c     if imain is 0, then dpsi is not normalized.
c     default value of imain is land mass with longest perimeter
c
      imain = min(2,nisle)
      do isle=1,nisle
        if (nippts(isle) .gt. nippts(imain)) then
          imain = isle
        end if
      end do
c
c     if any island perimeter is not calculated, imain must be one such
c
      do isle=1,nisle
        if (.not.imask(isle)) then
          imain = isle
        end if
      enddo
#ifdef symmetry
c
c     do not normalize dpsi when using symmetry about the equator
c
      imain = 0
#endif
c
      if (imain .gt. 0 .and. imain .le. nisle) then
        write (stdout,'(a,i4)')
     &   'dpsi normalized to zero on land mass',imain
      else if (imain .eq. 0) then
        write (stdout,*) 'no normalization on dpsi'
      else
        write (stdout,*) 'ERROR: illegal value for  normalization of '
     &,           'land mass, imain =', imain
      end if
      write (stdout,*) ' ( "imain" to any valid land mass number)'
c
c---------------------------------------------------------------------
c     compute checksum of density coefficients
c---------------------------------------------------------------------
c
      write (stdout,*) ' '
      call print_checksum (c(1,1), km, 9
     &,                   ' density coefficient checksum = ')
c
#if defined fourfil || defined firfil
      write (stdout,'(/a/)') 'Filtering details for latitudes'
c
c     compute sin and cos values for vector correction before filtering
c
      fx =  1.0e-10
      fxa = dxt(1)/radius
c
      do i=2,imtm1
        fxb = fxa*float(i-2)
        spsin(i) = sin(fxb)
        spcos(i) = cos(fxb)
        if (abs(spsin(i)) .lt. fx) spsin(i) = c0
        if (abs(spcos(i)) .lt. fx) spcos(i) = c0
      enddo
c
      spsin(1) = spsin(imt-1)
      spcos(1) = spcos(imt-1)
      spsin(imt) = spsin(2)
      spcos(imt) = spcos(2)
c
c     set up model indices for filtering polar latitudes
      call set_polar_filtering_indices

      call getunit (iodoc, 'document.dta','formatted sequential append')
      write (iodoc,'(2(a,e14.7),5(a,i6))')
     & 'filter_reflat_n=',filter_reflat_n
     &,'filter_reflat_s=',filter_reflat_s, 'jfrst=', jfrst
     &, 'jft1=',jft1, 'jft2=',jft2,  'jfu1=',jfu1, 'jfu2=',jfu2
      call relunit (iodoc)
      jskpt  = jft2-jft1
      jskpu  = jfu2-jfu1
      njtbft = (jft1-jfrst+1)+(jmtm1-jft2+1)
      njtbfu = (jfu1-jfrst+1)+(jmtm1-jfu2+1)
      if (njtbft .gt. jmtfil .or. njtbfu .gt. jmtfil) then
        write (stdout,9599) njtbft, njtbfu
        write (stderr,9599) njtbft, njtbfu
        call abort()
      endif
      if (jfrst .gt. jft1) then
        write (stdout,*) '=>Error: jfrst=',jfrst,' > jft1=',jft1
	call abort()
      endif
9551  format (/' ==== start and end indices for fourier filtering ====')
9552  format (' only 11 sets of indices fit across the page.',
     &       '  others will not be printed.'/)
9553  format (/,' == filtering indices for t,s ==')
9554  format (/,' == filtering indices for u,v ==')
9555  format (/,' == filtering indices for stream function ==')
9599  format (/,' error => jmtfil must be >= max(njtbft,njtbfu)',
     &        /,'          njtbft=',i8,' njtbfu=',i8)
c
# ifdef firfil
c
c     set "numflt" and "numflu" to filter polar latitudes
c
c     USER INPUT: set numflt and numflu only large enough to eliminate
c                 numerical instability due to convergence of meridians.
c                 some experimentation may have to be done to see what
c                 works for a given geometry.
c
      write(stdout,9501) filter_reflat_s, filter_reflat_n
      numfmx = imt * p25
      do jrow=1,jmt
        if ((jrow .le. jft1 .or. jrow .ge. jft2) .and. jrow .ge. jfrst)
     &    then
          if (phi(jrow) .gt. c0) then
            refcos = refcosn
          else
            refcos = refcoss
          endif
	  if (abs(yt(jrow)) .ge. 90.0) then
	    numflt(jrow) = 0
	  else
	    arg = refcos/cst(jrow)
            numflt(jrow) = min(max(1,int(arg)),numfmx)
	  endif
          write(stdout,9502) jrow, numflt(jrow), yt(jrow)
        else
          numflt(jrow) = 0
	endif
      enddo
      write(stdout,9503) 
      do jrow=1,jmt
        if ((jrow .le. jfu1 .or. jrow .ge. jfu2) .and. jrow .ge. jfrst)
     &    then
          if (phi(jrow) .gt. c0) then
            refcos = refcosn
          else
            refcos = refcoss
          endif
	  if (abs(yu(jrow)) .ge. 90.0) then
	    numflu(jrow) = 0
	  else
	    arg = refcos/csu(jrow)
            numflu(jrow) = min(max(1,int(arg)),numfmx)
	  endif
          write(stdout,9502) jrow, numflu(jrow), yu(jrow)
        else
	  numflu(jrow) = 0
        endif
      enddo
c
9501  format(/' southern reference latitude =',f7.3,' deg.  '
     &, ' northern reference latitude =',f7.3,' deg.'
     &, /' T-cell jrow  numflt  latitude')
9502  format(3x,2i8,6x,f7.3)
9503  format(//' U-cell jrow  numflu  latitude')
# endif
# ifdef fourfil
c
c     this section sets up tables for filter; it must be called once
c     per execution of ocean
c
c     note: lqmsum is the sum of (im-1)/2 for im=1,imt+1
c           lhsum is the sum of im-1 for im=1,imt+1
c
      imsave = im
c
c     assemble index array
c
      do i=1,imtx8
        ind(i) = i
      enddo
c
c     calculate and save all cosines which will be needed
c
      ibase = 0
      jbase = 0
c
      do im=1,imt+1
        fimr = c1/float(im)
        imm1 = im-1
        if (imm1 .ne. 0) then
          do i=1,imm1
            denmsv(ibase+i) = c1/(c1-cos(pi*float(i)*fimr))
          enddo
        endif
        idbase(im) = ibase
        ibase = ibase + imm1
        imqc = (im-1)/2
        if (imqc .ne. 0) then
          do i=1,imqc
            cossav(jbase+i) = cos(pi*float(i)*fimr)
          enddo
        endif
        icbase(im) = jbase
        jbase = jbase + imqc
      enddo
c
c     calculate adjustments for general fourier case if im=2*n
c
      circle(1) = 0.0
      circle(2) = -1.0
      circle(3) = 0.0
      circle(4) = 1.0
      do im=1,imt
        cosnpi(im) = circle(mod(im-1,4)+1)
      enddo
# endif
c
c---------------------------------------------------------------------
c     find and print start & end indices for filtering
c---------------------------------------------------------------------
c
      do j=1,jmtfil
        do l=1,lsegf
          iszf(j,l) = 0 
          iezf(j,l) = 0 
          do k=1,km
            istf(j,l,k) = 0
            ietf(j,l,k) = 0
            isuf(j,l,k) = 0
            ieuf(j,l,k) = 0
          enddo
        enddo
      enddo 
      write (stdout,9551)
      if (lsegf.gt.11) write (stdout,9552)
c
      write (stdout,9553)
      call findex (kmt, jmtfil, km, jft1, jft2, imt, istf, ietf)
c
      write (stdout,9554)
      call findex (kmu, jmtfil, km, jfu1, jfu2, imt, isuf, ieuf)
c
      write (stdout,9555)
# if defined rigid_lid_surface_pressure || defined free_surf_on
      call findex (kmu, jmtfil, 1, jfu1, jfu2, imt, iszf, iezf)
# endif
# if defined stream_function
c
c     compute an array to indicate "interior" stream function grid cells
c
      do jrow=1,jmt
        kmz(1,jrow)   = 0
        kmz(imt,jrow) = 0
      enddo
c
      do i=1,imt
        kmz(i,1)   = 0
        kmz(i,jmt) = 0
      enddo
c
      do jrow=2,jmt-1
        do i=2,imt
          kmz(i,jrow) = min(kmu(i-1,jrow-1), kmu(i,jrow-1)
     &,                     kmu(i-1,jrow), kmu(i,jrow))
        enddo
      enddo
#  ifdef cyclic
      do jrow=1,jmt
        kmz(1,jrow) = kmz(imtm1,jrow)
      enddo
#  endif
      call findex (kmz, jmtfil, 1, jft1, jft2, imt, iszf, iezf)
# endif
#endif
#ifdef meridional_tracer_budget
c
c-----------------------------------------------------------------------
c     set basins and initialize arrays
c-----------------------------------------------------------------------
c
      if (tmbint .ge. c0) then
        allocate (msktmb(imt,jmt))
c
c       set all points to signify basin # 1
c
        do jrow=1,jmt
          do i=2,imtm1
            msktmb(i,jrow) = 1
          enddo
          msktmb(1,jrow)   = 0
          msktmb(imt,jrow) = 0
        enddo
c
c       divide "msktmb" into other basins (2 .. ntmbb) here if desired.
c
c
c       set all land points to basin # 0
c
        do jrow=1,jmt
          do i=1,imt
            if (kmt(i,jrow) .eq. 0) msktmb(i,jrow) = 0
          enddo
        enddo
c
        if (ntmbb .gt. 1) then
          write (stdout,*) 
     &  ' Basin arrangement for Meridional Tracer Budget diagnostic'
          call iplot (msktmb, imt, imt, jmt)
        endif
c
        do jrow=1,jmt
          do i=1,imt
	    kmz(i,jrow) = msktmb(i,jrow)
	  enddo
        enddo
        deallocate (msktmb)
        allocate (msktmb(imt,jstask:jetask))
        do jrow=jstask,jetask
          do i=1,imt
            msktmb(i,jrow) = kmz(i,jrow)
          enddo
        enddo
      endif
#endif
c
c---------------------------------------------------------------------
c     print the timestep multipliers
c---------------------------------------------------------------------
c
      write (stdout,9601) (dtxcel(k),k=1,km)
c
c-----------------------------------------------------------------------
c     de-allocate existing 2-D arrays and re-allocate on each processor
c     for only that portion of the domain covered by the processor
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        do i=1,imt
	  kmz(i,jrow) = kmu(i,jrow)
	enddo
      enddo
      deallocate (kmu)
      allocate (kmu(imt,jstask:jetask))
      do jrow=jstask,jetask
        do i=1,imt
          kmu(i,jrow) = kmz(i,jrow)
        enddo
      enddo
c
      do jrow=1,jmt
        do i=1,imt
	  kmz(i,jrow) = map(i,jrow)
	enddo
      enddo
      deallocate (map)
      allocate (map(imt,jstask:jetask))
      do jrow=jstask,jetask
        do i=1,imt
          map(i,jrow) = kmz(i,jrow)
        enddo
      enddo
c
      do jrow=1,jmt
        do i=1,imt
	  kmz(i,jrow) = kmt(i,jrow)
	enddo
      enddo
      deallocate (kmt)
#ifdef partial_cell
      if(jetask+1 .GT. jmt) then
        allocate (kmt(imt,jstask:jetask))
        do jrow=jstask,jetask
          do i=1,imt
            kmt(i,jrow) = kmz(i,jrow)
          enddo
        enddo
      else
        allocate (kmt(imt,jstask:jetask+1))
        do jrow=jstask,jetask+1
          do i=1,imt
            kmt(i,jrow) = kmz(i,jrow)
          enddo
        enddo
      endif
#else
      allocate (kmt(imt,jstask:jetask))
      do jrow=jstask,jetask
        do i=1,imt
          kmt(i,jrow) = kmz(i,jrow)
        enddo
      enddo
#endif
c
#ifdef partial_cell
c
c     do not de-allocate because hup is used in read_restart
c
c      deallocate (hup)
#endif
c
      allocate (a2d(imt,jmt))
c
#ifdef partial_cell
      do jrow=1,jmt
        do i=1,imt
	  a2d(i,jrow) = htp(i,jrow)
	enddo
      enddo
      deallocate (htp)
      if(jetask+1 .GT. jmt) then
        allocate (htp(imt,jstask:jetask))
        do jrow=jstask,jetask
          do i=1,imt
            htp(i,jrow) = a2d(i,jrow)
          enddo
        enddo
      else
        allocate (htp(imt,jstask:jetask+1))
        do jrow=jstask,jetask+1
          do i=1,imt
            htp(i,jrow) = a2d(i,jrow)
          enddo
        enddo
      endif
#endif
c
      do jrow=1,jmt
        do i=1,imt
	  a2d(i,jrow) = hr(i,jrow)
	enddo
      enddo
      deallocate (hr)
      allocate (hr(imt,jstask:jetask))
      do jrow=jstask,jetask
        do i=1,imt
          hr(i,jrow) = a2d(i,jrow)
        enddo
      enddo
c
      do jrow=1,jmt
        do i=1,imt
	  a2d(i,jrow) = h(i,jrow)
	enddo
      enddo
      deallocate (h)
      allocate (h(imt,jstask:jetask))
      do jrow=jstask,jetask
        do i=1,imt
          h(i,jrow) = a2d(i,jrow)
        enddo
      enddo
c
      deallocate (a2d)
c
#if !defined explicit_free_surface 
      allocate (ptd(imt,jstask:jetask))
#endif
#ifdef stream_function
      allocate (psi(imt,jstask:jetask,2))
#endif
#ifdef neptune
      allocate (unep(imt,jstask:jetask,2))
      allocate (pnep(imt,jstask:jetask))
#endif
#ifdef kppvmix
      allocate (hblt(imt,jstask:jetask))
      allocate (hmxl(imt,jstask:jetask))
#endif
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
      allocate (ps(imt,jstask:jetask,2))
      allocate (uhat(imt,jstask:jetask,2))
      allocate (ubar(imt,jstask:jetask,2))
      allocate (ubarm1(imt,jstask:jetask,2))
      allocate (pguess(imt,jstask:jetask))
      allocate (divf(imt,jstask:jetask))
#endif
#if defined explicit_free_surface
      jbeg = jstask-1
      jend = jetask+1
      if(jstask .EQ. 1  ) jbeg = jstask
      if(jetask .EQ. jmt) jend = jetask
      allocate (etat(imt,jbeg:jend,3))
      allocate (etau(imt,jstask:jetask))
      allocate (ps(imt,jstask:jetask,2))
      allocate (depthu(imt,jstask:jetask))
      allocate (depthur(imt,jstask:jetask))
      allocate (thicktau(imt,jstask:jetask))
      allocate (emsk(imt,jstask:jetask))
      allocate (umsk(imt,jstask:jetask))
      allocate (ubar(imt,jstask:jetask,2))
      allocate (ubarm1(imt,jstask:jetask,2))
      allocate (convU(imt,jstask:jetask))
      allocate (rhosurf(imt,jstask:jetask))
# if defined time_step_monitor && defined explicit_free_surface_nonlinear      
      allocate (surftracer(imt,jstask:jetask,nt,2))
# endif
# if defined snapshots || defined time_averages
      allocate (psiu(imt,jmt))      
      allocate (psiv(imt,jstask:jetask))
# endif
# ifdef explicit_fresh_water_flux
      allocate (sfft(imt,jstask:jetask))
      allocate (sffu(imt,jstask:jetask))
#  ifdef need_sfft_expl
      allocate (sfft_expl(imt,jstask:jetask))
#  endif
#  ifndef simple_sbc 
      allocate (strf(imt,jstask:jetask,nt))
#  endif
# endif
# ifdef restnosp
      allocate (ptd(imt,jstask:jetask))
# endif
#endif
      
#if defined max_window && defined free_surf_on
      ubar(:,:,:)   = c0
      ubarm1(:,:,:) = c0
#endif

#if defined new_memory_window
c
c-----------------------------------------------------------------------
c     allocate arrays within memory window on each processor for
c     3-D baroclinic and tracer equations
c-----------------------------------------------------------------------
c
      allocate ( u(imt,km,jmw,2,-1:1) )
      allocate ( t(imt,km,jmw,nt,-1:1) )
      allocate ( adv_vet(imt,km,2:jmw) )
      allocate ( adv_vnt(imt,km,1:jmw) )
      allocate ( adv_veu(imt,km,2:jmw-1) )
      allocate ( adv_vnu(imt,km,1:jmw-1) )
      allocate ( adv_vbt(imt,0:km,2:jmw) )
      allocate ( grad_p(imt,km,1:jmw-1,2) )
# ifdef trajectories
      allocate ( adv_vbu(imt,0:km,1:jmw-1) )
# else
      allocate ( adv_vbu(imt,0:km,2:jmw-1) )
# endif
      allocate ( rho(imt,km,1:jmw) )
      allocate ( tmask(imt,km,1:jmw) )
      allocate ( umask(imt,km,1:jmw) )
      allocate ( adv_fe(imt,km,2:jmw-1) )
      allocate ( adv_fn(imt,km,1:jmw-1) )
      allocate ( adv_fb(imt,0:km,2:jmw-1) )
      allocate ( diff_fe(imt,km,2:jmw-1) )
      allocate ( diff_fn(imt,km,1:jmw-1) )
      allocate ( diff_fb(imt,0:km,2:jmw-1) )
# ifdef air_pressure
      allocate ( smf(imt,1:jmw,3) )
# else
      allocate ( smf(imt,1:jmw,2) )
# endif
      allocate ( bmf(imt,1:jmw,2) )      
      allocate ( stf(imt,1:jmw,nt) )
      allocate ( stf_turb(imt,1:jmw,nt) )
      allocate ( btf(imt,1:jmw,nt) )
      allocate ( cstdxtr(imt,jmw) )
      allocate ( cstdxur(imt,jmw) )
      allocate ( cstdxt2r(imt,jmw) )
      allocate ( ah_cstdxur(imt,jmw) )
      allocate ( csudxtr(imt,jmw) )
      allocate ( csudxur(imt,jmw) )
      allocate ( csudxu2r(imt,jmw) )
      allocate ( am_csudxtr(imt,jmw) )
      allocate ( advmet(jmt,2) )
# ifdef pressure_gradient_average
      allocate ( rhotaum1(imt,km,1:jmw) )
      allocate ( rhotaup1(imt,km,1:jmw) )
      allocate ( rhotilde(imt,km,1:jmw) )
# endif
# ifdef partial_cell
      allocate ( dht(imt,km,jmw+1) )
      allocate ( dhu(imt,km,jmw) )
      allocate ( dhwt(imt,0:km,jmw) )
      allocate ( pc_sink(imt,km,2:jmw-1) )
# endif
# if defined fourth_order_tracer_advection || defined quicker
      allocate ( adv_f4n(imt,km,1:jmw-1,nt) )
# endif
# if defined velocity_horz_biharmonic
      allocate ( del2_vel(imt,km,2:jmw-1,2) )
# endif
# if defined tracer_horz_biharmonic
      allocate ( del2_tracer(imt,km,2:jmw-1,nt) )
# endif
# if defined isoneutralmix
      allocate ( diff_fbiso(imt,0:km,2:jmw-1) )
# endif
# ifdef source_term
      allocate ( source(imt,km,2:jmw-1) )
# endif
# if defined implicitvmix || defined isoneutralmix
      allocate ( zzi(imt,km,2:jmw-1) )
# endif
# if defined linearized_advection || defined equatorial_thermocline
      allocate ( tbarz(km,nt) )
# endif
# ifdef tracer_advection_fct
      allocate ( anti_fe(imt,km,2:jmw-1,nt) )
      allocate ( anti_fn(imt,km,1:jmw-1,nt) )
      allocate ( anti_fb(imt,0:km,2:jmw-1,nt) )
      allocate ( R_plusY(imt,km,1:jmw-1,nt) )
      allocate ( R_minusY(imt,km,1:jmw-1,nt) )
#  ifdef fct_3d
      allocate ( R_plus3(imt,km,1:jmw-1,nt) )
      allocate ( R_minus3(imt,km,1:jmw-1,nt) )
#  endif
#  ifdef fct_dlm2
      allocate ( t_lo(imt,km,2:jmw-1,nt) )
#  endif
# endif
#endif

#ifdef pressure_depth
      rho(:,:,:) = c0
#endif
#ifdef explicit_free_surface
c
c     initialize some explicit free surface fields 
c
      call init_surface
c
#endif
c
c-----------------------------------------------------------------------
c     set prognostic quantities to either initial conditions or restart
c-----------------------------------------------------------------------
c
      if (init) then
#ifdef neptune
c
c       calculate neptune velocities
c
        call neptune (jstask, jetask)
#endif
c
c       initialize two dimensional fields on disk
c
#if !defined explicit_free_surface
        call surf (jstask, jetask)
#endif
c
c       initialize all latitude rows
c
        call rowi (jstask, jetask, kmz)
c
c       initialize time step counter = 0
c
        itt    = 0
	irstdy = 0
	msrsdy = 0
c
      else
c
c       initialize all prognostic quantities from the restart
c
#if !defined explicit_free_surface
        call surf (jstask, jetask)
#endif
        call read_restart (jstask, jetask)
c
#ifdef neptune
c
c       calculate neptune velocities
c
        call neptune (jstask, jetask)
#endif
c
c       compute a topography checksum
c
        cksum = 0.0
        do jrow=jstask+1,jetask-1
          do i=1,imt
	    cksum = cksum + i*jrow*kmt(i,jrow)
	  enddo
        enddo
#ifdef parallel_1d
        call mpp_sum (cksum, 1)
#endif
        write (stdout,*) ' "kmt" checksum = ', cksum
      endif
#if defined ramdrive && defined parallel_1d
      call mpp_update_domains (ramdisk3d, mom_ramdisk_domain)
      if (jstask .eq. 1  ) ramdisk3d(:,1,:) = 0.
      if (jetask .eq. jmt) ramdisk3d(:,jmt,:) = 0.
#endif
c
#if defined max_window && defined parallel_1d
      ptr_t = loc(t(1,1,1,1,tau))
      call mpp_update_domains( tmp_var_t, mom_prognostic_domain )
      ptr_t = loc(t(1,1,1,1,taup1))
      call mpp_update_domains( tmp_var_t, mom_prognostic_domain )
      ptr_u = loc(u(1,1,1,1,tau))
      call mpp_update_domains( tmp_var_u, mom_prognostic_domain )
      ptr_u = loc(u(1,1,1,1,taup1))
      call mpp_update_domains( tmp_var_u, mom_prognostic_domain )
#endif
c
#if defined stream_function
c
c-----------------------------------------------------------------------
c     checksum the starting stream function.
c-----------------------------------------------------------------------
c
      call print_checksum (psi(1,1,1), imt, jmt,' cksum psi(,,1) = ')
      call print_checksum (psi(1,1,2), imt, jmt,' cksum psi(,,2) = ')
#endif
c
c-----------------------------------------------------------------------
c     initialize the time manager with specified initial conditions
c     time, user reference time, model time, and how long to integrate.
c-----------------------------------------------------------------------
c
      call tmngri (year0, month0, day0, hour0, min0, sec0
     &,            ryear, rmonth, rday, rhour, rmin, rsec
     &,            irstdy, msrsdy
     &,            runlen, rununits, rundays, dtts)
c
c-----------------------------------------------------------------------
c     initialize various things
c-----------------------------------------------------------------------
c
#if defined sponges && !defined equatorial_thermocline
c
c     initialize sponge zones adjacent to northern and southern walls
c
      call spongei
#endif
#if defined obc
      call setobc
#endif
# ifdef posdef
        call posdefi
# endif
c
      allocate (zu(imt,jstask:jetask,2))
#ifdef stream_function
      allocate (ztd(imt,jstask:jetask))
#endif
      do jrow=jstask,jetask
        do i=1,imt
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
	  uhat(i,jrow,1) = c0
	  uhat(i,jrow,2) = c0
	  divf(i,jrow) = c0
#endif
#ifdef stream_function
          ztd(i,jrow) = c0
#endif
          zu(i,jrow,1)  = c0
          zu(i,jrow,2)  = c0
        enddo
      enddo
#if defined explicit_free_surface
# ifdef explhmix
c
c     calculate arrays used in the horizontal diffusion terms
c
      allocate (bhd(jstask:jetask))
      allocate (chd(jstask:jetask))
      allocate (dhd(jstask:jetask))
      allocate (ghd(jstask:jetask))
      allocate (hhd(jstask:jetask))
c
      do jrow = jstask,jetask 
         jp1       = min(jetask,jrow+1)
	 bhd(jrow) = csur(jrow)**2
	 chd(jrow) = dyur(jrow)*csur(jrow)*cst(jp1)
	 dhd(jrow) = dyur(jrow)*csur(jrow)*cst(jrow)
	 ghd(jrow) = (c1-tng(jrow)**2)/(radius**2)
	 hhd(jrow) = sine(jrow)/(radius*csu(jrow)**2)
      enddo
# endif
#endif
c
c     coriolis factors
c
#ifdef rot_grid
c
      allocate (cori(imt,jstask:jetask,2))
c
c     read geographic latitude from a file
c
      write (stdout,'(/,(1x,a))')
     & '==> Warning: when using "rot_grid", rotation of datasets      '
     &,'             must be done by the researcher following the     '
     &,'             examples in rotation.F. After this has been done,'
     &,'             remove the stop in setocn.F to proceed.          '
      if (imt .eq. imt) call abort()
c
      write (stdout,'(//a,a/)')
     & ' =>Reading velocity grid latitude from file "lat.dta" '
      call getunit (iolat, 'lat.dta'
     &,               'unformatted sequential rewind ieee')
      read (iolat) imax, jmax
      if (imax .ne. imt .or. jmax .ne. jmt) then
        write (stdout,*)
     &  '=>Error: imt from latitude file is ',imax, ' MOM expects ',imt
     &, '         jmt from latitude file is ',jmax, ' MOM expects ',jmt
        call abort()
      else
c
c       use array cori to temporarily hold latitude in order to save mem
c
        allocate (corbuf(imt,jmt))
	read (iolat) ((corbuf(i,jrow),i=1,imt),jrow=1,jmt)
	do jrow=jstask,jetask
	  do i=1,imt
            cori(i,jrow,1) = c2*omega*sin(corbuf(i,jrow)/radian)
	    cori(i,jrow,2) = -cori(i,jrow,1)
	  enddo
	enddo
        deallocate (corbuf)
      endif
      call relunit (iolat)
#else
c
      allocate (cori(jmt,2))
c
# if defined beta_plane || defined f_plane
c
c     beta plane with f = f0 + beta*y where f0 is at yu(1)
c     if f_plane then beta = 0
c
      degtcm = radius/radian
#  ifdef f_plane
      dfdy   = c0
#  else
      dfdy   = c2*omega*csu(1)/radius
#  endif
      do jrow=1,jmt
        y            = (yu(jrow)-yu(1))*degtcm
        cori(jrow,1) = c2*omega*sine(1)  + y*dfdy
        cori(jrow,2) = -c2*omega*sine(1) - y*dfdy
      enddo
#  ifdef f_plane
      write (stdout,'(//,a,1pg14.7,a,g14.7,a,g14.7//)')
     & ' "f plane" approximation set: f0 = ', cori(1,1)
     &, ' at latitude (deg) =', yu(1), ' with beta =',dfdy
#  else
      write (stdout,'(//,a,1pg14.7,a,g14.7,a,g14.7//)')
     & ' "beta plane" approximation set: f0 = ', cori(1,1)
     &, ' at latitude (deg) =', yu(1), ' with beta =',dfdy
#  endif
# else
      do jrow=1,jmt
        cori(jrow,1) = c2*omega*sine(jrow)
        cori(jrow,2) = -c2*omega*sine(jrow)
      enddo
# endif
#endif
c
c     metric advection factors
c
      do jrow=1,jmt
#if defined beta_plane || defined f_plane
c
c       set metric factors = 0 on beta and f plane
c
        advmet(jrow,1) = c0
	advmet(jrow,2) = c0
#else
        advmet(jrow,1) = tng(jrow)/radius
	advmet(jrow,2) = -advmet(jrow,1)
#endif
      enddo
c
c     diffusive flux through bottom of cells
c
      do j=2,jmw-1
        do k=0,km
	  do i=1,imt
	    diff_fb(i,k,j) = c0
	    adv_fb(i,k,j)  = c0
	  enddo
	enddo
      enddo	
c
#if defined bryan_lewis_vertical || defined bryan_lewis_horizontal
c
c-----------------------------------------------------------------------
c     initialize hybrid scheme:Bryan_Lewis tracer diffusion coefficients
c-----------------------------------------------------------------------
c
      call blmixi
#endif
c
c-----------------------------------------------------------------------
c     initialize horizontal mixing coefficients
c-----------------------------------------------------------------------
c
      call hmixci (error, cifdef, ifdmax, nifdef
     &,            htr_mixset, hvel_mixset)
c
#ifdef ppvmix
c
c-----------------------------------------------------------------------
c     initialize pacanowski-philander vertical mixing scheme
c-----------------------------------------------------------------------
c
      call ppmixi (error, cifdef, ifdmax, nifdef, vmixset)
#endif
#ifdef kppvmix
c
c-----------------------------------------------------------------------
c     initialize kpp vertical mixing scheme
c-----------------------------------------------------------------------
c
      call kppmixi (error, cifdef, ifdmax, nifdef, vmixset
     &,             jstask, jetask)
#endif
#ifdef isoneutralmix
c
c-----------------------------------------------------------------------
c    initialize hybrid scheme: isoneutral mixing
c-----------------------------------------------------------------------
c
      call isopi (error, am, ah)
c
# ifdef nonconst_diffusivity
c
c-----------------------------------------------------------------------
c     initialize nonconstant isoneutral diffusivity schemes 
c-----------------------------------------------------------------------
c
      call ncdiffi (error)
c
# endif
#endif
c
c-----------------------------------------------------------------------
c     initialize diagnostics
c-----------------------------------------------------------------------
c
#ifdef meridional_tracer_budget
      numtmb = 0
      do mask=0,ntmbb
        do jrow=1,jmt
          smdvol(jrow,mask)  = c0
        enddo
      enddo
      do mask=0,ntmbb
        do n=1,nt
          do jrow=1,jmt
            tstor(jrow,n,mask) = c0
            tdiv(jrow,n,mask)  = c0
            tflux(jrow,n,mask) = c0
            tdif(jrow,n,mask)  = c0
            tsorc(jrow,n,mask) = c0
          enddo
        enddo
      enddo
#endif
c
#ifdef trajectories
c
c-----------------------------------------------------------------------
c     initialize particle trajectories
c-----------------------------------------------------------------------
c
      if (initpt) call ptraji
        
#endif
#ifdef save_xbts
c
c-----------------------------------------------------------------------
c     initialize XBT locations and averaging arrays 
c-----------------------------------------------------------------------
c
      call xbti
c
#endif
#ifdef xlandmix
c
c-----------------------------------------------------------------------
c     initialize crossland mixing scheme
c-----------------------------------------------------------------------
c
      call xland1st(kmz)
c
      do nxl=1,nxland
        call xlandvchk (nxl)
      end do
#endif
c
c-----------------------------------------------------------------------
c     do all consistency checks last
c-----------------------------------------------------------------------
c
      call checks (error, vmixset, htr_mixset, hvel_mixset)
c
c-----------------------------------------------------------------------
c     list options enabled for this model run
c-----------------------------------------------------------------------
c
      call options (cifdef, ifdmax, nifdef)
c
#ifdef timing
      call toc ('driver', 'setocn')
#endif
c
c-----------------------------------------------------------------------
c     allocate storage for remaining 2D arrays 
c-----------------------------------------------------------------------
c
      allocate (res(imt,jstask:jetask))
#if !defined explicit_free_surface
      allocate (cf(imt,jstask:jetask,-1:1,-1:1))
#endif
c
      deallocate (kmz)
c
      return
c
9341  format (//,'  Global ocean statistics:'
     &,/,' the total ocean surface area (t cells)  =',1pe15.8,'cm**2'
     &,/,' the total resting ocean volume (t cells)=',1pe15.8,'cm**3'
     &,/,' the total ocean surface area (u cells)  =',1pe15.8,'cm**2'
     &,/,' the total resting ocean volume (u cells)=',1pe15.8,'cm**3',/)
9401  format(/t50,'number of levels on "u,v" grid')
9402  format(/t50,'number of levels on "t,s" grid')
9601  format(/,' "dtxcel(km)" tracer timestep multipliers:',/,10(f9.3))
      end



      subroutine calc_kmu (kmt, htp, imt, kmu, hup, js, je, jmt)
c
c=======================================================================
c     calculate depth arrays associated with U cells.
c
c     input:
c       kmt = number of oecan T cells from surface to bottom of ocean
#ifdef partial_cell
c       htp = depth to the base of partial T cells
#endif
c       imt = longitudinal dimension of arrays
c       jmt = latitudinal dimension of arrays
c       zw  = depth to bottom of T cells
c       km  = max number of depths
c
c     output:
c       kmu = number of ocean U cells from surface to bottom of ocean
#ifdef partial_cell
c       hup = depth to the base of partial U cells
#endif
c
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
      use mom3_mpp_mod
      dimension kmt(imt,js:je), kmu(imt,js:je)
#ifdef partial_cell
      dimension htp(imt,js:je), hup(imt,js:je)
#endif
      real, dimension(:,:), allocatable :: abuf
# include "taskrows.h"
c
c-----------------------------------------------------------------------
c     set some constants
c-----------------------------------------------------------------------
c
      c0 = 0.0
      c1 = 1.0
c
c-----------------------------------------------------------------------
c     compute number of vertical levels on the U grid
c-----------------------------------------------------------------------
c
      do jrow=js,je
        kmu(imt,jrow) = 0
#ifdef partial_cell
        hup(imt,jrow) = 0
#endif
      enddo
c
      do i=1,imt
        kmu(i,je) = 0
#ifdef partial_cell
        hup(i,je) = 0
#endif
      enddo
c
      do jrow=js,je-1
        do i=1,imt-1
           kmu(i,jrow) = min (kmt(i,jrow), kmt(i+1,jrow), kmt(i,jrow+1)
     &,                       kmt(i+1,jrow+1))
#ifdef partial_cell
           hup(i,jrow) = min (htp(i,jrow), htp(i+1,jrow), htp(i,jrow+1)
     &,                       htp(i+1,jrow+1))
#endif
        enddo
      enddo
#ifdef obc_north
      if (je .eq. jmt) then
        do i=1,imt
          kmu(i,jmt) = kmu(i,jmt-1)
# ifdef partial_cell
          hup(i,jmt) = hup(i,jmt-1)
# endif
        enddo
      endif
#endif
#ifdef obc_east
      do jrow=js,je
        kmu(imt,jrow) = kmu(imt-1,jrow)
# ifdef partial_cell
        hup(imt,jrow) = hup(imt-1,jrow)
# endif
      enddo
#endif
#ifdef cyclic
      do jrow=js,je
        kmu(imt,jrow) = kmu(2,jrow)
# ifdef partial_cell
        hup(imt,jrow) = hup(2,jrow)
# endif
      enddo
#endif
#ifdef symmetry
      if (je .eq. jmt) then
        do i=1,imt
          kmu(i,jmt) = kmu(i,jmt-2)
# ifdef partial_cell
          hup(i,jmt) = hup(i,jmt-2)
# endif
        enddo
      endif
#endif
#ifdef parallel_1d
      if (js .ne. 1 .or. je .ne. jmt) then
        allocate (abuf(imt,js:je))
	do j=js,je
	  do i=1,imt
	    abuf(i,j) = kmu(i,j)
	  enddo
	enddo
        call mpp_update_domains( abuf, mom_domain )
	do j=js,je
	  do i=1,imt
	    kmu(i,j) = nint(abuf(i,j))
	  enddo
	enddo
c
# ifdef partial_cell
	do j=js,je
	  do i=1,imt
	    abuf(i,j) = hup(i,j)
	  enddo
	enddo
        call mpp_update_domains( abuf, mom_domain )
	do j=js,je
	  do i=1,imt
	    hup(i,j) = nint(abuf(i,j))
	  enddo
	enddo
# endif
	deallocate (abuf)
      endif
#endif
c
      return
      end


      subroutine depth_u (kmu,hup,imt_local,zw,km_local,h,hr,js,je)
c
c=======================================================================
c     calculate depth arrays associated with U cells.
c
c     input:
c       kmu = number of oecan U cells from surface to bottom of ocean
c       imt_local = longitudinal dimension of arrays
c       jmt = latitudinal dimension of arrays
c       zw  = depth to bottom of T cells
c       km_local  = max number of depths
#ifdef partial_cell
c       hup = depth to bottom of partial U cells not including BBL
#endif
c
c     output:
c       h   = depth to ocean floor over U cells (cm)
#ifdef bbl_ag
c             including BBL thickness
#endif
c       hr  = reciprocal of "h"       
c
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
      dimension kmu(imt_local,js:je)
      dimension h(imt_local,js:je), hr(imt_local,js:je)
# ifdef partial_cell
      dimension hup(imt_local,js:je)
# endif
      dimension zw(km_local)
#include "stdunits.h"
#include "taskrows.h"
#if defined bbl_ag && defined bbl_mask
# include "size.h"
# include "bbl.h"
#endif
c
c-----------------------------------------------------------------------
c     set some constants
c-----------------------------------------------------------------------
c
      c0 = 0.0
      c1 = 1.0
c
c---------------------------------------------------------------------
c     compute depths and reciprocal depths over U cells
c---------------------------------------------------------------------
c
      do j=js,je
        do i=1,imt_local
	  if (kmu(i,j) .lt. 0 .or. kmu(i,j) .gt. km_local) then
            write (stdout,*)
     &      'depth_u pn=',pn,' i=',i,' j=',j,' kmu=',kmu(i,j)
            call abort()
	  endif
	enddo
      enddo
c
#ifdef bbl_ag
      dzbbl = zw(km_local) - zw(km_local-1)
#else 
      dzbbl = c0
#endif
c
      do jrow=js,je
        do i=1,imt_local
          hr(i,jrow) = c0
          h(i,jrow)  = c0
c
          if (kmu(i,jrow) .ne. 0) then
            h (i,jrow) =  dzbbl
#ifdef bbl_mask
     &           *mask_u_bbl(i,jrow)
#endif
#if defined partial_cell
     &                         + hup(i,jrow)
#else
     &                         + zw(kmu(i,jrow))
#endif
            hr(i,jrow) = c1/h(i,jrow)
          endif
        enddo
      enddo
      return
      end






      subroutine options (cifdef, ifdmax, nifdef)
      character*(*) cifdef(ifdmax)
#include "derived_options.h"
#include "stdunits.h"
c
c-----------------------------------------------------------------------
c     list "ifdef" options that are enabled (turned on)
c-----------------------------------------------------------------------
c
      logical error
c
      error = .false.
c
        write (stdout,'(/,10x,a,/)')
     &  ' For available options ---- consult the manual'
c
        write (stdout,'(/,a)')
     &   ' THE FOLLOWING "IFDEF" OPTIONS HAVE BEEN ENABLED (TURNED ON):'
#ifdef smooth_topo
        nifdef = nifdef + 1
        cifdef(nifdef) = 'smooth_topo'
#endif
#ifdef isotropic_grid
        nifdef = nifdef + 1
        cifdef(nifdef) = 'isotropic_grid'
#endif
#ifdef extend_isotropic_grid
        nifdef = nifdef + 1
        cifdef(nifdef) = 'extend_isotropic_grid'
#endif
#ifdef parallel_1d
        nifdef = nifdef + 1
        cifdef(nifdef) = 'parallel_1d'
#endif
#ifdef use_libSMA
        nifdef = nifdef + 1
        cifdef(nifdef) = 'use_libSMA'
#endif
#ifdef use_libMPI
        nifdef = nifdef + 1
        cifdef(nifdef) = 'use_libMPI'
#endif
#ifdef partial_cell
        nifdef = nifdef + 1
        cifdef(nifdef) = 'partial_cell'
#endif
#ifdef bbl_ag
        nifdef = nifdef + 1
        cifdef(nifdef) = 'bbl_ag'
#endif
#ifdef coupled
        nifdef = nifdef + 1
        cifdef(nifdef) = 'coupled     '
#endif
#ifdef simple_sbc
        nifdef = nifdef + 1
        cifdef(nifdef) = 'simple_sbc'
#endif
#ifdef no_sbc
        nifdef = nifdef + 1
        cifdef(nifdef) = 'no_sbc'
#endif
#ifdef time_mean_sbc_data
        nifdef = nifdef + 1
        cifdef(nifdef) = 'time_mean_sbc_data '
#endif
#ifdef time_varying_sbc_data
        nifdef = nifdef + 1
        cifdef(nifdef) = 'time_varying_sbc_data '
#endif
#ifdef minimize_sbc_memory
        nifdef = nifdef + 1
        cifdef(nifdef) = 'minimize_sbc_memory '
#endif
#ifdef trace_coupled_fluxes
        nifdef = nifdef + 1
        cifdef(nifdef) = 'trace_coupled_fluxes'
#endif
#ifdef linearized_advection
        nifdef = nifdef + 1
        cifdef(nifdef) = 'linearized_advection'
#endif
#ifdef equatorial_thermocline
        nifdef = nifdef + 1
        cifdef(nifdef) = 'equatorial_thermocline'
#endif
#ifdef constant_taux
        nifdef = nifdef + 1
        cifdef(nifdef) = 'constant_taux'
#endif
#ifdef constant_tauy
        nifdef = nifdef + 1
        cifdef(nifdef) = 'constant_tauy'
#endif
#ifdef linearized_density
        nifdef = nifdef + 1
        cifdef(nifdef) = 'linearized_density'
#endif
#ifdef linear_tstar
        nifdef = nifdef + 1
        cifdef(nifdef) = 'linear_tstar'
#endif
#ifdef analytic_zonal_winds
        nifdef = nifdef + 1
        cifdef(nifdef) = 'analytic_zonal_winds'
#endif
#ifdef timing
        nifdef = nifdef + 1
        cifdef(nifdef) = 'timing      '
#endif
#ifdef knudsen
        nifdef = nifdef + 1
        cifdef(nifdef) = 'knudsen'
#endif
#ifdef f_plane
        nifdef = nifdef + 1
        cifdef(nifdef) = 'f_plane'
#endif
#ifdef beta_plane
        nifdef = nifdef + 1
        cifdef(nifdef) = 'beta_plane'
#endif
#ifdef rectangular_box
        nifdef = nifdef + 1
        cifdef(nifdef) = 'rectangular_box'
#endif
#ifdef flat_bottom
        nifdef = nifdef + 1
        cifdef(nifdef) = 'flat_bottom '
#endif
#ifdef rot_grid
        nifdef = nifdef + 1
        cifdef(nifdef) = 'rot_grid'
#endif
#ifdef scripps_kmt
        nifdef = nifdef + 1
        cifdef(nifdef) = 'scripps_kmt '
#endif
#ifdef etopo_kmt
        nifdef = nifdef + 1
        cifdef(nifdef) = 'etopo_kmt '
#endif
#ifdef read_my_kmt
        nifdef = nifdef + 1
        cifdef(nifdef) = 'read_my_kmt '
#endif
#ifdef idealized_kmt
        nifdef = nifdef + 1
        cifdef(nifdef) = 'idealized_kmt '
#endif
#ifdef gaussian_kmt
        nifdef = nifdef + 1
        cifdef(nifdef) = 'gaussian_kmt '
#endif
#ifdef generate_a_grid
        nifdef = nifdef + 1
        cifdef(nifdef) = 'generate_a_grid'
#endif
#ifdef read_my_grid
        nifdef = nifdef + 1
        cifdef(nifdef) = 'read_my_grid'
#endif
#ifdef equivalence_mw
        nifdef = nifdef + 1
        cifdef(nifdef) = 'equivalence_mw'
#endif
#ifdef YMP
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cray_ymp'
#endif
#ifdef _CRAYC90
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cray_c90'
#endif
#ifdef _CRAYT90
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cray_t90'
#endif
#ifdef _CRAYT3E
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cray_t3e'
#endif
#ifdef sgi_mipspro
        nifdef = nifdef + 1
        cifdef(nifdef) = 'sgi'
#endif
#ifdef rs6000
        nifdef = nifdef + 1
        cifdef(nifdef) = 'rs6000'
#endif
#ifdef solid_walls
        nifdef = nifdef + 1
        cifdef(nifdef) = 'solid_walls'
#endif
#ifdef cyclic
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cyclic      '
#endif
#ifdef tracer_advection_fct
        nifdef = nifdef + 1
        cifdef(nifdef) = 'tracer_advection_fct'
# ifdef fct_3d
        nifdef = nifdef + 1
        cifdef(nifdef) = 'fct_3d'
# endif
# if defined fct_dlm1
        nifdef = nifdef + 1
        cifdef(nifdef) = 'fct_dlm1'
# endif
# if defined fct_dlm2
        nifdef = nifdef + 1
        cifdef(nifdef) = 'fct_dlm2'
# endif
# ifdef tst_fct_los
        nifdef = nifdef + 1
        cifdef(nifdef) = 'tst_fct_los '
# endif
# ifdef tst_fct_his
        nifdef = nifdef + 1
        cifdef(nifdef) = 'tst_fct_his '
# endif
# ifdef save_fct
        nifdef = nifdef + 1
        cifdef(nifdef) = 'save_fct'
# endif
#endif
#ifdef robert_time_filter
        nifdef = nifdef + 1
        cifdef(nifdef) = 'robert_time_filter'
#endif
#ifdef damp_inertial_oscillation
        nifdef = nifdef + 1
        cifdef(nifdef) = 'damp_inertial_oscillation'
#endif
#ifdef fourth_order_tracer_advection
        nifdef = nifdef + 1
        cifdef(nifdef) = 'fourth_order_tracer_advection'
#endif
#ifdef quicker
        nifdef = nifdef + 1
        cifdef(nifdef) = 'quicker'
#endif
#ifdef bottom_upwind
        nifdef = nifdef + 1
        cifdef(nifdef) = 'bottom_upwind'
#endif
#ifdef pressure_gradient_average
        nifdef = nifdef + 1
        cifdef(nifdef) = 'pressure_gradient_average'
#endif
#ifdef obc_south
        nifdef = nifdef + 1
        cifdef(nifdef) = 'obc_south '
#endif
#ifdef obc_north
        nifdef = nifdef + 1
        cifdef(nifdef) = 'obc_north '
#endif
#ifdef obc_west
        nifdef = nifdef + 1
        cifdef(nifdef) = 'obc_west '
#endif
#ifdef obc_east
        nifdef = nifdef + 1
        cifdef(nifdef) = 'obc_east '
#endif
#ifdef orlanski
        nifdef = nifdef + 1
        cifdef(nifdef) = 'orlanski '
#endif
#ifdef obc_relax_eta
        nifdef = nifdef + 1
        cifdef(nifdef) = 'obc_relax_eta '
#endif
#ifdef symmetry
        nifdef = nifdef + 1
        cifdef(nifdef) = 'symmetry    '
#endif
#ifdef source_term
        nifdef = nifdef + 1
        cifdef(nifdef) = 'source_term '
#endif
#ifdef firfil
        nifdef = nifdef + 1
        cifdef(nifdef) = 'firfil      '
#endif
#ifdef fourfil
        nifdef = nifdef + 1
        cifdef(nifdef) = 'fourfil     '
#endif
#if defined implicit_free_surface
        nifdef = nifdef + 1
        cifdef(nifdef) = 'implicit_free_surface'
#endif
#if defined explicit_free_surface
	nifdef = nifdef + 1
	cifdef(nifdef) = 'explicit_free_surface'
#endif
#if defined explicit_free_surface_linear
	nifdef = nifdef + 1
	cifdef(nifdef) = 'explicit_free_surface_linear'
#endif
#if defined explicit_free_surface_nonlinear
	nifdef = nifdef + 1
	cifdef(nifdef) = 'explicit_free_surface_nonlinear'
#endif
#if defined time_averaged_surface
	nifdef = nifdef + 1
	cifdef(nifdef) = 'time_averaged_surface'
#endif
#if defined salinity_psu
	nifdef = nifdef + 1
	cifdef(nifdef) = 'salinity_psu'
#endif
#if defined explicit_eb
	nifdef = nifdef + 1
	cifdef(nifdef) = 'explicit_eb'
#endif
#if defined explicit_efb
	nifdef = nifdef + 1
	cifdef(nifdef) = 'explicit_efb'
#endif
#if defined explicit_eta_dpdc
	nifdef = nifdef + 1
	cifdef(nifdef) = 'explicit_eta_dpdc'
#endif
#if defined explicit_eta_dpdc_all
	nifdef = nifdef + 1
	cifdef(nifdef) = 'explicit_eta_dpdc_all'
#endif
#if defined explicit_eta_dpdc_once
	nifdef = nifdef + 1
	cifdef(nifdef) = 'explicit_eta_dpdc_once'
#endif
#ifdef explicit_fresh_water_flux
        nifdef = nifdef + 1
        cifdef(nifdef) = 'explicit_fresh_water_flux'
#endif
#ifdef river_inflow
        nifdef = nifdef + 1
        cifdef(nifdef) = 'river_inflow'
#endif
#ifdef explicit_eta_laplacian
        nifdef = nifdef + 1
        cifdef(nifdef) = 'explicit_eta_laplacian'
#endif
#ifdef explicit_polar_filter
        nifdef = nifdef + 1
        cifdef(nifdef) = 'explicit_polar_filter'
#endif
#if defined explhmix
        nifdef = nifdef + 1
	cifdef(nifdef) = 'explhmix'
#endif
#ifdef air_pressure
        nifdef = nifdef + 1
        cifdef(nifdef) = 'air_pressure'
#endif
#if defined rigid_lid_surface_pressure
        nifdef = nifdef + 1
        cifdef(nifdef) = 'rigid_lid_surface_pressure'
#endif
#ifdef zero_island_flow
        nifdef = nifdef + 1
        cifdef(nifdef) = 'zero_island_flow'
#endif
#ifdef stream_function
        nifdef = nifdef + 1
        cifdef(nifdef) = 'stream_function'
#endif
#ifdef sf_5_point
        nifdef = nifdef + 1
        cifdef(nifdef) = 'sf_5_point'
#endif
#ifdef sf_9_point
        nifdef = nifdef + 1
        cifdef(nifdef) = 'sf_9_point'
#endif
#ifdef implicitvmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'implicitvmix'
#endif
#ifdef oldconvect
        nifdef = nifdef + 1
        cifdef(nifdef) = 'oldconvect'
#endif
#ifdef redi_diffusion
        nifdef = nifdef + 1
        cifdef(nifdef) = 'redi_diffusion   '
#endif
#ifdef small_tensor
        nifdef = nifdef + 1
        cifdef(nifdef) = 'small_tensor   '
#endif
#ifdef full_tensor
        nifdef = nifdef + 1
        cifdef(nifdef) = 'full_tensor   '
#endif
#ifdef dm_taper
        nifdef = nifdef + 1
        cifdef(nifdef) = 'dm_taper   '
#endif
#ifdef gkw_taper
        nifdef = nifdef + 1
        cifdef(nifdef) = 'gkw_taper   '
#endif
#ifdef gent_mcwilliams
        nifdef = nifdef + 1
        cifdef(nifdef) = 'gent_mcwilliams'
#endif
#ifdef gm_skew
        nifdef = nifdef + 1
        cifdef(nifdef) = 'gm_skew'
#endif
#ifdef gm_advect
        nifdef = nifdef + 1
        cifdef(nifdef) = 'gm_advect'
#endif
#ifdef isotropic_mixed
        nifdef = nifdef + 1
        cifdef(nifdef) = 'isotropic_mixed'
#endif
#ifdef biharmonic_rm
        nifdef = nifdef + 1
        cifdef(nifdef) = 'biharmonic_rm'
#endif
#ifdef save_density_terms
        nifdef = nifdef + 1
        cifdef(nifdef) = 'save_density_terms'
#endif
#ifdef nonconst_diffusivity
        nifdef = nifdef + 1
        cifdef(nifdef) = 'nonconst_diffusivity'
#endif
#ifdef save_ncdiff
        nifdef = nifdef + 1
        cifdef(nifdef) = 'save_ncdiff'
#endif
#ifdef hl_diffusivity
        nifdef = nifdef + 1
        cifdef(nifdef) = 'hl_diffusivity'
#endif
#ifdef vmhs_diffusivity
        nifdef = nifdef + 1
        cifdef(nifdef) = 'vmhs_diffusivity'
#endif
#ifdef beta_term
        nifdef = nifdef + 1
        cifdef(nifdef) = 'beta_term'
#endif
#ifdef save_convection
        nifdef = nifdef + 1
        cifdef(nifdef) = 'save_convection'
#endif
#ifdef save_mixing_coeff
        nifdef = nifdef + 1
        cifdef(nifdef) = 'save_mixing_coeff'
#endif
#ifdef cross_flow
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cross_flow'
#endif
#ifdef save_pressure
        nifdef = nifdef + 1
        cifdef(nifdef) = 'save_pressure'
#endif
#ifdef neptune
        nifdef = nifdef + 1
        cifdef(nifdef) = 'neptune'
#endif
#ifdef velocity_horz_laplacian
        nifdef = nifdef + 1
        cifdef(nifdef) = 'velocity_horz_laplacian  '
#endif
#ifdef velocity_horz_biharmonic
        nifdef = nifdef + 1
        cifdef(nifdef) = 'velocity_horz_biharmonic  '
#endif
#ifdef velocity_horz_mix_const
        nifdef = nifdef + 1
        cifdef(nifdef) = 'velocity_horz_mix_const  '
#endif
#ifdef velocity_horz_mix_var
        nifdef = nifdef + 1
        cifdef(nifdef) = 'velocity_horz_mix_var  '
#endif
#ifdef velocity_horz_mix_smag
        nifdef = nifdef + 1
        cifdef(nifdef) = 'velocity_horz_mix_smag  '
#endif
#if defined posdef
        nifdef = nifdef + 1
        cifdef(nifdef) = 'posdef  '
#endif
#ifdef am_cosine
        nifdef = nifdef + 1
        cifdef(nifdef) = 'am_cosine '
#endif
#ifdef am_cosinep1
        nifdef = nifdef + 1
        cifdef(nifdef) = 'am_cosinep1 '
#endif
#ifdef am_taper_highlats
        nifdef = nifdef + 1
        cifdef(nifdef) = 'am_taper_highlats '
#endif
#ifdef tracer_horz_laplacian
        nifdef = nifdef + 1
        cifdef(nifdef) = 'tracer_horz_laplacian  '
#endif
#ifdef tracer_horz_biharmonic
        nifdef = nifdef + 1
        cifdef(nifdef) = 'tracer_horz_biharmonic  '
#endif
#ifdef tracer_horz_mix_const
        nifdef = nifdef + 1
        cifdef(nifdef) = 'tracer_horz_mix_const  '
#endif
#ifdef tracer_horz_mix_var
        nifdef = nifdef + 1
        cifdef(nifdef) = 'tracer_horz_mix_var  '
#endif
#ifdef tracer_horz_mix_smag
        nifdef = nifdef + 1
        cifdef(nifdef) = 'tracer_horz_mix_smag  '
#endif
#ifdef bryan_lewis_horizontal
        nifdef = nifdef + 1
        cifdef(nifdef) = 'bryan_lewis_horizontal'
#endif
#ifdef bryan_lewis_vertical
        nifdef = nifdef + 1
        cifdef(nifdef) = 'bryan_lewis_vertical'
#endif
#ifdef constvmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'constvmix   '
#endif
#ifdef shortwave
        nifdef = nifdef + 1
        cifdef(nifdef) = 'shortwave   '
#endif
#ifdef tcvmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'tcvmix      '
        write (stdout,*) ' Error ==> tcvmix not implemented yet'
	error = .true.	
#endif
#ifdef leq
        nifdef = nifdef + 1
        cifdef(nifdef) = 'leq         '
        write (stdout,*) ' Error ==> leq not fully checked out yet'
	error = .true.	
#endif
#ifdef lalg
        nifdef = nifdef + 1
        cifdef(nifdef) = 'lalg        '
        write (stdout,*) ' Error ==> lalg not fully checked out yet'
	error = .true.	
#endif
#ifdef readrmsk
        nifdef = nifdef + 1
        cifdef(nifdef) = 'readrmsk  '
#endif
#ifdef restorst
        nifdef = nifdef + 1
        cifdef(nifdef) = 'restorst    '
#endif
#ifdef sponges
# ifdef sponges_old
        nifdef = nifdef + 1
        cifdef(nifdef) = 'sponges_old'
# endif
        nifdef = nifdef + 1
        cifdef(nifdef) = 'sponges    '
#endif
#ifdef constvmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'constvmix   '
#endif
#ifdef diagnostics_lomem
        nifdef = nifdef + 1
        cifdef(nifdef) = 'diagnostics_lomem   '
#endif
#ifdef diagnostics_himem
        nifdef = nifdef + 1
        cifdef(nifdef) = 'diagnostics_himem   '
#endif
#ifdef snapshots
        nifdef = nifdef + 1
        cifdef(nifdef) = 'snapshots   '
#endif
#ifdef snapshots_regions
        nifdef = nifdef + 1
        cifdef(nifdef) = 'snapshots_regions   '
#endif
#ifdef trajectories
        nifdef = nifdef + 1
        cifdef(nifdef) = 'trajectories'
#endif
#ifdef lyapunov
        nifdef = nifdef + 1
        cifdef(nifdef) = 'lyapunov    '
#endif
#ifdef time_averages
        nifdef = nifdef + 1
        cifdef(nifdef) = 'time_averages   '
#endif
#ifdef monthly_averages
        nifdef = nifdef + 1
        cifdef(nifdef) = 'monthly_averages'
#endif
#ifdef tracer_yz
        nifdef = nifdef + 1
        cifdef(nifdef) = 'tracer_yz'
#endif
#ifdef meridional_overturning
        nifdef = nifdef + 1
        cifdef(nifdef) = 'meridional_overturning'
#endif
#ifdef meridional_overturning_density
        nifdef = nifdef + 1
        cifdef(nifdef) = 'meridional_overturning_density'
#endif
#ifdef meridional_overturning_basin
        nifdef = nifdef + 1
        cifdef(nifdef) = 'meridional_overturning_basin'
#endif
#ifdef gyre_components
        nifdef = nifdef + 1
        cifdef(nifdef) = 'gyre_components   '
#endif
#ifdef tracer_averages
        nifdef = nifdef + 1
        cifdef(nifdef) = 'tracer_averages   '
#endif
#ifdef show_zonal_mean_of_sbc
        nifdef = nifdef + 1
        cifdef(nifdef) = 'show_zonal_mean_of_sbc'
#endif
#ifdef show_external_mode
        nifdef = nifdef + 1
        cifdef(nifdef) = 'show_external_mode'
#endif
#ifdef diagnostic_surf_height
        nifdef = nifdef + 1
        cifdef(nifdef) = 'diagnostic_surf_height'
#endif
#ifdef term_balances
        nifdef = nifdef + 1
        cifdef(nifdef) = 'term_balances'
#endif
#ifdef meridional_tracer_budget
        nifdef = nifdef + 1
        cifdef(nifdef) = 'meridional_tracer_budget'
#endif
#ifdef time_step_monitor
        nifdef = nifdef + 1
        cifdef(nifdef) = 'time_step_monitor   '
#endif
#ifdef energy_analysis
        nifdef = nifdef + 1
        cifdef(nifdef) = 'energy_analysis'
#endif
#ifdef stability_tests
        nifdef = nifdef + 1
        cifdef(nifdef) = 'stability_tests'
#endif
#ifdef matrix_sections
        nifdef = nifdef + 1
        cifdef(nifdef) = 'matrix_sections'
#endif
#ifdef save_xbts
        nifdef = nifdef + 1
        cifdef(nifdef) = 'save_xbts '
#endif
#ifdef read_ic
        nifdef = nifdef + 1
        cifdef(nifdef) = 'read_ic   '
#endif
#ifdef levitus_ic
        nifdef = nifdef + 1
        cifdef(nifdef) = 'levitus_ic   '
#endif
#ifdef idealized_ic
        nifdef = nifdef + 1
        cifdef(nifdef) = 'idealized_ic'
#endif
#ifdef constant_initial_temp 
        nifdef = nifdef + 1
        cifdef(nifdef) = 'constant_initial_temp'
#endif
#ifdef ramdrive
        nifdef = nifdef + 1
        cifdef(nifdef) = 'ramdrive    '
#endif
#ifdef  crayio
        nifdef = nifdef + 1
        cifdef(nifdef) = 'crayio      '
#endif
#ifdef  ssread_sswrite
        nifdef = nifdef + 1
        cifdef(nifdef) = 'ssread_sswrite'
#endif
#ifdef fio
        nifdef = nifdef + 1
        cifdef(nifdef) = 'fio         '
#endif
#ifdef show_details
        nifdef = nifdef + 1
        cifdef(nifdef) = 'show_details'
#endif
#ifdef use_netCDF
        nifdef = nifdef + 1
        cifdef(nifdef) = 'use_netCDF'
#endif
#ifdef netcdf
        nifdef = nifdef + 1
        cifdef(nifdef) = 'netcdf'
#endif
#ifdef xlandmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'xlandmix'
#endif
#ifdef max_window 
        nifdef = nifdef + 1
        cifdef(nifdef) = 'max_window'
#endif
#ifdef pressure_depth 
        nifdef = nifdef + 1
        cifdef(nifdef) = 'pressure_depth'
#endif
c
#ifdef coupled
        write (stdout,'(/a/a/)')
     &  ' Note ==> MOM S.B.C. will come from an atmospheric model'
     &, '          with spatial interpolations for all S.B.C.'
#else
# if !defined simple_sbc
        write (stdout,'(/a/a/)')
     & ' Note ==> MOM S.B.C. will come from a prepared dataset'
     &,'          (requiring no spatial interpolations)'
# else
        write (stdout,'(/a,a/)')
     & ' Note ==> MOM S.B.C. are internally generated idealized zonal'
     &,' means.'
# endif
#endif
      if (error) then
        write (stdout,'(a)')
     & ' Regrets: sorry but the above options are not implemented yet.'
	call abort()
      endif
      nacros = 4
      nn = nifdef/float(nacros) + 1
      do i=1,nn
        ns = 1+nacros*(i-1)
        ne = ns + nacros-1
        if (ne .gt. ifdmax) then
          ne = ifdmax
          write (stdout,'(/a/)')
     &'=>Error: option list is truncated. increase parameter "ifdmax"'
          call abort()
	endif
        if (ns .le. nifdef) write (stdout,9101) (cifdef(n),n=ns,ne)
      enddo
      write (stdout,'(/a/)')  ' END OF ENABLED OPTIONS LIST'
c
      return
9101  format(4(1x,a30))
      end

#if !defined explicit_free_surface
      subroutine surf (js, je)
c
c=======================================================================
c
c     initialize prognositc barotropic quantities at "tau-1" and "tau"
c
c=======================================================================
c
      use mom3_mpp_mod
      use arrays_2d
#include "param.h"
#include "scalar.h"
#include "iounit.h"
#include "emode.h"
# if defined rigid_lid_surface_pressure || defined free_surf_on
c
c     initialize barotropic fields in memory
c
      do jrow=js,je
        do i=1,imt
          ps(i,jrow,1)     = c0
	  ps(i,jrow,2)     = c0
	  ubar(i,jrow,1)   = c0
	  ubar(i,jrow,2)   = c0
	  ubarm1(i,jrow,1) = c0
	  ubarm1(i,jrow,2) = c0
          pguess(i,jrow)   = c0
        enddo
      enddo
# endif
# ifdef stream_function
c
c     initialize stream function fields in memory
c
      do n=1,2
        do jrow=js,je
          do i=1,imt
#  ifdef neptune
c
c           initialize to neptune streamfunction
c
            psi(i,jrow,n) = pnep(i,jrow)
            ptd(i,jrow)   = c0
#  else
            psi(i,jrow,n) = c0
	    ptd(i,jrow)   = c0
#  endif
          enddo
        enddo
      enddo
c
c     initialize stream function guess fields on disk
c     block`s 1 & 2 are for the stream function guess field on disk
c
      do n=1,nkflds
        call oput (kflds, nwds, n, ptd(1,js))
      enddo
# endif

      return
      end
#endif /* !defined explicit_free_surface */


      subroutine rowi (js, je, kmt)
c
c=======================================================================
c
c     initialize prognositc quantities at "tau-1" and "tau"
c
c=======================================================================
c
      use mom3_mpp_mod
      use memory_window
#if defined sponges && !defined sponges_old
      use  xsponge_mod
#endif
#if defined read_ic
      use mpp_io_mod
      use arrays_2d, kmtz => kmt
#else
# if defined partial_cell
      use arrays_2d, only: htp
# endif
#endif
#include "param.h"
#include "coord.h"
#include "iounit.h"
#include "index.h"
#include "mw.h"
#include "taskrows.h"
#if defined firfil
      dimension tmsk(imt,km)
#endif
#if defined bbl_ag
# include "bbl.h"
#endif
      dimension kmt(imt,jmt)
#ifdef read_ic
# include "readic.h"
# if defined sgi_mipspro && defined use_libSMA
      real, dimension(imt, kbot)      :: tbuf
      real, dimension(imt, jmt )      :: hbuf
      integer(LONG_KIND)              :: error_8
      pointer(ptr1,tbuf)
      integer, save :: len1=0
      save :: ptr1

      pointer(ptr2,hbuf)
      integer, save :: len2=0
      save :: ptr2
# else
      real, allocatable               :: hbuf(:,:)
      real, allocatable               :: tbuf(:,:)
# endif  
#else
# ifdef levitus_ic
! mjh:fixes for 64-bit conversion on sgi in order to ensure compatability
! with CRAY files.  Also syymetric allocation for SHMEM on SGI 

      character*32 stampl, seqrecl
      dimension xtlev(imt), ztlev(kbot)
# endif
# ifdef obctest
      real tempsouth(km),tempnorth(km),saltsouth(km),saltnorth(km)
      data tempsouth / 10.30,  7.32,  4.74,  2.50,  0.52, -0.65 /
      data tempnorth /  4.60, -0.54, -1.27, -1.10, -1.15, -1.20 /
      data saltsouth / 35.053, 35.113, 34.985, 34.914, 34.915, 34.912 /
      data saltnorth / 34.500, 34.782, 34.873, 34.894, 34.984, 34.889 /
# endif

      dimension tbuf(imt,kbot),sbuf(imt,kbot)
      
      real*4 xtlev4(imt),ztlev4(kbot),dpml4,alat4,tbuf4(imt,kbot)
      
# if defined sgi_mipspro && defined use_libSMA

      pointer(ptr1,tbuf)
      pointer(ptr2,sbuf)
      integer, save :: len1=0, len2=0
      save :: ptr1,ptr2

      call mpp_malloc( ptr1, size(tbuf), len1 )
      call mpp_malloc( ptr2, size(sbuf), len2 )
# endif  
#endif
c
c-----------------------------------------------------------------------
c     update pointers to tau-1, tau, & tau+1 data on disk.
c     for latitude rows they point to latdisk(1) or latdisk(2)
c     for 2D fields they point to records on kflds
c-----------------------------------------------------------------------
c
      itt   = 0
      taum1disk = mod(itt+1,2) + 1
      taudisk   = mod(itt  ,2) + 1
      taup1disk = taum1disk
c
c-----------------------------------------------------------------------
c     update pointers to tau-1, tau, & tau+1 data in the MW
c-----------------------------------------------------------------------
c
#if defined max_window
c
c     rotate time levels to prevent moving data
c
      taum1 = mod(itt+0,3) - 1
      tau   = mod(itt+1,3) - 1
      taup1 = mod(itt+2,3) - 1
#endif
c
      ucksum = 0.0
      vcksum = 0.0
      tcksum = 0.0
      scksum = 0.0
#ifdef free_surf_on
      ecksum = c0
      pcksum = c0
      ubksum = c0
      vbksum = c0
#endif
c
c-----------------------------------------------------------------------
c     initialize every latitude jrow either in the MW (when wide opened)
c     or on disk (when jmw < jmt)
c-----------------------------------------------------------------------
c
#ifdef read_ic
!-----------------------------------------------------------------------
!     read initial model fields from a netCDF file prepared by user     
!-----------------------------------------------------------------------
c
c     uses read_ic_mod from readic_mod.F90           
c
c     * data are read and transferred slab by slab with buffer tbuf(imt,km)
c     * 2dim arrays are read first (if given) since code is independent of
c       max_window, and surface elevation is required for depth calculation
c     * the barotropic mode will be subtracted from input 3dim velocities
c       and be stored into 2dim ubar (overwriting former input possibly)
c            
c     authors: Martin Schmidt e-mail=> mschmidt@jule.io-warnemuende.de 
c              Torsten Seifert e-mail=> torsten.seifert@.io-warnemuende.de 
c            
c     using "NetCDF User's Guide for FORTRAN" by Russ Rew et al., 
c     Unidata Program Center, Version 3, June 1997       
c            
c-----------------------------------------------------------------------
c      PE 0 reads initial values from the netCDF ic_file 
c-----------------------------------------------------------------------
c
      call getunit 
     &        (ioic, 'namelist.icdata','formatted sequential rewind')
      read  (ioic, icdata)
      write (stdout,icdata)
      call relunit (ioic)
      call getunit 
     &        (ioic, 'document.dta','formatted sequential append')
      write (ioic,icdata)
      call relunit (ioic)
      if(nic3d.gt.nt+2) then
        write(stdout,*) '=> Error! Cannot read more than ',nt+2,
     &       'variables! Check "readic.h" and "namelist.icdata"!'
        call abort()	
      endif
      
      if(idim_ic.lt.imtm2.or.jdim_ic.lt.jmtm2.or.kdim_ic.lt.kbot) then      
        write(stdout,*) '=> Error! The model grid cannot be filled'//
     &       ' from this file. '
        write(stdout,*) 'The following file grid has been defined in '//
     &       '"readic.h" and "namelist.icdata":'
        write(stdout,*) 'idim_ic = ', idim_ic
        write(stdout,*) 'jdim_ic = ', jdim_ic
        write(stdout,*) 'kdim_ic = ', kdim_ic
        call abort()	
      endif
            
      allocate (csum(nic2d+nic3d))
      csum(:) = c0
      
#if defined use_netCDF
      ind_nc = index(ic_file,'.nc') -1 
      filename = ic_file(1:ind_nc)
#endif      
      call mpp_open(ioic, filename, action=MPP_RDONLY, 
     &       form=INPUT_FORM, threading=MPP_SINGLE, fileset=MPP_SINGLE,
     &       iospec= '-F cachea:256:100:2')
      call mpp_read_meta(ioic,filemap)
!        call open_ncfile(ic_file, filemap)
c-----------------------------------------------------------------------
c       check for variables in file
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     get 2dim ic data
c-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!     get the id of the time axis
!-----------------------------------------------------------------------
!
      if (pe .eq. 0) then 
        id = getvarid(filemap,'time')
        if (id.gt.0) then
          call mpp_read
     &       (filemap,id,varindex=ictimestep,value=timeic)
	  write(stdout,*) 'Reading timestep ', ictimestep 
        else
     	  write(stdout,*) '=> no ic time found in '
     &                    ,trim(filemap%filename)
     	  ictimestep=1
        endif
c
c-----------------------------------------------------------------------
c     get 2dim data from netCDF variables (time on axis_3)
c     do the job slab by slab using the same buffers as for 3dim fields
c-----------------------------------------------------------------------
c
        error  = .false.
        write(stdout,*) 'Searching '//trim(filemap%filename)
     &                //' for the following 2D-variables:'
        write(stdout,*) (trim(icmap2d(1,n))//' ',n=1,nic2d)
        do n=1,nic2d
          id = getvarid(filemap,icmap2d(1,n))
          if(id.lt.0) error = .true.
        enddo
        if (error) then
          write(stdout,*) '=> Error! A 2-D variable is not found in '//
     &     filemap%filename
          write(stdout,*) 'The following variables are found: ',
     &    (trim(filemap%ncvar(n)%name)//' ',n=1,filemap%nvars)
          write(stdout,*) 'Check "nic2d" in "readic.h" or namelist '//
     &    '"namelist.icdata"!'
          call abort()
        endif ! pe = 0
        allocate(buff8(idim_ic,jdim_ic,1,1))
      endif
# if defined sgi_mipspro && defined use_libSMA
      if (nic2d.gt.0) call mpp_malloc( ptr2, size(hbuf), len2 )
# else
      if (nic2d.gt.0) allocate(hbuf(imt,jmt))
# endif  
      do n=1, nic2d
        if (pe. eq. 0) then
	  hbuf(:,:) = 0.
          id = getvarid(filemap,icmap2d(1,n))   
          count=(/idim_ic,jdim_ic,1,1/)
          start=(/1,1,ictimestep,1/)
	  jstart = max(1+joff_ic ,1)
	  jend   = min(jmt-joff_ic,jmt)
	  istart = max(1+ioff_ic ,1)
	  iend   = min(imt-ioff_ic,imt)
          call mpp_read
     &       (filemap,varname=icmap2d(1,n),start=start,count=count
     &       ,buff=buff8)
          do jrow=jstart,jend
	    do i=istart,iend  
	      hbuf(i,jrow) =  buff8(i-ioff_ic,jrow-joff_ic,1,1)
	    enddo
	  enddo
c
c     zero out missing values
c
          where (hbuf(:,:) .le. missing) hbuf(:,:) = c0
	endif
c
#  if defined parallel_1d
c
c         send data to the other PEs
c
        call mpp_transmit
     &    (hbuf(1,1), imt*jmt, ALL_PES,
     &     hbuf(1,1), imt*jmt, 0)
#  endif
c
c         transfer the data
c
	error = .true.
# ifdef free_surf_on
#  ifdef explicit_free_surface
        if (icmap2d(2,n).eq. 'etat') then
	  error = .false.
	  jbeg = max(jstask-1,1)
          jend = min(jetask+1,jmt)
	  do jrow=jbeg, jend
	    do i=1,imt
              etat(i,jrow,1) = hbuf(i,jrow) * min(1,kmt(i,jrow))
	    enddo
	  enddo
	  do jrow=jscomp, jecomp
	    csum(n) = csum(n) + checksum (etat(1,jrow,1), imt, 1)
	  enddo
	  do jrow=jstask, jetask
	    do i=1,imt
              ps(i,jrow,1) = etat(i,jrow,1) * grav*rho0
	    enddo
	  enddo
	endif
#  else
        if (icmap2d(1,n).eq. 'ps') then
	  error = .false.
	  do jrow=jstask, jetask
	    do i=1,imt
              ps(i,jrow,1) = 
     &             grav*rho0 * hbuf(i,jrow) * min(1,kmt(i,jrow))
	    enddo
	  enddo
	  do jrow=jscomp, jecomp
	    csum(n) = csum(n) + checksum (ps(1,jrow,1), imt, 1)
	  enddo
	endif
#  endif
        ps(:,:,2)   = ps(:,:,1)
#  ifdef explicit_free_surface
        etat(:,:,2) = etat(:,:,1)
        etat(:,:,3) = etat(:,:,1)
#  else
        pguess(:,:) = ps(:,:,1)
#  endif
	if(icmap2d(2,n).eq. 'ubar'.or.icmap2d(2,n).eq. 'vbar') then
	  error = .false.
	  if(icmap2d(2,n).eq. 'ubar') nu = 1	  
	  if(icmap2d(2,n).eq. 'vbar') nu = 2
	  do jrow=jstask, jetask
	    do i=1,imt
              ubar(i,jrow,nu) = hbuf(i,jrow) * min(1,kmu(i,jrow))
            enddo
	  enddo
	  do jrow=jscomp, jecomp
	    csum(n) = csum(n) + checksum (ubar(1,jrow,nu), imt, 1)
	  enddo
          ubarm1(:,:,nu) = ubar(:,:,nu)
	endif
# else /* stream function */
        if (icmap2d(2,n).eq. 'psi') then
	  error = .false.
	  do jrow=jstask, jetask
	    do i=1,imt
              psi(i,jrow,1) = hbuf(i,jrow) * min(1,kmt(i,jrow))
 	    enddo
	  enddo
	  do jrow=jscomp, jecomp
	    csum(n) = csum(n) + checksum (psi(1,jrow,1), imt, 1)
	  enddo
	  
	  psi(:,:,2) = psi(:,:,1)
	  ptd(:,:)   = psi(:,:,1)
	  num_wrds = (jetask-jstask+1)*imt
          call oput (kflds, num_wrds, nkflds-1, ptd(1,js))
	  call oput (kflds, num_wrds, nkflds  , ptd(1,js))
        endif
# endif /* free_surf_on */	  
        if(error) then
	  write(stdout,*) '=>Error! Variable '//trim(icmap2d(1,n))//
     &     ' is read from file but not assigned to any model variable!'
	  write(stdout,*) 'Check "readic.h" and "namelist.icdata"!'
          call abort()
        endif
      enddo ! n=1,nic2d
      if (allocated(buff8)) deallocate(buff8)
# if defined sgi_mipspro && defined use_libSMA
      call SHPDEALLC( ptr2, error_8, -1 )       
# else  
      if (allocated(hbuf)) deallocate(hbuf)
# endif  
#endif /* read_ic */
c
!-----------------------------------------------------------------------
!     Now initialize the 3D fields    
!-----------------------------------------------------------------------
#ifdef read_ic
!
c-----------------------------------------------------------------------
c       check for variables in file and allocate buffers
c-----------------------------------------------------------------------
      if (pe.eq.0) then
        error  = .false.
        write(stdout,*) 'Searching '//trim(filemap%filename)
     &                  //' for the following 3D-variables:'
        write(stdout,*) (trim(icmap3d(1,n))//' ',n=1,nic3d)
        do n=1,nic3d
          id = getvarid(filemap,icmap3d(1,n))
          if(id.lt.0) error = .true.
        enddo
        if (error) then
          write(stdout,*) '=> Error! A 3-D variable is not found in '//
     &     filemap%filename
          write(stdout,*) 'The following variables are found: ',
     &    (trim(filemap%ncvar(n)%name)//' ',n=1,filemap%nvars)
        endif 
        allocate(buff8(idim_ic,1,kdim_ic,1))
!     Check, where to put the data
      endif
      trmap(:) = -1
      idmap(:) = -1
      error = .false.
# if defined sgi_mipspro && defined use_libSMA
      if (nic3d.gt.0) call mpp_malloc( ptr1, size(tbuf), len1 )
# else
      if (nic3d.gt.0) allocate(tbuf(imt,kbot))
# endif  
      do m=1, nic3d
!-----------------------------------------------------------------------
!       set trmap to point to the corresponding model tracer index 
!-----------------------------------------------------------------------
	if(icmap3d(2,m)(1:7).eq.'tracer_') then
          len = len_trim(icmap3d(2,m))
	  read(icmap3d(2,m)(8:len),*) trmap(m)
!         Error stop for n > nt
	  if (trmap(m).gt.nt) then
	    write(stdout,*) '=> Error! '//
     &  		  'index('//trim(icmap3d(2,m))//') > ',nt
	    call abort()
	  endif
	endif
	if(icmap3d(2,m).eq.'u') then
	  trmap(m) = nt + 1
	endif
	if(icmap3d(2,m).eq.'v') then
	  trmap(m) = nt + 2
	endif
	do n=1, nt
	  if (trim(icmap3d(2,m)).eq.trim(trname(n))) then
	    if (trmap(m).lt.0) then
	      trmap(m) = n
	    else
!           Error if more than one file variable is assigned to a model variable
	      write(stdout,*) '=> Error! Tracer ',n,' has the name'//
     &            trim(trname(n))//' This confuses readic!'//
     &           ' A tracer name must not start with "tracer_"!'
              error = .true.
	    endif
	  endif
	enddo 
!       Error stop, if a file variable cannot be assigned to a model variable.	
	if(trmap(m).lt.0) then
	  write(stdout,*) '=> Error! Reading "'//trim(icmap3d(1,m))//
     &  		  '" but there is no variable "'//
     &  		  trim(icmap3d(2,m))//'" to put the data in!'
	  write(stdout,*) 'Check "readic.h" and "namelist.icdata"!'
	  error = .true.
        endif
	do m1=1,m-1
	  if (trmap(m1).eq.trmap(m)) then
            write(stdout,*) 
     &       '=> Error! File variable '//trim(icmap3d(1,m))//' and '
     &       //trim(icmap3d(1,m1))//' initialize the same tracer '//
     &       'variable '//trim(trname(trmap(m1)))//'!'
	    write(stdout,*) 'Check "readic.h" and "namelist.icdata"!'
            error = .true.
	  endif
	enddo
	if (pe.eq.0) then
	  idmap(m)=getvarid(filemap,icmap3d(1,m))
	  do m1=1,m-1
	    if (idmap(m1).eq.idmap(m)) then
              write(stdout,*) 
     &         '=> Error! Tracer '//trim(trname(trmap(m)))//
     &         ' is initialized from the same file variable as tracer '
     &         //trim(trname(trmap(m1)))//'!'
	      write(stdout,*) 'Check "readic.h" and "namelist.icdata"!'
              write(stdout,*) 'Remove the error stop if this is really'
     &             //' desired!'
              error = .true.
	    endif
	  enddo
        endif
      enddo
      if (error) call abort()
!-----------------------------------------------------------------------
!     Check, if all model tracers are present in the file. Stop, if one 
!     tracer is absent. Remove the error stop if this is desired.       
!-----------------------------------------------------------------------
      if (pe.eq.0) then
        do n=1,nt
          error = .true.
	  do m=1,nic3d
	    if(trmap(m).eq.n) error = .false.
  	  enddo 
	  if (error) then 
            write(stdout,*) '=> Error! Tracer '//trim(trname(n))//
     &         ' is not in file '//trim(filemap%filename)//'!'
            write(stdout,*) 'To continue without initializing '
     &                    //trim(trname(n))//' remove this error stop.'
           call abort()
          endif
        enddo
      endif	
      do m=1, nic3d
        trmap_check : select case(trmap(m))
	  case(1:nt) 
            write(stdout,*) trim(icmap3d(1,m))//' --> '//
     &	          'tracer ',trmap(m),' ('//trim(trname(trmap(m)))//')'    
	  case(nt+1)
            write(stdout,*) trim(icmap3d(1,m))//' --> '//
     &	          'zonal velocity (u)'    
	  case(nt+2)
            write(stdout,*) trim(icmap3d(1,m))//' --> '//
     &	          'meridional velocity (v)'    
	end select trmap_check
      enddo
# endif
#ifdef max_window
c
c-----------------------------------------------------------------------
c     begin section of code for max_window
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
	joff = jscomp - jbuf - 1
	j = jrow - joff
c
c-----------------------------------------------------------------------
c       zero out all variables. velocities are internal modes only
c-----------------------------------------------------------------------
c
        if (pe .eq. pe_of_lat(jrow)) then
          u(:,:,j,:,taup1) = c0
	  t(:,:,j,:,taup1) = c0
c       Initialize passive tracers to 1.0 (k=1) and ~0.0 (k>1)
c       (~0.0 is used to distinguish ocean from land for NetCDF)
          if (nt .gt. 2) then
            do n=3,nt
              do k=1,km
                do i=1,imt
                  t(i,k,j,n,taup1) = passive0 (yt(jrow), zt(k), k)
                enddo
              enddo
            enddo
	  endif
        endif ! (pe .eq. pe_of_lat(jrow))
# ifdef read_ic
c
c-----------------------------------------------------------------------
c     read ic data slab by slab from netCDF variables
c-----------------------------------------------------------------------
c
        if (jrow.ge.2 .and. jrow.le.jmtm1) then
c	
          if (pe .eq. 0) then
            count=(/idim_ic,1,kdim_ic,1/)
            start=(/1,jrow-joff_ic,1,ictimestep/)
          endif
          do n=1,nic3d
            if (pe .eq. 0) then
	      id=idmap(n)
	      tbuf(:,:) = c0
              call mpp_read
     &             (filemap,id,start=start,count=count,buff=buff8)
	      do i=2,imtm1
                tbuf(i,:) = buff8(i-ioff_ic,1,:,1)
	      enddo
            endif
c
#  if defined parallel_1d
c
c     send data to the other PEs
c
            call mpp_sync()
            if (pe_of_lat(jrow) .ne. 0) then
              if (pe .eq. 0) call mpp_transmit (tbuf(1,1), imt*kbot,
     &            pe_of_lat(jrow), tbuf(1,1), imt*kbot, NULL_PE)
              if (pe .eq. pe_of_lat(jrow)) call mpp_transmit 
     &           (tbuf(1,1),imt*kbot, NULL_PE, tbuf(1,1), imt*kbot, 0)
            endif
            call mpp_sync()
#  endif
c
            if (pe .eq. pe_of_lat(jrow)) then
!             remove missing values for correct checksums
              where (tbuf .le. missing) tbuf = c0
!
c-----------------------------------------------------------------------
!             velocities
c-----------------------------------------------------------------------
!
	      if (trmap(n).gt.nt) then
c     zero out missing velocities
c
!               u
                if (trmap(n).eq.nt+1) then
		  nu = 1
                  do k=1,kbot
	            do i=1,imt
	              u(i,k,j,1,taup1) = tbuf(i,k)
	            enddo
	          enddo
		  csum(n+nic2d) = csum(n+nic2d) 
     &                      + checksum (u(1,1,j,1,taup1), imt, kbot)
	        else
!               v
		  nu = 2
                  do k=1,kbot
  	            do i=1,imt
	              u(i,k,j,2,taup1) = tbuf(i,k)
	            enddo
	          enddo
		  csum(n+nic2d) = csum(n+nic2d) 
     &                      + checksum (u(1,1,j,2,taup1), imt, kbot)
		endif
c
c-----------------------------------------------------------------------
c               subtract barotropic mode
c-----------------------------------------------------------------------
c
#  ifdef partial_cell
c
c               construct thickness of T cells and U cells
c               (an aberration of the canonical indexing is used here)
c  
                jp1 = min(j,jmt)
                do k=1,km
	          dht(:,k,j) = dzt(k)
	          dht(:,k,jp1) = dzt(k)
                enddo
#   if defined explicit_free_surface_nonlinear 
                dht(:,1,j) = max(epsln,dzt(1)+etat(:,jrow,1))
#   endif
                do i=2,imt-1
	          kb = kmt(i,jrow)
	          if (kb .gt. 1) dht(i,kb,j) = htp(i,jrow) - zw(kb-1)
                enddo
                minj = min(jrow+1,jmt)
#   if defined explicit_free_surface_nonlinear
                dht(:,1,jp1) = max(epsln,dzt(1)+etat(:,minj,1))
#   endif
                do i=2,imt-1
	          kb = kmt(i,minj)
	          if (kb .gt. 1) dht(i,kb,jp1) = htp(i,minj) - zw(kb-1)
                enddo
#   ifdef cyclic
                call setbcx (dht(1,1,j), imt, km)
                call setbcx (dht(1,1,jp1), imt, km)
#   else
c       do not set zero boundaries on dht
#   endif
c
                do k=1,km
                  do i=2,imt-1
	            dhu(i,k,j) = min(dht(i,k,j)  , dht(i+1,k,j)
     &,         		     dht(i,k,jp1), dht(i+1,k,jp1))
	          enddo
                  dhu(1,k,j)   = dht(1,k,j)
                  dhu(imt,k,j) = dht(imt,k,j)
                enddo
#   if defined explicit_free_surface_nonlinear 
                umsk(:,jrow) = min(kmu(:,jrow),1)
                etau(:,jrow)   = (dhu(:,1,j)-dzt(1))*umsk(:,jrow)
                depthu(:,jrow) = h(:,jrow) + etau(:,jrow)
                depthur(:,jrow) = c0
                where(umsk(:,jrow).eq.c1)
     &                            depthur(:,jrow) = c1/depthu(:,jrow)
#   endif
#   ifdef cyclic
                dhu(1,:,j)	  = dhu(imt-1,:,j)
                dhu(imt,:,j)	  = dhu(2,:,j)
#    if defined explicit_free_surface_nonlinear 
                etau(1,jrow)	  = etau(imt-1,jrow)
                etau(imt,jrow)    = etau(2,jrow)
                depthu(1,jrow)    = depthu(imt-1,jrow)
                depthu(imt,jrow)  = depthu(2,jrow)
                depthur(1,jrow)   = depthur(imt-1,jrow)
                depthur(imt,jrow) = depthur(2,jrow)
#    endif
#   else
                dhu(1,:,j)	  = dht(1,:,j)
                dhu(imt,:,j)	  = dht(imt,:,j)
#   endif   
#  endif
c
	     	do i=2,imtm1
             	  baru = c0
             	  do k=1,kbot
#  ifdef partial_cell
                    baru = baru + u(i,k,j,nu,taup1)*dhu(i,k,j)
#  else
                    baru = baru + u(i,k,j,nu,taup1)*dzt(k)
#  endif
                  enddo
#  if defined explicit_free_surface_nonlinear && defined partial_cell 
                  baru = baru*depthur(i,jrow)
#  else
                  baru = baru*hr(i,jrow)
#  endif
                  do k=1,kbot
		    if(kmu(i,jrow) .ge. k) then
                      u(i,k,j,nu,taup1) = u(i,k,j,nu,taup1) - baru
                    else
	              u(i,k,j,nu,taup1) = c0
		    endif
		  enddo
#  ifdef bbl_ag
		  if(kmu(i,jrow) .ge. km) then
                    u(i,km,j,nu,taup1) = u(i,km,j,nu,taup1) - baru
                  else
	            u(i,km,j,nu,taup1) = c0
		  endif
#  endif
#  if defined rigid_lid_surface_pressure || defined implicit_free_surface || defined explicit_free_surface
c
c
c     if necessary, save barotropic mode to ubar
c
                  if(nu.eq.1) then
		    do m=1,nic2d
	              if(icmap2d(2,m).ne. 'ubar') then
		        ubar(i,jrow,nu) = baru
                        ubarm1(i,jrow,nu) = baru
		      endif
		    enddo
 		  endif
                  if(nu.eq.2) then
		    do m=1,nic2d
	              if(icmap2d(2,m).ne. 'vbar') then
                        ubar(i,jrow,nu) = baru
                        ubarm1(i,jrow,nu) = baru
		      endif
		    enddo
		  endif
#  endif
                enddo
                call setbcx (u(1,1,j,nu,taup1), imt, km)
#  if defined rigid_lid_surface_pressure || defined implicit_free_surface || defined explicit_free_surface
                call setbcx (ubar(1,jrow,nu), imt,1 )
                call setbcx (ubarm1(1,jrow,nu), imt, 1)
#  endif
c
c-----------------------------------------------------------------------
!               tracers
c-----------------------------------------------------------------------
!
              else
c               check where to put the data                
		do k=1,kbot
	          do i=1,imt
	            t(i,k,j,trmap(n),taup1) = tbuf(i,k)
	          enddo
	        enddo
		csum(n+nic2d) = csum(n+nic2d) 
     &                   + checksum (t(1,1,j,trmap(n),taup1), imt, kbot)
                call setbcx (t(1,1,j,n,taup1), imt, km)
              endif ! (n)
c
            endif ! (pe .eq. pe_of_lat(jrow))
c
          enddo ! n=1,nic3d
c
        endif ! (jrow.ge.2 .and. jrow.le.jmtm1)
c
        if (jrow .eq. jmt) then
          if (pe .eq. 0) then
	    call mpp_close(ioic)
            deallocate(filemap%ncvar) 
            if (allocated(buff8)) deallocate(buff8)
	  endif
# if defined sgi_mipspro && defined use_libSMA
          call SHPDEALLC( ptr1, error_8, -1 )       
# else
         if (allocated(tbuf)) deallocate(tbuf)
# endif  
!         print the checksums
# ifdef parallel_1d
	  do n=1,nic2d+nic3d 
            call mpp_sum (csum(n), 1)
	  enddo          
# endif  
          write(stdout,'(/a/)') 'Checksums after reading:'
	  do n=1,nic2d 
	    write(stdout,'(a10,$)') trim(icmap2d(2,n))//': '
	    write(stdout,*) csum(n)
	  enddo          
	  do n=1,nic3d 
	    write(stdout,'(a10,$)') trim(icmap3d(2,n))//': '
	    write(stdout,*) csum(n+nic2d)
	  enddo          
	endif
# endif /* read_ic */

# ifdef levitus_ic
c
c       set T & S using climatological monthly Levitus data prepared
c       (by user) from the MOM database using subroutine "ic.F"
c       note: initial condition time in "tmngr.F" should be set to the
c       middle of the starting month where this data is defined. 
c
        if (pn .eq. 1) then
          if (jrow .eq. 1) then
            write (seqrecl,'(a,i12)') ' words=',imt*km + imt + km + 6+8
            call getunit (iolevt, 'temp.mom.ic'
     &,                'unformatted sequential rewind ieee'//seqrecl)
            call getunit (iolevs, 'salt.mom.ic'
     &,                'unformatted sequential rewind ieee'//seqrecl)
          endif
c          
          read (iolevt)
#  if defined _CRAY || defined rs6000          
          read (iolevt) stampl, dpml, im, kk, jro, alat, mm, xtlev
     &         ,                 ztlev,tbuf          
#  else          
          read (iolevt) stampl, dpml4, im, kk, jro, alat4, mm, xtlev4
     &,                 ztlev4,tbuf4
          tbuf(:,:)=tbuf4(:,:)
#  endif
          read (iolevs)
#  if defined _CRAY || defined rs6000  
          read (iolevs) stampl, dpml, im, kk, jro, alat, mm, xtlev
     &         ,                 ztlev,sbuf
#  else          
          read (iolevs) stampl, dpml4, im, kk, jro, alat4, mm, xtlev4
     &         ,                 ztlev4,tbuf4
          dpml=dpml4
          alat=alat4
          xtlev(:)=xtlev4(:)
          ztlev(:)=ztlev4(:)
          sbuf(:,:)=tbuf4(:,:)
#  endif          
          
#ifdef salinity_psu            
          do i=1,imt
            do k=1,kbot
              if (sbuf(i,k).ne.0.0) sbuf(i,k) = sbuf(i,k)*1.e3 + 35.0
            end do
          end do
#endif
c
          if (im .ne. imt .or. kk .ne. kbot .or. jro .ne. jrow) then
            write (stdout,'(/a,a,/,6(a,i5)/)')
     &     ' ==> Error reading levitus data for ', stampl
     &,    ' jrow=',jrow, ', jro=',jro, ', imt=', imt, ', im=',im
     &,    ', km=',km, ', kk=',kk
            call abort()
          endif
	endif
#  if defined parallel_1d
        call mpp_sync()
        if (pe_of_lat(jrow) .ne. 0)then
            if (pe .eq. 0) call mpp_transmit (tbuf(1,1), imt*kbot,
     &           pe_of_lat(jrow), tbuf(1,1), imt*kbot, NULL_PE)
            if (pe .eq. pe_of_lat(jrow)) call mpp_transmit (tbuf(1,1),
     &           imt*kbot, NULL_PE, tbuf(1,1), imt*kbot, 0)
            call mpp_sync()
            if (pe .eq. 0 ) call mpp_transmit (sbuf(1,1), imt*kbot,
     &           pe_of_lat(jrow), sbuf(1,1), imt*kbot, NULL_PE)
            if (pe .eq. pe_of_lat(jrow)) call mpp_transmit (sbuf(1,1),
     &           imt*kbot, NULL_PE, sbuf(1,1), imt*kbot, 0)
        endif
        call mpp_sync()
#  endif
c
c       transfer the T and S data
c
        if (pe .eq. pe_of_lat(jrow)) then
          do k=1,kbot
	    do i=1,imt
	      t(i,k,j,1,taup1) = tbuf(i,k)
	      t(i,k,j,2,taup1) = sbuf(i,k)
	    enddo
	  enddo
	endif ! (pe .eq. pe_of_lat(jrow))
c
        if (jrow .eq. jmt .and. pn .eq. 1) then
          call relunit (iolevt)
	  call relunit (iolevs)
	  write (stdout,'(/a,a,/a,g14.7,a/)')
     &    '=> Initialized T & S from Levitus data for ', stampl
     &,   '   averaged over previous ', dpml, ' days'
        endif
# endif /*levitus_ic*/
# ifdef idealized_ic
c
c       construct highly idealized initial density profile
c       theta0---defined function at bottom of setocn.F
c       salt0 ---defined function at bottom of setocn.F
c       passive0 ---defined function at bottom of setocn.F
c
        if (pe .eq. pe_of_lat(jrow)) then
        do i=1,imt
          do k=1,km
            t(i,k,j,1,taup1) = theta0 (yt(jrow), zt(k))
            t(i,k,j,2,taup1) = salt0 (yt(jrow), zt(k))
          enddo
        enddo
        if (jrow .eq. je) then
	  write (stdout,'(a,a)')
     & '=> Note: initialized T & S to internally generated idealized'
     &,                    ' values for a test case.'
        endif
        endif ! (pe .eq. pe_of_lat(jrow))
# endif
# if defined linearized_advection || defined equatorial_thermocline
c
c       initialize density profile with tbarz
c
        if (pe .eq. pe_of_lat(jrow)) then
        do n=1,nt
          do i=1,imt
            do k=1,km
              t(i,k,j,n,taup1) = tbarz(k,n)
            enddo
          enddo
        enddo
        if (jrow .eq. je) then
	  write (stdout,'(a,a)')
     & '=> Note: initialized T & S with idealized profile T(z), S=35ppt'
        endif
        endif ! (pe .eq. pe_of_lat(jrow))
# endif
c
# if defined firfil
c
c       remove small scales from filtering latitudes
c       USER INPUT: change number of filter applications to suit needs.
c
        if (pe .eq. pe_of_lat(jrow)) then
        if (numflt(jrow) .ne. 0) then
c        if (jrow .le. 8 .or. jrow .ge. 59) then
          do k=1,km
	    do i=1,imt
	      tmsk(i,k) = min(1,kmt(i,jrow))
	    enddo
	  enddo
	  num_fil = 4
	  write (stdout,*) '=> Note: filtering I.C. on jrow =',jrow
     &,' (lat=',yt(jrow),') using ',num_fil,' passes of fast_fir.'
	  do n=1,nt
	    call fast_fir (t(1,1,j,n,taup1), imt, km, 1, tmsk, num_fil)
	  enddo
	endif
        endif ! (pe .eq. pe_of_lat(jrow))
# endif
# ifdef obctest
c
c       construct idealized initial density profile a la Stevens (1990)
c
        if (pe .eq. pe_of_lat(jrow)) then
        do i=1,imt
          do k=1,km
            if ( jrow .le. (jmt/2) ) then
               t(i,k,j,1,taup1) = tempsouth(k)
#  if defined salinity_psu
               t(i,k,j,2,taup1) = saltsouth(k)
#  else
               t(i,k,j,2,taup1) = (saltsouth(k)-35.0)*0.001
#  endif
            else
               t(i,k,j,1,taup1) = tempnorth(k)
#  if defined salinity_psu
               t(i,k,j,2,taup1) = saltnorth(k)
#  else
               t(i,k,j,2,taup1) = (saltnorth(k)-35.0)*0.001
#  endif
            endif
          enddo
        enddo
        if (jrow .eq. je) then
	  write (stdout,'(a,a)')
     & '=> Note: initialized T & S to Stevens (1990) first test case'
     &,                    '                        '
        endif
        endif ! (pe .eq. pe_of_lat(jrow))
# endif
# ifdef obctest2
c
c       construct idealized two basins a la Stevens (1990), case two
c
        if (pe .eq. pe_of_lat(jrow)) then
        do i=1,imt
          do k=1,km
            t(i,k,j,1,taup1) = 6.
#  if defined salinity_psu
            t(i,k,j,2,taup1) = 34.9
#  else
            t(i,k,j,2,taup1) = (34.9-35.0)*0.001
#  endif
          enddo
        enddo
c
        if (jrow .eq. je) then
	  write (stdout,'(a,a)')
     & '=> Note: initialized T & S to Stevens (1990) second test case'
     &,                    '                        '
        endif
        endif ! (pe .eq. pe_of_lat(jrow))
# endif
# ifdef partial_cell
c
c       linearly interpolate T & S to partial bottom cells
c
        if (pe .eq. pe_of_lat(jrow)) then
	if (jrow .ge. js .and. jrow .le. je) then
          do i=1,imt
            kb = kmt(i,jrow)
	    if (kb .gt. 1) then
	      dhtkb = htp(i,jrow) - zw(kb-1)
	      fx = (zt(kb)-zw(kb-1))/dzt(kb)
	      ztkb = zw(kb-1) + fx*dhtkb
              dhwkbm1 = ztkb - zt(kb-1)
	      fact = dhwkbm1/dzw(kb-1)
              do n=1,nt
                t(i,kb,j,n,taup1) = t(i,kb-1,j,n,taup1)
     &                 - fact*(t(i,kb-1,j,n,taup1) - t(i,kb,j,n,taup1))
              enddo
            endif
          enddo
        endif
        endif ! (pe .eq. pe_of_lat(jrow))
# endif
c
# ifdef bbl_ag
c
c       linearly interpolate T & S to bbl cells
c
        if (pe .eq. pe_of_lat(jrow)) then
	if (jrow .ge. js .and. jrow .le. je) then
          do i=1,imt
            kb = kmt(i,jrow)
	    if (kb .gt. 1
#  ifdef bbl_mask
     &           .and.mask_t_bbl(i,jrow).eq.1
#  endif
     &           ) then
#  ifdef partial_cell
	      dhtkb   = htp(i,jrow) - zw(kb-1)
	      fx      = (zt(kb)-zw(kb-1))/dzt(kb)
	      ztkb    = zw(kb-1) + fx*dhtkb
              dz1_bbl = ztkb-zt(kb-1)
              dz2_bbl = htp(i,jrow)-ztkb+dzt(km)/2.
#  else
              dz1_bbl = zt(kb)-zt(kb-1)
              dz2_bbl = zw(kb)-zt(kb)+dzt(km)/2.
#  endif
              fact_bbl=(dz2_bbl+dz1_bbl)/dz1_bbl
              do n=1,nt
                t(i,km,j,n,taup1) = t(i,kb-1,j,n,taup1) - fact_bbl*
     &               (t(i,kb-1,j,n,taup1) - t(i,kb,j,n,taup1))
              enddo
            endif
          enddo
        endif
        endif ! (pe .eq. pe_of_lat(jrow))
# endif
c
c-----------------------------------------------------------------------
c       zero out Temperature and Salinity in land points
c-----------------------------------------------------------------------
c
        if (pe .eq. pe_of_lat(jrow)) then
        do i=1,imt
          kz = kmt(i,jrow)
          do k=1,kbot
            if (k .gt. kz) then
	      do n=1,nt
                t(i,k,j,n,taup1) = c0
	      enddo
            endif
          enddo
# ifdef bbl_ag
          if (kz .eq. 0
#  ifdef bbl_mask
     &           .or.mask_t_bbl(i,jrow).eq.0
#  endif
     &         ) then 
	    do n=1,nt
	      t(i,km,j,n,taup1) = c0
            enddo
          endif
# endif
        enddo
        endif ! (pe .eq. pe_of_lat(jrow))
c
c-----------------------------------------------------------------------
c       checksum the initial conditions
c-----------------------------------------------------------------------
c
        if (pe .eq. pe_of_lat(jrow)) then 
          if (jrow .ge. jscomp .and. jrow .le. jecomp) then
            ucksum = ucksum + checksum (u(1,1,j,1,taup1), imt, km)
            vcksum = vcksum + checksum (u(1,1,j,2,taup1), imt, km)
            tcksum = tcksum + checksum (t(1,1,j,1,taup1), imt, km)
            scksum = scksum + checksum (t(1,1,j,2,taup1), imt, km)
# if defined read_ic && defined free_surf_on
#  ifdef explicit_free_surface
            ecksum = ecksum + checksum (etat(1,jscomp,1), imt, 1)
#  endif
            pcksum = pcksum + checksum (ps  (1,jscomp,1), imt, 1)
            ubksum = ubksum + checksum (ubar(1,jscomp,1), imt, 1)
            vbksum = vbksum + checksum (ubar(1,jscomp,2), imt, 1)
# endif
          endif
        endif ! (pe .eq. pe_of_lat(jrow))
c
c-----------------------------------------------------------------------
c       initialize every latitude jrow either on disk (when jmw < jmt)
c       or in the MW (when the last jrow is complete and jmw = jmt)
c-----------------------------------------------------------------------
c
        if (pe .eq. pe_of_lat(jrow)) then
        if (jrow .eq. 1) then
c
c         when jrow=1, initalize (also do buffer rows < jrow=1)
c
          do jj=1,jbuf
            do k=1,km
	      do i=1,imt
	        tmask(i,k,jj) = c0
	        umask(i,k,jj) = c0
	        do n=1,nt
	          t(i,k,jj,n,taum1) = c0
	          t(i,k,jj,n,tau)   = c0
	          t(i,k,jj,n,taup1) = c0
	        enddo
	        do n=1,2
	          u(i,k,jj,n,taum1) = c0
	          u(i,k,jj,n,tau)   = c0
	          u(i,k,jj,n,taup1) = c0
	        enddo
	      enddo
	    enddo
	  enddo
        endif
        endif ! pe .eq. pe_of_lat(jrow))

        if (pe .eq. pe_of_lat(jrow)) then
	if (jrow .eq. jmt) then
c
c         when jrow=jmt, initialize (also do buffer rows > jrow=jmt)
c
          do jj=jmw-jbuf+1,jmw
            do k=1,km
	      do i=1,imt
	        tmask(i,k,jj) = c0
	        umask(i,k,jj) = c0
	        do n=1,nt
	          t(i,k,jj,n,taum1) = c0
	          t(i,k,jj,n,tau)   = c0
	          t(i,k,jj,n,taup1) = c0
	        enddo
	        do n=1,2
	          u(i,k,jj,n,taum1) = c0
	          u(i,k,jj,n,tau)   = c0
	          u(i,k,jj,n,taup1) = c0
	        enddo
	      enddo
	    enddo
	  enddo
	endif
        endif ! (pe .eq. pe_of_lat(jrow))

        if (jrow .eq. je) then
	  call copy_all_rows (taup1, tau)
	  call copy_all_rows (tau, taum1)
	endif
c	
      enddo ! jrow=1,jmt
c
c-----------------------------------------------------------------------
c     end section of code for max_window
c-----------------------------------------------------------------------
c
#else
c
c-----------------------------------------------------------------------
c     begin section of code for NOT max_window
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     get 3dim ic data
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        j = jmw
c
c-----------------------------------------------------------------------
c       zero out all variables. velocities are internal modes only
c-----------------------------------------------------------------------

        u(:,:,j,:,taup1) = c0
	t(:,:,j,:,taup1) = c0
c
c       Initialize passive tracers to 1.0 (k=1) and ~0.0 (k>1)
c       (~0.0 is used to distinguish ocean from land for NetCDF)
c
        if (nt .gt. 2) then
          do n=3,nt
            do k=1,km
              do i=1,imt
                t(i,k,j,n,taup1) = passive0 (yt(jrow), zt(k), k)
              enddo
            enddo
          enddo
	endif
c
# ifdef read_ic
c
c-----------------------------------------------------------------------
c     read ic data slab by slab from netCDF variables
c-----------------------------------------------------------------------
c
        if (jrow.ge.2 .and. jrow.le.jmtm1) then
c	
          if (pe .eq. 0) then
            count=(/idim_ic,1,kdim_ic,1/)
            start=(/1,jrow-joff_ic,1,ictimestep/)
	  endif
          do n=1,nic3d
            if (pe .eq. 0) then
	      id=idmap(n)
	      tbuf(:,:) = c0
              call mpp_read
     &              (filemap,id,start=start,count=count,buff=buff8)
              do i=2,imtm1
		tbuf(i,:) = buff8(i-ioff_ic,1,:,1)
	      enddo
            endif
c
#  if defined parallel_1d
c
c     send data to the other PEs
c
            call mpp_sync()
            if (pe_of_lat(jrow) .ne. 0) then
              if (pe .eq. 0) call mpp_transmit (tbuf(1,1), imt*kbot,
     &            pe_of_lat(jrow), tbuf(1,1), imt*kbot, NULL_PE)
              if (pe .eq. pe_of_lat(jrow)) call mpp_transmit 
     &           (tbuf(1,1),imt*kbot, NULL_PE, tbuf(1,1), imt*kbot, 0)
            endif
            call mpp_sync()
#  endif
c
c     transfer the data
c
            if (pe .eq. pe_of_lat(jrow)) then
!             remove missing values for correct checksums
              where (tbuf .le. missing) tbuf = c0
!
c-----------------------------------------------------------------------
!             velocities
c-----------------------------------------------------------------------
!
	      if (trmap(n).gt.nt) then
c     zero out missing velocities
c
!               u
                if (trmap(n).eq.nt+1) then
                  nu = 1
		  do k=1,kbot
	            do i=1,imt
	              u(i,k,j,1,taup1) = tbuf(i,k)
	            enddo
	          enddo
		  csum(n+nic2d) = csum(n+nic2d) 
     &                    + checksum (u(1,1,j,1,taup1), imt, kbot)
	        else
!               v
                  nu = 2
                  do k=1,kbot
  	            do i=1,imt
	              u(i,k,j,2,taup1) = tbuf(i,k)
	            enddo
	          enddo
		  csum(n+nic2d) = csum(n+nic2d) 
     &                    + checksum (u(1,1,j,2,taup1), imt, kbot)
		endif
c
c-----------------------------------------------------------------------
c               subtract barotropic mode
c-----------------------------------------------------------------------
c
#  ifdef partial_cell
c
c               construct thickness of T cells and U cells
c               (an aberration of the canonical indexing is used here)
c
                do k=1,km
	          dht(:,k,j) = dzt(k)
	          dht(:,k,j+1) = dzt(k)
                enddo
#   if defined explicit_free_surface_nonlinear 
                dht(:,1,j) = max(epsln,dzt(1)+etat(:,jrow,1))
#   endif
                do i=2,imt-1
	          kb = kmt(i,jrow)
	          if (kb .gt. 1) dht(i,kb,j) = htp(i,jrow) - zw(kb-1)
                enddo
                minj = min(jrow+1,jmt)
#   if defined explicit_free_surface_nonlinear 
                dht(:,1,j+1) = max(epsln,dzt(1)+etat(:,minj,1))
#   endif
                do i=2,imt-1
	          kb = kmt(i,minj)
	          if (kb .gt. 1) dht(i,kb,j+1) = htp(i,minj) - zw(kb-1)
                enddo
#   ifdef cyclic
                call setbcx (dht(1,1,j), imt, km)
                call setbcx (dht(1,1,j+1), imt, km)
#   else
c       do not set zero boundaries on dht
#   endif
c
                do k=1,km
                  do i=2,imt-1
	            dhu(i,k,j) = min(dht(i,k,j)  , dht(i+1,k,j)
     &,         		     dht(i,k,j+1), dht(i+1,k,j+1))
	          enddo
                  dhu(1,k,j)   = dht(1,k,j)
                  dhu(imt,k,j) = dht(imt,k,j)
                enddo
#   if defined explicit_free_surface_nonlinear 
                umsk(:,jrow) = min(kmu(:,jrow),1)
                etau(:,jrow)   = (dhu(:,1,j)-dzt(1))*umsk(:,jrow)
                depthu(:,jrow) = h(:,jrow) + etau(:,jrow)
                depthur(:,jrow) = c0
                where(umsk(:,jrow).eq.c1)
     &                            depthur(:,jrow) = c1/depthu(:,jrow)
#   endif
#   ifdef cyclic
                dhu(1,:,j)	  = dhu(imt-1,:,j)
                dhu(imt,:,j)	  = dhu(2,:,j)
#    if defined explicit_free_surface_nonlinear 
                etau(1,jrow)	  = etau(imt-1,jrow)
                etau(imt,jrow)    = etau(2,jrow)
                depthu(1,jrow)    = depthu(imt-1,jrow)
                depthu(imt,jrow)  = depthu(2,jrow)
                depthur(1,jrow)   = depthur(imt-1,jrow)
                depthur(imt,jrow) = depthur(2,jrow)
#    endif
#   else
                dhu(1,:,j)	  = dht(1,:,j)
                dhu(imt,:,j)	  = dht(imt,:,j)
#   endif   
#  endif

	     	do i=2,imtm1
             	  baru = c0
             	  do k=1,kbot
#  ifdef partial_cell
                    baru = baru + u(i,k,j,nu,taup1)*dhu(i,k,j)
#  else
                    baru = baru + u(i,k,j,nu,taup1)*dzt(k)
#  endif
                  enddo
#  if defined explicit_free_surface_nonlinear && defined partial_cell 
                  baru = baru*depthur(i,jrow)
#  else
                  baru = baru*hr(i,jrow)
#  endif
                  do k=1,kbot
		    if(kmu(i,jrow) .ge. k) then
                      u(i,k,j,nu,taup1) = u(i,k,j,nu,taup1) - baru
                    else
	              u(i,k,j,nu,taup1) = c0
		    endif
		  enddo
#  ifdef bbl_ag
		  if(kmu(i,jrow) .ge. km) then
                    u(i,km,j,nu,taup1) = u(i,km,j,nu,taup1) - baru
                  else
	            u(i,km,j,nu,taup1) = c0
		  endif
#  endif
#  if defined rigid_lid_surface_pressure || defined implicit_free_surface || defined explicit_free_surface
c
c     if necessary, save barotropic mode to ubar
c
                  if(nu.eq.1) then
		    do m=1,nic2d
	              if(icmap2d(2,m).ne. 'ubar') then
		        ubar(i,jrow,nu) = baru
                        ubarm1(i,jrow,nu) = baru
		      endif
		    enddo
 		  endif
                  if(nu.eq.2) then
		    do m=1,nic2d
	              if(icmap2d(2,m).ne. 'vbar') then
                        ubar(i,jrow,nu) = baru
                        ubarm1(i,jrow,nu) = baru
		      endif
		    enddo
		  endif
#  endif
                enddo
                call setbcx (u(1,1,j,nu,taup1), imt, km)
#  if defined rigid_lid_surface_pressure || defined implicit_free_surface || defined explicit_free_surface
                call setbcx (ubar(1,jrow,nu), imt,1 )
                call setbcx (ubarm1(1,jrow,nu), imt, 1)
#  endif
!
c-----------------------------------------------------------------------
!               tracers
c-----------------------------------------------------------------------
!
              else
c               check where to put the data                
		do k=1,kbot
	          do i=1,imt
	            t(i,k,j,trmap(n),taup1) = tbuf(i,k)
	          enddo
	        enddo
		csum(n+nic2d) = csum(n+nic2d) 
     &                   + checksum (t(1,1,j,trmap(n),taup1), imt, kbot)
                call setbcx (t(1,1,j,n,taup1), imt, km)
	      endif ! (n)
c
            endif ! (pe .eq. pe_of_lat(jrow))
c
          enddo ! n=1,nic3d
c
        endif ! (jrow.ge.2 .and. jrow.le.jmtm1)
c
        if (jrow .eq. jmt) then
	  if(pe .eq. 0) then
	    call mpp_close(ioic)  
	    deallocate(filemap%ncvar) 
            if (allocated(buff8)) deallocate(buff8)
	  endif
# if defined sgi_mipspro && defined use_libSMA
          call SHPDEALLC( ptr1, error_8, -1 )       
# else
         if (allocated(tbuf)) deallocate(tbuf)
# endif  
!         print the checksums
          write(stdout,'(/a/)') 'Checksums after reading:'
# ifdef parallel_1d
	  do n=1,nic2d+nic3d 
            call mpp_sum (csum(n), 1)
	  enddo          
# endif  
	  do n=1,nic2d 
	    write(stdout,'(a10,$)') trim(icmap2d(2,n))//': '
	    write(stdout,*) csum(n)
	  enddo          
	  do n=1,nic3d 
	    write(stdout,'(a10,$)') trim(icmap3d(2,n))//': '
	    write(stdout,*) csum(n+nic2d)
	  enddo          
	endif
#endif /* read_ic */
#ifdef levitus_ic
c
c       set T & S using climatological monthly Levitus data prepared
c       (by user) from the MOM database using subroutine "ic.F"
c       note: initial condition time in "tmngr.F" should be set to the
c       middle of the starting month where this data is defined. 
c
        if (pn .eq. 1) then
          if (jrow .eq. 1) then
            write (seqrecl,'(a,i12)') ' words=',imt*km + imt + km + 6+8
            call getunit (iolevt, 'temp.mom.ic'
     &,                'unformatted sequential rewind ieee'//seqrecl)
            call getunit (iolevs, 'salt.mom.ic'
     &,                'unformatted sequential rewind ieee'//seqrecl)
          endif
c
          read (iolevt)
#if defined _CRAY || defined rs6000       
          read (iolevt) stampl, dpml, im, kk, jro, alat, mm, xtlev
     &         ,                 ztlev,tbuf
#else          
          read (iolevt) stampl, dpml4, im, kk, jro, alat4, mm, xtlev4
     &,                 ztlev4,tbuf4
          tbuf(:,:)=tbuf4(:,:)
#endif
          read (iolevs)
#if defined _CRAY || defined rs6000           
          read (iolevs) stampl, dpml, im, kk, jro, alat, mm, xtlev
     &         ,                 ztlev,sbuf
#else
          read (iolevs) stampl, dpml4, im, kk, jro, alat4, mm, xtlev4
     &         ,                 ztlev4,tbuf4          
          dpml=dpml4
          alat=alat4
          xtlev(:)=xtlev4(:)
          ztlev(:)=ztlev4(:)
          sbuf(:,:)=tbuf4(:,:)
#endif
          
#ifdef salinity_psu            
          do i=1,imt
            do k=1,kbot
              if (sbuf(i,k).ne.0.0) sbuf(i,k) = sbuf(i,k)*1.e3 + 35.0
            end do
          end do
#endif
c
          if (im .ne. imt .or. kk .ne. kbot .or. jro .ne. jrow) then
            write (stdout,'(/a,a,/,6(a,i5)/)')
     &     ' ==> Error reading levitus data for ', stampl
     &,    ' jrow=',jrow, ', jro=',jro, ', imt=', imt, ', im=',im
     &,    ', km=',km, ', kk=',kk
            call abort()
          endif
	endif
# if defined parallel_1d
        call mpp_sync()  
        if (pe_of_lat(jrow) .ne. 0)then
            if (pe .eq. 0) call mpp_transmit (tbuf(1,1), imt*kbot,
     &           pe_of_lat(jrow), tbuf(1,1), imt*kbot, NULL_PE)
            if (pe .eq. pe_of_lat(jrow)) call mpp_transmit (tbuf(1,1),
     &           imt*kbot, NULL_PE, tbuf(1,1), imt*kbot, 0)
            call mpp_sync()
            if (pe .eq. 0 ) call mpp_transmit (sbuf(1,1), imt*kbot,
     &           pe_of_lat(jrow), sbuf(1,1), imt*kbot, NULL_PE)
            if (pe .eq. pe_of_lat(jrow)) call mpp_transmit (sbuf(1,1),
     &           imt*kbot, NULL_PE, sbuf(1,1), imt*kbot, 0)
        endif
        call mpp_sync()
# endif
c
c       transfer the T and S data
c
        if (pe.eq.pe_of_lat(jrow)) then
          do k=1,kbot
	    do i=1,imt
	      t(i,k,j,1,taup1) = tbuf(i,k)
	      t(i,k,j,2,taup1) = sbuf(i,k)
	    enddo
	  enddo
        endif ! (pe.eq.pe_of_lat(jrow)
c
        if (jrow .eq. jmt .and. pn .eq. 1) then
          call relunit (iolevt)
	  call relunit (iolevs)
	  write (stdout,'(/a,a,/a,g14.7,a/)')
     &    '=> Initialized T & S from Levitus data for ', stampl
     &,   '   averaged over previous ', dpml, ' days'
        endif
#endif
#ifdef idealized_ic
c
c       construct highly idealized initial density profile
c       theta0---defined function at bottom of setocn.F
c       salt0 ---defined function at bottom of setocn.F
c       passive0 ---defined function at bottom of setocn.F
c
        do i=1,imt
          do k=1,km
            t(i,k,j,1,taup1) = theta0 (yt(jrow), zt(k))
            t(i,k,j,2,taup1) = salt0 (yt(jrow), zt(k))
          enddo
        enddo
        if (jrow .eq. je) then
	  write (stdout,'(a,a)')
     & '=> Note: initialized T & S to internally generated idealized'
     &,                    ' values for a test case.'
        endif
#endif
#if defined linearized_advection || defined equatorial_thermocline
c
c       initialize density profile with tbarz
c
        do n=1,nt
          do i=1,imt
            do k=1,km
              t(i,k,j,n,taup1) = tbarz(k,n)
            enddo
          enddo
        enddo
        if (jrow .eq. je) then
	  write (stdout,'(a,a)')
     & '=> Note: initialized T & S with idealized profile T(z), S=35ppt'
        endif
#endif
c
#if defined firfil
c
c       remove small scales from filtering latitudes
c       USER INPUT: change number of filter applications to suit needs.
c
        if (numflt(jrow) .ne. 0) then
c        if (jrow .le. 8 .or. jrow .ge. 59) then
          do k=1,km
	    do i=1,imt
	      tmsk(i,k) = min(1,kmt(i,jrow))
	    enddo
	  enddo
	  num_fil = 4
	  write (stdout,*) '=> Note: filtering I.C. on jrow =',jrow
     &,' (lat=',yt(jrow),') using ',num_fil,' passes of fast_fir.'
	  do n=1,nt
	    call fast_fir (t(1,1,j,n,taup1), imt, km, 1, tmsk, num_fil)
	  enddo
	endif
#endif
#ifdef obctest
c
c       construct idealized initial density profile a la Stevens (1990)
c
        do i=1,imt
          do k=1,km
            if ( jrow .le. (jmt/2) ) then
               t(i,k,j,1,taup1) = tempsouth(k)
# if defined salinity_psu
               t(i,k,j,2,taup1) = saltsouth(k)
# else
               t(i,k,j,2,taup1) = (saltsouth(k)-35.0)*0.001
# endif
            else
               t(i,k,j,1,taup1) = tempnorth(k)
# if defined salinity_psu
               t(i,k,j,2,taup1) = saltnorth(k)
# else
               t(i,k,j,2,taup1) = (saltnorth(k)-35.0)*0.001
# endif
            endif
          enddo
        enddo
        if (jrow .eq. je) then
	  write (stdout,'(a,a)')
     & '=> Note: initialized T & S to Stevens (1990) first test case'
     &,                    '                        '
        endif
#endif
#ifdef obctest2
c
c       construct idealized two basins a la Stevens (1990), case two
c
        do i=1,imt
          do k=1,km
            t(i,k,j,1,taup1) = 6.
# if defined salinity_psu
            t(i,k,j,2,taup1) = 34.9
# else
            t(i,k,j,2,taup1) = (34.9-35.0)*0.001
# endif
          enddo
        enddo
c
        if (jrow .eq. je) then
	  write (stdout,'(a,a)')
     & '=> Note: initialized T & S to Stevens (1990) second test case'
     &,                    '                        '
        endif
#endif
#ifdef partial_cell
c
c       linearly interpolate T & S to partial bottom cells
c
	if (jrow .ge. jscomp .and. jrow .le. jecomp) then
          do i=1,imt
            kb = kmt(i,jrow)
	    if (kb .gt. 1) then
	      dhtkb = htp(i,jrow) - zw(kb-1)
	      fx = (zt(kb)-zw(kb-1))/dzt(kb)
	      ztkb = zw(kb-1) + fx*dhtkb
              dhwkbm1 = ztkb - zt(kb-1)
	      fact = dhwkbm1/dzw(kb-1)
              do n=1,nt
                t(i,kb,j,n,taup1) = t(i,kb-1,j,n,taup1)
     &                 - fact*(t(i,kb-1,j,n,taup1) - t(i,kb,j,n,taup1))
              enddo
            endif
          enddo
        endif
#endif
c
#ifdef bbl_ag
c
c       linearly interpolate T & S to bbl bottom cells
c
	if (jrow .ge. jscomp .and. jrow .le. jecomp) then
          do i=1,imt
            kb = kmt(i,jrow)
	    if (kb .gt. 1
# ifdef bbl_mask
     &           .and.mask_t_bbl(i,jrow).eq.1
# endif
     &           ) then
# ifdef partial_cell
	      dhtkb   = htp(i,jrow) - zw(kb-1)
	      fx      = (zt(kb)-zw(kb-1))/dzt(kb)
	      ztkb    = zw(kb-1) + fx*dhtkb
              dz1_bbl = ztkb-zt(kb-1)
              dz2_bbl = htp(i,jrow)-ztkb+dzt(km)/2.
# else
              dz1_bbl = zt(kb)-zt(kb-1)
              dz2_bbl = zw(kb)-zt(kb)+dzt(km)/2.
# endif
              fact_bbl=(dz2_bbl+dz1_bbl)/dz1_bbl
              do n=1,nt
                t(i,km,j,n,taup1) = t(i,kb-1,j,n,taup1)
     &               - fact_bbl*
     &               (t(i,kb-1,j,n,taup1) - t(i,kb,j,n,taup1))
              enddo
            endif
          enddo
        endif
#endif
c
c-----------------------------------------------------------------------
c       zero out tracer variables in land points
c-----------------------------------------------------------------------
c
        do i=1,imt
          kz = kmt(i,jrow)
          do k=1,kbot
            if (k .gt. kz) then
	      do n=1,nt
                t(i,k,j,n,taup1) = c0
	      enddo
            endif
          enddo
#ifdef bbl_ag
          if (kz .eq. 0
# ifdef bbl_mask
     &           .or.mask_t_bbl(i,jrow).eq.0
# endif
     &         ) then 
	    do n=1,nt
	      t(i,km,j,n,taup1) = c0
            enddo
          endif
#endif
        enddo
c
c-----------------------------------------------------------------------
c       checksum the initial conditions
c-----------------------------------------------------------------------
c
        
        if (jrow .ge. jscomp .and. jrow .le. jecomp) then
          ucksum = ucksum + checksum (u(1,1,j,1,taup1), imt, km)
          vcksum = vcksum + checksum (u(1,1,j,2,taup1), imt, km)
          tcksum = tcksum + checksum (t(1,1,j,1,taup1), imt, km)
          scksum = scksum + checksum (t(1,1,j,2,taup1), imt, km)
#if defined read_ic && defined free_surf_on
# ifdef explicit_free_surface
          ecksum = ecksum + checksum (etat(1,jrow,1), imt, 1)
# endif
          pcksum = pcksum + checksum (ps  (1,jrow,1), imt, 1)
          ubksum = ubksum + checksum (ubar(1,jrow,1), imt, 1)
          vbksum = vbksum + checksum (ubar(1,jrow,2), imt, 1)
#endif
        endif
c
c-----------------------------------------------------------------------
c       initialize every latitude jrow either on disk (when jmw < jmt)
c       or in the MW (when the last jrow is complete and jmw = jmt)
c-----------------------------------------------------------------------
	if ((jrow .ge. jscomp .and. jrow .le. jecomp)
     &       .or. (jrow .eq. 1 .and. pe .eq. 0)
     &       .or. (jrow .eq. jmt .and. pe .eq. npes)) then
          call putrow (latdisk(taudisk),  nslab, jrow
     &,                u(1,1,j,1,taup1), t(1,1,j,1,taup1)) 
          call putrow (latdisk(taup1disk), nslab, jrow
     &,                u(1,1,j,1,taup1), t(1,1,j,1,taup1)) 
        endif
c
      enddo !  jrow=1,jmt
c
# if defined read_ic && defined free_surf_on      
#   ifdef parallel_1d
      call mpp_update_domains (ps, mom_domain)
      call mpp_update_domains (ubar, mom_domain)
      call mpp_update_domains (ubarm1, mom_domain)
#   endif
# endif
c
c-----------------------------------------------------------------------
c     end section of code for NOT max_window
c-----------------------------------------------------------------------
c
#endif
c      
#ifdef parallel_1d
      call mpp_sum (tcksum, 1)
      call mpp_sum (scksum, 1)
      call mpp_sum (ucksum, 1)  
      call mpp_sum (vcksum, 1)
# ifdef free_surf_on
#  ifdef explicit_free_surface
      call mpp_sum (ecksum, 1)
#  endif
      call mpp_sum (pcksum, 1)
      call mpp_sum (ubksum, 1)
      call mpp_sum (vbksum, 1)
# endif
#endif
#if defined debug_setocn
      write (stdout,*) ' I.C. checksum for t =',tcksum
      write (stdout,*) ' I.C. checksum for s =',scksum
      write (stdout,*) ' I.C. checksum for u =',ucksum
      write (stdout,*) ' I.C. checksum for v =',vcksum
# ifdef free_surf_on
      write (stdout,*) ' I.C. checksum for etat =',ecksum
      write (stdout,*) ' I.C. checksum for ps	=',pcksum
      write (stdout,*) ' I.C. checksum for ubar =',ubksum
      write (stdout,*) ' I.C. checksum for vbar =',vbksum
# endif
#endif
c
#ifdef max_window
c
c-----------------------------------------------------------------------
c     begin section of code for max_window
c-----------------------------------------------------------------------
c
#if defined sponges && !defined sponges_old
c
c     reset i.c. within sponges
c
      if (j_north_spng .ne. jmt) then
        do jrow=j_north_spng+1,jmt
          if (pe .eq. pe_of_lat(jrow)) then
	  joff = jscomp - jbuf - 1
	  j = jrow - joff
          do i=1,imt
            do k=1,km
              u(i,k,j,1,tau)   = u(i,k,j_north_spng,1,tau)
              u(i,k,j,2,tau)   = u(i,k,j_north_spng,2,tau)
              u(i,k,j,1,taup1) = u(i,k,j_north_spng,1,taup1)
              u(i,k,j,2,taup1) = u(i,k,j_north_spng,2,taup1)
              do n=1,nt
                t(i,k,j,n,tau)   = t(i,k,j_north_spng,n,tau)
		t(i,k,j,n,taup1) = t(i,k,j_north_spng,n,taup1)
              enddo
            enddo
          enddo
          write (stdout,*)
     &  '=>Note: reset i.c. on row',jrow,' to i.c. on row', j_north_spng
          endif ! (pe .eq. pe_of_lat(jrow))     
        enddo
      endif
c
      if (j_south_spng .ne. 1) then
        do jrow=j_north_spng+1,jmt
          if (pe .eq. pe_of_lat(jrow)) then
	  joff = jscomp - jbuf - 1
	  j = jrow - joff
          do i=1,imt
            do k=1,km
              u(i,k,j,1,tau)   = u(i,k,j_south_spng,1,tau)
              u(i,k,j,2,tau)   = u(i,k,j_south_spng,2,tau)
              u(i,k,j,1,taup1) = u(i,k,j_south_spng,1,taup1)
              u(i,k,j,2,taup1) = u(i,k,j_south_spng,2,taup1)
              do n=1,nt
                t(i,k,j,n,tau)   = t(i,k,j_south_spng,n,tau)
                t(i,k,j,n,taup1) = t(i,k,j_south_spng,n,taup1)
              enddo
            enddo
          enddo
          write (stdout,*)
     &  '=>Note: reset i.c. on row',jrow,' to i.c. on row', j_south_spng
          endif ! (pe .eq. pe_of_lat(jrow))
        enddo
      endif
#endif
c
c-----------------------------------------------------------------------
c     end section of code for max_window
c-----------------------------------------------------------------------
c
#else
c
c-----------------------------------------------------------------------
c     begin section of code for NOT max_window
c-----------------------------------------------------------------------
c
#if defined sponges && !defined sponges_old
c
c     reset i.c. within sponges
c
      if (j_north_spng .ne. jmt) then
        call getrow (latdisk(taudisk),  nslab, j_north_spng
     &       ,            u(1,1,j,1,taup1), t(1,1,j,1,taup1)) 
        do jrow=j_north_spng+1,jmt
          j = jmw
          write (stdout,*)
     &  '=>Note: reset i.c. on row',jrow,' to i.c. on row', j_north_spng
          call putrow (latdisk(taudisk),  nslab, jrow
     &         ,            u(1,1,j,1,taup1), t(1,1,j,1,taup1)) 
          call putrow (latdisk(taup1disk), nslab, jrow
     &         ,            u(1,1,j,1,taup1), t(1,1,j,1,taup1)) 
        enddo
      endif
c
      if (j_south_spng .ne. 1) then 
        call getrow (latdisk(taudisk),  nslab, j_south_spng
     &       ,            u(1,1,j,1,taup1), t(1,1,j,1,taup1)) 
        do jrow=2,j_south_spng-1
          j = jmw
          write (stdout,*)
     &  '=>Note: reset i.c. on row',jrow,' to i.c. on row', j_south_spng
          call putrow (latdisk(taudisk),  nslab, jrow
     &         ,            u(1,1,j,1,taup1), t(1,1,j,1,taup1)) 
          call putrow (latdisk(taup1disk), nslab, jrow
     &         ,            u(1,1,j,1,taup1), t(1,1,j,1,taup1)) 
        enddo
      endif
#endif
c
c-----------------------------------------------------------------------
c     end section of code for NOT max_window
c-----------------------------------------------------------------------
c
#endif
      return
      end


#if defined term_balances || defined gyre_components
# define balances_gyres_on
#endif
#if defined tracer_averages || defined balances_gyres_on
      subroutine sethr (nr, xstart, xend, ystart, yend)
c
c=======================================================================
c     discretizes the horizontal region to nearset model grid points
c
c     nr     = the horizontal region number
c     xstart = starting longitude at edge of T cell region
c     xend   = ending longitude at edge of  T cell region
c     ystart = starting latitude at edge of  T cell region
c     yend   = ending latitude at edge of T cell region
c=======================================================================
c
      use arrays_2d
#include "param.h"
#include "coord.h"
#include "diag_region.h"
c
c     bail out if start and stop latitudes are outside of domain
c
      if (ystart .lt. yu(1) .and. yend .lt. yu(1)) return
      if (ystart .gt. yu(jmt) .and. yend .gt. yu(jmt)) return
c
c     find the nearest T cell indicies within the region
c
      jsr = min(indp (ystart, yu, jmt)+1, jmt)
      jer = indp (yend, yu, jmt)
c
      isr = min(indp (xstart, xu, imt)+1, imt)
      ier = indp (xend, xu, imt)
c
c     define "edges" of T cell region
c
      if (isr .eq. 1) then
        xsrl = xu(1) - dxudeg(1)
      else
        xsrl = xu(isr-1)
      endif
      xerl = xu(ier)
      if (jsr .eq. 1) then
        ysrl = yu(1) - dyudeg(1)
      else
        ysrl = yu(jsr-1)
      endif
      yerl = yu(jer)
c
      write (hregnm(nr),9000) xsrl, xerl, ysrl, yerl
      write (stdout,*) ' Defining horizontal region # ',nr
     &, ' as T cells within ', hregnm(nr)
      if (isr .gt. ier) then
        write (stdout,*) ' Error: isr=',isr,' >  ier=',ier,' in sethr'
        call abort()
      endif
      if (jsr .gt. jer) then
        write (stdout,*) ' Error: jsr=',jsr,' >  jer=',jer, 'in sethr'
        call abort()
      endif
      do j=jsr,jer
        do i=isr,ier
          if (kmt(i,j) .gt. 0)  then
            mskhr(i,j) = nr
          endif
        enddo
      enddo
c
      nr = nr + 1
c
9000  format ('lon: ',f5.1,' => ',f5.1,'  lat: ',f5.1,' => ',f5.1)
      return
      end


      subroutine setvr (nr, zstart, zend)
c
c=======================================================================
c     discretizes the vertical region to nearset model grid points
c
c     nr     = the vertical region number
c     zstart = starting depth at edge of T cell region in cm.
c     zend   = ending depth at edge of T cell region in cm.
c=======================================================================
c
#include "param.h"
#include "coord.h"
#include "diag_region.h"
c
c     find the nearest T cell indicies within the region
c
      if (zstart .lt. p5*zw(1)) then
        ksr = 1
        ztopb = 0.0
      else
        ksr = min(indp (zstart, zw, km)+1, km)
        ztopb = zw(ksr-1)*0.01
      endif
      ker = indp (zend, zw, km)
      llvreg(nr,1) = ksr
      llvreg(nr,2) = ker
c
      write (vregnm(nr),9000) ztopb, zw(ker)*0.01
      write (stdout,*) ' Defining vertical region # ',nr
     &,  ' as T cells within ',  vregnm(nr)
      if (ksr .gt. ker) then
        write (stdout,*) ' Error: ksr=',ksr,' >  ker=',ker, 'in setvr'
        call abort()
      endif
c
      nr = nr + 1
c
9000  format (' dpt:',f6.1, '=>',f6.1, 'm')
      return
      end
#endif


#if defined sponges
      subroutine spongei
      use timeinterp_mod
      use xsponge_mod
# include "size.h"
# include "stdunits.h"
# include "pconst.h"
# include "coord.h"
# include "iounit.h"
# include "tmngr.h"
      
      real(kind=4) spdpm4,spng4(jmt)
     *,            bufmb14(imt,km,max_rows,2)

      type(time_type) :: scratch_time
c
c-----------------------------------------------------------------------
c     set latitude functions for newtonian damping term in sponge &
c     sponge zones in limited domain models.
c     
c     schematically: damping = -spng(jrow) * (t(i,k,j) - levitus(i,k,j))
c     all related data is assumed to have been prepared using the
c     "sponge" routines in PREP_DATA.
c     disk resource factor of 13 is for 12 months + 1 annual mean
c-----------------------------------------------------------------------
c
      lrec = 8 + 6 + 2*jmt + max_rows*2*imt*km
      write (opt_sponge,'(a,i12)') ' words=',lrec
      call getunit (ionew1, 'sponge.mom'
     &,             'unformatted sequential rewind ieee'//opt_sponge)
      write (stdout,'(/a,1pg10.3,a)')
     & ' Sequential access disk resource for file "sponge.mom" = '
     &,lrec*13*1.e-6,' (MW)'
# if defined equatorial_thermocline
c
c     read "spng" for use with equatorial_thermocline
c
      read (ionew1)
#   if defined _CRAY || defined rs6000
      read (ionew1) sstamp, spdpm, im, kk, jm, num_row, mm
     &     , spng, sponge_row, bufmb1
#   else      
      read (ionew1) sstamp, spdpm4, im, kk, jm, num_row, mm
     &     , spng4, sponge_row, bufmb14
      spdpm=spdpm4
      spng=spng4
      bufmb1=bufmb14
#  endif
        
# else
      lrec = max_rows*2*imt*km
#  ifdef ssread_sswrite
      call ostart (ionew2, 'sponges',13*lrec, lrec, 1)
#  else
      write (stdout,'(/a,1pg10.3,a)')
     & ' Direct access disk resource for file "sponges" = '
     &,lrec*13*1.e-6,' (MW)'
      write (opt_sponge,'(a,1x,i8)') 'unformatted direct words =',lrec
      call getunit (ionew2, 'sponges', opt_sponge) 
#  endif
c
      sum = 0.0
      cksum = 0.0
      do m=1,13
        read (ionew1)
#   if defined _CRAY || defined rs6000
        read (ionew1) sstamp, spdpm, im, kk, jm, num_row, mm
     &       , spng, sponge_row, bufmb1
#   else        
        read (ionew1) sstamp, spdpm4, im, kk, jm, num_row, mm
     &       , spng4, sponge_row, bufmb14
        spdpm=spdpm4
        spng=spng4
        bufmb1=bufmb14
#  endif
        
#  ifdef ssread_sswrite
        call oput (ionew2, lrec, m, bufmb1)
#  else
        write (ionew2, rec=m) bufmb1
#  endif
        if (m .eq. 1) then
c
c         construct time at beginning boundary of dataset
c
          call get_stamp (sstamp, iyr, imon, iday, ihr, imin, isec)
          if (iyr .eq. 0 .and. imon .eq. 1 .and. iday .eq. 1) then
	    iyr = iyr + 2
	  else
	    iyr = iyr + 1
	  endif
          scratch_time = set_date (iyr, imon-1, iday, ihr, imin, isec)
          ndim  = days_in_month(scratch_time)
          sp_time = set_date (iyr, imon, iday, ihr, imin, isec)
	  sp_start_time = decrement_date (sp_time, 0,0,ndim, 0,0,0)
          call get_date (sp_start_time,iyr, imon, iday,ihr,imin,isec)
          call set_stamp (sstamp, iyr, imon, iday, ihr, imin, isec)
	endif
	if (m .eq. 12) then
c
c         construct time at end boundary of dataset
c
          call get_stamp (sstamp, iyr, imon, iday, ihr, imin, isec) 
          if (iyr .eq. 0 .and. imon .eq. 1 .and. iday .eq. 1) then
	    iyr = iyr + 2
	  else
	    iyr = iyr + 1
	  endif
          sp_end_time = set_date (iyr, imon, iday, ihr, imin, isec)
          call set_stamp (sstamp, iyr, imon, iday, ihr, imin, isec)
	endif
c
        if (m .le. 12) then
          spgdpm(m) = spdpm
	  sum = sum + spdpm
	  tspng(m) = sum - 0.5*spdpm
	endif
	do n=1,max_rows
	  cksum = cksum + checksum (bufmb1(1,1,n,1), imt, km)
	  cksum = cksum + checksum (bufmb1(1,1,n,2), imt, km)
	enddo
      enddo
c
      write (stdout,*) ' checksum for sponge data = ',cksum
      if (annlev) then
        write (stdout,'(/a)')
     &  ' => Annual mean levitus data will be used for sponge'
        do m=1,12
#  ifdef ssread_sswrite
          call oput (ionew2, lrec, m, bufmb1)
#  else
	  write (ionew2, rec=m) bufmb1
#  endif
	enddo
      else
        write (stdout,'(/a)')
     &  ' => Monthly levitus data will be used for sponge'
      endif
      write (stdout,'(1a/)')  ' Sponge zone setup:'
      do jrow=1,jmt
        if (spng(jrow) .ne. c0) then
          write (stdout,'(a,i3,a,f7.2,a,e10.3,a,e10.3,a,i3)') ' jrow='
     &,   jrow,', yt(jrow)=',yt(jrow), ', spng(jrow) =',spng(jrow)
     &,   ', Damping time scale (days)=',secday/spng(jrow)
     &,   ', using sponge_row=', sponge_row(jrow)
        endif
      enddo
      write (stdout,*) ' '
c
c     express model time as number of days into dataset
c
      begtim = position_within_data (model_time, sp_start_time
     &,                              sp_end_time, .true.)
c
      if (.not.eqyear) then
        write (stdout,*) '=>Warning leap year being used with '
     &, ' climatological sponges?'
      endif
      iprev  = 1
      inext  = 2
      indxsp = -1
      method = 3
      call timeinterp(begtim, indxsp, tspng, spgdpm, 12, .true., method
     &,               inextd, iprevd, wprev, readsp, inext, iprev)
c
      stprev = "(stamp not available) "
#  ifdef ssread_sswrite
      call oget (ionew2, lrec, iprevd, bufmb(1,1,1,1,iprev))
#  else
      read (ionew2, rec=iprevd) bufmb1
      do n=1,max_rows
        do k=1,km
          do i=1,imt
            bufmb(i,k,n,1,iprev) = bufmb1(i,k,n,1)
            bufmb(i,k,n,2,iprev) = bufmb1(i,k,n,2)
          enddo
        enddo
      enddo
#  endif
      stnext = "(stamp not available) "
#  ifdef ssread_sswrite
      call oget (ionew2, lrec, inextd, bufmb(1,1,1,1,inext))
#  else
      read (ionew2, rec=inextd) bufmb1
      do n=1,max_rows
        do k=1,km
          do i=1,imt
            bufmb(i,k,n,1,inext) = bufmb1(i,k,n,1)
            bufmb(i,k,n,2,inext) = bufmb1(i,k,n,2)
          enddo
        enddo
      enddo
#  endif
c
      write (stdout,'(2(/a,i3,1x,a,a,i2)/a,g14.7,1x,a,a,g14.7/)')
     & ' reading sponge record ', iprevd, stprev, ' into buf ', iprev
     &,' reading sponge record ', inextd, stnext, ' into buf ', inext
     &,' for day =', begtim, stamp, ' weight =',wprev
      write (stdout,'(/a,i2,a,i2/)') 'sponge is dataset index ',indxsp
     &,' for time interpolation using method #',method
      write (stdout,*)
     & ' 1st sponge cksum ='
     &, checksum (bufmb(1,1,1,1,iprev), imt*km, max_rows*2)
     &,' 2nd sponge cksum ='
     &, checksum (bufmb(1,1,1,1,inext), imt*km, max_rows*2)
      call relunit (ionew1)
#  ifndef ssread_sswrite
      call relunit (ionew2)
#  endif
# endif
      return
      end
#endif


      function theta0 (ydeg, depth)
c
c=======================================================================
c     this subroutine returns estimates of global mean potential
c     temperature for model initialization as a function of depth.
c     it is used to produce a reference thermal stratification for the
c     upper 2000m of the MOM`s test case.  below 2000m, the
c     potential temperature returned is 2.0 degrees C.  surface
c     values are set slightly above 18.4 degrees C at the reference
c     latitude "reflat".
c     the estimates are produced from a 7th order ploynomial fit to
c     the annual mean world ocean potential temperature observations
c     of Levitus (1982).
c
c     input [units]:
c       a latitdue (ydeg): [degrees]
c       a zt value (depth): [centimeters]
c     output [units]:
c       potential temperature estimate (est): [degrees centigrade]
c
c     variables:
c       coeft     = coefficients for the polynomial fit of potential
c                   temperature vs. depth
c       reflat    = reference latitude at which observed surface
c                   temperatures approximately equal coeft(1)
c       factor    = the ratio of the cosine of the latitude requested
c                   ("ydeg") to the reference latitude ("reflat")
c                   used to scale the upper 2000 meters of the vertical
c                   temperature profile
c       tmin,tmax = the minumum and maximum potential temperatures
c                   allowed at the time of model initialization
c
c     reference:
c       Levitus, S., Climatological atlas of the world ocean, NOAA
c     Prof. Paper 13, US Gov`t printing Office, Washington, DC, 1982.
c
#include "param.h"
#include "coord.h"
      parameter (ndeg=7)
      dimension coeft(ndeg+1)
      save coeft, tmin, tmax, reflat
      data coeft / 0.184231944E+02,-0.430306621E-01, 0.607121504E-04
     &           ,-0.523806281E-07, 0.272989082E-10,-0.833224666E-14
     &           , 0.136974583E-17,-0.935923382E-22/
      data tmin, tmax, reflat /2.0, 25.0, 34.0/
c
c=======================================================================
c
      pi = atan(1.0) * 4.0
      refcos = abs(cos(pi*reflat/180.))
c
      coslat = abs(cos(pi*ydeg/180.))
      factor = coslat/refcos
      z = depth * 0.01
c
      if (z .gt. 2000.) then
        est = 2.0
      else
        est = c0
        bb = 1.0
        do nn=1,ndeg+1
c          if (nn.gt.1) bb = z**(nn-1)
          est = est + coeft(nn)*bb
	  bb = bb*z
        enddo
        est = est * factor
      endif
c
      if (est .gt. tmax) est = tmax
      if (est .lt. tmin) est = tmin
c
      theta0 = est
c
#ifdef ideal_thermocline
      T0 = 7.5
      T1 = 10.0
      H  = 10.0
      Z0 = 2.0
      theta = T0*(1.-tanh((z-H)/Z0))+T1*(1.-100.*z/zt(km))      
c
c     for cosine weighting in lat
c      theta0 = (coslat**2)*theta
c     for no latitudinal dependence
c
      theta0 = theta
#endif
c
#ifdef constant_initial_temp 
      theta0 = 12.0
#endif
      return
      end

      function salt0 (ydeg, depth)
c
c=======================================================================
c     this subroutine returns a salt profile
c
c     input [units]:
c       a latitdue (ydeg): [degrees]
c       a zt value (depth): [centimeters]
c     output [units]:
c       salt: [model units]
c
c=======================================================================
c
#include "size.h"
#include "coord.h"
c
#if defined salinity_psu
      salt0 = 34.72
#else
      salt0 = 0.03472-0.035
#endif

#ifdef ideal_pycnocline

      c0     = 0.0
      pi     = atan(1.0) * 4.0
      reflat = 34.0
      refcos = abs(cos(pi*reflat/180.))
c
      coslat = abs(cos(pi*ydeg/180.))
      factor = coslat/refcos
c
c     depth in metres     
c
      z = depth * 0.01

      S0 = 1.
      S1 = 35.2
      H1 = 75.0
      H2 = 150.
      H3 = 1./200.
      salt = S1*(1-H3*100.*z/zt(km)) 
# if defined salinity_psu
      salt0 = salt
# else
      salt0 = -(salt-35.0)*.001
# endif

#endif

      return
      end


      function passive0 (ydeg, depth, k)
c
c=======================================================================
c     this subroutine returns a passive tracer profile
c
c     input [units]:
c       a latitdue (ydeg): [degrees]
c       a zt value (depth): [centimeters]
c     output [units]:
c       passive units
c
c=======================================================================
#include "pconst.h"
c
      if (k .eq. 1 ) then
        passive0 = c1
      else
        passive0 = epsln
      endif
      return
      end

      subroutine set_polar_filtering_indices
c=======================================================================
c     set up model indices for filtering polar latitudes
c     called by: setocn.F, topog.F
c=======================================================================
#include "size.h"
#include "index.h"
#include "coord.h"
#include "ndcon.h"

      jfrst  = indp (rjfrst, yt, jmt)
      refcosn = cos(filter_reflat_n/radian)
      refcoss = cos(filter_reflat_s/radian)
      jft1   = indp (filter_reflat_s, yt, jmt)
      if (yt(jft1) .gt. filter_reflat_s) jft1 = jft1 - 1
      jft2   = indp (filter_reflat_n, yt, jmt)
      if (yt(jft2) .lt. filter_reflat_n) jft2 = jft2 + 1
      jfu1   = jft1 - 1
      jfu2   = jft2

      return
      end
