#if !defined coupler_main
      program driver
c
c=======================================================================
c
c     This is the main driver for the Modular Ocean Model (MOM). The  
c     structure allows for surface boundary conditions (SBC) to
c     come from a dataset or an atmospheric model. Integration time
c     is divided into a number of equal time segments and SBC are held 
c     fixed for each time segment. When coupling to an atmospheric
c     model, SBC are supplied each time segment (the coupling period)
c     and held fixed for that period. When supplying SBC to MOM from a
c     dataset, the length of each time segment is one ocean time step.
c     
c     author:        r. c. pacanowski      e-mail=> rcp@gfdl.gov
c                    a. rosati             e-mail=>  ar@gfdl.gov
c    (mpp items)     v. balaji             e-mail=> vb@gfdl.gov
c    (fresh water)   m. schmidt e-mail=>martin.schmidt@io-warnemuende.de
c=======================================================================
c

      use mom3_mpp_mod
      use mpp_io_mod
      use arrays_2d
      use memory_window
      use sbc_info
      use sbc_ocean
#if defined river_inflow
      use river_mod
#endif
#if defined coupled
      use sbc_atmos
#endif
      character*128 file_name, file_stamp
#include "param.h"
#include "coord.h"
#include "iounit.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "taskrows.h"
#include "tmngr.h"
      namelist /mbcin/  ocean_sbc, atmos_sbc, numpas, bwidth
#ifdef new_timing
      common / timers_vb / iitk(99)
      
      do i = 1,99
         iitk(i) = 0
      enddo
      call system_clock(iit0)
#endif
!     there cannot be a parallel_1d option around the following call
      call mpp_init()
      call mpp_domains_init(halosize=imt*km*16)
      call mpp_io_init()
      num_processors = mpp_npes()
      pe = mpp_pe()
      pn = mpp_pe() + 1
      npes = mpp_npes()
c  
c-----------------------------------------------------------------------
c     define standard output as 6 (pn = 1) or 7 (pn > 1)
c-----------------------------------------------------------------------
c
      if (pn .eq. 1) then
        stdout = 6
      else
        stdout = pn + 5
        open
     &   (stdout,FORM='formatted',ACCESS='sequential',POSITION='rewind')
      endif
# ifdef parallel_1d
      call mpp_sync()
# endif
c
      write (stdout,'(/,48x,a,/,/)') '====  Driving MOM  ===='
      write (stdout,'(a,i5,a,i5,a)')
     &' (Executing on ',num_processors
     &,' processor(s): this processor number {pn} =',pn,' )'
c
c-----------------------------------------------------------------------
c     define allowable i/o units for MOM
c-----------------------------------------------------------------------
c
      call ioinit
c
c-----------------------------------------------------------------------
c     set up tasks for each processor
c-----------------------------------------------------------------------
c
      call define_tasks
c
c-----------------------------------------------------------------------
c     initialize timing analysis for various parts of MOM
c-----------------------------------------------------------------------
c
      call inittimers
c
c-----------------------------------------------------------------------
c     choose surface boundary conditions from a list of possibilities
c     (note: index of S.B.C. is determined by order of declaration)
c-----------------------------------------------------------------------
c
c     default ordering for ocean(1=taux, 2=tauy, 3=heatflux, 4=saltflux)
c
      n = 0
      call setup_sbc ('taux', n)
      call setup_sbc ('tauy', n)
#if defined air_pressure
      call setup_sbc ('air press', n)
#endif
#if defined shortwave
      call setup_sbc ('short wave', n)
#endif
c
c     stack all surface tracer fluxes in order (1..nt)
c
      call setup_sbc ('heatflux', n)
      call setup_sbc ('saltflux', n)
#ifdef explicit_fresh_water_flux
      call setup_sbc ('fresh wtr', n)
#endif
      numosbc = n
c
c     default ordering for the atmos (5=sst, 6=sss)
c
      call setup_sbc ('sst', n)
      call setup_sbc ('sss', n)
#ifdef coupled
      call setup_sbc ('surface u', n)
      call setup_sbc ('surface v', n)
#endif
      numasbc = n - numosbc
      numsbc  = numasbc + numosbc
c
c     provide for S.B.C. selection and ordering via namelist
c     (note: if specifying through namelist, the entire list of ocean
c            and atmos names must be specified)
c
      do n=1,maxsbc
        ocean_sbc(n) = 'none'
        atmos_sbc(n) = 'none'
      enddo
c
      bwidth   = 0.0
      numpas   = 100
      call getunit (io, 'namelist.mbcin', 'formatted sequential rewind')
      read  (io, mbcin)
      call relunit (io)
c
      n_ocean = 0
      n_atmos = 0
      do n=1,maxsbc
        if (ocean_sbc(n) .ne. 'none') n_ocean = n
        if (atmos_sbc(n) .ne. 'none') n_atmos = n
      enddo
      if (n_ocean .ne. 0 .and. n_atmos .ne. 0) then
        numosbc = n_ocean
	numasbc = n_atmos
        numsbc  = numosbc + numasbc
        m = 0
        do n=1,numosbc
          call setup_sbc(trim(ocean_sbc(n)), m)
        enddo
        do n=1,numasbc
          call setup_sbc(trim(atmos_sbc(n)), m)
        enddo
      elseif ((n_ocean .eq. 0 .and. n_atmos .ne. 0) .or. (n_ocean .ne. 0
     &        .and. n_atmos .eq. 0)) then
        write (stdout,*)
     & '=>Error: specify ALL atmos AND ocean S.B.C. through namelist.' 
        call abort()
      endif
c
c-----------------------------------------------------------------------
c     printout the setup for surface boundary conditions
c-----------------------------------------------------------------------
c
      write (stdout,*) ' '
      write (stdout,*)
     & ' Surface Boundary Conditions are defined as follows:'
      write (stdout,*) ' '
      write (stdout,*) ' for the ocean (comming from the atmos):'
      do m=1,numsbc
        if (m .le. numosbc) then
          write (stdout,*) ' S.B.C. # ',m,' is ',sbcname(m)
     &,', MOM expects units of ', dunits(m)
#ifdef coupled
     &,', conversion factor from atmos units to MOM is ',coabc(m)
     &,', extrap crit =',crits(m)
#endif
        endif
        if (m .eq. numosbc+1) then
          write (stdout,*) ' '
          write (stdout,*) ' for the atmosphere (comming from MOM):'
        endif
        if (m .gt. numosbc) then
          write (stdout,*) ' S.B.C. # ',m,' is ',sbcname(m)
     &,', atmos expects units of ', dunits(m)
#ifdef coupled
     &,', conversion factor from MOM to atmos is ',coabc(m)
     &,', extrap crit =',crits(m)
#endif
        endif
      enddo
c
#if !defined simple_sbc
c
c     define space for S.B.C on ocean grid and atmos grid
c
      allocate (sbcocn(imt,jmt,numsbc))
#endif        
c
c-----------------------------------------------------------------------
c     do the introductory ocean and river setup once per run
c-----------------------------------------------------------------------
c
      call setocn
c
#ifdef river_inflow
!
      call setriver
!
#endif
#ifdef sea_ice
c-----------------------------------------------------------------------
c     do the introductory sea ice setup once per run
c-----------------------------------------------------------------------
!
      call setice(dtice)
!
#endif
c-----------------------------------------------------------------------
c     do the introductory atmosphere setup once per run
c-----------------------------------------------------------------------
c
      write (stdout,'(/a36/)') ' ==> Note: the atmos setup follows:'
c
#ifdef simple_sbc
      write (stdout,'(/,a,/)')
     &'==> Note: Simple S.B.C. are implemented in setvbc.F'
      dtatm = 0.0
#else
c
c     "setatm" must do the following:
c       1) set up the atmospheric S.B.C. grid definition
c       2) define the atmosphere land/sea mask
c       3) set the atmosphere time step "dtatm" {seconds}
c
      call setatm (dtatm)
#endif
c
#if !defined coupled
c
c     when the MOM S.B.C. come from a dataset, force the segment time
c     and atmospheric time step to one MOM time step. This will force
c     the number of segments to one and the number of time steps per
c     segment to represent the length of the run in days. 
c
      dtatm = dtts
      segtim = dtts*secday
#endif
c
c-----------------------------------------------------------------------
c     compute the number of ocean time steps "numots" for this run and 
c     the number of ocean time steps per ocean segment "ntspos".
c     compute the number of atmos time steps "numats" for this run and
c     the number of atmos time steps per atmos segment "ntspas".
c     divide the integration time "days" into "numseg" segments.
c     each will be length "segtim" days. Surface boundary conditions
c     are supplied every "segtim" days.
c-----------------------------------------------------------------------
c
      numots = nint(rundays/(dtts*secday))
      ntspos = nint(segtim/(dtts*secday))
      numats = nint(rundays/(dtatm*secday))
      ntspas = nint(segtim/(dtatm*secday))
      numseg = numots/ntspos
c
#ifdef coupled
      write (stdout,8800) rundays, numseg, segtim, ntspos, ntspas, dtts
     &,                   dtatm
#else
      write (stdout,8800) rundays, numseg, segtim, ntspos, dtts
#endif
c
c-----------------------------------------------------------------------
c     prepare initial condition S.B.C. for the models
c-----------------------------------------------------------------------
c
#if !defined simple_sbc
      write (stdout,'(/,1x,a)')
     &    '==> Initializing all S.B.C. for MOM now.'
c
c     initialize S.B.C. for the atmosphere (on the ocean S.B.C. grid)
c
      do n=1,numsbc
        do jrow=1,jmt
	  do i=1,imt
	    sbcocn(i,jrow,n) = c0
          enddo
	enddo
      enddo
c
c     load the tracers (SST & SSS) for each row "j". (zero on land)
c     load from the MW if fully opened otherwise load from disk
c
      isst  = index_of_sbc ('sst')
      isss  = index_of_sbc ('sss')
      do jrow=jscomp,jecomp
# if defined max_window
	joff = jscomp - jbuf - 1
	j = jrow - joff
	do i=1,imt
	  if (isst .ne. 0) sbcocn(i,jrow,isst) = t(i,1,j,1,tau)        
	  if (isss .ne. 0) sbcocn(i,jrow,isss) = t(i,1,j,2,tau)
	enddo
# else
	if (isst .ne. 0) call getst (jrow, sbcocn(1,1,isst), 1)         
	if (isss .ne. 0) call getst (jrow, sbcocn(1,1,isss), 2)
# endif
      enddo
             
# ifdef parallel_1d
      call mpp_update_domains (sbcocn, mom_global_domain)
# endif
#endif
c
#ifdef coupled
      if (init .and. numasbc .gt. 0) then
c
c       do nothing. S.B.C. are initialized in setatm.F
c
      else
c
c       when restarting ... read in the necessary S.B.C. data
c       to preserve values in land for extrapolations.
c       (use the generic form for the restart file because the time
c        stamp is unknown)
c
        file_name = 'restrt_sbc.dta'
        call getunit (iosbc, trim(file_name)
     &,               'unformatted sequential rewind')
        read (iosbc) stamp, iotext
        read (iosbc) iimt, ijmt, iima, ijma, inbc
        read (iosbc) stamp, iotext
        read (iosbc) sbcatm
        read (iosbc) stamp, iotext
        read (iosbc) sbcocn
        write (stdout,9100) stamp
	call relunit (iosbc)
      endif
#endif
c
c-----------------------------------------------------------------------
c     check for consistency in the S.B.C. setup
c-----------------------------------------------------------------------
c
      call chkcpl (dtatm, dtts)
c
#ifdef new_timing
      call system_clock(iit1)
      iitk(1) = iitk(1) + iit1 - iit0
      iit0 = iit1
#endif
c
c
c
c-----------------------------------------------------------------------
c     S T A R T    S E G M E N T    L O O P
c-----------------------------------------------------------------------
c
      do n=1,numseg
#ifdef coupled
        write (stdout,9000) 'A T M O S', n, stamp
c
c-----------------------------------------------------------------------
c       get the atmospheric S.B.C. from MOM
c-----------------------------------------------------------------------
c
        call gasbc (n)
#endif
c
#if !defined simple_sbc && !defined minimize_sbc_memory
c
c-----------------------------------------------------------------------
c       call the atmospheric model once for each time step until one
c       segment of "segtim" days is complete. hold atmos S.B.C. fixed
c       during each segment and predict average S.B.C. for MOM
c-----------------------------------------------------------------------
c
        do loop=1,ntspas
# ifdef timing
          call tic ('driver', 'atmos model')
# endif
          asegs   = (loop .eq. 1)
          asege   = (loop .eq. ntspas)
          afirst  = (loop*n .eq. 1)
          alast   = (loop*n .eq. numats)
c
          call atmos
c
# ifdef timing
          call toc ('driver', 'atmos model')
# endif
# ifdef sea_ice
!       update ice with the ice-atmosphere interaction
#  ifdef timing
          call tic ('driver', 'ice_fast model')
#  endif
!
          call ice_fast
c
#  ifdef timing
          call toc ('driver', 'ice_fast model')
#  endif
# endif
        enddo
#endif
c
#ifdef river_inflow
!       reads the fresh water and tracer flux of rivers from a database
        call river
#endif
#ifdef coupled
        write (stdout,9000) 'O C E A N', n, stamp
c
c-----------------------------------------------------------------------
c       get MOM S.B.C.s from the atmosphere
c-----------------------------------------------------------------------
c
        call gosbc (n)
c
#endif
#if !defined sea_ice && defined explicit_fresh_water_flux
!-----------------------------------------------------------------------
!       collect all surface fluxes
!-----------------------------------------------------------------------
!       Without sea ice it is sufficient to call call get_ocean_bc once
!      
        call get_ocean_bc        
#endif
c-----------------------------------------------------------------------
c       call the ocean model once for each time step until one 
c       segment of "segtim" days is complete. hold MOM S.B.C. fixed
c       during each segment and predict average S.B.C. for atmos
c-----------------------------------------------------------------------
c
        do loop=1,ntspos
          osegs   = (loop .eq. 1)
          osege   = (loop .eq. ntspos)
#ifdef sea_ice
!         update ice with the ice-ocean interaction
# ifdef timing
          call tic ('driver', 'ice_slow model')
# endif
!
          call ice_slow
c
# ifdef timing
          call toc ('driver', 'ice_slow model')
# endif
# ifdef explicit_fresh_water_flux
!         With sea ice get_ocean_bc must be called every time step.
          call get_ocean_bc        ! collect all surface fluxes
# endif
#endif
          call mom 
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     E N D    S E G M E N T    L O O P
c-----------------------------------------------------------------------
c
#ifdef new_timing
      call system_clock(iit1)
      iitk(2) = iitk(2) + iit1 - iit0
      iit0 = iit1
#endif
c
c
c
#ifdef coupled
c
c-----------------------------------------------------------------------
c     save all S.B.C. for restarting coupled models
c-----------------------------------------------------------------------
c
c
c     reduce size of old restrt_sbc.dta to save space
c
      file_name = 'restrt_sbc.dta'
      call getunit (iosbc, trim(file_name)
     &,             'unformatted sequential rewind')
      write (iosbc) 'restrt_sbc.dta removed'
      call relunit (iosbc)
c
c     open new restrt_sbc.dta file with date stanp on filename
c
      file_name = file_stamp ('restrt_sbc', xstamp, '.dta')
      call getunit (iosbc, trim(file_name)
     &,             'unformatted sequential rewind')
c
      iotext =  ' read (iosbc) imt, jmt, ima, jma, numsbc'
      write (iosbc) stamp, iotext
      write (iosbc) imt, jmt, ima, jma, numsbc
c
      iotext =
     &  ' read (iosbc) (((sbcatm(i,j,n),i=1,imap2),j=1,jma),n=1,numsbc)'
      write (iosbc) stamp, iotext
      write (iosbc) sbcatm
c
      iotext =
     &  ' read (iosbc) (((sbcocn(i,j,n),i=1,imt),j=1,jmt),n=1,numsbc)'
      write (iosbc) stamp, iotext
      write (iosbc) sbcocn
      write (stdout,9200) stamp
c
      call relunit (iosbc)
#endif
c
      write (stdout,9300)
c
c     write suffix for appending to printout file for identification
c
      call getunit (iostop, 'date.suffix','formatted sequential rewind')
      file_name = file_stamp ('', stop_stamp, '')
      write (iostop,'(a13)') trim(file_name)
      call relunit (iostop)
c
c     show timing analysis
c
#ifdef timing
      write (stdout,'(/,/,/,/a/a/,/)')
     & ' WARNING: timing analysis significantly increases cpu time'
     &,'          so turn it off when making long integrations'
      call showtimers
#endif
c
      call mpp_io_exit()
#ifdef new_timing
      call system_clock(iit1)
      iitk(1) = iitk(1) + iit1 - iit0
      if (pn .eq. 1) then
        write(stdout,*) '(6x,5a10,3(i3,1hx),i2)', 'overhead', 'segloop',
     &         'momloop', 'freesurf', 'diag',
     &         num_processors, imt, jmt, km
      endif
      call flush(stdout)             !flush all buffered writes to stdout
      call mpp_sync()
c system_clock is an ANSI f90 intrinsic: we retrieve inverse clock freq
      call system_clock (count_rate=ncps)
      write( stdout,'(x,a2,i2.2,x,5f10.5,3(i3,1hx),i2)' )'PE', pn-1,
     &     (iitk(i)/float(ncps),i=1,5), num_processors, imt, jmt, km
#endif
#ifdef parallel_1d
      call mpp_exit()
#endif
c
#ifdef coupled
8800  format (/,/,1x,'The model has been configured to run for'
     &,g14.7,' days in ',i6,' segments of ',g14.7,' days each.'
     &,/1x,'There will be ',i6,' ocean time steps per segment'
     &,' and ',i6,' atmosphere time steps per segment.'/
     &,/1x,' The ocean "dtts" =', g14.7, ' seconds'
     &,', and the atmosphere "dtatm" =', g14.7,' seconds'
     &,/,/)
9000  format(/,/,1x,'==== E N T E R I N G    ',a9,'    S E G M E N T  #'
     &,       i4, ',  time = ',a32,' ===='/,/)
9100  format (10x,'==> Finished reading MOM S.B.C. data from file'
     &,    ' restrt_sbc.dta at time = ',a32)
9200  format (10x,'==> Finished writing MOM S.B.C. data to file '
     &,    'restrt_sbc.dta at time = ',a32)
#else
8800  format (/,/,1x,'The model has been configured to run for'
     &,g14.7,' days in ',i6,' segments of ',g14.7,' days each.'
     &,/1x,'There will be ',i6,' ocean time steps per segment.'
     &,/1x,' The ocean "dtts" =', g14.7, ' seconds'
     &,/,/)
#endif
9300  format (/,10x,' ==>  M.O.M. integration is complete.')
      end
#else
      subroutine driver
      return
      end
#endif

      subroutine chkcpl (dtatm, dtts)
      use sbc_info
      use sbc_ocean
#if defined coupled
      use sbc_atmos
#endif
      logical errorc
#include "param.h"
#include "switch.h"
c
c-----------------------------------------------------------------------
c     do consistency checks before allowing model to continue
c-----------------------------------------------------------------------
c
      errorc = .false.
      write (stdout,*) ' '
      write (stdout,*) '    (checking MOM S.B.C. setup)'
c
#ifdef coupled
      if (jma .gt. jmt) then
          write (stdout,9000)
     & '==> Warning: the atmospheric "jma" is greater than "jmt".      '
     &,'    Verify that the atmosphere resolution is less than the     '
     &,'    ocean resolution.                                          '
      endif
      if (imap2 .gt. imt) then
          write (stdout,9000)
     & '==> Warning: the atmospheric "imap2" is greater than "imt".    '
     &,'    Verify that the atmosphere resolution is less than the     '
     &,'    ocean resolution.                                          '
      endif
      if (bwidth .eq. c0) then
        write (stdout,9000)
     & '==> Note: no blending region was set. Verify that the ocean   '
     &,'    domain is global.                                         '
      else
        write (stdout,9000)
     & '==> Warning: blending region for SST set. verify the blending '
     &,'    zone is correctly set by inspecting the "bzone" printout  '
     &,'    this assumes a limited domain ocean and global atmosphere '
        write (stdout,'(1x,a27,g10.3,a8)') ' Note: "bwidth" is set to '
     &,                                  bwidth, ' degrees'
      endif
#endif
      if (dtatm .eq. c0) then
          write (stdout,9000)
     & '==> Error: the atmospheric time step must be set in "setatm"  '
          errorc = .true.
          dtatm = 1.e-6
      endif
c      critv = 1.e-6
      critv = 1.e-4
      if (segtim .ne. c0) then
        r1 = rundays/segtim
      else
        r1 = 0.5
      endif
      r2 = segtim/(dtts*secday)
      r3 = segtim/(dtatm*secday)
      if (segtim .eq. c0) then
          write (stdout,9000)
     & '==> Error: coupling period "segtim" must be specified when    '
     &,'   the "coupled" option is enabled                            '
          errorc = .true.
      elseif (abs(r1-nint(r1)) .gt. critv) then
          write (stdout,9000)
     & '==> Error: there must be an integral number of segments       '
     &,'    "segtim"  within "rundays" (the length of the run)        '
          errorc = .true.
      elseif (abs(r2-nint(r2)) .gt. critv) then
          write (stdout,9000)
     & '==> Error: there must be an integral number of density time   '
     &,'    steps "dtts"  within "segtim" (the segment time)          '
          errorc = .true.
      elseif (abs(r3-nint(r3)) .gt. critv) then
          write (stdout,9000)
     & '==> Error: there must be an integral number of atmos time     '
     &,'    steps "dtatm"  within "segtim" (the segment time)         '
          errorc = .true.
      endif
#if defined coupled && defined restorst
        write (stdout,9000)
     & '==> Warning: restoring to surface tracers ("restorst" enabled)'
     &, '   when "coupled" is also enabled                            '
#endif
#if defined coupled && defined simple_sbc
        write (stdout,9000)
     & '==> Error: ifdef "simple_sbc" is not compatible               '
     &, '   with "coupled"                                            '
        errorc = .true.
#endif
#ifdef coupled
        if (numasbc .eq. 0) then
          write (stdout,9000) 
     & '==> Warning: "coupled" enabled but number of atmosphere       '
     &,'    boundary conditions "numasbc" is zero!                    '
        endif
#endif
      if (numsbc .gt. maxsbc) then
        write (stdout,9000) 
     & '==> Error: "numsbc" is greater than "maxsbc"                  '
        errorc = .true.
      endif
      write (stdout,*) '    (End of MOM S.B.C. checks) '
      write (stdout,*) ' '
      if (errorc) call abort()
c
9000  format (/,(1x,a80))
      return
      end


      subroutine getst (jrow, ocnout, ntabc)
c
c-----------------------------------------------------------------------
c     read surface tracers from disk row "jrow"
c-----------------------------------------------------------------------
c
      use memory_window
#include "param.h"
#include "iounit.h"
#include "mw.h"
#include "tmngr.h"
c
      dimension ocnout(imt,jmt), uu(imt,km,jmw,2), tt(imt,km,jmw,nt)
c
      call getrow (latdisk(taudisk), nslab, jrow, uu, tt)
      do i=1,imt
        ocnout(i,jrow) = tt(i,1,1,ntabc)
      enddo
      return
      end

