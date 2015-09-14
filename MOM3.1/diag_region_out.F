
!=======================================================================
! Write out snapshot in model subregions either netCDF or IEEE. 
! This is a modification of file snapshots_out.F 
!
! There are 2 main options:
!
! -Ddiagnostics_himem allocates buffer areas dimensioned as in 
! max_window (jstask:jetask).
! This runs considerably faster, but at the expense of
! much extra memory. -Ddiagnostics_himem is default.
!
! -Ddiagnostics_lomem allocates buffer areas dimensioned (js:je),
! and does smaller, more frequent writes. 
!
! For himem: call region_out_out whenever within_window_loop=.TRUE. 
! This is done from within the window loop (diag_insidemw.F), and 
! if within_window_loop=.FALSE. from outside the window 
! (in diag_out_window.F).
! For lomem: do same, and write to disk each window as it comes.
!
! Buffering is done when within_window_loop=.TRUE. 
! In lomem case, data is also written if within_window_loop=.TRUE.
! THIS IS NOT READY FOR region_out but can be taken mostly from snapshots_out.F!!!
!
! Each data record is preceded by a character record describing 
! the data record. 
!
! How to use region_out?
! Enable option "snapshots_regions"
! The number of regions must be defined in file diag_region_out.h
! The setup of regions must be defined in file diag_snapregions.h.
! The output filenames are build from the regions name as defined in 
! diag_snapregions.h.
! The types needed to organise data in region_out are defined in module
! snapregions.F90. Compile snapregions.F90 prior region_out.F.
!
! Each region has it's own switch to trigger the output. 
! The switches are defined in the same manner as the interval for 
! snapshots. Add a variable "snapint_reg" to namelist "namelist.diagn"    
! as follows:
! snapint_reg=$snap_region_1, $snap_region_2, ...
! The number of switches must be the same as the number of regions. 
! The unit is the same as that of "snapint".
!
!     authors:   R. C. Pacanowski  e-mail=> rcp@gfdl.gov
!                A. Rosati         e-mail=> ar@gfdl.gov
!                M. Harrison       e-mail=> mjh@gfdl.gov
!
!                Updated with mpp_io routines Summer 1999 by  
!                C. Kerr           e-mail=> ck@gfdl.gov 
!                Balaji            e-mail=> vb@gfdl.gov
!                S. M. Griffies    e-mail=> smg@gfdl.gov 
!     modified as  region_out
!                M. Schmidt        e-mail=> martin.schmidt@io-warnmeunede.de
!=======================================================================
      subroutine region_out  
     & (is, ie, joff, js, je, within_window_loop, nr)
#if defined snapshots_regions 
!
!-----------------------------------------------------------------------
!     save snapshot data
!-----------------------------------------------------------------------
!
      use arrays_2d
      use memory_window
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
      use diag_mod

# include "param.h"
# include "coord.h"
# include "derived_options.h"
# include "diag.h"
# include "docnam.h"
# include "emode.h"
# include "iounit.h"
# include "mw.h"
# include "switch.h"
# include "taskrows.h"
# include "tmngr.h"
# include "vers.h"
# include "diag_snapregions.h"

# ifdef isoneutralmix
#  include "isopyc.h"
#  include "grdvar.h"
# endif

# ifdef kppvmix
#   include "vmixc.h"
# endif

# if defined bbl_ag
#  include "bbl.h"
# endif

# if defined rigid_lid_surface_pressure || defined free_surf_on
#  include "scalar.h"
# endif

      character*128 file_name(nregions), file_stamp, gvtitle(nregions)
      character(len=16) :: cvarnam(nregions)

# if defined salinity_psu
      character*3 saldim
# else 
      character*9 saldim
# endif

      logical :: within_window_loop
      dimension bufzt(km), bufzw(km)
      integer :: num_snaps(nregions), io(nregions)
      logical :: first_out
      data (num_snaps(n),n=1, nregions) /nregions*0/
      data first_out /.true./
      save num_snaps
      save io, first_out

      real, dimension(0:imt), save :: edges_xt, edges_xu
      real, dimension(0:jmt), save :: edges_yt, edges_yu
      real, dimension(0:km ), save :: edges_zt, edges_zw

      type(region2D), allocatable, dimension(:), target:: diag_regions
      save diag_regions

# ifdef diagnostics_himem
      type(region_output), allocatable, target :: output(:)
      type(region_output), pointer :: ou
      save output
# else
      write(stdout,*) 'region_out does not work with diagnostics_lomem'
      call abort()
# endif
!-----------------------------------------------------------------------
!     axis edge arrays
!-----------------------------------------------------------------------

# ifdef diagnostics_himem
      type(domain2D), pointer :: domain, region_domain
# endif

# ifdef diagnostics_lomem
      type(domain2D) :: domain
!
!     dimension within mw and set domain to a domain2D type
!
      real :: trcr(imt,js+joff:je+joff,km,nt)
      real :: u(imt,js+joff:je+joff,km)
      real :: v(imt,js+joff:je+joff,km)
      real :: w(imt,js+joff:je+joff,km)
#  ifdef gent_mcwilliams
      real :: u_gm(imt,js+joff:je+joff,km)
      real :: v_gm(imt,js+joff:je+joff,km)
      real :: w_gm(imt,js+joff:je+joff,km)
#  endif
      real :: stf(imt,js+joff:je+joff,nt)
#  ifdef air_pressure
      real :: smf(imt,js+joff:je+joff,3 )
#  else
      real :: smf(imt,js+joff:je+joff,2 )
#  endif
#  ifdef explicit_free_surface && defined time_averaged_surface
      real :: psiu(imt,js+joff:je+joff)
      real :: psiv(imt,js+joff:je+joff)
      real :: convU(imt,js+joff:je+joff)
#  endif
#  ifdef explicit_fresh_water_flux
      real :: sfft(imt,js+joff:je+joff)
#  endif
#  ifdef stream_function
      real :: psi(imt,js+joff:je+joff)
#  endif
#  if defined(rigid_lid_surface_pressure) || defined(free_surf_on)
      real :: eta(imt,js+joff:je+joff)
#  endif
#  ifdef kppvmix
      real :: hmxl(imt,js+joff:je+joff)
      real :: hblt(imt,js+joff:je+joff)
#  endif
#  if defined bbl_ag
      real :: oppsignx(imt,js+joff:je+joff)
      real :: oppsigny(imt,js+joff:je+joff)
#  endif
# endif


# ifdef timing
      call tic ('diagnostic', 'snapshots ')
# endif
!-----------------------------------------------------------------------
!     Allocate a domain for each region
!     This is needed to concatinate the netcdf output files 
!-----------------------------------------------------------------------
      if (first_out) then
        allocate(diag_regions(nregions) )
        do n=1,nregions       
          allocate (diag_regions(n)%domains(0:npes-1) )
          diag_regions(n)%domain => diag_regions(n)%domains(pe)
          diag_regions(n)%name = regionname(n)
	  call mpp_define_domains( 
     &      (/isreg(n),iereg(n),jsreg(n),jereg(n)/)
     &      , diag_regions(n)%domains, xflags=GLOBAL_COMPUTE_DOMAIN 
     &      , xhalo=0, yhalo=0 )
!-----------------------------------------------------------------------
!         mpp_define_domains distributes a model region over all PEs
!         So the data and compute index are wrong now andneed correction.
!         This is not clean code and could be changed by someone who is
!         familiar with mpp modules and the netcdf stuff! But it works.
!         (m.s.)
!-----------------------------------------------------------------------
          ibeg = isreg(n)
          iend = iereg(n)
          jbeg = max(jsreg(n),jscomp)
          jend = min(jereg(n),jecomp) 
          diag_regions(n)%domain%y%data%start_index   = jbeg
          diag_regions(n)%domain%y%data%end_index     = jend
          diag_regions(n)%domain%y%data%size          = jend - jbeg +1
          diag_regions(n)%domain%y%compute%start_index= jbeg
          diag_regions(n)%domain%y%compute%end_index  = jend
          diag_regions(n)%domain%y%compute%size	      = jend - jbeg +1
        enddo
        allocate(output(nregions))
	first_out = .false.
      endif
!-----------------------------------------------------------------------

      region_domain => diag_regions(nr)%domain
      ou            => output(nr)

!-----------------------------------------------------------------------
!     set averaging period, time, title, ocean depth, salinity units 
!-----------------------------------------------------------------------

      file_name(nr) = file_stamp (diag_regions(nr)%name, xstamp, '.dta')
      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle(nr) = 'Snapshots' // momver
# if defined salinity_psu
      saldim = 'psu' 
# else 
      saldim = 'ppm-0.035' 
# endif
      bufzt = zt(:)*0.01
      bufzw = zw(:)*0.01

!-----------------------------------------------------------------------
!     open diagnostic file and write meta information 
!-----------------------------------------------------------------------
      if ( num_snaps(nr) .eq. 0) then
        call mpp_open(io(nr), file_name(nr), action=MPP_WRONLY,  
     &        form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,  
     &        iospec= '-F cachea')

        num_snaps(nr) = num_snaps(nr) + 1

!-----------------------------------------------------------------------
!       axes (on global domain)
!-----------------------------------------------------------------------

	ibeg = region_domain%x%global%start_index
	iend = region_domain%x%global%end_index
	jbeg = region_domain%y%global%start_index
	jend = region_domain%y%global%end_index
	if (jbeg.gt.jecomp.or.jend.lt.jscomp)     return
        if (iend.lt.ibeg.or.jend.lt.jbeg)         return 

        call mpp_write_meta(  
     &        io(nr), ou%fields%axis_xt, 
     &        'xt_i', 'degrees_E', 'Longitude of T points',  
     &        cartesian='X', domain=region_domain%x, 
     &        data=xt(ibeg:iend))

        call mpp_write_meta(  
     &          io(nr), ou%fields%axis_yt, 
     &          'yt_j', 'degrees_N', 'Latitude of T points',  
     &          cartesian='Y', domain=region_domain%y,  
     &          data=yt(jbeg:jend))

        call mpp_write_meta(  
     &          io(nr), ou%fields%axis_zt, 'zt_k', 'm', 
     &          'Depth of T grid point',  
     &          cartesian='Z', sense=-1, data=bufzt)
 
        call mpp_write_meta(  
     &        io(nr), ou%fields%axis_xu, 
     &        'xu_i', 'degrees_E', 'Longitude of U points',  
     &        cartesian='X', domain=region_domain%x,  
     &        data=xu(ibeg:iend))

        call mpp_write_meta(  
     &          io(nr), ou%fields%axis_yu, 
     &          'yu_j', 'degrees_N', 'Latitude of U points',  
     &          cartesian='Y', domain=region_domain%y,  
     &          data=yu(jbeg:jend))

        call mpp_write_meta(  
     &          io(nr), ou%fields%axis_zw, 
     &          'zw_k', 'm', 'Depth of T cell bottom',  
     &          cartesian='Z', sense=-1, data=bufzw)
 
        call mpp_write_meta(  
     &          io(nr), ou%fields%axis_surf, 'Surface', 'm', 
     &          'Depth of surface',  
     &          cartesian='Z', sense=-1, data=(/0./))

        call mpp_write_meta(  
     &          io(nr), ou%fields%axis_t, 'Time', timeunits,  
     &          'Time since initial condition', cartesian='T')


!-----------------------------------------------------------------------
!       define auxiliary "edges" axes for xt, yt, zt coordinates
!-----------------------------------------------------------------------

        call edge_maker (1, edges_xt, xt, xu, imt)
        call mpp_write_meta(  
     &          io(nr), axis_xt_edges, 'xt_i_edges', 'degrees_E',  
     &          'Longitude of T edges', cartesian='X'  
     & ,        data=edges_xt)

        call edge_maker (1, edges_yt, yt, yu, jmt)
        call mpp_write_meta(  
     &          io(nr), axis_yt_edges, 'yt_j_edges', 'degrees_N',  
     &          'Latitude of T edges', cartesian='Y'  
     & ,        data=edges_yt)
     
        call edge_maker (1, edges_zt, zt, zw, km)
        call mpp_write_meta(  
     &          io(nr), axis_zt_edges, 'zt_k_edges', 'cm',  
     &          'Depth of T edges', cartesian='Z', sense=-1  
     & ,        data=edges_zt)

!-----------------------------------------------------------------------
!       define auxiliary "edges" axes for xu, yu, zw coordinates
!-----------------------------------------------------------------------

        call edge_maker (2, edges_xu, xt, xu, imt)
        call mpp_write_meta(  
     &          io(nr), axis_xu_edges, 'xu_i_edges', 'degrees_E',  
     &          'Longitude of U edges', cartesian='X'  
     & ,         data=edges_xu)

        call edge_maker (2, edges_yu, yt, yu, jmt)
        call mpp_write_meta(  
     &          io(nr), axis_yu_edges, 'yu_j_edges', 'degrees_N',  
     &          'Latitude of U edges', cartesian='Y'  
     & ,         data=edges_yu)
    
        call edge_maker (2, edges_zw, zt, zw, km)
        call mpp_write_meta(  
     &          io(nr), axis_zw_edges, 'zw_k_edges', 'cm',  
     &          'Depth of W edges', cartesian='Z'  
     & ,         sense=-1, data=edges_zw)

!-----------------------------------------------------------------------
!       define fields (on data domain)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!       tracers
!-----------------------------------------------------------------------

        allocate (ou%fields%trcr(nt))
        allocate (ou%fields%stf(nt))
        do n=1,nt 
          if (n .EQ. 1) then
            call mpp_write_meta(  
     &         io(nr), ou%fields%trcr(n), 
     &         (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_zt
     &         ,ou%fields%axis_t/),  
     &         'temp', 'deg C', 'potential temperature',   
     &         -5.0, 50.0, -1.0E+34, +1.0E+34)

          elseif (n .EQ. 2) then
            call mpp_write_meta(  
     &         io(nr), ou%fields%trcr(n), 
     &         (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_zt
     &         ,ou%fields%axis_t/),  
     &         'salinity', saldim, 'salinity',  
     &         -100.00, +100.00, -1.0E+34, +1.0E+34)
          else
            write( cvarnam(nr), '(a,i2.2)' ) 'tracer_', n
            call mpp_write_meta(  
     &         io(nr), ou%fields%trcr(n),
     &         (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_zt
     &         ,ou%fields%axis_t/),  
     &         cvarnam(nr), 'tracer units', trname(n),  
     &         -1.0E+10, +1.0E+10, -1.0E+34, +1.0E+34)
          endif
        enddo

!-----------------------------------------------------------------------
!       velocities
!-----------------------------------------------------------------------

        call mpp_write_meta(  
     &       io(nr), ou%fields%u, 
     &       (/ou%fields%axis_xu,ou%fields%axis_yu,ou%fields%axis_zt
     &       ,ou%fields%axis_t/),  
     &       'u', 'cm/s', 'Zonal velocity',  
     &       -3000.0, +3000.0, -1.0E+34, +1.0E+34)
 
        call mpp_write_meta(  
     &       io(nr), ou%fields%v, 
     &       (/ou%fields%axis_xu,ou%fields%axis_yu,ou%fields%axis_zt
     &       ,ou%fields%axis_t/),  
     &       'v', 'cm/s', 'Meridional velocity',  
     &       -3000.0, +3000.0, -1.0E+34, +1.0E+34)
 
        call mpp_write_meta(  
     &       io(nr), ou%fields%w, 
     &       (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_zw
     &       ,ou%fields%axis_t/),  
     &       'w', 'cm/s', 'W at T cell bottom',  
     &       -100.0, +100.0, -1.0E+34, +1.0E+34)
 
!-----------------------------------------------------------------------
!       additional three-dimensional variables
!-----------------------------------------------------------------------

#  if defined gent_mcwilliams
        call mpp_write_meta(  
     &       io(nr), ou%fields%u_gm, 
     &       (/ou%fields%axis_xu,ou%fields%axis_yu,ou%fields%axis_zt
     &       ,ou%fields%axis_t/),  
     &       'u_gm', 'cm/s', 'Zonal GM velocity on U-cell',  
     &       -1.0E+06, +1.0E+06, -1.0E+34, +1.0E+34)

        call mpp_write_meta(  
     &       io(nr), ou%fields%v_gm, 
     &       (/ou%fields%axis_xu,ou%fields%axis_yu,ou%fields%axis_zt
     &       ,ou%fields%axis_t/),  
     &       'v_gm', 'cm/s', 'Meridional GM velocity on U-cell',  
     &       -1.0E+06, +1.0E+06, -1.0E+34, +1.0E+34)

        call mpp_write_meta(  
     &       io(nr), ou%fields%w_gm, 
     &       (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_zw
     &       ,ou%fields%axis_t/),  
     &       'w_gm', 'cm/s', 'Vertical GM velocity at T-cell bottom',  
     &       -1.0E+06, +1.0E+06, -1.0E+34, +1.0E+34)
#  endif

!-----------------------------------------------------------------------
!       two-dimensional variables
!-----------------------------------------------------------------------

	do n=1,nt
          if(n .eq. 1) then
            call mpp_write_meta(  
     &       io(nr), ou%fields%stf(1), 
     &       (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_surf
     &       ,ou%fields%axis_t/),  
     &       'hflx', 'cal/cm**2/s', 'Surface heat flux',  
     &       -1.0E+00, +1.0E+00, -1.0E+34, +1.0E+34)

          elseif(n.eq.2) then
            call mpp_write_meta(  
     &       io(nr), ou%fields%stf(2), 
     &       (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_surf
     &       ,ou%fields%axis_t/),  
     &       'sflx', 'gm/cm**2/s', 'Surface salt flux',  
     &       -1.0E+00, +1.0E+00, -1.0E+34, +1.0E+34)

          else 
            write( cvarnam(nr), '(a,i2.2)' ) 'stf_', n
            call mpp_write_meta(  
     &         io(nr), ou%fields%stf(n), 
     &       (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_surf
     &       ,ou%fields%axis_t/),  
     &         cvarnam(nr), 'tracer units', 'Surface tracer flux',  
     &         -1.0E+06, +1.0E+06, -1.0E+34, +1.0E+34)
          endif
        enddo

        call mpp_write_meta(  
     &       io(nr), ou%fields%taux, 
     &       (/ou%fields%axis_xu,ou%fields%axis_yu,ou%fields%axis_surf
     &       ,ou%fields%axis_t/),  
     &       'taux', 'dyn/cm**2', 'Zonal windstress',  
     &       -20.0, +20.0, -1.0E+34, +1.0E+34)

        call mpp_write_meta(  
     &       io(nr), ou%fields%tauy, 
     &       (/ou%fields%axis_xu,ou%fields%axis_yu,ou%fields%axis_surf
     &       ,ou%fields%axis_t/),  
     &       'tauy', 'dyn/cm**2', 'Meridional windstress',  
     &       -20.0, +20.0, -1.0E+34, +1.0E+34)

#  ifdef air_pressure
        call mpp_write_meta(  
     &       io(nr), ou%fields%pair, 
     &       (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_surf
     &       ,ou%fields%axis_t/),  
     &       'pair', 'hPa', 'Air pressure',  
     &       0.0, +2000.0, -1.0E+34, +1.0E+34)
#  endif

#  if defined explicit_free_surface && defined time_averaged_surface
       call mpp_write_meta(  
     &       io(nr), ou%fields%convU, 
     &       (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_surf
     &       ,ou%fields%axis_t/),  
     &       'convU', 'cm/s', 'convergence of H*ubar',  
     &       -300000.0, +300000.0, -1.0E+34, +1.0E+34)
       call mpp_write_meta(  
     &       io(nr), ou%fields%psiu, 
     &       (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_surf
     &       ,ou%fields%axis_t/),  
     &       'psiu', 'Sv', 'Streamfunction for H*ubar',  
     &       -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)

       call mpp_write_meta(  
     &       io(nr), ou%fields%psiv, 
     &       (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_surf
     &       ,ou%fields%axis_t/),  
     &       'psiv', 'Sv', 'Streamfunction for H*vbar',  
     &       -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)
#  endif

#  if defined explicit_fresh_water_flux
        call mpp_write_meta(  
     &       io(nr), ou%fields%sfft, 
     &       (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_surf
     &       ,ou%fields%axis_t/),  
     &       'sff', 'cm/s', 'Fresh water flux',  
     &       -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)
#  endif

#  if defined stream_function
        call mpp_write_meta(  
     &       io(nr), ou%fields%psi, 
     &       (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_t/),  
     &       'psi', 'cm**3/s', 'Stream function',  
     &       -300.0, +300.0, -1.0E+34, +1.0E+34)
#  endif

#  if defined rigid_lid_surface_pressure || defined free_surf_on
        call mpp_write_meta(  
     &       io(nr), ou%fields%eta, 
     &       (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_t/),  
     &       'eta', 'cm', 'Surface height',  
     &       -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)
#  endif

#  if defined kppvmix
        call mpp_write_meta(  
     &       io(nr), ou%fields%hmxl, 
     &       (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_t/),  
     &       'hmxl', 'cm', 'Mixed layer depth',  
     &       -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)

        call mpp_write_meta(  
     &       io(nr), ou%fields%hblt, 
     &       (/ou%fields%axis_xt,ou%fields%axis_yt,ou%fields%axis_t/),  
     &       'hblt', 'cm', 'Boundary layer depth',  
     &       -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)
#  endif

#  if defined bbl_ag
        call mpp_write_meta(  
     &       io(nr), ou%fields%oppsignx, 
     &       (/ou%fields%axis_xu,ou%fields%axis_yu,ou%fields%axis_t/),  
     &       'oppsignx', 'cm', 'x-detrainment indicator in bbl',  
     &       -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)

        call mpp_write_meta(  
     &       io(nr), ou%fields%oppsigny, 
     &       (/ou%fields%axis_xu,ou%fields%axis_yu,ou%fields%axis_t/),  
     &       'oppsigny', 'cm', 'y-detrainment indicator in bbl',  
     &       -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)
#  endif
!
!       averaging period 
!
        call mpp_write_meta(  
     &       io(nr), ou%fields%period, (/ou%fields%axis_t/),  
     &       'period', 'day', 'averaging period',  
     &       +0.0E+00, +1.0E+20, -1.0E+34, +1.0E+34)

!
!       figure title 
!
        call mpp_write_meta(io(nr), 'title', cval=gvtitle(nr) )
!
!-----------------------------------------------------------------------
!     write axis now that all meta data has been written
!-----------------------------------------------------------------------

        call mpp_write (io(nr),ou%fields%axis_surf)
        call mpp_write (io(nr),ou%fields%axis_xt)
        call mpp_write (io(nr),ou%fields%axis_xu)
        call mpp_write (io(nr),ou%fields%axis_yt)
        call mpp_write (io(nr),ou%fields%axis_yu)
        call mpp_write (io(nr),ou%fields%axis_zt)
        call mpp_write (io(nr),ou%fields%axis_zw)

        call mpp_write (io(nr),axis_xt_edges)
        call mpp_write (io(nr),axis_xu_edges)
        call mpp_write (io(nr),axis_yt_edges)
        call mpp_write (io(nr),axis_yu_edges)
        call mpp_write (io(nr),axis_zt_edges)
        call mpp_write (io(nr),axis_zw_edges)
!-----------------------------------------------------------------------
!       allocate himem arrays over full 3d domain on processor
!-----------------------------------------------------------------------

        ibeg = region_domain%x%data%start_index
        iend = region_domain%x%data%end_index
        jbeg = region_domain%y%data%start_index
        jend = region_domain%y%data%end_index
	if (jbeg.gt.jecomp.or.jend.lt.jscomp)     return
        if (iend.lt.ibeg.or.jend.lt.jbeg)         return 


# ifdef diagnostics_himem
        allocate( ou%buffs%trcr(ibeg:iend,jbeg:jend,km,nt) )
        allocate( ou%buffs%u(ibeg:iend,jbeg:jend,km) )
        allocate( ou%buffs%v(ibeg:iend,jbeg:jend,km) )
        allocate( ou%buffs%w(ibeg:iend,jbeg:jend,km) )
#  ifdef gent_mcwilliams
        allocate( ou%buffs%u_gm(ibeg:iend,jbeg:jend,km) )
        allocate( ou%buffs%v_gm(ibeg:iend,jbeg:jend,km) )
        allocate( ou%buffs%w_gm(ibeg:iend,jbeg:jend,km) )
#  endif
        allocate( ou%buffs%stf(ibeg:iend,jbeg:jend,nt) )
#  ifdef air_pressure
        allocate( ou%buffs%smf(ibeg:iend,jbeg:jend,3 ) )
#  else
        allocate( ou%buffs%smf(ibeg:iend,jbeg:jend,2 ) )
#  endif
#  if defined explicit_free_surface && defined time_averaged_surface
        allocate( ou%buffs%convU(ibeg:iend,jbeg:jend) )
        allocate( ou%buffs%psiu(ibeg:iend,jbeg:jend) )
        allocate( ou%buffs%psiv(ibeg:iend,jbeg:jend) )
#  endif
#  ifdef explicit_fresh_water_flux
        allocate( ou%buffs%sfft(ibeg:iend,jbeg:jend) )
#  endif
#  ifdef stream_function
        allocate( ou%buffs%psi(ibeg:iend,jbeg:jend) )
#  endif
#  if defined(rigid_lid_surface_pressure) || defined(free_surf_on)
        allocate( ou%buffs%eta(ibeg:iend,jbeg:jend) )
#  endif
#  ifdef kppvmix
        allocate( ou%buffs%hmxl(ibeg:iend,jbeg:jend) )
        allocate( ou%buffs%hblt(ibeg:iend,jbeg:jend) )
#  endif
#  if defined bbl_ag
        allocate( ou%buffs%oppsignx(ibeg:iend,jbeg:jend) )
        allocate( ou%buffs%oppsigny(ibeg:iend,jbeg:jend) )
#  endif
# endif /* diagnostics_himem */


      endif
!-----------------------------------------------------------------------
!     endif for opening of diagnostic file
!-----------------------------------------------------------------------

# ifdef diagnostics_lomem
      if( .NOT.within_window_loop ) return
# endif /* diagnostics_lomem*/


!-----------------------------------------------------------------------
!     point "domain" to the x-y-domain from "region_domain" 
!-----------------------------------------------------------------------

# ifdef diagnostics_himem
!
!     domain was declared pointer because it needs to inherit linked list
!
      ibeg = region_domain%x%data%start_index
      iend = region_domain%x%data%end_index
      jbeg = region_domain%y%data%start_index
      jend = region_domain%y%data%end_index
      if (jbeg.gt.jecomp.or.jend.lt.jscomp)     return
      if (iend.lt.ibeg.or.jend.lt.jbeg)         return 

      domain => region_domain
# endif
# ifdef diagnostics_lomem
!
!     domain was not declared pointer
!
      domain                       = region_domain
      domain%y%global              = domain%y%compute
      domain%y%compute%start_index = js+joff
      domain%y%compute%end_index   = je+joff
      domain%y%data                = domain%y%compute
# endif

!-----------------------------------------------------------------------
!     data indexed by (i,k,j) must be re-arranged to match (i,j,k). 
!-----------------------------------------------------------------------

      if( within_window_loop )then

!     test whether the memory window contains rows which must be stored
        js_region = max(js+joff, jbeg) - joff 
        je_region = min(je+joff, jend) - joff  
        if(js_region.gt.je_region) return 
!-----------------------------------------------------------------------
!     define variables (on data domain)
!
!     incoming data is already masked to 0, but will still add 
!     (1-mask)*missing_val to be compatible with with Ferret. 
!-----------------------------------------------------------------------
        
!-----------------------------------------------------------------------
!     tracers
!-----------------------------------------------------------------------
          do n=1,nt
            do k=1,km
              do j=js_region,je_region
                do i=ibeg,iend
                   ou%buffs%trcr(i,j+joff,k,n) = t(i,k,j,n,tau)  
     &              + (1-tmask(i,k,j))*ou%fields%trcr(n)%missing
                enddo
              enddo
            enddo
          enddo
!-----------------------------------------------------------------------
!     velocities
!-----------------------------------------------------------------------
          do k=1,km
            do j=js_region,je_region
              do i=ibeg,iend
                ou%buffs%u     (i,j+joff,k) =	 u(i,k,j,1,tau)  
     &                + (1-umask(i,k,j))*ou%fields%u%missing
                ou%buffs%v     (i,j+joff,k) =	 u(i,k,j,2,tau)  
     &                + (1-umask(i,k,j))*ou%fields%v%missing
                ou%buffs%w(i,j+joff,k) = adv_vbt(i,k,j)  
     &                + (1-tmask(i,k,j))*ou%fields%w%missing
              enddo
            enddo
          enddo
!-----------------------------------------------------------------------
!        additional three-dimensional variables
!-----------------------------------------------------------------------

#  if defined gent_mcwilliams
!-----------------------------------------------------------------------
!        define the horizontal GM velocities on the velocity grid.
!
!        NOTE:  Do not require longitudinal boundary conditions at 
!        i=1,i=imt for buffer fields u_gm and v_gm
!        since these points do not get written; they are not
!        part of the computational domain.  
!-----------------------------------------------------------------------

          do k=1,km
            do j=js_region,je_region
              do i=ibeg,iend
                ou%buffs%u_gm(i,j+joff,k) = umask(i,k,j)  
     &                *(adv_vetiso(i,k,j)+adv_vetiso(i,k,j+1))/2.0  
     &                +  (1-umask(i,k,j))*ou%fields%u_gm%missing
                ou%buffs%v_gm(i,j+joff,k) = umask(i,k,j)  
     &                *(adv_vntiso(i,k,j)+adv_vntiso(i+1,k,j))/(2.0*  
     &                csu(j+joff))  
     &                +  (1-umask(i,k,j))*ou%fields%v_gm%missing
              enddo
            enddo
          enddo
          do k=1,km
            do j=js_region,je_region
              do i=ibeg,iend
                ou%buffs%w_gm(i,j+joff,k) = adv_vbtiso(i,k,j)  
     &               + (1-tmask(i,k,j))*ou%fields%w_gm%missing
              enddo
            enddo
          enddo
#  endif

!-----------------------------------------------------------------------
!     two-dimensional variables
!-----------------------------------------------------------------------
          do n=1,nt
            do j=js_region,je_region
              do i=ibeg,iend
                ou%buffs%stf(i,j+joff,n) = stf(i,j,n)  
     &                + (1-tmask(i,1,j))*ou%fields%stf(n)%missing
              enddo
            enddo
          enddo

          do n=1,2
            do j=js_region,je_region
              do i=ibeg,iend
            	ou%buffs%smf(i,j+joff,n) = smf(i,j,n)  
     &      	   + (1-umask(i,1,j))*ou%fields%taux%missing
              enddo
            enddo
          enddo
#  ifdef air_pressure
! output air presure in hPa (mbar) = 1000 dyn/cm**2
          do j=js_region,je_region
            do i=ibeg,iend
              ou%buffs%smf(i,j+joff,3) = 1.0e-3*smf(i,j,3)  
     &     	   + (1-tmask(i,1,j))*ou%fields%pair%missing
            enddo
          enddo
#  endif

#  if defined explicit_free_surface && defined time_averaged_surface 
!-----------------------------------------------------------------------
!     convU is convergence of H*ubar on t-cells. 
!-----------------------------------------------------------------------
          do j = js_region,je_region
            do i =ibeg,iend
              ou%buffs%convU(i,j+joff) = convU(i,j+joff)  
     &      	     + (1-tmask(i,1,j))*ou%fields%convU%missing
            enddo
          enddo

!-----------------------------------------------------------------------
!     psiu and psiv are streamfunctions for ubar*H and vbar*H
!     They are purposefully not masked so that they may be normalized
!     to a value over land (e.g., America) using Ferret.  
!-----------------------------------------------------------------------

          do j = js_region,je_region
            do i =ibeg,iend
              ou%buffs%psiu(i,j+joff) = psiu(i,j+joff)
!     &               + (1-tmask(i,1,j))*psiu%missing
            enddo
          enddo
          do j = js_region,je_region
            do i =ibeg,iend
              ou%buffs%psiv(i,j+joff) = psiv(i,j+joff)
!     &               + (1-tmask(i,1,j))*psiv%missing
            enddo
          enddo
#  endif

#  ifdef explicit_fresh_water_flux
!-----------------------------------------------------------------------
!         sfft is the fresh water flux on t-cells
!-----------------------------------------------------------------------
          do j = js_region,je_region
            do i =ibeg,iend
              ou%buffs%sfft(i,j+joff) = sfft(i,j+joff)  
     &              + (1-tmask(i,1,j))*ou%fields%sfft%missing
            enddo
          enddo
#  endif

#  ifdef stream_function

!-----------------------------------------------------------------------
!     psi is at "tau" since the external mode has not been updated yet
!-----------------------------------------------------------------------
          do j = js_region,je_region
            do i =ibeg,iend
              ou%buffs%psi(i,j+joff) = psi(i,j+joff,1)  
     &                + (1-tmask(i,1,j))*ou%fields%psi%missing
            enddo
          enddo
#  endif

#  if defined rigid_lid_surface_pressure || defined free_surf_on

!-----------------------------------------------------------------------
!     ps is at "tau" since the external mode has not been updated yet.
!     map the free surface height in cm.  
!-----------------------------------------------------------------------
          do j=js_region,je_region
            do i=ibeg,iend
#   if defined time_averaged_surface
              ou%buffs%eta(i,j+joff) = etat(i,j+joff,1)
#   else
              ou%buffs%eta(i,j+joff) = rho0r*ps(i,j+joff,1)/grav
#   endif  
     &              + (1-tmask(i,1,j))*ou%fields%eta%missing
            enddo
          enddo
#  endif

#  ifdef kppvmix

!-----------------------------------------------------------------------
!     hmxl is depth of surface mixed layer on t-cells. hblt is depth of 
!     planetary ocean boundary layer on t-cells
!-----------------------------------------------------------------------
          do j = js_region,je_region
            do i =ibeg,iend
              ou%buffs%hmxl(i,j+joff) = hmxl(i,j+joff)  
     &              + (1-tmask(i,1,j))*ou%fields%hmxl%missing
            enddo
          enddo
          do j = js_region,je_region
            do i =ibeg,iend
              ou%buffs%hblt(i,j+joff) = hblt(i,j+joff)  
     &              + (1-tmask(i,1,j))*ou%fields%hblt%missing
            enddo
          enddo
#  endif

#  if defined bbl_ag
          do j=js_region,je_region
            do i=ibeg,iend
              ou%buffs%oppsignx(i,j+joff) = 
     &                 oppsignx(i,j)*umask(i,1,j)  
     &              +  (1-umask(i,1,j))*ou%fields%oppsignx%missing
              ou%buffs%oppsigny(i,j+joff) = 
     &                 oppsigny(i,j)*umask(i,1,j)  
     &              +  (1-umask(i,1,j))*ou%fields%oppsigny%missing
            enddo
          enddo
#  endif
             
      endif  !end of within_window_loop if-test
!    
!     Data is written only outside window loop for himem
!     but within window loop for lomem
!
# ifdef diagnostics_himem
      if (.not. within_window_loop) then
# endif
!-----------------------------------------------------------------------
!         tracers
!-----------------------------------------------------------------------
          do n=1,nt
            call mpp_write(io(nr),ou%fields%trcr(n), domain, 
     &	                   ou%buffs%trcr(:,:,:,n),  timrec)
          enddo

!-----------------------------------------------------------------------
!         velocities
!-----------------------------------------------------------------------
          call mpp_write(io(nr), ou%fields%u, domain, 
     &	                   ou%buffs%u     , timrec)
          call mpp_write(io(nr), ou%fields%v, domain, 
     &	                   ou%buffs%v     , timrec)
          call mpp_write(io(nr), ou%fields%w, domain, 
     &	                   ou%buffs%w, timrec) 

!-----------------------------------------------------------------------
!         additional three-dimensional variables
!-----------------------------------------------------------------------

#  if defined gent_mcwilliams
!-----------------------------------------------------------------------
!         define the horizontal GM velocities on the velocity grid.
!-----------------------------------------------------------------------
          call mpp_write(io(nr), ou%fields%u_gm, domain, 
     &	                   ou%buffs%u_gm, timrec)
          call mpp_write(io(nr), ou%fields%v_gm, domain, 
     &	                   ou%buffs%v_gm, timrec)
          call mpp_write(io(nr), ou%fields%w_gm, domain, 
     &	                   ou%buffs%w_gm,timrec)
#  endif

!-----------------------------------------------------------------------
!         two-dimensional variables
!-----------------------------------------------------------------------
          do n=1,nt
            call mpp_write(io(nr), ou%fields%stf(n), domain, 
     &	                   ou%buffs%stf(:,:,n),  timrec)
          enddo

          call mpp_write(io(nr),ou%fields%taux, domain, 
     &	                   ou%buffs%smf(:,:,1),timrec)
          call mpp_write(io(nr),ou%fields%tauy, domain, 
     &	                   ou%buffs%smf(:,:,2),timrec)
#  ifdef air_pressure
          call mpp_write(io(nr),ou%fields%pair, domain, 
     &	                   ou%buffs%smf(:,:,3),timrec)
#  endif
#  if defined explicit_free_surface && defined time_averaged_surface

!-----------------------------------------------------------------------
!         convU is convergence of H*ubar on t-cells 
!-----------------------------------------------------------------------
          call mpp_write(io(nr), ou%fields%convU,domain,
     &	                 ou%buffs%convU, timrec)

!-----------------------------------------------------------------------
!         psiu and psiv are streamfunctions for ubar*H and vbar*H
!-----------------------------------------------------------------------
          call mpp_write(io(nr), ou%fields%psiu, domain, 
     &	                 ou%buffs%psiu, timrec)
          call mpp_write(io(nr), ou%fields%psiv, domain, 
     &	                 ou%buffs%psiv, timrec)
#  endif

#  ifdef explicit_fresh_water_flux
!-----------------------------------------------------------------------
!         sfft is the fresh water flux on t-cells
!-----------------------------------------------------------------------
          call mpp_write(io(nr), ou%fields%sfft, domain, 
     &	                 ou%buffs%sfft, timrec)
#  endif

#  ifdef stream_function

!-----------------------------------------------------------------------
!         psi is at "tau" as the external mode has not been updated yet
!-----------------------------------------------------------------------
          call mpp_write(io(nr), ou%fields%psi, domain, 
     &	                 ou%buffs%psi, timrec)
#  endif

#  if defined rigid_lid_surface_pressure || defined free_surf_on

!-----------------------------------------------------------------------
!         ps is at "tau" as the external mode has not been updated yet.
!-----------------------------------------------------------------------
          call mpp_write(io(nr), ou%fields%eta, domain, 
     &	                 ou%buffs%eta, timrec)
#  endif

#  ifdef kppvmix

!-----------------------------------------------------------------------
!         hmxl is depth of surface mixed layer on t-cells. hblt is
!         depth of planetary ocean boundary layer on t-cells
!-----------------------------------------------------------------------
          call mpp_write(io(nr), ou%fields%hmxl, domain, 
     &	                 ou%buffs%hmxl, timrec)
          call mpp_write(io(nr), ou%fields%hblt, domain, 
     &	                 ou%buffs%hblt, timrec)
#  endif

#  if defined bbl_ag 
          call mpp_write(io(nr), ou%fields%oppsignx, domain,  
     &	                 ou%buffs%oppsignx, timrec)
          call mpp_write(io(nr), ou%fields%oppsigny, domain,  
     &	                 ou%buffs%oppsigny, timrec)
#  endif

!-----------------------------------------------------------------------
!       write out averaging period 
!-----------------------------------------------------------------------
        call mpp_write(io(nr), ou%fields%period, period, timrec)

        if(js+joff .eq. jscomp) then
          write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8,a,i8)')  
     & '    ==>  Regions written to file ', file_name(nr)  
     & ,'   on ts = ', itt,' ', stamp  
     & ,'   between latitudes:',jscomp,' and ', jecomp  
     & ,'   and from the surface down to depth level ',km
        endif

# ifdef diagnostics_himem
      endif
# endif

# ifdef timing
      call toc ('diagnostic', 'snapshots')
# endif

      return
      end

#else
      return
      end
#endif


