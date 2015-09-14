!=======================================================================
! Write out the snapshot in either netCDF or IEEE
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
! It runs many times slower. 

! For himem: call snapshots_out whenever within_window_loop=.TRUE. 
! This is done from within the window loop (diag_insidemw.F), and 
! if within_window_loop=.FALSE. from outside the window 
! (in diag_out_window.F).
! For lomem: do same, and write to disk each window as it comes.
!
! Buffering is done when within_window_loop=.TRUE. 
! In lomem case, data is also written if within_window_loop=.TRUE.
!
! Each data record is preceded by a character record describing 
! the data record. The latitude and depth range for saving data is 
! in "diag.h" and is initialized via "blkdta.F"
!
!     authors:   R. C. Pacanowski  e-mail=> rcp@gfdl.gov
!                A. Rosati         e-mail=> ar@gfdl.gov
!                M. Harrison       e-mail=> mjh@gfdl.gov
!
!                Updated with mpp_io routines Summer 1999 by  
!                C. Kerr           e-mail=> ck@gfdl.gov 
!                Balaji            e-mail=> vb@gfdl.gov
!                S. M. Griffies    e-mail=> smg@gfdl.gov  
!=======================================================================
      subroutine snapshots_out (joff, js, je, within_window_loop)
#if defined snapshots 
c
c-----------------------------------------------------------------------
c     save snapshot data
c-----------------------------------------------------------------------
c
      use arrays_2d
      use memory_window
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod

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

      character*128 file_name, file_stamp, gvtitle
      character(len=16) :: cvarnam

# if defined salinity_psu
      character*3 saldim
# else 
      character*9 saldim
# endif

      logical :: within_window_loop
      dimension bufzt(km), bufzw(km)

      data num_snaps /0/
      save num_snaps
      save io

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------

      type(fieldtype), save :: field_tracer(nt)
      type(fieldtype), save :: field_stf(nt)
      type(fieldtype), save :: field_u
      type(fieldtype), save :: field_v
      type(fieldtype), save :: field_w
      type(fieldtype), save :: field_convU

      type(fieldtype), save :: field_ustar
      type(fieldtype), save :: field_vstar
      type(fieldtype), save :: field_wstar

      type(fieldtype), save :: field_sfft
      type(fieldtype), save :: field_taux
      type(fieldtype), save :: field_tauy
#  ifdef air_pressure
      type(fieldtype), save :: field_pair
#  endif

      type(fieldtype), save :: field_eta
      type(fieldtype), save :: field_psiu
      type(fieldtype), save :: field_psiv
      type(fieldtype), save :: field_psi

      type(fieldtype), save :: field_hblt
      type(fieldtype), save :: field_hmxl
      type(fieldtype), save :: field_oppsignx
      type(fieldtype), save :: field_oppsigny

      type(fieldtype), save :: field_period

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_xt
      type (axistype), save :: axis_xt_edges
      type (axistype), save :: axis_xu
      type (axistype), save :: axis_xu_edges
      type (axistype), save :: axis_yt
      type (axistype), save :: axis_yt_edges
      type (axistype), save :: axis_yu
      type (axistype), save :: axis_yu_edges
      type (axistype), save :: axis_zt
      type (axistype), save :: axis_zt_edges
      type (axistype), save :: axis_zw
      type (axistype), save :: axis_zw_edges
      type (axistype), save :: axis_surf
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     axis edge arrays
c-----------------------------------------------------------------------

      real, dimension(0:imt), save :: edges_xt, edges_xu
      real, dimension(0:jmt), save :: edges_yt, edges_yu
      real, dimension(0:km ), save :: edges_zt, edges_zw

# ifdef diagnostics_himem
      type(domain2D), pointer :: domain

      real, dimension(:,:,:,:), save, allocatable :: buf_t

      real, dimension(:,:,:), save, allocatable :: buf_u1
      real, dimension(:,:,:), save, allocatable :: buf_u2
      real, dimension(:,:,:), save, allocatable :: buf_adv_vbt
      real, dimension(:,:,:), save, allocatable :: buf_u_gm
      real, dimension(:,:,:), save, allocatable :: buf_v_gm
      real, dimension(:,:,:), save, allocatable :: buf_adv_vbtiso

      real, dimension(:,:,:), save, allocatable :: buf_stf
      real, dimension(:,:,:), save, allocatable :: buf_smf

      real, dimension(:,:), save, allocatable :: buf_surf_ht
      real, dimension(:,:), save, allocatable :: buf_oppsignx
      real, dimension(:,:), save, allocatable :: buf_oppsigny
      real, dimension(:,:), save, allocatable :: buf_convU
      real, dimension(:,:), save, allocatable :: buf_psiu
      real, dimension(:,:), save, allocatable :: buf_psiv
      real, dimension(:,:), save, allocatable :: buf_psi
      real, dimension(:,:), save, allocatable :: buf_sfft
      real, dimension(:,:), save, allocatable :: buf_hmxl
      real, dimension(:,:), save, allocatable :: buf_hblt
# endif

# ifdef diagnostics_lomem
      type(domain2D) :: domain
c
c     dimension within mw and set domain to a domain2D type
c
      real :: buf_t(1:imt,js+joff:je+joff,km,nt)
      real :: buf_u1(imt,js+joff:je+joff,km)
      real :: buf_u2(imt,js+joff:je+joff,km)
      real :: buf_adv_vbt(imt,js+joff:je+joff,km)
#  ifdef gent_mcwilliams
      real :: buf_u_gm(imt,js+joff:je+joff,km)
      real :: buf_v_gm(imt,js+joff:je+joff,km)
      real :: buf_adv_vbtiso(imt,js+joff:je+joff,km)
#  endif
      real :: buf_stf(imt,js+joff:je+joff,nt)
#  ifdef air_pressure
      real :: buf_smf(imt,js+joff:je+joff,3 )
#  else
      real :: buf_smf(imt,js+joff:je+joff,2 )
#  endif
#  ifdef explicit_free_surface && defined time_averaged_surface
      real :: buf_psiu(imt,js+joff:je+joff)
      real :: buf_psiv(imt,js+joff:je+joff)
      real :: buf_convU(imt,js+joff:je+joff)
#  endif
#  ifdef explicit_fresh_water_flux
      real :: buf_sfft(imt,js+joff:je+joff)
#  endif
#  ifdef stream_function
      real :: buf_psi(imt,js+joff:je+joff)
#  endif
#  if defined(rigid_lid_surface_pressure) || defined(free_surf_on)
      real :: buf_surf_ht(imt,js+joff:je+joff)
#  endif
#  ifdef kppvmix
      real :: buf_hmxl(imt,js+joff:je+joff)
      real :: buf_hblt(imt,js+joff:je+joff)
#  endif
#  if defined bbl_ag
      real :: buf_oppsignx(imt,js+joff:je+joff)
      real :: buf_oppsigny(imt,js+joff:je+joff)
#  endif
# endif


# ifdef timing
      call tic ('diagnostic', 'snapshots ')
# endif

c-----------------------------------------------------------------------
c     set averaging period, time, title, ocean depth, salinity units 
c-----------------------------------------------------------------------

      file_name = file_stamp ('snapshots', xstamp, '.dta')
      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Snapshots' // momver
# if defined salinity_psu
      saldim = 'psu' 
# else 
      saldim = 'ppm-0.035' 
# endif
      bufzt = zt(:)*0.01
      bufzw = zw(:)*0.01

c-----------------------------------------------------------------------
c     open diagnostic file and write meta information 
c-----------------------------------------------------------------------
 
      if ( num_snaps .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_snaps = num_snaps + 1

c-----------------------------------------------------------------------
c       axes (on global domain)
c-----------------------------------------------------------------------

        ibeg = mom_domain%x%global%start_index
        iend = mom_domain%x%global%end_index
        jbeg = mom_domain%y%global%start_index
        jend = mom_domain%y%global%end_index

        call mpp_write_meta(
     &       io, axis_xt, 'xt_i', 'degrees_E', 'Longitude of T points',
     &       cartesian='X', domain=mom_domain%x, 
     &       data=xt(ibeg:iend))

        call mpp_write_meta(
     &         io, axis_yt, 'yt_j', 'degrees_N', 'Latitude of T points',
     &         cartesian='Y', domain=mom_domain%y, 
     &         data=yt(jbeg:jend))

        call mpp_write_meta(
     &         io, axis_zt, 'zt_k', 'm', 'Depth of T grid point',
     &         cartesian='Z', sense=-1, data=bufzt)
 
        call mpp_write_meta(
     &       io, axis_xu, 'xu_i', 'degrees_E', 'Longitude of U points',
     &       cartesian='X', domain=mom_domain%x, 
     &       data=xu(ibeg:iend))

        call mpp_write_meta(
     &         io, axis_yu, 'yu_j', 'degrees_N', 'Latitude of U points',
     &         cartesian='Y', domain=mom_domain%y, 
     &         data=yu(jbeg:jend))

        call mpp_write_meta( 
     &         io, axis_zw, 'zw_k', 'm', 'Depth of T cell bottom',
     &         cartesian='Z', sense=-1, data=bufzw)
 
        call mpp_write_meta(
     &         io, axis_surf, 'Surface', 'm', 'Depth of surface',
     &         cartesian='Z', sense=-1, data=(/0./))

        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')


c-----------------------------------------------------------------------
c       define auxiliary "edges" axes for xt, yt, zt coordinates
c-----------------------------------------------------------------------

        call edge_maker (1, edges_xt, xt, xu, imt)
        call mpp_write_meta(
     &         io, axis_xt_edges, 'xt_i_edges', 'degrees_E', 
     &         'Longitude of T edges', cartesian='X'
     &,        data=edges_xt)

        call edge_maker (1, edges_yt, yt, yu, jmt)
        call mpp_write_meta(
     &         io, axis_yt_edges, 'yt_j_edges', 'degrees_N', 
     &         'Latitude of T edges', cartesian='Y'
     &,        data=edges_yt)
     
        call edge_maker (1, edges_zt, zt, zw, km)
        call mpp_write_meta(
     &         io, axis_zt_edges, 'zt_k_edges', 'cm', 
     &         'Depth of T edges', cartesian='Z', sense=-1
     &,        data=edges_zt)

c-----------------------------------------------------------------------
c       define auxiliary "edges" axes for xu, yu, zw coordinates
c-----------------------------------------------------------------------

        call edge_maker (2, edges_xu, xt, xu, imt)
        call mpp_write_meta(
     &         io, axis_xu_edges, 'xu_i_edges', 'degrees_E', 
     &         'Longitude of U edges', cartesian='X'
     &,         data=edges_xu)

        call edge_maker (2, edges_yu, yt, yu, jmt)
        call mpp_write_meta(
     &         io, axis_yu_edges, 'yu_j_edges', 'degrees_N', 
     &         'Latitude of U edges', cartesian='Y'
     &,         data=edges_yu)
    
        call edge_maker (2, edges_zw, zt, zw, km)
        call mpp_write_meta(
     &         io, axis_zw_edges, 'zw_k_edges', 'cm', 
     &         'Depth of W edges', cartesian='Z'
     &,         sense=-1, data=edges_zw)

c-----------------------------------------------------------------------
c       define fields (on data domain)
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c       tracers
c-----------------------------------------------------------------------

        do n=1,nt 
          if (n .EQ. 1) then
            call mpp_write_meta( 
     &        io, field_tracer(n), (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &        'temp', 'deg C', 'potential temperature',  
     &        -5.0, 50.0, -1.0E+34, +1.0E+34)

          elseif (n .EQ. 2) then
            call mpp_write_meta(
     &        io, field_tracer(n), (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &        'salinity', saldim, 'salinity', 
     &        -100.00, +100.00, -1.0E+34, +1.0E+34)
          else
            write( cvarnam, '(a,i2.2)' ) 'tracer_', n
            call mpp_write_meta(
     &        io, field_tracer(n),(/axis_xt,axis_yt,axis_zt,axis_t/), 
     &        cvarnam, 'tracer units', trname(n), 
     &        -1.0E+10, +1.0E+10, -1.0E+34, +1.0E+34)
          endif
        enddo

c-----------------------------------------------------------------------
c       velocities
c-----------------------------------------------------------------------

        call mpp_write_meta(
     &      io, field_u, (/axis_xu,axis_yu,axis_zt,axis_t/), 
     &      'u', 'cm/s', 'Zonal velocity',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
 
        call mpp_write_meta(
     &      io, field_v, (/axis_xu,axis_yu,axis_zt,axis_t/), 
     &      'v', 'cm/s', 'Meridional velocity', 
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
 
        call mpp_write_meta(
     &      io, field_w, (/axis_xt,axis_yt,axis_zw,axis_t/), 
     &      'w', 'cm/s', 'W at T cell bottom', 
     &      -100.0, +100.0, -1.0E+34, +1.0E+34)
 
c-----------------------------------------------------------------------
c       additional three-dimensional variables
c-----------------------------------------------------------------------

#  if defined gent_mcwilliams
        call mpp_write_meta(
     &      io, field_ustar, (/axis_xu,axis_yu,axis_zt,axis_t/), 
     &      'ustar', 'cm/s', 'Zonal GM velocity on U-cell',
     &      -1.0E+06, +1.0E+06, -1.0E+34, +1.0E+34)

        call mpp_write_meta(
     &      io, field_vstar, (/axis_xu,axis_yu,axis_zt,axis_t/), 
     &      'vstar', 'cm/s', 'Meridional GM velocity on U-cell',
     &      -1.0E+06, +1.0E+06, -1.0E+34, +1.0E+34)

        call mpp_write_meta(
     &      io, field_wstar, (/axis_xt,axis_yt,axis_zw,axis_t/), 
     &      'wstar', 'cm/s', 'Vertical GM velocity at T-cell bottom',
     &      -1.0E+06, +1.0E+06, -1.0E+34, +1.0E+34)
#  endif

c-----------------------------------------------------------------------
c       two-dimensional variables
c-----------------------------------------------------------------------

        do n=1,nt
          if(n .eq. 1) then
            call mpp_write_meta(
     &      io, field_stf(1), (/axis_xt,axis_yt,axis_surf,axis_t/), 
     &      'hflx', 'cal/cm**2/s', 'Surface heat flux',
     &      -1.0E+00, +1.0E+00, -1.0E+34, +1.0E+34)

          elseif(n.eq.2) then
            call mpp_write_meta(
     &      io, field_stf(2), (/axis_xt,axis_yt,axis_surf,axis_t/), 
     &      'sflx', 'gm/cm**2/s', 'Surface salt flux',
     &      -1.0E+00, +1.0E+00, -1.0E+34, +1.0E+34)

          else 
            write( cvarnam, '(a,i2.2)' ) 'stf_', n
            call mpp_write_meta(
     &        io, field_stf(n), (/axis_xt,axis_yt,axis_surf,axis_t/), 
     &        cvarnam, 'tracer units', 'Surface tracer flux',
     &        -1.0E+06, +1.0E+06, -1.0E+34, +1.0E+34)
          endif
        enddo

        call mpp_write_meta(
     &      io, field_taux, (/axis_xu,axis_yu,axis_surf,axis_t/), 
     &      'taux', 'dyn/cm**2', 'Zonal windstress',
     &      -20.0, +20.0, -1.0E+34, +1.0E+34)

        call mpp_write_meta(
     &      io, field_tauy, (/axis_xu,axis_yu,axis_surf,axis_t/), 
     &      'tauy', 'dyn/cm**2', 'Meridional windstress',
     &      -20.0, +20.0, -1.0E+34, +1.0E+34)

#  ifdef air_pressure
        call mpp_write_meta(
     &      io, field_pair, (/axis_xt,axis_yt,axis_surf,axis_t/), 
     &      'pair', 'hPa', 'Air pressure',
     &      0.0, +2000.0, -1.0E+34, +1.0E+34)
#  endif

#  if defined explicit_free_surface && defined time_averaged_surface
       call mpp_write_meta(
     &      io, field_convU, (/axis_xt,axis_yt,axis_surf,axis_t/), 
     &      'convU', 'cm/s', 'convergence of H*ubar',
     &      -300000.0, +300000.0, -1.0E+34, +1.0E+34)
       call mpp_write_meta(
     &      io, field_psiu, (/axis_xt,axis_yt,axis_surf,axis_t/), 
     &      'psiu', 'Sv', 'Streamfunction for H*ubar',
     &      -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)

       call mpp_write_meta(
     &      io, field_psiv, (/axis_xt,axis_yt,axis_surf,axis_t/), 
     &      'psiv', 'Sv', 'Streamfunction for H*vbar',
     &      -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)
#  endif

#  if defined explicit_fresh_water_flux
        call mpp_write_meta(
     &      io, field_sfft, (/axis_xt,axis_yt,axis_surf,axis_t/), 
     &      'sff', 'cm/s', 'Fresh water flux',
     &      -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)
#  endif

#  if defined stream_function
        call mpp_write_meta(
     &      io, field_psi, (/axis_xt,axis_yt,axis_t/), 
     &      'psi', 'cm**3/s', 'Stream function',
     &      -300.0, +300.0, -1.0E+34, +1.0E+34)
#  endif

#  if defined rigid_lid_surface_pressure || defined free_surf_on
        call mpp_write_meta(
     &      io, field_eta, (/axis_xt,axis_yt,axis_t/), 
     &      'eta', 'cm', 'Surface height',
     &      -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)
#  endif

#  if defined kppvmix
        call mpp_write_meta(
     &      io, field_hmxl, (/axis_xt,axis_yt,axis_t/), 
     &      'hmxl', 'cm', 'Mixed layer depth',
     &      -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)

        call mpp_write_meta(
     &      io, field_hblt, (/axis_xt,axis_yt,axis_t/), 
     &      'hblt', 'cm', 'Boundary layer depth',
     &      -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)
#  endif

#  if defined bbl_ag
        call mpp_write_meta(
     &      io, field_oppsignx, (/axis_xu,axis_yu,axis_t/), 
     &      'oppsignx', 'cm', 'x-detrainment indicator in bbl',
     &      -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)

        call mpp_write_meta(
     &      io, field_oppsigny, (/axis_xu,axis_yu,axis_t/), 
     &      'oppsigny', 'cm', 'y-detrainment indicator in bbl',
     &      -3.0E+06, +3.0E+06, -1.0E+34, +1.0E+34)
#  endif
c
c       averaging period 
c
        call mpp_write_meta(
     &      io, field_period, (/axis_t/), 
     &      'period', 'day', 'averaging period',
     &      +0.0E+00, +1.0E+20, -1.0E+34, +1.0E+34)

c
c       figure title 
c
        call mpp_write_meta(io, 'title', cval=gvtitle )
c
c-----------------------------------------------------------------------
c     write axis now that all meta data has been written
c-----------------------------------------------------------------------

        call mpp_write (io,axis_surf)
        call mpp_write (io,axis_xt)
        call mpp_write (io,axis_xu)
        call mpp_write (io,axis_yt)
        call mpp_write (io,axis_yu)
        call mpp_write (io,axis_zt)
        call mpp_write (io,axis_zw)

        call mpp_write (io,axis_xt_edges)
        call mpp_write (io,axis_xu_edges)
        call mpp_write (io,axis_yt_edges)
        call mpp_write (io,axis_yu_edges)
        call mpp_write (io,axis_zt_edges)
        call mpp_write (io,axis_zw_edges)

c-----------------------------------------------------------------------
c       allocate himem arrays over full 3d domain on processor
c-----------------------------------------------------------------------

# ifdef diagnostics_himem
        allocate( buf_t(1:imt,jstask:jetask,km,nt) )
        allocate( buf_u1(imt,jstask:jetask,km) )
        allocate( buf_u2(imt,jstask:jetask,km) )
        allocate( buf_adv_vbt(imt,jstask:jetask,km) )
#  ifdef gent_mcwilliams
        allocate( buf_u_gm(imt,jstask:jetask,km) )
        allocate( buf_v_gm(imt,jstask:jetask,km) )
        allocate( buf_adv_vbtiso(imt,jstask:jetask,km) )
#  endif
        allocate( buf_stf(imt,jstask:jetask,nt) )
#  ifdef air_pressure
        allocate( buf_smf(imt,jstask:jetask,3 ) )
#  else
        allocate( buf_smf(imt,jstask:jetask,2 ) )
#  endif
#  if defined explicit_free_surface && defined time_averaged_surface
        allocate( buf_convU(imt,jstask:jetask) )
        allocate( buf_psiu(imt,jstask:jetask) )
        allocate( buf_psiv(imt,jstask:jetask) )
#  endif
#  ifdef explicit_fresh_water_flux
        allocate( buf_sfft(imt,jstask:jetask) )
#  endif
#  ifdef stream_function
        allocate( buf_psi(imt,jstask:jetask) )
#  endif
#  if defined(rigid_lid_surface_pressure) || defined(free_surf_on)
        allocate( buf_surf_ht(imt,jstask:jetask) )
#  endif
#  ifdef kppvmix
        allocate( buf_hmxl(imt,jstask:jetask) )
        allocate( buf_hblt(imt,jstask:jetask) )
#  endif
#  if defined bbl_ag
        allocate( buf_oppsignx(imt,jstask:jetask) )
        allocate( buf_oppsigny(imt,jstask:jetask) )
#  endif
# endif diagnostics_himem


      endif
c-----------------------------------------------------------------------
c     endif for opening of diagnostic file
c-----------------------------------------------------------------------

# ifdef diagnostics_lomem
      if( .NOT.within_window_loop ) return
# endif diagnostics_lomem


c-----------------------------------------------------------------------
c     point "domain" to the x-y-domain from "mom_domain" 
c-----------------------------------------------------------------------

# ifdef diagnostics_himem
c
c     domain was declared pointer because it needs to inherit linked list
c
      domain => mom_domain
# endif
# ifdef diagnostics_lomem
c
c     domain was not declared pointer
c
      domain                       = mom_domain
      domain%y%global              = domain%y%compute
      domain%y%compute%start_index = js+joff
      domain%y%compute%end_index   = je+joff
      domain%y%data                = domain%y%compute
# endif

c-----------------------------------------------------------------------
c     data indexed by (i,k,j) must be re-arranged to match (i,j,k). 
c-----------------------------------------------------------------------

      if( within_window_loop )then

c-----------------------------------------------------------------------
c     define variables (on data domain)
c
c     incoming data is already masked to 0, but will still add 
c     (1-mask)*missing_val to be compatible with with Ferret. 
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     tracers
c-----------------------------------------------------------------------
          do n=1,nt
             do k=1,km
                do j=js,je
                   do i=1,imt
                      buf_t(i,j+joff,k,n) = t(i,k,j,n,tau)
     &                     + (1-tmask(i,k,j))*field_tracer(n)%missing
                   enddo
                enddo
             enddo
          enddo
c-----------------------------------------------------------------------
c     velocities
c-----------------------------------------------------------------------
          do k=1,km
             do j=js,je
                do i=1,imt
                   buf_u1     (i,j+joff,k) =       u(i,k,j,1,tau)
     &                  + (1-umask(i,k,j))*field_u%missing
                   buf_u2     (i,j+joff,k) =       u(i,k,j,2,tau)
     &                  + (1-umask(i,k,j))*field_v%missing
                   buf_adv_vbt(i,j+joff,k) = adv_vbt(i,k,j)
     &                  + (1-tmask(i,k,j))*field_w%missing
                enddo
             enddo
          enddo
c-----------------------------------------------------------------------
c        additional three-dimensional variables
c-----------------------------------------------------------------------

#  if defined gent_mcwilliams
c-----------------------------------------------------------------------
c        define the horizontal GM velocities on the velocity grid.
c
c        NOTE:  Do not require longitudinal boundary conditions at 
c        i=1,i=imt for buffer fields buf_u_gm and buf_v_gm
c        since these points do not get written; they are not
c        part of the computational domain.  
c-----------------------------------------------------------------------

          do k=1,km
             do j=js,je
                do i=2,imt-1
                   buf_u_gm(i,j+joff,k) = umask(i,k,j)
     &                  *(adv_vetiso(i,k,j)+adv_vetiso(i,k,j+1))/2.0
     &                  +  (1-umask(i,k,j))*field_ustar%missing
                   buf_v_gm(i,j+joff,k) = umask(i,k,j)
     &                  *(adv_vntiso(i,k,j)+adv_vntiso(i+1,k,j))/(2.0*
     &                  csu(j+joff))
     &                  +  (1-umask(i,k,j))*field_vstar%missing
                enddo
             enddo
          enddo
          do k=1,km
             do j=js,je
                do i=1,imt
                   buf_adv_vbtiso(i,j+joff,k) = adv_vbtiso(i,k,j)
     &                  + (1-tmask(i,k,j))*field_wstar%missing
                enddo
             enddo
          enddo
#  endif

c-----------------------------------------------------------------------
c     two-dimensional variables
c-----------------------------------------------------------------------
          do n=1,nt
             do j=js,je
                do i=1,imt
                   buf_stf(i,j+joff,n) = stf(i,j,n)
     &                  + (1-tmask(i,1,j))*field_stf(n)%missing
                enddo
             enddo
          enddo

          do n=1,2
             do j=js,je
                do i=1,imt
                   buf_smf(i,j+joff,n) = smf(i,j,n)
     &                  + (1-umask(i,1,j))*field_taux%missing
                enddo
             enddo
          enddo
#  ifdef air_pressure
! output air presure in hPa (mbar) = 1000 dyn/cm**2
          do j=js,je
             do i=1,imt
          	buf_smf(i,j+joff,3) = 1.0e-3*smf(i,j,3)
     &    	     + (1-tmask(i,1,j))*field_pair%missing
             enddo
          enddo
#  endif

#  if defined explicit_free_surface && defined time_averaged_surface 
c-----------------------------------------------------------------------
c     convU is convergence of H*ubar on t-cells. 
c-----------------------------------------------------------------------
          do j = js,je
             do i =1,imt
                buf_convU(i,j+joff) = convU(i,j+joff)
     &               + (1-tmask(i,1,j))*field_convU%missing
             enddo
          enddo

c-----------------------------------------------------------------------
c     psiu and psiv are streamfunctions for ubar*H and vbar*H
c     They are purposefully not masked so that they may be normalized
c     to a value over land (e.g., America) using Ferret.  
c-----------------------------------------------------------------------

          do j = js,je
             do i =1,imt
                buf_psiu(i,j+joff) = psiu(i,j+joff)
c     &               + (1-tmask(i,1,j))*field_psiu%missing
             enddo
          enddo
          do j = js,je
             do i =1,imt
                buf_psiv(i,j+joff) = psiv(i,j+joff)
c     &               + (1-tmask(i,1,j))*field_psiv%missing
             enddo
          enddo
#  endif

#  ifdef explicit_fresh_water_flux
c-----------------------------------------------------------------------
c         sfft is the fresh water flux on t-cells
c-----------------------------------------------------------------------
          do j = js,je
             do i =1,imt
                buf_sfft(i,j+joff) = sfft(i,j+joff)
     &               + (1-tmask(i,1,j))*field_sfft%missing
             enddo
          enddo
#  endif

#  ifdef stream_function

c-----------------------------------------------------------------------
c     psi is at "tau" since the external mode has not been updated yet
c-----------------------------------------------------------------------
          do j = js,je
             do i =1,imt
                buf_psi(i,j+joff) = psi(i,j+joff,1)
     &               + (1-tmask(i,1,j))*field_psi%missing
             enddo
          enddo
#  endif

#  if defined rigid_lid_surface_pressure || defined free_surf_on

c-----------------------------------------------------------------------
c     ps is at "tau" since the external mode has not been updated yet.
c     map the free surface height in cm.  
c-----------------------------------------------------------------------
          do j=js,je
             do i=1,imt
#   if defined time_averaged_surface
                buf_surf_ht(i,j+joff) = etat(i,j+joff,1)
#   else
                buf_surf_ht(i,j+joff) = rho0r*ps(i,j+joff,1)/grav
#   endif
     &               + (1-tmask(i,1,j))*field_eta%missing
             enddo
          enddo
#  endif

#  ifdef kppvmix

c-----------------------------------------------------------------------
c     hmxl is depth of surface mixed layer on t-cells. hblt is depth of 
c     planetary ocean boundary layer on t-cells
c-----------------------------------------------------------------------
          do j = js,je
             do i =1,imt
                buf_hmxl(i,j+joff) = hmxl(i,j+joff)
     &               + (1-tmask(i,1,j))*field_hmxl%missing
             enddo
          enddo
          do j = js,je
             do i =1,imt
                buf_hblt(i,j+joff) = hblt(i,j+joff)
     &               + (1-tmask(i,1,j))*field_hblt%missing
             enddo
          enddo
#  endif

#  if defined bbl_ag
          do j=js,je
             do i=1,imt
                buf_oppsignx(i,j+joff) = oppsignx(i,j)*umask(i,1,j)
     &               +  (1-umask(i,1,j))*field_oppsignx%missing
                buf_oppsigny(i,j+joff) = oppsigny(i,j)*umask(i,1,j)
     &               +  (1-umask(i,1,j))*field_oppsigny%missing
             enddo
          enddo
#  endif
             
      endif  !end of within_window_loop if-test
c    
c     Data is written only outside window loop for himem
c     but within window loop for lomem
c
# ifdef diagnostics_himem
      if (.not. within_window_loop) then
# endif
c-----------------------------------------------------------------------
c         tracers
c-----------------------------------------------------------------------
          do n=1,nt
            call mpp_write(io,field_tracer(n), domain, buf_t(:,:,:,n), 
     &                     timrec)
          enddo

c-----------------------------------------------------------------------
c         velocities
c-----------------------------------------------------------------------
          call mpp_write(io, field_u, domain, buf_u1     , timrec)
          call mpp_write(io, field_v, domain, buf_u2     , timrec)
          call mpp_write(io, field_w, domain, buf_adv_vbt, timrec) 

c-----------------------------------------------------------------------
c         additional three-dimensional variables
c-----------------------------------------------------------------------

#  if defined gent_mcwilliams
c-----------------------------------------------------------------------
c         define the horizontal GM velocities on the velocity grid.
c-----------------------------------------------------------------------
          call mpp_write(io, field_ustar, domain, buf_u_gm, timrec)
          call mpp_write(io, field_vstar, domain, buf_v_gm, timrec)
          call mpp_write(io,field_wstar,domain,buf_adv_vbtiso,timrec)
#  endif

c-----------------------------------------------------------------------
c         two-dimensional variables
c-----------------------------------------------------------------------
          do n=1,nt
            call mpp_write(io, field_stf(n), domain, buf_stf(:,:,n),
     &                     timrec)
          enddo

          call mpp_write(io,field_taux, domain, buf_smf(:,:,1),timrec)
          call mpp_write(io,field_tauy, domain, buf_smf(:,:,2),timrec)
#  ifdef air_pressure
          call mpp_write(io,field_pair, domain, buf_smf(:,:,3),timrec)
#  endif
#  if defined explicit_free_surface && defined time_averaged_surface

c-----------------------------------------------------------------------
c         convU is convergence of H*ubar on t-cells 
c-----------------------------------------------------------------------
          call mpp_write(io, field_convU,domain,buf_convU, timrec)

c-----------------------------------------------------------------------
c         psiu and psiv are streamfunctions for ubar*H and vbar*H
c-----------------------------------------------------------------------
          call mpp_write(io, field_psiu, domain, buf_psiu, timrec)
          call mpp_write(io, field_psiv, domain, buf_psiv, timrec)
#  endif

#  ifdef explicit_fresh_water_flux
c-----------------------------------------------------------------------
c         sfft is the fresh water flux on t-cells
c-----------------------------------------------------------------------
          call mpp_write(io, field_sfft, domain, buf_sfft, timrec)
#  endif

#  ifdef stream_function

c-----------------------------------------------------------------------
c         psi is at "tau" as the external mode has not been updated yet
c-----------------------------------------------------------------------
          call mpp_write(io, field_psi, domain, buf_psi, timrec)
#  endif

#  if defined rigid_lid_surface_pressure || defined free_surf_on

c-----------------------------------------------------------------------
c         ps is at "tau" as the external mode has not been updated yet.
c-----------------------------------------------------------------------
          call mpp_write(io, field_eta, domain, buf_surf_ht, timrec)
#  endif

#  ifdef kppvmix

c-----------------------------------------------------------------------
c         hmxl is depth of surface mixed layer on t-cells. hblt is
c         depth of planetary ocean boundary layer on t-cells
c-----------------------------------------------------------------------
          call mpp_write(io, field_hmxl, domain, buf_hmxl, timrec)
          call mpp_write(io, field_hblt, domain, buf_hblt, timrec)
#  endif

#  if defined bbl_ag 
          call mpp_write(io, field_oppsignx, domain, 
     &         buf_oppsignx, timrec)
          call mpp_write(io, field_oppsigny, domain, 
     &         buf_oppsigny, timrec)
#  endif

c-----------------------------------------------------------------------
c       write out averaging period 
c-----------------------------------------------------------------------
        call mpp_write(io, field_period, period, timrec)

        if(js+joff .eq. jscomp) then
          write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8,a,i8)')
     &'    ==>  Snapshots written to file ', trim(file_name)
     &,'   on ts = ', itt,' ', stamp
     &,'   between latitudes:',jscomp,' and ', jecomp
     &,'   and from the surface down to depth level ',km
        endif
!-----------------------------------------------------------------------
!       flush the buffer
!-----------------------------------------------------------------------

        call mpp_flush(io)       

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








