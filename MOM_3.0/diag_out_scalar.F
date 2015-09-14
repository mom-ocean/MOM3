      subroutine time_step_out
#if defined time_step_monitor

c-----------------------------------------------------------------------
c     save time step monitor data
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c     updated with mpp_io routines July 1999 by 
c                  s. m. griffies        e-mail=> smg@gfdl.gov  
c-----------------------------------------------------------------------
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "taskrows.h"
# include "docnam.h"
# include "diag.h"
# include "emode.h"
# include "tmngr.h"
# include "vers.h"
# include "switch.h"
      character*128 file_name, file_stamp, gvtitle
      character(len=60) :: cvarnam
# if defined salinity_psu
      character*3 saldim
# else 
      character*9 saldim
# endif
      data num_tsint /0/
      save num_tsint
      save io

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_ke
      type(fieldtype), save :: field_tracer_tot(nt)
      type(fieldtype), save :: field_tracer_bar(nt)
      type(fieldtype), save :: field_tracer_var(nt)
      type(fieldtype), save :: field_tracer_tend(nt)
# ifdef explicit_free_surface
      type(fieldtype), save :: field_volume_tot
      type(fieldtype), save :: field_volrate
      type(fieldtype), save :: field_etabar
      type(fieldtype), save :: field_kinext
      type(fieldtype), save :: field_potext
# endif
# if !defined explicit_free_surface
      type(fieldtype), save :: field_scans
# endif
      type(fieldtype), save :: field_period

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where energy analysis is written
c-----------------------------------------------------------------------

      file_name = file_stamp ('ts_intgrls', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, and salinity units 
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Time Step Monitor' // momver
# if defined salinity_psu
      saldim = 'psu' 
# else 
      saldim = 'ppm-0.035' 
# endif

c-----------------------------------------------------------------------
c     open diagnostic file and write meta information 
c-----------------------------------------------------------------------
  
      if (num_tsint .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_tsint = num_tsint + 1

c-----------------------------------------------------------------------
c       define time axis
c-----------------------------------------------------------------------

        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       define fields
c-----------------------------------------------------------------------

        call mpp_write_meta(
     &      io, field_ke, (/axis_t/), 
     &      'ke', 'g/(cm s**2)', 'kinetic energy per unit volume',
     &      -1e10, +1e10, -1.0E+34, +1.0E+34)

        do n=1,nt 

          if (n .eq. 1) then
            call mpp_write_meta(
     &      io, field_tracer_tot(n), (/axis_t/), 
     &      'heat', '1e26 Joule', 'heat',
     &      -1e10, +1e10, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tracer_bar(n), (/axis_t/), 
     &      'tbar', 'deg K', 'global average temp',
     &      -1e10, +1e10, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tracer_var(n), (/axis_t/), 
     &      'tvar', '(deg C)**2', 'variance of temp',
     &      -1e10, +1e10, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tracer_tend(n), (/axis_t/), 
     &      'tend_t', 'deg C/s', '|dT/dt|',
     &      -1e10, +1e10, -1.0E+34, +1.0E+34)

          elseif(n .eq. 2) then

            call mpp_write_meta(
     &      io, field_tracer_tot(n), (/axis_t/), 
     &      'salt', '1e26 g', 'salt deviation from 35',
     &      -1e10, +1e10, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tracer_bar(n), (/axis_t/), 
     &      'sbar', saldim, 'global average salinity - 35',
     &      -1e10, +1e10, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tracer_var(n), (/axis_t/), 
     &      'svar', saldim//'**2', 'variance of salinity',
     &      -1e10, +1e10, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tracer_tend(n), (/axis_t/), 
     &      'tend_s', saldim//'/s', '|dS/dt|',
     &      -1e10, +1e10, -1.0E+34, +1.0E+34)

          else

            write( cvarnam, '(a,i3.3)' ) 'tot_', n
            call mpp_write_meta(
     &      io, field_tracer_tot(n), (/axis_t/), 
     &      cvarnam, '1e26 g', 'total tracer mass + ocn_vol * rho0',
     &      -1e20, +1e20, -1.0E+34, +1.0E+34)

            write( cvarnam, '(a,i3.3)' ) 'avg_', n
            call mpp_write_meta(
     &      io, field_tracer_bar(n), (/axis_t/), 
     &      cvarnam, 'tracer unit', 'avg tracer concentration + 1',
     &      -1e10, +1e10, -1.0E+34, +1.0E+34)

            write( cvarnam, '(a,i3.3)' ) 'var_', n
            call mpp_write_meta(
     &      io, field_tracer_var(n), (/axis_t/), 
     &      cvarnam, 'tracer unit**2', 'variance of tracer',
     &      -1e10, +1e10, -1.0E+34, +1.0E+34)

            write( cvarnam, '(a,i3.3)' ) 'tend_', n
            call mpp_write_meta(
     &      io, field_tracer_tend(n), (/axis_t/), 
     &      cvarnam, 'tracer unit/s', '|Dtracer/Dt|',
     &      -1e10, +1e10, -1.0E+34, +1.0E+34)

          endif
        enddo

# ifdef explicit_free_surface
        call mpp_write_meta(
     &      io, field_volume_tot, (/axis_t/), 
     &      'vol_tot', 'cm**3', 'accumulated volume of sff',
     &      -1e20, +1e20, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_volrate, (/axis_t/), 
     &      'volrate', 'Sv', 'global transport from sff',
     &      -1e20, +1e20, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_etabar, (/axis_t/), 
     &      'etabar', 'cm', 'volume in free surface per unit area',
     &      -1e20, +1e20, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_kinext, (/axis_t/), 
     &      'kinext','g/(cm s**2)','ext kinetic energy per volume',
     &      -1e20, +1e20, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_potext, (/axis_t/), 
     &      'potext','g/(cm s**2)','ext potential energy per volume',
     &      -1e20, +1e20, -1.0E+34, +1.0E+34)
# endif
# if !defined explicit_free_surface
        call mpp_write_meta(
     &      io, field_scans, (/axis_t/), 
     &      'scans','dimensionless','elliptic solver iterations',
     &      0e0, +1e20, -1.0E+34, +1.0E+34)
# endif
c
c       averaging period 
c
        call mpp_write_meta(
     &      io, field_period, (/axis_t/), 
     &      'period', 'year', 'averaging period',
     &      +0.0E+00, +1.0E+20, -1.0E+34, +1.0E+34)
c
c       figure title 
c
        call mpp_write_meta(io, 'title', cval=gvtitle )

c-----------------------------------------------------------------------
c     endif for opening of diagnostic file
c-----------------------------------------------------------------------
      endif

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

      call mpp_write(io, field_ke, ektot(0,1), timrec)

      do n=1,nt 
        call mpp_write(io, field_tracer_tot(n), tratot(n),timrec)
        call mpp_write(io, field_tracer_bar(n), tbar(0,n,1), timrec)
        call mpp_write(io, field_tracer_var(n), travar(0,n,1), timrec)
        call mpp_write(io, field_tracer_tend(n),dtabs(0,n,1),timrec)
      enddo
# ifdef explicit_free_surface
      call mpp_write(io, field_volume_tot, volume_tot, timrec)
      call mpp_write(io, field_volrate, etaana(1), timrec)
      call mpp_write(io, field_etabar, etaana(2), timrec)
      call mpp_write(io, field_kinext, etaana(3), timrec)
      call mpp_write(io, field_potext, etaana(4), timrec)
# endif
# if !defined explicit_free_surface
      call mpp_write(io, field_scans, float(mscan), timrec)
# endif
      call mpp_write(io, field_period, period, timrec)

      write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8)')
     &'  ==>  Time step integrals written to file '
     &,  file_name,' on ts = ', itt,' ', stamp
     &,' between latitudes:',jscomp,' and ', jecomp

      return
      end
#else
      return
      end
#endif


#if defined energy_analysis
      subroutine energy_out
c
c-----------------------------------------------------------------------
c     save energy analysis data
c
c     author:  r. c. pacanowski   e-mail=> rcp@gfdl.gov
c
c     updated with mpp_io routines July 1999 by 
c              s. m. griffies     e-mail=> smg@gfdl.gov             
c-----------------------------------------------------------------------
c
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "taskrows.h"
# include "diag.h"
# include "tmngr.h"
# include "switch.h"
# include "vers.h"
      character*128 file_name, file_stamp, gvtitle

      data num_energy /0/
      save num_energy
      save io
c
c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
c     
      type(fieldtype), save :: field_ddt_int
      type(fieldtype), save :: field_hadv_int
      type(fieldtype), save :: field_vadv_int
      type(fieldtype), save :: field_hfric_int
      type(fieldtype), save :: field_vfric_int
      type(fieldtype), save :: field_press_int
      type(fieldtype), save :: field_imbal_int
      type(fieldtype), save :: field_wind_int
      type(fieldtype), save :: field_drag_int
c
      type(fieldtype), save :: field_ddt_ext
      type(fieldtype), save :: field_hadv_ext
      type(fieldtype), save :: field_vadv_ext
      type(fieldtype), save :: field_hfric_ext
      type(fieldtype), save :: field_vfric_ext
      type(fieldtype), save :: field_press_ext
      type(fieldtype), save :: field_imbal_ext
      type(fieldtype), save :: field_wind_ext
      type(fieldtype), save :: field_drag_ext
c
      type(fieldtype), save :: field_buoy
      type(fieldtype), save :: field_buoy_err
      type(fieldtype), save :: field_energy_leak
c
      type(fieldtype), save :: field_period
c
c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------
c
      type (axistype), save :: axis_t
c
c-----------------------------------------------------------------------
c     specify file where energy analysis is written
c-----------------------------------------------------------------------
c
      file_name = file_stamp ('energy_int', xstamp, '.dta')
c
c-----------------------------------------------------------------------
c     set averaging period, time, and title
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Energy Analysis' // momver

c-----------------------------------------------------------------------
c     open diagnostic file and write meta information 
c-----------------------------------------------------------------------
  
      if (num_energy .eq. 0) then
c-----------------------------------------------------------------------
c       open io unit
c
c       Note: 
c       MPP_NETCDF saves 32 bit netCDF snapshots 
c       MPP_IEEE32 saves 32 bit IEEE snapshots 
c       MPP_NATIVE saves snapshots data in same form as computed.
c-----------------------------------------------------------------------
c
        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_energy = num_energy + 1

c-----------------------------------------------------------------------
c       define time axis
c-----------------------------------------------------------------------

        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       define fields
c-----------------------------------------------------------------------
c
c       internal mode terms
c 
        call mpp_write_meta(
     &      io, field_ddt_int, (/axis_t/), 
     &      'ddt_int', 'g/cm/s**3', 'int mode energy change',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_hadv_int, (/axis_t/), 
     &      'hadv_int', 'g/cm/s**3','int mode work by horz adv',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vadv_int, (/axis_t/), 
     &      'vadv_int', 'g/cm/s**3','int mode work by vert adv',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_hfric_int, (/axis_t/), 
     &      'hfric_int', 'g/cm/s**3','int mode work by horz frict',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vfric_int, (/axis_t/), 
     &      'vfric_int', 'g/cm/s**3','int mode work by vert frict',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_press_int, (/axis_t/), 
     &      'press_int', 'g/cm/s**3','int mode work by pressure',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_imbal_int, (/axis_t/), 
     &      'imbal_int', 'g/cm/s**3','int mode imbalance',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_wind_int, (/axis_t/), 
     &      'wind_int', 'g/cm/s**3','int mode work by wind',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_drag_int, (/axis_t/), 
     &      'drag_int', 'g/cm/s**3','int mode work by bottom drag',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
c
c       external mode terms
c 
        call mpp_write_meta(
     &      io, field_ddt_ext, (/axis_t/), 
     &      'ddt_ext', 'g/cm/s**3', 'ext mode energy change',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_hadv_ext, (/axis_t/), 
     &      'hadv_ext', 'g/cm/s**3','ext mode work by horz adv',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vadv_ext, (/axis_t/), 
     &      'vadv_ext', 'g/cm/s**3','ext mode work by vert adv',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_hfric_ext, (/axis_t/), 
     &      'hfric_ext', 'g/cm/s**3','ext mode work by horz frict',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vfric_ext, (/axis_t/), 
     &      'vfric_ext', 'g/cm/s**3','ext mode work by vert frict',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_press_ext, (/axis_t/), 
     &      'press_ext', 'g/cm/s**3','ext mode work by press',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_imbal_ext, (/axis_t/), 
     &      'imbal_ext', 'g/cm/s**3','ext mode imbalance',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_wind_ext, (/axis_t/), 
     &      'wind_ext', 'g/cm/s**3','ext mode work by wind',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_drag_ext, (/axis_t/), 
     &      'drag_ext', 'g/cm/s**3','ext mode work by bottom drag',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
c
c       misc terms 
c
        call mpp_write_meta(
     &      io, field_buoy, (/axis_t/), 
     &      'buoy', 'g/cm/s**3','buoyancy',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_buoy_err, (/axis_t/), 
     &      'buoy_err', 'g/cm/s**3','buoyancy error',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_energy_leak, (/axis_t/), 
     &      'energy_leak', 'g/cm/s**3','non-linear error',
     &      -3000.0, +3000.0, -1.0E+34, +1.0E+34)
c
c       averaging period 
c
        call mpp_write_meta(
     &      io, field_period, (/axis_t/), 
     &      'period', 'year', 'averaging period',
     &      +0.0E+00, +1.0E+20, -1.0E+34, +1.0E+34)
c
c       figure title 
c
        call mpp_write_meta(io, 'title', cval=gvtitle )

c-----------------------------------------------------------------------
c     endif for opening of diagnostic file(s)
c-----------------------------------------------------------------------
      endif

      call mpp_write(io, field_ddt_int,  engint(0,1,1), timrec)
      call mpp_write(io, field_hadv_int, engint(0,2,1), timrec)
      call mpp_write(io, field_vadv_int, engint(0,3,1), timrec)
      call mpp_write(io, field_hfric_int,engint(0,4,1), timrec)
      call mpp_write(io, field_vfric_int,engint(0,5,1), timrec)
      call mpp_write(io, field_press_int,engint(0,6,1), timrec)
      call mpp_write(io, field_imbal_int,plicin       , timrec)
      call mpp_write(io, field_wind_int, engint(0,7,1), timrec)
      call mpp_write(io, field_drag_int, engint(0,8,1), timrec)
c
      call mpp_write(io, field_ddt_ext,  engext(1,1)  , timrec)
      call mpp_write(io, field_hadv_ext, engext(2,1)  , timrec)
      call mpp_write(io, field_vadv_ext, engext(3,1)  , timrec)
      call mpp_write(io, field_hfric_ext,engext(4,1)  , timrec)
      call mpp_write(io, field_vfric_ext,engext(5,1)  , timrec)
      call mpp_write(io, field_press_ext,engext(6,1)  , timrec)
      call mpp_write(io, field_imbal_ext,plicex       , timrec)
      call mpp_write(io, field_wind_ext, engext(7,1)  , timrec)
      call mpp_write(io, field_drag_ext, engext(8,1)  , timrec)

      call mpp_write(io, field_buoy,     buoy(0,1)    , timrec)
      call mpp_write(io, field_buoy_err, buoerr       , timrec)
      call mpp_write(io, field_energy_leak, enleak    , timrec)

      call mpp_write(io, field_period, period         , timrec)

      write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8)')
     &'  ==>  Global energy analysis written to file '
     &,  file_name,' on ts = ', itt,' ', stamp
     &,' between latitudes:',jscomp,' and ', jecomp

      return
      end
#endif


#if defined tracer_averages
      subroutine traceravg_out

c-----------------------------------------------------------------------
c     save tracer averages over specified domain 
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c     updated with mpp_io routines July 1999 by 
c                  s. m. griffies        e-mail=> smg@gfdl.gov  
c-----------------------------------------------------------------------
      use arrays_2d, only : mskhr
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "taskrows.h"
# include "coord.h"
# include "diag.h"
# include "diag_region.h"
# include "tmngr.h"
# include "switch.h"
# include "vers.h"
      character*128 file_name, file_stamp, gvtitle
      character*9 basin_num
      character(len=50) :: cvarnam
      character(len=50) :: cvarlnam
# if defined salinity_psu
      character*3 saldim
# else 
      character*9 saldim
# endif
      dimension bufbk(km)
      dimension bufmskhr(2:imt-1,jscomp:jecomp)
      data num_travg /0/
      save num_travg
      save io

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_stf(nt)
      type(fieldtype), save :: field_stf_basin(nhreg,nt)
      type(fieldtype), save :: field_tracer_basin(nhreg,nt)
      type(fieldtype), save :: field_tracer_level(nt)
      type(fieldtype), save :: field_tracer_level_basin(nhreg,nt)
      type(fieldtype), save :: field_mskhr
      type(fieldtype), save :: field_period

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_xt
      type (axistype), save :: axis_yt
      type (axistype), save :: axis_zt
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where energy analysis is written
c-----------------------------------------------------------------------

      file_name = file_stamp ('tracer_avg', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, and salinity units  
c-----------------------------------------------------------------------

      period  = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Tracer Averages' // momver
# if defined salinity_psu
      saldim = 'psu' 
# else 
      saldim = 'ppm-0.035' 
# endif
c-----------------------------------------------------------------------
c     open diagnostic file and write meta information
c-----------------------------------------------------------------------
  
      if (num_travg .eq. 0) then
c-----------------------------------------------------------------------
c       open io unit
c
c       Note: 
c       MPP_NETCDF saves 32 bit netCDF snapshots 
c       MPP_IEEE32 saves 32 bit IEEE snapshots 
c       MPP_NATIVE saves snapshots data in same form as computed
c-----------------------------------------------------------------------

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_travg = num_travg + 1

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
     &         io, axis_zt, 'zt_k', 'cm', 'Depth of T grid point',
     &         cartesian='Z', sense=-1, data=zt)
        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c      define fields (on data domain)
c-----------------------------------------------------------------------
c
c      global averaged surface tracer flux 
c
       do n=1,nt 

          if (n .eq. 1) then
            call mpp_write_meta(
     &        io, field_stf(n), (/axis_t/),
     &        'stf_1', 'cal/cm**2/s', 'global surface heat flux',  
     &        -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          elseif (n .eq. 2) then
            call mpp_write_meta(
     &        io, field_stf(n), (/axis_t/),
     &        'stf_2', 'g/cm**2/s', 'global surface salt flux',  
     &        -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          else
            write( cvarnam, '(a,i3.3)' ) 'stf_', n
            call mpp_write_meta(
     &        io, field_stf(n), (/axis_t/),
     &        cvarnam, 'tracer/cm**2/s', 'global surface tracer flux',  
     &        -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          endif

        enddo
c
c       surface tracer flux (basin)
c
        do n=1,nt
	  do m=1,nhreg
            write (basin_num,'(a,i3.3)') '_basin',m
            if (n .eq. 1) then
              call mpp_write_meta(
     &        io, field_stf_basin(m,n), (/axis_t/)
     &        ,'stf_1'//basin_num, 'cal/cm**2/s'
     &        ,'surface heat flux'//basin_num  
     &        ,-1.e6, 1.e6, -1.0E+34, +1.0E+34)
            elseif (n .eq. 2) then
              call mpp_write_meta(
     &        io, field_stf_basin(m,n), (/axis_t/)
     &       ,'stf_2'//basin_num, 'g/cm**2/s'
     &       ,'surface salt flux'//basin_num  
     &       ,-1.e6, 1.e6, -1.0E+34, +1.0E+34)
            else
              write (cvarnam,'(a,i3.3,a9)') 'stf_',n, basin_num
              write (cvarlnam,'(a,i3.3,a9)')
     &              'surface tracer flux',n, basin_num
              call mpp_write_meta(
     &        io, field_stf_basin(m,n), (/axis_t/)
     &       ,cvarnam, 'tracer/cm**2/s'
     &       ,cvarlnam  
     &       ,-1.e6, 1.e6, -1.0E+34, +1.0E+34)
	    endif

          enddo
        enddo
c
c       average tracer (basin)
c
        do n=1,nt
	  do m=1,nhreg
            write (basin_num,'(a,i3.3)') '_basin',m
            if (n .eq. 1) then
              call mpp_write_meta(
     &        io, field_tracer_basin(m,n), (/axis_t/)
     &        ,'temp'//basin_num, 'deg C'
     &        ,'temp'//basin_num  
     &        ,-1.e6, 1.e6, -1.0E+34, +1.0E+34)
            elseif (n .eq. 2) then
             call mpp_write_meta(
     &        io, field_tracer_basin(m,n), (/axis_t/)
     &        ,'salinity'//basin_num, saldim  
     &        ,'salinity'//basin_num  
     &        ,-1.e6, 1.e6, -1.0E+34, +1.0E+34)
            else
              write (cvarnam,'(a,i3.3,a)') 'P_',n, basin_num
              write (cvarlnam,'(a,i3.3,a9)')
     &              'passive tracer ',n, basin_num
              call mpp_write_meta(
     &        io, field_tracer_basin(m,n), (/axis_t/)
     &       ,cvarnam, 'basin tracer'
     &       ,cvarlnam  
     &       ,-1.e6, 1.e6, -1.0E+34, +1.0E+34)
	    endif

          enddo
        enddo
c
c       global average tracer (level)
c
        do n=1,nt 

          if (n .eq. 1) then
            call mpp_write_meta(
     &        io, field_tracer_level(n), (/axis_zt,axis_t/)
     &       ,'T_k', 'cal/cm**2/s', 'global potential temp(k)'  
     &       ,-1.e6, 1.e6, -1.0E+34, +1.0E+34)
          elseif (n .eq. 2) then
            call mpp_write_meta(
     &        io, field_tracer_level(n), (/axis_zt,axis_t/)
     &       ,'S_k', saldim, 'global salinity(k)'  
     &       ,-1.e6, 1.e6, -1.0E+34, +1.0E+34)
	  else
            write (cvarnam,'(a,i3.3)') 'P_k_',n
            write (cvarlnam,'(a,i3.3)')
     &      'global passive tracer(k)_',n
            call mpp_write_meta(
     &        io, field_tracer_level(n), (/axis_zt,axis_t/)
     &        ,cvarnam, 'global tracer(k)',cvarlnam    
     &        ,-1.e6, 1.e6, -1.0E+34, +1.0E+34)
          endif

        enddo
c
c       basin average tracer (level)
c
        do n=1,nt
	  do m=1,nhreg
            write (basin_num,'(a,i3.3)') '_basin',m
            if (n .eq. 1) then
              call mpp_write_meta(
     &        io, field_tracer_level_basin(m,n), (/axis_zt,axis_t/)
     &       ,'T_k'//basin_num,'deg C'
     &       ,'temp(k)'//basin_num 
     &       ,-1.e6, 1.e6, -1.0E+34, +1.0E+34)

            elseif (n .eq. 2) then
              call mpp_write_meta(
     &        io, field_tracer_level_basin(m,n), (/axis_zt,axis_t/)
     &       ,'S_k'//basin_num, saldim, 'salinity(k)'//basin_num
     &       ,-1.e6, 1.e6, -1.0E+34, +1.0E+34)
	    else
              write (cvarnam,'(a,i3.3,a9)') 'P_k_',n, basin_num
              write (cvarlnam,'(a,i3.3,a9)')
     &        'tracer(k)_',n,basin_num
              call mpp_write_meta(
     &        io, field_tracer_level_basin(m,n), (/axis_zt,axis_t/)
     &       ,cvarnam,'tracer units',cvarlnam             
     &       ,-1.e6, 1.e6, -1.0E+34, +1.0E+34)
	    endif

          enddo
        enddo
c
c       averaging period 
c
        call mpp_write_meta(
     &      io, field_period, (/axis_t/), 
     &      'period', 'years', 'averaging period',
     &      +0.0E+00, +1.0E+20, -1.0E+34, +1.0E+34)
c
c       basin masking 
c
        call mpp_write_meta(
     &      io, field_mskhr, (/axis_xt,axis_yt,axis_t/), 
     &      'mskhr', 'basin numbers', 'numbers',
     &      0.0, 100.0, -1.0E+34, +1.0E+34)
c
c       figure title 
c
        call mpp_write_meta(io, 'title', cval=gvtitle )

c-----------------------------------------------------------------------
c     write axis after all meta data has been written
c-----------------------------------------------------------------------

        call mpp_write (io,axis_xt)
        call mpp_write (io,axis_yt)
        call mpp_write (io,axis_zt)
      
c-----------------------------------------------------------------------
c     endif for opening of diagnostic file
c-----------------------------------------------------------------------
      endif

c-----------------------------------------------------------------------
c     write the output 
c-----------------------------------------------------------------------
c
c     global average tracer flux  
c
      do n=1,nt
        call mpp_write(io, field_stf(n), avggf(n), timrec)
      enddo
c
c     basin average tracer flux  
c
      do n=1,nt
	do m=1,nhreg
          call mpp_write(io, field_stf_basin(m,n), avgbf(m,n), timrec)
        enddo
      enddo
c
c     basin average tracer 
c
      do n=1,nt
	do m=1,nhreg
          call mpp_write(io, field_tracer_basin(m,n), avgbt(m,n),timrec)
        enddo
      enddo
c
c     global average tracer on k-levels
c
      do n=1,nt
        call mpp_write(io, field_tracer_level(n), avggk(:,n), timrec)
      enddo
c
c     basin average tracer on k-levels
c
      do n=1,nt
	do m=1,nhreg
	  do k=1,km
	    bufbk(k) = avgbk(m,k,n)
	  enddo
          call mpp_write(io, field_tracer_level_basin(m,n)
     &,                  bufbk, timrec)
        enddo
      enddo
c
c     averaging period 
c
      call mpp_write(io, field_period, period, timrec)
c
c     to save space, only write basin mask if "itavg" is true
c 
      if (itavg .and. num_travg .eq. 1) then
        do j=jscomp,jecomp
          do i=2,imt-1
            bufmskhr(i,j) = float(mskhr(i,j))
            if(bufmskhr(i,j) .eq. 0.0) then
              bufmskhr(i,j) = field_mskhr%missing
            endif
	  enddo
	enddo
        call mpp_write(io, field_mskhr, bufmskhr, timrec)
      endif

      write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8)')
     &'  ==>  Tracer averages written to file '
     &,  file_name,' on ts = ', itt,' ', stamp
     &,' between latitudes:',jscomp,' and ', jecomp

      return
      end
#endif



#ifdef term_balances
      subroutine termbal_out
c
c-----------------------------------------------------------------------
c     save term balance data
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c     updated with mpp_io routines July 1999 by 
c                  s. m. griffies        e-mail=> smg@gfdl.gov  
c-----------------------------------------------------------------------
c
      use arrays_2d, only : mskhr
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "taskrows.h"
# include "coord.h"
# include "diag.h"
# include "diag_region.h"
# include "switch.h"
# include "tmngr.h"
# include "vers.h"
      dimension bufreg(numreg)
      dimension bufregp1(0:numreg)
      dimension bufhorz(nhreg)
      dimension bufvert(nvreg)
      dimension bufmskvr(km)
      dimension bufmskhr(2:imt-1,jscomp:jecomp)
      dimension bufzt(km)
      dimension regnum(nhreg,nvreg)
      character*128 file_name, file_stamp, gvtitle
      character*10 tname
# if defined salinity_psu
      character*3 saldim
# else 
      character*9 saldim
# endif
      data num_termbal /0/
      save num_termbal
      save io

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_regnum
      type(fieldtype), save :: field_mskhr
      type(fieldtype), save :: field_mskvr

      type(fieldtype), save :: field_volt
      type(fieldtype), save :: field_volu
      type(fieldtype), save :: field_areat
      type(fieldtype), save :: field_areau

      type(fieldtype), save :: field_smf(2)
      type(fieldtype), save :: field_stf(nt)
      type(fieldtype), save :: field_sstracer(nt)

      type(fieldtype), save :: field_uterm(17)
      type(fieldtype), save :: field_vterm(17)
      type(fieldtype), save :: field_w
      type(fieldtype), save :: field_tterm(15,nt)

      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_numreg
      type (axistype), save :: axis_numregp1
      type (axistype), save :: axis_nhreg
      type (axistype), save :: axis_nvreg

      type (axistype), save :: axis_xt
      type (axistype), save :: axis_yt
      type (axistype), save :: axis_zt
      type (axistype), save :: axis_surf
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where term balances are written
c-----------------------------------------------------------------------

      file_name = file_stamp ('term_bal', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, salinity units, and depth 
c-----------------------------------------------------------------------

      period  = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Term Balances' // momver
# if defined salinity_psu
      saldim = 'psu' 
# else 
      saldim = 'ppm-0.035' 
# endif
      bufzt(:) = zt(:)*0.01

c-----------------------------------------------------------------------
c     open diagnostic file and write meta information 
c-----------------------------------------------------------------------

      if (num_termbal .eq. 0) then

c-----------------------------------------------------------------------
c       open io unit
c
c       Note: 
c       MPP_NETCDF saves 32 bit netCDF snapshots 
c       MPP_IEEE32 saves 32 bit IEEE snapshots 
c       MPP_NATIVE saves snapshots data in same form as computed
c-----------------------------------------------------------------------
        
        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_termbal = num_termbal + 1

c-----------------------------------------------------------------------
c       axes (on global domain)
c-----------------------------------------------------------------------

        ibeg = mom_domain%x%global%start_index
        iend = mom_domain%x%global%end_index
        jbeg = mom_domain%y%global%start_index
        jend = mom_domain%y%global%end_index
c
c       axes which specify the region number 
c
        do n=1,numreg
          bufreg(n)   = float(n)
        enddo
        do n=0,numreg
          bufregp1(n) = float(n)
        enddo
        do n=1,nhreg
          bufhorz(n)  = float(n)
        enddo
        do n=1,nvreg
          bufvert(n)  = float(n)
        enddo        

        call mpp_write_meta(
     &       io, axis_numreg, 'numreg', 'N', 'region #',data=bufreg)
        call mpp_write_meta(
     &       io, axis_numregp1,'numregp1','N', 'region #',data=bufregp1)
        call mpp_write_meta(
     &       io, axis_nhreg, 'nhreg', 'N', 'region #',data=bufhorz)
        call mpp_write_meta(
     &       io, axis_nvreg, 'nvreg', 'N', 'region #',data=bufvert)
c
c       axes which specify space-time 
c
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
     &         io, axis_surf, 'surface', 'm', 'Depth of surface',
     &         cartesian='Z', sense=-1, data=(/0./))
        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')
c
c       horizontal and vertical region distribution
c
        call mpp_write_meta(
     &      io, field_mskhr, (/axis_xt,axis_yt,axis_t/), 
     &      'mskhr', 'numbers', 'horizontal regions', 
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_mskvr, (/axis_zt,axis_t/), 
     &      'mskvr', 'numbers', 'vertical regions',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_regnum, (/axis_nhreg,axis_nvreg,axis_t/), 
     &      'regnum', 'numbers', 'region numbers',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c
c       volumes and areas  
c
        call mpp_write_meta(
     &      io, field_volt, (/axis_numregp1,axis_t/), 
     &      'volt', 'cm**3', 'regional T volume',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_volu, (/axis_numregp1,axis_t/), 
     &      'volu', 'cm**3', 'regional U volume',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_areat, (/axis_numregp1,axis_t/), 
     &      'areat', 'cm**3', 'regional T surface area',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_areau, (/axis_numregp1,axis_t/), 
     &      'areau', 'cm**3', 'regional U surface area',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c
c       surface wind stress 
c
        call mpp_write_meta(
     &      io, field_smf(1), (/axis_nhreg,axis_t/), 
     &      'taux', 'dynes/cm**2', 'zonal wind stress',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_smf(2), (/axis_nhreg,axis_t/), 
     &      'tauy', 'dynes/cm**2', 'merid wind stress',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c
c       surface tracer fluxes and values
c
        do n=1,nt  
          if(n .eq. 1) then 
            call mpp_write_meta(
     &      io, field_stf(n), (/axis_nhreg,axis_t/), 
     &      'hflx', 'cal/cm**2/s', 'heat flux',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          elseif(n .eq. 2) then
            call mpp_write_meta(
     &      io, field_stf(n), (/axis_nhreg,axis_t/), 
     &      'sflx', 'gram/cm**2/s', 'salt flux',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          else
  	    write (tname,'(a,i3.3,a)') 'stf_',n
            call mpp_write_meta(
     &      io, field_stf(n), (/axis_nhreg,axis_t/), 
     &      tname, 'unknown', 'tracer flux',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          endif
        enddo
        do n=1,nt 
          if(n .eq. 1) then 
            call mpp_write_meta(
     &      io, field_sstracer(n), (/axis_nhreg,axis_t/), 
     &      'sst', 'deg C', 'sea surface temp',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          elseif(n .eq. 2) then
            call mpp_write_meta(
     &      io, field_sstracer(n), (/axis_nhreg,axis_t/), 
     &      'sss', saldim, 'sea surface salinity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          else
  	    write (tname,'(a,i3.3)') 'ssp_',n
            call mpp_write_meta(
     &      io, field_sstracer(n), (/axis_nhreg,axis_t/), 
     &      tname, 'unknown', 'sea surface passive tracer',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          endif
        enddo
c
c       zonal velocity terms 
c
        call mpp_write_meta(
     &      io, field_uterm(1), (/axis_numreg,axis_t/), 
     &      'dudt', 'cm/s**2', 'zonal velocity tendency',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(2), (/axis_numreg,axis_t/), 
     &      'neg_px', 'cm/s**2', '-zonal press gradient',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(3), (/axis_numreg,axis_t/), 
     &      'neg_uux', 'cm/s**2', '-d(uu)/dx',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(4), (/axis_numreg,axis_t/), 
     &      'neg_vuy', 'cm/s**2', '-d(vu)/dy',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(5), (/axis_numreg,axis_t/), 
     &      'neg_wuz', 'cm/s**2', '-d(wu)/dz',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(6), (/axis_numreg,axis_t/), 
     &      'zon_diff_u', 'cm/s**2', 'zonal diff of u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(7), (/axis_numreg,axis_t/), 
     &      'mer_diff_u', 'cm/s**2', 'merid diff of u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(8), (/axis_numreg,axis_t/), 
     &      'ver_diff_u', 'cm/s**2', 'vert diff of u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(9), (/axis_numreg,axis_t/), 
     &      'metric_diff_u', 'cm/s**2', 'metric diff of u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(10), (/axis_numreg,axis_t/), 
     &      'fv', 'cm/s**2', 'u coriolis term',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(11), (/axis_numreg,axis_t/), 
     &      'source_u', 'cm/s**2', 'u source term',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(12), (/axis_numreg,axis_t/), 
     &      'neg_surfpx', 'cm/s**2', '-zonal surf press grad',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(13), (/axis_numreg,axis_t/), 
     &      'adv_metric_u', 'cm/s**2', 'advect metric for u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(14), (/axis_numreg,axis_t/), 
     &      'neg_zonadv_u', 'cm/s**2', '-zonal advect of u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(15), (/axis_numreg,axis_t/), 
     &      'neg_meradv_u', 'cm/s**2', '-merid advect of u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(16), (/axis_numreg,axis_t/), 
     &      'neg_veradv_u', 'cm/s**2', '-vert advect of u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uterm(17), (/axis_numreg,axis_t/), 
     &      'zon_vel', 'cm/s', 'zonal velocity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c
c       meridional velocity terms 
c
        call mpp_write_meta(
     &      io, field_vterm(1), (/axis_numreg,axis_t/), 
     &      'dvdt', 'cm/s**2', 'merid velocity tendency',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(2), (/axis_numreg,axis_t/), 
     &      'neg_py', 'cm/s**2', '-merid press gradient',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(3), (/axis_numreg,axis_t/), 
     &      'neg_uvx', 'cm/s**2', '-d(uv)/dx',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(4), (/axis_numreg,axis_t/), 
     &      'neg_vvy', 'cm/s**2', '-d(vv)/dy',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(5), (/axis_numreg,axis_t/), 
     &      'neg_wvz', 'cm/s**2', '-d(wv)/dz',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(6), (/axis_numreg,axis_t/), 
     &      'zon_diff_v', 'cm/s**2', 'zonal diff of v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(7), (/axis_numreg,axis_t/), 
     &      'mer_diff_v', 'cm/s**2', 'merid diff of v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(8), (/axis_numreg,axis_t/), 
     &      'ver_diff_v', 'cm/s**2', 'vert diff of v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(9), (/axis_numreg,axis_t/), 
     &      'metric_diff_v', 'cm/s**2', 'metric diff of v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(10), (/axis_numreg,axis_t/), 
     &      'neg_fu', 'cm/s**2', 'v coriolis term',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(11), (/axis_numreg,axis_t/), 
     &      'source_v', 'cm/s**2', 'v source term',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(12), (/axis_numreg,axis_t/), 
     &      'neg_surfpy', 'cm/s**2', '-merd surface press grad',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(13), (/axis_numreg,axis_t/), 
     &      'adv_metric_v', 'cm/s**2', 'advect metric for v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(14), (/axis_numreg,axis_t/), 
     &      'neg_zonadv_v', 'cm/s**2', '-zonal advect of v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(15), (/axis_numreg,axis_t/), 
     &      'neg_meradv_v', 'cm/s**2', '-merid advect of v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(16), (/axis_numreg,axis_t/), 
     &      'neg_veradv_v', 'cm/s**2', '-vert advect of v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vterm(17), (/axis_numreg,axis_t/), 
     &      'mer_vel', 'cm/s', 'merd velocity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c
c       vertical velocity 
c
        call mpp_write_meta(
     &      io, field_w, (/axis_numreg,axis_t/), 
     &      'ver_vel', 'cm/s', 'vertical velocity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c
c       tracer terms 
c
        do n=1,nt

          if(n .eq. 1) then 
            call mpp_write_meta(
     &      io, field_tterm(1,n), (/axis_numreg,axis_t/), 
     &      'dTdt', 'deg C/s', 'd(temp)/dt',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(2,n), (/axis_numreg,axis_t/), 
     &      'neg_utx', 'deg C/s', '-d(u*temp)/dx',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(3,n), (/axis_numreg,axis_t/), 
     &      'neg_vty', 'deg C/s', '-d(v*temp)/dy',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(4,n), (/axis_numreg,axis_t/), 
     &      'neg_wtz', 'deg C/s', '-d(w*temp)/dz',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(5,n), (/axis_numreg,axis_t/), 
     &      'diff_tx', 'deg C/s', 'zonal diff of temp',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(6,n), (/axis_numreg,axis_t/), 
     &      'diff_ty', 'deg C/s', 'merid diff of temp',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(7,n), (/axis_numreg,axis_t/), 
     &      'diff_tz', 'deg C/s', 'vert diff of temp',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(8,n), (/axis_numreg,axis_t/), 
     &      'source_T', 'deg C/s', 'temp change by sources',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(9,n), (/axis_numreg,axis_t/), 
     &      'convect_T', 'deg C/s', 'temp change by convect',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(10,n), (/axis_numreg,axis_t/), 
     &      'filter_T', 'deg C/s', 'temp change by filter',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(11,n), (/axis_numreg,axis_t/), 
     &      'neg_zonadv_T', 'deg C/s', '-zonal temp advect',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(12,n), (/axis_numreg,axis_t/), 
     &      'neg_meradv_T', 'deg C/s', '-merid temp advect',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(13,n), (/axis_numreg,axis_t/), 
     &      'neg_veradv_T', 'deg C/s', '-vert temp advect',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(14,n), (/axis_numreg,axis_t/), 
     &      'var_T', '(deg C)**2', 'variance of temp',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(15,n), (/axis_numreg,axis_t/), 
     &      'temp', 'deg C', 'potential temp',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

          elseif(n .eq. 2) then 
            call mpp_write_meta(
     &      io, field_tterm(1,n), (/axis_numreg,axis_t/), 
     &      'dSdt', saldim//'/s', 'd(salinity)/dt',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(2,n), (/axis_numreg,axis_t/), 
     &      'neg_usx', saldim//'/s', '-d(u*salinity)/dx',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(3,n), (/axis_numreg,axis_t/), 
     &      'neg_vsy',saldim//'/s', '-d(v*salinity)/dy',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(4,n), (/axis_numreg,axis_t/), 
     &      'neg_wsz',saldim//'/s', '-d(w*salinity)/dz',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(5,n), (/axis_numreg,axis_t/), 
     &      'diff_sx',saldim//'/s', 'zonal diff of salinity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(6,n), (/axis_numreg,axis_t/), 
     &      'diff_sy',saldim//'/s', 'merid diff of salinity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(7,n), (/axis_numreg,axis_t/), 
     &      'diff_sz',saldim//'/s', 'vert diff of salinity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(8,n), (/axis_numreg,axis_t/), 
     &      'source_S',saldim//'/s', 'salinity change by source',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(9,n), (/axis_numreg,axis_t/), 
     &      'convect_S',saldim//'/s', 'salinity change by convect',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(10,n), (/axis_numreg,axis_t/), 
     &      'filter_S',saldim//'/s', 'salinity change by filter',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(11,n), (/axis_numreg,axis_t/), 
     &      'neg_zonadv_S',saldim//'/s', '-zonal salinity advect',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(12,n), (/axis_numreg,axis_t/), 
     &      'neg_meradv_S',saldim//'/s', '-merid  salinity advect',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(13,n), (/axis_numreg,axis_t/), 
     &      'neg_veradv_S',saldim//'/s', '-vert salinity advect',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(14,n), (/axis_numreg,axis_t/), 
     &      'var_S', saldim//'**2', 'variance of salinity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(15,n), (/axis_numreg,axis_t/), 
     &      'salinity', saldim, 'salinity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

          else
	  write (tname,'(a,i1)') '_',n

            call mpp_write_meta(
     &      io, field_tterm(1,n), (/axis_numreg,axis_t/), 
     &      'dTdt'//tname, 'tracer'//tname//'/s',
     &      'd(tracer)/dt'//tname,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(2,n), (/axis_numreg,axis_t/), 
     &      'neg_utx'//tname, 'tracer'//tname//'/s',
     &      '-d(u*tracer)/dx'//tname,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(3,n), (/axis_numreg,axis_t/), 
     &      'neg_vty'//tname,'tracer'//tname//'/s',
     &      '-d(v*tracer)/dy',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(4,n), (/axis_numreg,axis_t/), 
     &      'neg_wtz'//tname,'tracer'//tname//'/s',
     &      '-d(w*tracer)/dz',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(5,n), (/axis_numreg,axis_t/), 
     &      'diff_tx'//tname,'tracer'//tname//'/s',
     &      'zonal diff of tracer',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(6,n), (/axis_numreg,axis_t/), 
     &      'diff_ty'//tname,'tracer'//tname//'/s',
     &      'merid diff of tracer',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(7,n), (/axis_numreg,axis_t/), 
     &      'diff_tz'//tname,'tracer'//tname//'/s',
     &      'vert diff of tracer',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(8,n), (/axis_numreg,axis_t/), 
     &      'source_T'//tname,'tracer'//tname//'/s',
     &      'tracer change by source',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(9,n), (/axis_numreg,axis_t/), 
     &      'convect_T'//tname,'tracer'//tname//'/s',
     &      'tracer change by convect',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(10,n), (/axis_numreg,axis_t/), 
     &      'filter_T'//tname,'tracer'//tname//'/s',
     &      'tracer change by filter',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(11,n), (/axis_numreg,axis_t/), 
     &      'neg_zonadv_T'//tname,'tracer'//tname//'/s',
     &      '-zonal tracer advect',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(12,n), (/axis_numreg,axis_t/), 
     &      'neg_meradv_T'//tname,'tracer'//tname//'/s',
     &      '-merid tracer advect',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(13,n), (/axis_numreg,axis_t/), 
     &      'neg_veradv_T'//tname,'tracer'//tname//'/s',
     &      '-vert tracer advect',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(14,n), (/axis_numreg,axis_t/), 
     &      'var_T'//tname, 'tracer'//tname//'**2',
     &      'variance of tracer',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_tterm(15,n), (/axis_numreg,axis_t/), 
     &      'tracer'//tname,'passive tracer'//tname//' units',
     &      'passive tracer',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

          endif
        enddo
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
c     write axes after all meta data has been written
c-----------------------------------------------------------------------

        call mpp_write (io,axis_numreg)
        call mpp_write (io,axis_numregp1)
        call mpp_write (io,axis_nhreg)
        call mpp_write (io,axis_nvreg)
      
        call mpp_write (io,axis_xt)
        call mpp_write (io,axis_yt)
        call mpp_write (io,axis_zt)    
        call mpp_write (io,axis_surf)

c-----------------------------------------------------------------------
c     endif for opening of diagnostic file
c-----------------------------------------------------------------------
      endif

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

c     only save region distributions, areas and volumes only when
c     "itrmb" is true to save space.
c
      if (itrmb .and. num_termbal .eq. 1) then
        do j=jscomp,jecomp
	  do i=2,imt-1
	    bufmskhr(i,j) = float(mskhr(i,j))
            if(bufmskhr(i,j) .eq. 0.0) then
              bufmskhr(i,j) = field_mskhr%missing
            endif
	  enddo
	enddo
        call mpp_write(io, field_mskhr, bufmskhr, timrec)
        do k=1,km
	  bufmskvr(k) = float(mskvr(k))
	enddo
        call mpp_write(io, field_mskvr, bufmskvr, timrec)
        do j=1,nvreg
          do i=1,nhreg
	    regnum(i,j) = (j-1)*nhreg + i
	  enddo
        enddo
        call mpp_write(io, field_regnum, regnum, timrec)
        call mpp_write(io, field_volt, volt, timrec)
        call mpp_write(io, field_volu, volu, timrec)
        call mpp_write(io, field_areat, areat, timrec)
        call mpp_write(io, field_areau, areau, timrec)
      endif
c
c     surface momentum and tracer fluxes
c
      do i=1,2
        do n=1,nhreg
          bufhorz(n) = smflx(i,n)
        enddo
        call mpp_write(io, field_smf(i), bufhorz, timrec)
      enddo
      do i=1,nt
        do n=1,nhreg
          bufhorz(n) = stflx(i,n)
        enddo
        call mpp_write(io, field_stf(i), bufhorz, timrec)
      enddo
c
c     surface tracers
c
      do i=1,nt
        do n=1,nhreg
          bufhorz(n) = asst(i,n)
        enddo
        call mpp_write(io, field_sstracer(i), bufhorz, timrec)
      enddo
c
c     terms for u equation
c
      do m=1,nuterms      
        do n=1,numreg
	  bufreg(n) = termbm(0,m,1,n)
        enddo
        call mpp_write(io, field_uterm(m), bufreg, timrec)
      enddo
c
c     terms for v equation
c
      do m=1,nuterms      
        do n=1,numreg
	  bufreg(n) = termbm(0,m,2,n)
        enddo
        call mpp_write(io, field_vterm(m), bufreg, timrec)
      enddo
c
c     vertical velocity field
c
      call mpp_write(io, field_w, avgw, timrec)
c
c     terms for tracer equation
c
      do i=1,nt
        do m=1,ntterms      
          do n=1,numreg
	    bufreg(n) = termbt(0,m,i,n)
          enddo
          call mpp_write(io, field_tterm(m,i), bufreg, timrec)
        enddo
      enddo
c
c     averaging period 
c
      call mpp_write(io, field_period, period, timrec)
c
      write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8)')
     &'  ==>  Term balances written to file '
     &,  file_name,' on ts = ', itt,' ', stamp
     &,' between latitudes:',jscomp,' and ', jecomp

      return
      end
#endif
