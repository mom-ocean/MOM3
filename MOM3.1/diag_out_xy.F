      subroutine dsp_out (dspper, dsp)
#if defined diagnostic_surf_height
c
c-----------------------------------------------------------------------
c     save diagnostic surface height data
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c     updated with mpp_io routines July 1999 by 
c                  s. m. griffies        e-mail=> smg@gfdl.gov  
c-----------------------------------------------------------------------
      use arrays_2d
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "coord.h"
# include "diag.h"
# include "switch.h"
# include "tmngr.h"
# include "taskrows.h"
# include "vers.h"
      character*128 file_name, file_stamp, gvtitle
      dimension dsp(imt,jmt)
      data num_dsp /0/
      save num_dsp
      save io

c-----------------------------------------------------------------------
c     domain type for specifying x-y model domain 
c-----------------------------------------------------------------------

      type(domain2D), pointer :: domain

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_eta
      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_xt
      type (axistype), save :: axis_yt
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where zonal mean tracer components are written
c-----------------------------------------------------------------------

      file_name = file_stamp ('diag_surf', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, and title
c-----------------------------------------------------------------------

      period = dspper
      if (rununits.eq.'days') then
        timrec = prelyr - 0.5*dspper
      else
        timrec = prelyr - 0.5*dspper/yrlen
      endif
      gvtitle = 'Diagnostic Surface Height' // momver

c-----------------------------------------------------------------------
c     check whether diagnostic file(s) should be opened
c-----------------------------------------------------------------------

      if (num_dsp .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_dsp = num_dsp + 1

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
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       attributes of data 
c-----------------------------------------------------------------------

        call mpp_write_meta(
     &      io, field_eta, (/axis_xt,axis_yt,axis_t/), 
     &      'eta', 'cm', 'Diagnosed surface height',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

        call mpp_write_meta(
     &      io, field_period, (/axis_t/), 
     &      'dspper', 'day', 'averaging period',
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

c-----------------------------------------------------------------------
c     endif for opening of diagnostic file
c-----------------------------------------------------------------------

      endif

c-----------------------------------------------------------------------
c     associate domain with x-y-domain from mom_domain 
c-----------------------------------------------------------------------

      domain  => mom_domain

c-----------------------------------------------------------------------
c     mask the field
c-----------------------------------------------------------------------

      dsp(:,:) =  dsp(:,:) + (1-min(1,kmt(:,:)))*field_eta%missing

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

      call mpp_write(io, field_eta, domain, dsp, timrec)
      call mpp_write(io, field_period, period, timrec)

       write (stdout,'(a,a,f8.2,a,a,a,a,i10,a/a,a,i3,a,i3)')
     & ' ==> Diagnostic surface height (averaged '
     &,' over ',dspper,' days) written '
     &,' to file ',file_name, ' on ts=', itt, ' ',stamp
     &,' between latitudes:',jscomp,' and ', jecomp

      return
      end
#else
      return
      end
#endif



#if defined save_ncdiff
      subroutine ncdiff_out (joff, js, je, within_window_loop)
     
c=======================================================================
c     Save two-dimensional fields developed from 
c     nonconstant diffusivity calculation.
c     
c     author:  S. M. Griffies  e-mail=> smg@gfdl.gov
c
c     updates to mpp_io made July 1999  
c=======================================================================

      use arrays_2d
      use memory_window
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod    
#  include "param.h"
#  include "coord.h"
#  include "diag.h"
#  include "iounit.h"
#  include "isopyc.h"
#  include "mw.h"
#  include "switch.h"
#  include "taskrows.h"
#  include "tmngr.h"
#  include "vers.h"
      character*128 file_name, file_stamp, gvtitle
      logical within_window_loop
      dimension bufzt(km)
      data num_ncdiff /0/
      save num_ncdiff
      save io

      real, dimension(:,:), save, allocatable :: buf_growth2a
      real, dimension(:,:), save, allocatable :: buf_agm
#  if defined vmhs_diffusivity
      real, dimension(:,:), save, allocatable :: buf_betaeff
      real, dimension(:,:), save, allocatable :: buf_hllength
      real, dimension(:,:), save, allocatable :: buf_hldiff
#  endif
#  if defined vmhs_diffusivity
      real, dimension(:,:), save, allocatable :: buf_vmhslength
      real, dimension(:,:), save, allocatable :: buf_vmhsdiff
#  endif

c-----------------------------------------------------------------------
c     domain pointer 
c-----------------------------------------------------------------------
     
      type(domain2D), pointer :: domain

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_growth2a
      type(fieldtype), save :: field_agm
# if defined vmhs_diffusivity
      type(fieldtype), save :: field_betaeff
      type(fieldtype), save :: field_hllength
      type(fieldtype), save :: field_hldiff
# endif
#  if defined vmhs_diffusivity
      type(fieldtype), save :: field_vmhslength
      type(fieldtype), save :: field_vmhsdiff
#  endif
      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_xt
      type (axistype), save :: axis_yt
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where mixing coefficients are written
c-----------------------------------------------------------------------

      file_name = file_stamp ('ncdiff', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, and ocean depth in meters
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Nonconstant diffusivity terms' // momver

c-----------------------------------------------------------------------
c     check whether diagnostic file(s) should be opened
c-----------------------------------------------------------------------

      if (num_ncdiff .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_ncdiff = num_ncdiff + 1

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
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       attributes of data 
c-----------------------------------------------------------------------

        call mpp_write_meta(
     &      io, field_growth2a, (/axis_xt,axis_yt,axis_t/), 
     &      'growth2', 's^-2', 'Squared Eady growth rate',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_agm, (/axis_xt,axis_yt,axis_t/), 
     &      'agm', 'cm**2/s', 'GM diffusivity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
#  if defined hl_diffusivity
        call mpp_write_meta(
     &      io, field_betaeff, (/axis_xt,axis_yt,axis_t/), 
     &      'hllength', '1/(cm*s)', 'Effective beta',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_hllength, (/axis_xt,axis_yt,axis_t/), 
     &      'hllength', 'cm', 'HL length scale',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_hldiff, (/axis_xt,axis_yt,axis_t/), 
     &      'hldiff', 'cm**2/s', 'HL diffusivity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
#  endif
#  if defined vmhs_diffusivity
        call mpp_write_meta(
     &      io, field_vmhslength, (/axis_xt,axis_yt,axis_t/), 
     &      'vmhslength', 'cm', 'VMHS length scale',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vmhsdiff, (/axis_xt,axis_yt,axis_t/), 
     &      'vmhsdiff', 'cm**2/s', 'VMHS diffusivity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
#  endif
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
c       write axes now that all meta data has been written
c-----------------------------------------------------------------------
      
        call mpp_write (io,axis_xt)
        call mpp_write (io,axis_yt)

c-----------------------------------------------------------------------
c       allocate buffer data 
c-----------------------------------------------------------------------

        allocate( buf_growth2a(1:imt,jstask:jetask) )
        allocate( buf_agm(1:imt,jstask:jetask) )
#  if  defined hl_diffusivity
        allocate( buf_betaeff(1:imt,jstask:jetask) )
        allocate( buf_hldiff(1:imt,jstask:jetask) )
#  endif
#  if defined vmhs_diffusivity
        allocate( buf_vmhslength(1:imt,jstask:jetask) )
        allocate( buf_vmhsdiff(1:imt,jstask:jetask) )
#  endif

c-----------------------------------------------------------------------
c     endif for opening the diagnostic file
c-----------------------------------------------------------------------
      endif


c-----------------------------------------------------------------------
c     point "domain" to the x-y-domain from "mom_domain" 
c-----------------------------------------------------------------------

      domain => mom_domain

c-----------------------------------------------------------------------
c     data indexed by (i,k,j) must be re-arranged to match (i,j,k). 
c-----------------------------------------------------------------------

      if( within_window_loop )then

c-----------------------------------------------------------------------
c     define variables (on data domain)
c     incoming data is already masked to 0.
c     add (1-mask)*missing_val for Ferret. 
c-----------------------------------------------------------------------

        k=1
        do j=js,je
          jrow=j+joff
          do i=1,imt
           buf_growth2a(i,jrow) = growth2a(i,jrow)
     &          + (1-tmask(i,k,j))*field_growth2a%missing
            buf_agm(i,jrow) = agm(i,jrow)
     &          + (1-tmask(i,k,j))*field_agm%missing
#  if defined hl_diffusivity
            buf_betaeff(i,jrow) = betaeff(i,jrow)
     &          + (1-tmask(i,k,j))*field_betaeff%missing
            buf_hllength(i,jrow) = hllength(i,j)
     &          + (1-tmask(i,k,j))*field_hllength%missing
            buf_hldiff(i,jrow) = hldiff(i,j)
     &          + (1-tmask(i,k,j))*field_hldiff%missing
#  endif
#  if defined vmhs_diffusivity   
            buf_vmhslength(i,jrow) = vmhslength(i,j)
     &          + (1-tmask(i,k,j))*field_vmhslength%missing
            buf_vmhsdiff(i,jrow) = vmhsdiff(i,j)
     &          + (1-tmask(i,k,j))*field_vmhsdiff%missing
#  endif
          enddo
	enddo

      endif  !end of within_window_loop if-test

c-----------------------------------------------------------------------
c     write out the diagnostics after completed all latitudes 
c-----------------------------------------------------------------------

      if( .NOT. within_window_loop )then

        call mpp_write(io, field_growth2a, domain, buf_growth2a, timrec)
        call mpp_write(io, field_agm, domain, buf_agm, timrec)
#  if defined hl_diffusivity
        call mpp_write(io, field_betaeff, domain, buf_betaeff, timrec)
        call mpp_write(io, field_hllength, domain, buf_hllength, timrec)
        call mpp_write(io, field_hldiff, domain, buf_hldiff, timrec)
#  endif
#  if defined vmhs_diffusivity 
        call mpp_write(io,field_vmhslength,domain,buf_vmhslength,timrec)
        call mpp_write(io,field_vmhsdiff,domain,buf_vmhsdiff,timrec)
#  endif
        call mpp_write(io, field_period, period, timrec)

        write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8,a,i8)')
     &'  ==>  Nonconst diffusivity fields written to', file_name
     &,' on ts = ', itt,' ', stamp
     &,' between latitudes:',jscomp,' and ', jecomp
     &,' and from the surface down to depth level ',km

      endif

      return
      end
#endif

