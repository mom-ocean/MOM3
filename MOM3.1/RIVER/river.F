      subroutine river
!-----------------------------------------------------------------------
!     This subroutine interpolates data from the time dependent dataset 
!     'river_data' to a time 'model_time + 0.5*(segtim-dtts)'.  
!     author:      M. Schmidt   e-mail=> martin.schmidt@io-warnemuende.de
!-----------------------------------------------------------------------
#ifdef river_inflow
      use time_manager_mod
      use river_mod 
      use timeinterp_mod 
      integer nb, indriv
      logical rdriv
      type (time_type) :: a_time
# include "stdunits.h"
# include "scalar.h"
# include "river.h"
# if defined show_river_details
      write(stdout,'(a)') '==> This is river'
# endif      
      a_time = model_time + river_time_offset
      do nb=1, rivbox_in_pe
        indriv = river_data(nb)%index
        dayriv = position_within_data (a_time
     &           , river_data(nb)%start_time
     &           , river_data(nb)%end_time, river_data(nb)%perriv)
# ifdef debug_river
      write (stdout,*) 'Time interpolation'
      call get_time(model_time,isec, iday)
      write (stdout,*) 'Model time days    ', iday
      write (stdout,*) 'Model time seconds ', isec
      call get_time(a_time,isec, iday)
      write (stdout,*) 'river_interpolation_time days ', iday
      write (stdout,*) 'river_interpolation_time seconds ', isec
      call get_time(river_data(nb)%start_time,isec, iday)
      write (stdout,*) 'river_data_start_time days ', iday
      write (stdout,*) 'river_data_start_time seconds ', isec
      call get_time(river_data(nb)%end_time,isec, iday)
      write (stdout,*) 'river_data_end_time days ', iday
      write (stdout,*) 'river_data_end_time seconds ', isec
# endif
# if defined show_river_details
        write (stdout,*) ' river: ',indriv,' dayriv=',dayriv
# endif
!-----------------------------------------------------------------------
!       determine the disk pointers, time weight interpolation factor,
!       and whether or not it is time to new data records.
!-----------------------------------------------------------------------
!
!       initialize the pointers (indices) to memory buffers and choose
!       interpolation method
        call timeinterp (dayriv, river_data(nb)%ninterp
     &      , river_data(nb)%tdrec, river_data(nb)%aprec
     &      , river_data(nb)%mrecriv, river_data(nb)%perriv
     &      , methodriv, inextdriv(nb), iprevdriv(nb), wprevriv(nb)
     &      , rdriv, inextmriv(nb), iprevmriv(nb))
!
	in = inextdriv(nb)
	ip = iprevdriv(nb)
# ifdef debug_river
!       Write out the geometry, name and index of the river        
        write (stdout,*) 'Time interpolation'
        write (stdout,*) 'Pointer to data records: '
        write (stdout,*) 'iprevdriv ',iprevdriv(nb)
        write (stdout,*) 'iprevmriv ',iprevmriv(nb)
        write (stdout,*) 'inextdriv ',inextdriv(nb)
        write (stdout,*) 'inextmriv ',inextmriv(nb)
        write (stdout,*) 'wprevriv  ',wprevriv(nb)
        write (stdout,*) 'river_data(ip)  ',river_data(nb)%rff(ip)
        write (stdout,*) 'river_data(in)  ',river_data(nb)%rff(in)
# endif
        river_rec(nb)%rff    = wprevriv(nb)*river_data(nb)%rff(ip)
     &                     + (1-wprevriv(nb))*river_data(nb)%rff(in)
# ifndef simple_sbc
        do n=1,ntriv
          river_rec(nb)%rtf(n) = wprevriv(nb)*river_data(nb)%rtf(n,ip)
     &                     + (1-wprevriv(nb))*river_data(nb)%rtf(n,in)
        enddo
# endif      
# if defined show_river_details
	write(stdout,'(a,i2,a,i2,a,a,a,1pe12.5,a,3(a,1pe12.5,a))') 
     &      'Setting the flux of riverbox ',nb
     &     ,' for river ',indriv,': ',trim(river_rec(nb)%name)
     &	   ,', surface area: ', river_rec(nb)%area,' cm**2'
     &	   ,', fresh water flux: ',river_rec(nb)%rff,' cm/s'
#  ifndef simple_sbc 
        write(stdout,'(a,100(a12,1x))') 
     &      'Tracer name	  : ',(rivertrname(n), n=1,ntriv)
!!	write(stdout,'(a,100(1pe12.5,1x))') 
!!     &      'Tracer concentration : ',rtc(indriv,n), n=1,ntriv
!!	write(stdout,'(a,100(1pe12.5,1x))') 
!!     &      'Tracer flux          : ',river_rec(nb)%rtf(n), n=1,ntriv
#  endif      
# endif      
      enddo
# ifdef diag_river
      call diag_out_river
# endif
      return
      end


# ifdef diag_river
      subroutine diag_out_river
c-----------------------------------------------------------------------
c     save river inflow data
c
c     author:      M. Schmidt  e-mail=> martin.schmidt@io-warnemuende.de
c-----------------------------------------------------------------------
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
      use datamod
      use time_manager_mod
      use river_mod 
# include "param.h"
# include "taskrows.h"
# include "docnam.h"
# include "diag.h"
# include "vers.h"
# include "tmngr.h"
# include "switch.h"
# include "river.h"
      character*128, dimension(nriv) :: file_name, gvtitle
      character*128 file_stamp
      character(len=60) :: cvarnam
# if defined salinity_psu
      character*3 saldim
# else 
      character*9 saldim
# endif
      data num_rec /0/
      save num_rec
      integer nb, indriv
      integer yriv, mriv, driv, hriv, miriv, sriv, psriv, pdriv
      real                  :: period_seconds
      integer, dimension(nriv) :: io
      save io
      real, dimension(nriv) :: sfft_riv, sfft_acc
      save sfft_acc
      real, dimension(nriv,ntriv) :: strf_riv, strf_acc
      save strf_acc
      type (time_type) :: a_time
c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_sfft_riv(nriv)
      type(fieldtype), save :: field_sfft_acc(nriv)
      type(fieldtype), save :: field_strf_riv(nriv,ntriv)
      type(fieldtype), save :: field_strf_acc(nriv,ntriv)
      type(fieldtype), save :: field_period

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_t

!-----------------------------------------------------------------------
!     prepare the time stamp
!-----------------------------------------------------------------------
      a_time = model_time + river_time_offset - start_time + dt_time
      call get_date (a_time, yriv, mriv, driv, hriv, miriv, sriv)
      call set_stamp (xstamp, yriv, mriv, driv, hriv, miriv, sriv )
c-----------------------------------------------------------------------
c     specify file where energy analysis is written
c-----------------------------------------------------------------------
      do n=1,nriv
        file_name(n) = file_stamp (trim(rivername(n)), xstamp, '.dta')
        gvtitle(n) = 'River Monitor for '//trim(rivername(n)) // momver
      enddo
c-----------------------------------------------------------------------
c     set averaging period, time, title, and salinity units 
c-----------------------------------------------------------------------

      call get_time 
     & (river_time_offset + river_time_offset + dt_time, psriv, pdriv)
      period_seconds = float(pdriv) * 86400. + float(psriv)
      call get_time (a_time, sriv, driv)
      if (rununits.eq.'days') then
        timrec = float(driv) + float(sriv)/86400.
        period = float(pdriv) + float(psriv)/86400.
      else
        timrec = (float(driv) + float(sriv)/86400.)/yrlen
        period = (float(pdriv) + float(psriv)/86400.)/yrlen
      endif
# if defined salinity_psu
      saldim = 'psu' 
# else 
      saldim = 'ppm-0.035' 
# endif

!-----------------------------------------------------------------------
!     calculate the diagnostic quantities
!-----------------------------------------------------------------------
      sfft_riv(:)   = c0
      strf_riv(:,:) = c0
      if (num_rec .eq. 0) then
        sfft_acc(:)   = c0
        strf_acc(:,:) = c0
      endif
      do nb=1,rivbox_in_pe
        indriv = river_rec(nb)%index
        sfft_riv(indriv) = sfft_riv(indriv)
     &                   + river_rec(nb)%rff * river_rec(nb)%area 
     &                     /nboxriv(indriv) 
        do n=1,ntriv
          strf_riv(indriv,n) = strf_riv(indriv,n)
     &                   + river_rec(nb)%rtf(n) * river_rec(nb)%area   
     &                     /nboxriv(indriv) 
        enddo
      enddo
!     Accumulate the fluxes
      sfft_acc(:) = sfft_acc(:) + sfft_riv(:) * period_seconds
      strf_acc(:,:) = strf_acc(:,:) + strf_riv(:,:) * period_seconds
!
# if defined parallel_1d
c
c       add partial sums from all processors and place in processor #1
      write(stdout,*) 'diag_river',pe,sfft_riv, nriv
      call mpp_sum( sfft_riv, nriv )
      call mpp_sum( strf_riv(1,1), nriv*ntriv )
      call mpp_sum( sfft_acc, nriv )
      call mpp_sum( strf_acc(1,1), nriv*ntriv )
# endif
c-----------------------------------------------------------------------
c     open diagnostic file and write meta information 
c-----------------------------------------------------------------------
  
      if (pn.eq.1 .and. num_rec.eq.0) then
        do nb=1,nriv
          call mpp_open(io(nb), file_name(nb), action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

c-----------------------------------------------------------------------
c       define time axis
c-----------------------------------------------------------------------

          call mpp_write_meta(
     &         io(nb), axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       define fields
c-----------------------------------------------------------------------

          call mpp_write_meta(
     &      io(nb), field_sfft_riv(nb), (/axis_t/), 
     &      'Q_r', 'cm**3/s', 'volume flux',
     &      0., +1e20, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &      io(nb), field_sfft_acc(nb), (/axis_t/), 
     &      'Q_I', 'cm**3', 'accumulated volume flux',
     &      0., +1e20, -1.0E+34, +1.0E+34)

          do n=1,ntriv 

            if (n .eq. 1) then
              call mpp_write_meta(
     &        io(nb), field_strf_riv(nb,n), (/axis_t/), 
     &        'enthalpy_flux', 'cal/s'
     &        ,'enthalpy flux',
     &        -1e10, +1e20, -1.0E+34, +1.0E+34)
              call mpp_write_meta(
     &        io(nb), field_strf_acc(nb,n), (/axis_t/), 
     &        'accumulated_enthalpy', 'cal'
     &        ,'total enthalpy inflow',
     &        -1e10, +1e20, -1.0E+34, +1.0E+34)
            endif
            if (n .ge. 3) then
              call mpp_write_meta(
     &        io(nb), field_strf_riv(nb,n), (/axis_t/), 
     &        trim(rivertrname(n)), 'tracerunits', 'other tracer',
     &        -1e10, +1e20, -1.0E+34, +1.0E+34)
              call mpp_write_meta(
     &        io(nb), field_strf_acc(nb,n), (/axis_t/), 
     &        trim(rivertrname(n)), 'tracerunits*time', 'other tracer',
     &        -1e10, +1e20, -1.0E+34, +1.0E+34)
            endif
	  enddo  
          call mpp_write_meta(
     &      io(nb), field_period, (/axis_t/), 
     &      'period', rununits, 'averaging period',
     &      +0.0E+00, +1.0E+20, -1.0E+34, +1.0E+34)
c
c       figure title 
c
          call mpp_write_meta(io(nb), 'title', cval=gvtitle(nb) )

        enddo
        num_rec = num_rec + 1
c-----------------------------------------------------------------------
c     endif for opening of diagnostic file
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

      endif
      if (pn.eq.1) then
        do nb=1,nriv
          call mpp_write(io(nb),field_sfft_riv(nb),sfft_riv(nb),timrec)
          call mpp_write(io(nb),field_sfft_acc(nb),sfft_acc(nb),timrec)
          do n=1,ntriv
            if(n.ne.2) then
              call mpp_write
     &            (io(nb), field_strf_riv(nb,n), strf_riv(nb,n), timrec)
              call mpp_write
     &            (io(nb), field_strf_acc(nb,n), strf_acc(nb,n), timrec)
            endif
	  enddo
          call mpp_write(io(nb), field_period, period, timrec)

          write (stdout,'(a,a)')
     &    '  ==>  River data written to file '//trim(file_name(nb))
     &    ,' for interpolation time '//xstamp
        enddo
      endif
      
      return
      end
# endif
#else
      return
      end
#endif      

