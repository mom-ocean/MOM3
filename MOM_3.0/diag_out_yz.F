      subroutine gyre_comp_out
#if defined gyre_components

c-----------------------------------------------------------------------
c     save meridional gyre transport components
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c     updated with mpp_io routines July 1999
c                  s. m. griffies        e-mail=> smg@gfdl.gov  
c-----------------------------------------------------------------------

      use arrays_2d
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "coord.h"
# include "diag.h"
# include "diag_region.h"
# include "switch.h"
# include "taskrows.h"
# include "tmngr.h"
# include "vers.h"
      character*128 file_name, file_stamp, gvtitle
      character*5 tname
      character*15 units
      character*5 basin_num
      character*20 id
      dimension bufmskhr(imt,jstask:jetask)
      dimension buf(jstask:jetask)
   
      data num_gyre /0/
      save num_gyre
      save io

c-----------------------------------------------------------------------
c     domain pointer 
c-----------------------------------------------------------------------

      type(domain1D), pointer :: domainy
      type(domain2D), pointer :: domainxy

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_ttn(8,nt)
      type(fieldtype), save :: field_ttn2(3,nt,nhreg)
      type(fieldtype), save :: field_period 
      type(fieldtype), save :: field_mskhr

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_xt
      type (axistype), save :: axis_yt
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where gyre components are written
c-----------------------------------------------------------------------

      file_name = file_stamp ('gyre_comp', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, and buffers 
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Gyre components' // momver

c-----------------------------------------------------------------------
c     open diagnostic file and write meta information 
c-----------------------------------------------------------------------

      if (num_gyre .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_gyre = num_gyre + 1

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

	do n=1,nt
	  if (n .eq. 1) then
	    id = "_T"
	    units = "PW"
	  elseif (n .eq. 2) then
	    id = "_S"
	    units = "giga cm**3/s"
	  else
            write (id,'(a1,i3.3)') "_",n
	    units = "tracer/s"
	  endif

          call mpp_write_meta(
     &      io, field_ttn(1,n), (/axis_yt,axis_t/), 
     &      'zonal_mean'//id, units, 'x mean'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &      io, field_ttn(2,n), (/axis_yt,axis_t/), 
     &      'zonal_eddy'//id, units, 'x eddy'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &      io, field_ttn(3,n), (/axis_yt,axis_t/), 
     &      'vert_mean'//id, units, 'z mean'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &      io, field_ttn(4,n), (/axis_yt,axis_t/), 
     &      'vert_eddy'//id, units, 'z eddy'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &      io, field_ttn(5,n), (/axis_yt,axis_t/), 
     &      'ekman'//id, units, 'ekman'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &      io, field_ttn(6,n), (/axis_yt,axis_t/), 
     &      'adv'//id, units, 'total advection'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &      io, field_ttn(7,n), (/axis_yt,axis_t/), 
     &      'diff'//id, units, 'diffusion'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &      io, field_ttn(8,n), (/axis_yt,axis_t/), 
     &      'total'//id, units, 'total transport'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        enddo
c
	do n=1,nt
          write (tname,'(a2,i3.3)') "_T",n
	  if (n .eq. 1) then
	    units = "PW"
	  elseif (n .eq. 2) then
	    units = "giga cm**3/s"
	  else
            units = "tracer/s"
	  endif
	  do m=1,nhreg
	    if (m .eq. 1) then
	      basin_num = " "
	    elseif (m .gt. 1) then
              write (basin_num,'(a4,i1)') '_bas',m
	    endif
	    id = tname//basin_num
            call mpp_write_meta(
     &      io, field_ttn2(1,n,m), (/axis_yt,axis_t/), 
     &      'tot_adv'//id, units, 'total advection'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_ttn2(2,n,m), (/axis_yt,axis_t/), 
     &      'tot_diff'//id, units, 'total diffusion'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
            call mpp_write_meta(
     &      io, field_ttn2(3,n,m), (/axis_yt,axis_t/), 
     &      'tot_xsport'//id, units, 'total transport'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
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
c       write axes after all meta data has been written
c-----------------------------------------------------------------------
      
        call mpp_write (io,axis_xt)
        call mpp_write (io,axis_yt)

c-----------------------------------------------------------------------
c     endif for opening of diagnostic file
c-----------------------------------------------------------------------
      endif

c-----------------------------------------------------------------------
c     point domainY to y-domain from mom_domain 
c-----------------------------------------------------------------------

      domainy   => mom_domain%y
      domainxy  => mom_domain

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

      do n=1,nt
        do l=1,8
	  do j=jscomp,jecomp
	    buf(j) = ttn(l,j,n)
	  enddo
          call mpp_write(io, field_ttn(l,n), domainy, buf, timrec)
        enddo
      enddo
      do n=1,nt
	do m=1,nhreg
          do l=6,8
	    do j=jscomp,jecomp
	      buf(j) = ttn2(l,j,n,m)
	    enddo
            call mpp_write(io, field_ttn2(l-5,n,m), domainy, buf,timrec)
          enddo
        enddo
      enddo
      call mpp_write(io, field_period, period, timrec)
      if (nhreg .gt. 1 .and. igyre .and. num_gyre .eq. 1) then 
        do j=jscomp,jecomp
          do i=1,imt
            bufmskhr(i,j) = float(mskhr(i,j))
            if(bufmskhr(i,j) .eq. 0.0) then
              bufmskhr(i,j) = field_mskhr%missing
            endif
          enddo
        enddo
        call mpp_write(io, field_mskhr, domainxy, bufmskhr, timrec)
      endif

      write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8)')
     &'  ==>  Gyre transports written to file '
     &,  file_name,' on ts = ', itt,' ', stamp
     &,' between latitudes:',jscomp,' and ', jecomp

      return
      end
#else
      return
      end
#endif



#if defined tracer_yz
      subroutine tyz_out 
c-----------------------------------------------------------------------
c     save zonal mean tracer components
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c     updated with mpp_io routines July 1999
c                  s. m. griffies        e-mail=> smg@gfdl.gov  
c-----------------------------------------------------------------------
c
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "taskrows.h"
# include "coord.h"
# include "diag.h"
# include "switch.h"      
# include "tmngr.h"
# include "vers.h"
      character*128 file_name, file_stamp, gvtitle
      character*64 string
      character*15 units
      dimension bufzt(km)
      data num_tyz /0/
      save num_tyz
      save io

c-----------------------------------------------------------------------
c     domain pointer 
c-----------------------------------------------------------------------

      type(domain1D), pointer :: domainy

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_tyz(nt,5)
      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_yt
      type (axistype), save :: axis_zt
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where zonal mean tracer components are written
c-----------------------------------------------------------------------

      file_name = file_stamp ('tracer_yz', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, and ocean depth in meters
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Zonally Averaged Tracer Components' // momver
      bufzt(:) = zt(:)*0.01

c-----------------------------------------------------------------------
c     open diagnostic file and write meta information 
c-----------------------------------------------------------------------

      if (num_tyz .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_tyz = num_tyz + 1

c-----------------------------------------------------------------------
c       axes (on global domain)
c-----------------------------------------------------------------------

        jbeg = mom_domain%y%global%start_index
        jend = mom_domain%y%global%end_index

        call mpp_write_meta(
     &         io, axis_yt, 'yt_j', 'degrees_N', 'Latitude of T points',
     &         cartesian='Y', domain=mom_domain%y, 
     &         data=yt(jbeg:jend))
        call mpp_write_meta(
     &         io, axis_zt, 'zt_k', 'm', 'Depth of T grid point',
     &         cartesian='Z', sense=-1, data=bufzt)
        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       attributes of data 
c-----------------------------------------------------------------------

        do n=1,nt
	  if (n .eq. 1) then
	    units = "deg C"
	  elseif (n .eq. 2) then
# if defined salinity_psu
	    units = "psu"
# else 
	    units = "ppm-0.035"
# endif
	  else
	    units = "tracer units"
	  endif
          write(string,'(a1,i3.3)') '_',n
          units = trim(units)

          call mpp_write_meta(
     &      io, field_tyz(n,1), (/axis_yt,axis_zt,axis_t/), 
     &      'tbarx'//string, units, 'tracer_yz'//string,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &      io, field_tyz(n,2), (/axis_yt,axis_zt,axis_t/), 
     &      'dTdt'//string, units//'/s', 'd_tracer_yz_dt'//string,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &      io, field_tyz(n,3), (/axis_yt,axis_zt,axis_t/), 
     &      'advT'//string, units//'/s', 'advection_yz'//string,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &      io, field_tyz(n,4), (/axis_yt,axis_zt,axis_t/), 
     &      'diffT'//string, units//'/s', 'diffusion_yz'//string,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &      io, field_tyz(n,5), (/axis_yt,axis_zt,axis_t/), 
     &      'source'//string, units//'/s', 'source_yz'//string,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

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
c       write axes now that all meta data has been written
c-----------------------------------------------------------------------
      
        call mpp_write (io,axis_yt)
        call mpp_write (io,axis_zt)

c-----------------------------------------------------------------------
c     endif for opening of diagnostic file
c-----------------------------------------------------------------------
      endif

c-----------------------------------------------------------------------
c     point domain to y-domain from mom_domain 
c-----------------------------------------------------------------------

      domainy => mom_domain%y

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

      do n=1,nt
        do m=1,5
          call mpp_write(io, field_tyz(n,m), domainy, tyz(:,:,n,m)
     &,                  timrec)
        enddo
      enddo
      call mpp_write(io, field_period, period, timrec)

      write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8)')
     &'  ==>  Zonal mean of tracer written to file '
     &,  file_name,' on ts = ', itt,' ', stamp
     &,' between latitudes:',jscomp,' and ', jecomp

      return
      end
#endif



#if defined meridional_overturning
# ifdef gent_mcwilliams
      subroutine vms_out (vms, vms_gm)
# else
      subroutine vms_out (vms)
# endif

c-----------------------------------------------------------------------
c     write vertical-meridional streamfunction
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c     updated with mpp_io routines July 1999
c                  s. m. griffies        e-mail=> smg@gfdl.gov  
c-----------------------------------------------------------------------
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "diag.h"      
# include "coord.h"
# include "switch.h"      
# include "tmngr.h"
# include "taskrows.h"
# include "vers.h"
      character*128 file_name, file_stamp, gvtitle
      dimension vms(jmt,km+1), zw0(0:km)
# ifdef gent_mcwilliams
      dimension vms_gm(jmt,km+1) 
# endif
      data num_vms /0/
      save num_vms
      save io

c-----------------------------------------------------------------------
c     domain pointer 
c-----------------------------------------------------------------------

      type(domain1D), pointer :: domainy

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_vms
      type(fieldtype), save :: field_vms_gm
      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_yu
      type (axistype), save :: axis_zw
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where vertical-meridional streamfunction is written
c-----------------------------------------------------------------------

      file_name = file_stamp ('overturn', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, and ocean depth in meters
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Vertical-meridional overturning streamfunction'
     &          // momver
      zw0(0) = c0
      do k=1,km
       zw0(k) = zw(k)*0.01
      enddo

c-----------------------------------------------------------------------
c     open diagnostic file and write meta information 
c-----------------------------------------------------------------------

      if (num_vms .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_vms = num_vms + 1

c-----------------------------------------------------------------------
c       axes (on global domain)
c-----------------------------------------------------------------------

        jbeg = mom_domain%y%global%start_index
        jend = mom_domain%y%global%end_index

        call mpp_write_meta(
     &         io, axis_yu, 'yu_j', 'degrees_N', 'Latitude of U points',
     &         cartesian='Y', domain=mom_domain%y, 
     &         data=yu(jbeg:jend))
        call mpp_write_meta(
     &         io, axis_zw, 'zw_k', 'm', 'Depth of U cell bottom',
     &         cartesian='Z', sense=-1, data=zw0)
        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       attributes of data 
c-----------------------------------------------------------------------

# ifdef gent_mcwilliams
        call mpp_write_meta(
     &    io, field_vms, (/axis_yu,axis_zw,axis_t/), 
     &    'vms', 'Sv', 'Global strmfct from v+vgm',
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
	call mpp_write_meta(
     &    io, field_vms_gm, (/axis_yu,axis_zw,axis_t/), 
     &    'vms_gm', 'Sv', 'Global strmfct from vgm',
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
# else 
        call mpp_write_meta(
     &    io, field_vms, (/axis_yu,axis_zw,axis_t/), 
     &    'vms', 'Sv', 'Global strmfct from v',
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
# endif
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
c       write axes now that all meta data has been written
c-----------------------------------------------------------------------
      
        call mpp_write (io,axis_yu)
        call mpp_write (io,axis_zw)

c-----------------------------------------------------------------------
c     endif for opening of diagnostic file
c-----------------------------------------------------------------------

      endif

c-----------------------------------------------------------------------
c     point domain to y-domain from mom_domain 
c-----------------------------------------------------------------------

      domainy => mom_domain%y


c-----------------------------------------------------------------------
c     set missing values 
c-----------------------------------------------------------------------

      do k=1,km+1
        do j=jstask,jetask
          if(vms(j,k) .eq. 0.0) vms(j,k) = field_vms%missing
# ifdef gent_mcwilliams
           if(vms_gm(j,k) .eq. 0.0) vms_gm(j,k) = field_vms_gm%missing           
# endif
        enddo
      enddo
      
c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------
      
      call mpp_write(io, field_vms, domainy, vms, timrec)
# ifdef gent_mcwilliams
      call mpp_write(io, field_vms_gm, domainy, vms_gm, timrec)
# endif
      call mpp_write(io, field_period, period, timrec)

      write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8)')
     &'  ==>  Vertical-meridional streamfcn written to file '
     &,  file_name,' on ts = ', itt,' ', stamp
     &,' between latitudes:',jscomp,' and ', jecomp

      return
      end
#endif


#if defined meridional_overturning_density 
      subroutine vms_sigma_outA(sigmsf, sigvert)

c-----------------------------------------------------------------------
c     save sigma-meridional streamfunction to netcdf file
c
c     authors:  r. c. pacanowski      e-mail=> rcp@gfdl.gov
c               s. m. griffies        e-mail=> smg@gfdl.gov
c
c     updated with mpp_io routines July 1999
c               s. m. griffies        e-mail=> smg@gfdl.gov
c
c     NOTE: This diagnostic has not been exercised much at GFDL
c-----------------------------------------------------------------------
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
      dimension sigmsf(jmt,km), sigvert(km)
      data num_vmssig /0/
      save num_vmssig
      save io

c-----------------------------------------------------------------------
c     domain pointer 
c-----------------------------------------------------------------------

      type(domain1D), pointer :: domainy

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_vmssig
      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_yu
      type (axistype), save :: axis_sig
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where vertical-meridional streamfunction is written
c-----------------------------------------------------------------------

      file_name = file_stamp ('overturn_sigmaA', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Sigma-Meridional Streamfunction' // momver

c-----------------------------------------------------------------------
c     open diagnostic file and write meta information 
c-----------------------------------------------------------------------

      if (num_vmssig .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_vmssig = num_vmssig + 1

c-----------------------------------------------------------------------
c       axes (on global domain)
c-----------------------------------------------------------------------

        jbeg = mom_domain%y%global%start_index
        jend = mom_domain%y%global%end_index

        call mpp_write_meta(
     &         io, axis_yu, 'yu_j', 'degrees_N', 'Latitude of U points',
     &         cartesian='Y', domain=mom_domain%y, 
     &         data=yu(jbeg:jend))
        call mpp_write_meta(
     &         io, axis_sig, 'sigma', 'm', 'Potential density level',
     &         cartesian='Z', sense=-1, data=sigvert)
        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       attributes of data 
c-----------------------------------------------------------------------

        call mpp_write_meta(
     &    io, field_vmssig, (/axis_yu,axis_sig,axis_t/), 
     &    'sigmsf', 'Sv', 'Global strmfct on sigma0 surfaces',
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
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
c       write axes now that all meta data has been written
c-----------------------------------------------------------------------
      
        call mpp_write (io,axis_yu)
        call mpp_write (io,axis_sig)

c-----------------------------------------------------------------------
c     endif for opening of diagnostic file
c-----------------------------------------------------------------------

      endif

c-----------------------------------------------------------------------
c     point domain to y-domain from mom_domain 
c-----------------------------------------------------------------------

      domainy => mom_domain%y

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

      call mpp_write(io, field_vmssig, domainy, sigmsf, timrec)
      call mpp_write(io, field_period, period, timrec)

      write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8)')
     &'  ==>  Density-meridional streamfcnA written to file '
     &,  file_name,' on ts = ', itt,' ', stamp
     &,' between latitudes:',jscomp,' and ', jecomp

      return
      end
#endif

#if defined meridional_overturning_density 
      subroutine vms_sigma_outB(sigmsf, sigvert)

c-----------------------------------------------------------------------
c     save sigma-meridional streamfunction to netcdf file
c
c     authors:  r. c. pacanowski      e-mail=> rcp@gfdl.gov
c               s. m. griffies        e-mail=> smg@gfdl.gov
c
c     updated with mpp_io routines July 1999
c               s. m. griffies        e-mail=> smg@gfdl.gov
c-----------------------------------------------------------------------
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
      dimension sigmsf(jmt,km), sigvert(km)
      data num_vmssig /0/
      save num_vmssig
      save io

c-----------------------------------------------------------------------
c     domain pointer 
c-----------------------------------------------------------------------

      type(domain1D), pointer :: domainy

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_vmssig
      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_yu
      type (axistype), save :: axis_sig
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where vertical-meridional streamfunction is written
c-----------------------------------------------------------------------

      file_name = file_stamp ('overturn_sigmaB', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Sigma-Meridional Streamfunction' // momver

c-----------------------------------------------------------------------
c     open diagnostic file and write meta information 
c-----------------------------------------------------------------------

      if (num_vmssig .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_vmssig = num_vmssig + 1

c-----------------------------------------------------------------------
c       axes (on global domain)
c-----------------------------------------------------------------------

        jbeg = mom_domain%y%global%start_index
        jend = mom_domain%y%global%end_index

        call mpp_write_meta(
     &         io, axis_yu, 'yu_j', 'degrees_N', 'Latitude of U points',
     &         cartesian='Y', domain=mom_domain%y, 
     &         data=yu(jbeg:jend))
        call mpp_write_meta(
     &         io, axis_sig, 'sigma', 'm', 'Potential density level',
     &         cartesian='Z', sense=-1, data=sigvert)
        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       attributes of data 
c-----------------------------------------------------------------------

        call mpp_write_meta(
     &    io, field_vmssig, (/axis_yu,axis_sig,axis_t/), 
     &    'sigmsf', 'Sv', 'Global strmfct on sigma2 surfaces',
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
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
c       write axes now that all meta data has been written
c-----------------------------------------------------------------------
      
        call mpp_write (io,axis_yu)
        call mpp_write (io,axis_sig)

c-----------------------------------------------------------------------
c     endif for opening of diagnostic file
c-----------------------------------------------------------------------

      endif

c-----------------------------------------------------------------------
c     point domain to y-domain from mom_domain 
c-----------------------------------------------------------------------

      domainy => mom_domain%y

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

      call mpp_write(io, field_vmssig, domainy, sigmsf, timrec)
      call mpp_write(io, field_period, period, timrec)

      write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8)')
     &'  ==>  Density-meridional streamfcnB written to file '
     &,  file_name,' on ts = ', itt,' ', stamp
     &,' between latitudes:',jscomp,' and ', jecomp

      return
      end
#endif


#if defined meridional_overturning_basin
# ifdef gent_mcwilliams
      subroutine vms_basin_out (vms, vms_gm)
# else
      subroutine vms_basin_out (vms)
# endif
c
c--------------------------------------------------------------------------
c     save vertical-meridional streamfunction by basin to netcdf file
c     author:      k. dixon      e-mail=> kd@gfdl.gov
c
c     updated with mpp_io routines July 1999
c                  s. m. griffies        e-mail=> smg@gfdl.gov 
c--------------------------------------------------------------------------
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
      dimension vms(jmt,km+1,0:nbasin), zw0(0:km)
# ifdef gent_mcwilliams
      dimension vms_gm(jmt,km+1,0:nbasin) 
# endif
      dimension vmstmp(jmt,km+1)
      data num_vms /0/
      save num_vms
      save io

c-----------------------------------------------------------------------
c     domain pointer 
c-----------------------------------------------------------------------

      type(domain1D), pointer :: domainy

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_vms(0:nbasin)
      type(fieldtype), save :: field_vms_gm(0:nbasin)
      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_yu
      type (axistype), save :: axis_zw
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where vertical-meridional streamfunction is written
c-----------------------------------------------------------------------

      file_name = file_stamp ('overturn', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, and ocean depth in meters
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Vertical-meridional overturning streamfunction'
     &          // momver
      zw0(0) = c0
      do k=1,km
       zw0(k) = zw(k)*0.01
      enddo

c-----------------------------------------------------------------------
c     open diagnostic file and write meta information 
c-----------------------------------------------------------------------

      if (num_vms .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_vms = num_vms + 1

c-----------------------------------------------------------------------
c       axes (on global domain)
c-----------------------------------------------------------------------

        jbeg = mom_domain%y%global%start_index
        jend = mom_domain%y%global%end_index

        call mpp_write_meta(
     &         io, axis_yu, 'yu_j', 'degrees_N', 'Latitude of U points',
     &         cartesian='Y', domain=mom_domain%y, 
     &         data=yu(jbeg:jend))
        call mpp_write_meta(
     &         io, axis_zw, 'zw_k', 'm', 'Depth of U cell bottom',
     &         cartesian='Z', sense=-1, data=zw0)
        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       Attributes of data.  Note that Test Case 0 
c       only computes the global streamfunction.
c       Examples are given of Atlantic+Arctic plus Indian+Pacific,
c       whose values must be computed in subroutine merid1st.
c-----------------------------------------------------------------------

# ifdef gent_mcwilliams
        call mpp_write_meta(
     &    io, field_vms(0), (/axis_yu,axis_zw,axis_t/), 
     &    'vms', 'Sv', 'Global strmfct from v+vgm',
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &    io, field_vms(1), (/axis_yu,axis_zw,axis_t/), 
     &    'vms2', 'Sv', 'Global strmfct from v+vgm',
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)

c       call mpp_write_meta(
c     &    io, field_vms(1), (/axis_yu,axis_zw,axis_t/), 
c     &    'vmsa', 'Sv', 'Atlantic+Arctic strmfct from v+vgm',
c     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c        call mpp_write_meta(
c     &    io, field_vms(2), (/axis_yu,axis_zw,axis_t/), 
c     &    'vmspi', 'Sv', 'Indian+Pacific strmfct from v+vgm',
c     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c
        call mpp_write_meta(
     &    io, field_vms_gm(0), (/axis_yu,axis_zw,axis_t/), 
     &    'vms_gm', 'Sv', 'Global strmfct from vgm',
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &    io, field_vms_gm(1), (/axis_yu,axis_zw,axis_t/), 
     &    'vms_gm2', 'Sv', 'Global strmfct from vgm',
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)

c       call mpp_write_meta(
c     &    io, field_vms_gm(1), (/axis_yu,axis_zw,axis_t/), 
c     &    'vmsa_gm', 'Sv', 'Atlantic+Arctic strmfct from vgm',
c     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c        call mpp_write_meta(
c     &    io, field_vms_gm(2), (/axis_yu,axis_zw,axis_t/), 
c     &    'vmspi_gm', 'Sv', 'Indian+Pacific strmfct from vgm',
c     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)

# else  
        call mpp_write_meta(
     &    io, field_vms(0), (/axis_yu,axis_zw,axis_t/), 
     &    'vms', 'Sv', 'Global strmfct from v',
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &    io, field_vms(1), (/axis_yu,axis_zw,axis_t/), 
     &    'vms2', 'Sv', 'Global strmfct from v',
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)

c       call mpp_write_meta(
c     &    io, field_vms(1), (/axis_yu,axis_zw,axis_t/), 
c     &    'vmsa', 'Sv', 'Atlantic+Arctic strmfct from v',
c     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c        call mpp_write_meta(
c     &    io, field_vms(2), (/axis_yu,axis_zw,axis_t/), 
c     &    'vmspi', 'Sv', 'Indian+Pacific strmfct from v',
c     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)

# endif 
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
c       write axes now that all meta data has been written
c-----------------------------------------------------------------------
      
        call mpp_write (io,axis_yu)
        call mpp_write (io,axis_zw)

c-----------------------------------------------------------------------
c     endif for opening of diagnostic file
c-----------------------------------------------------------------------

      endif

c-----------------------------------------------------------------------
c     point domain to y-domain from mom_domain 
c-----------------------------------------------------------------------

      domainy => mom_domain%y

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

      do l=0,nbasin
        do k=1,km+1
          do j=jstask,jetask
            if(vms(j,k,l) .eq. 0.0) then
              vmstmp(j,k) = field_vms(l)%missing
            else
              vmstmp(j,k) = vms(j,k,l)
            endif
          enddo
        enddo
        call mpp_write(io, field_vms(l), domainy, vmstmp, timrec)
      enddo
# ifdef gent_mcwilliams
      do l=0,nbasin
        do k=1,km+1
          do j=jstask,jetask
            if(vms_gm(j,k,l) .eq. 0.0) then
              vmstmp(j,k) = field_vms_gm(l)%missing
            else
              vmstmp(j,k) = vms_gm(j,k,l)
            endif
          enddo
        enddo
        call mpp_write(io, field_vms_gm(l), domainy, vmstmp, timrec)
      enddo
# endif

      write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8)')
     &'  ==>  Vertical-meridional streamfcn written to '
     &,  file_name,' on ts = ', itt,' ', stamp
     &,' between latitudes:',jscomp,' and ', jecomp

      return
      end
#endif

#if defined show_zonal_mean_of_sbc
      subroutine zmsbc_out
c-----------------------------------------------------------------------
c     save zonal means of surface boundary conditions data
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c     updated with mpp_io routines July 1999
c                  s. m. griffies        e-mail=> smg@gfdl.gov 
c-----------------------------------------------------------------------
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "taskrows.h"
# include "coord.h"
# include "diag.h"
# include "switch.h"
# include "tmngr.h"
# include "vers.h"
      character*128 file_name, file_stamp, gvtitle
      character*15 var, unit
      character*30 namelong

# if defined salinity_psu
      character*3 saldim
# else 
      character*9 saldim
# endif
      data num_zmbc /0/
      save num_zmbc
      save io

c-----------------------------------------------------------------------
c     domain pointer 
c-----------------------------------------------------------------------

      type(domain1D), pointer :: domainy

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_zmsmf(2)
      type(fieldtype), save :: field_zmsm(2)
      type(fieldtype), save :: field_zmstf(nt)
      type(fieldtype), save :: field_zmsff
      type(fieldtype), save :: field_zmst(nt)
      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_yu
      type (axistype), save :: axis_yt
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where zonal mean sbc are written
c-----------------------------------------------------------------------

      file_name = file_stamp ('zmean_sbc', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, salinity units, 
c     and ocean depth in meters
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Zonal Mean SBCs' // momver
# if defined salinity_psu
      saldim = 'psu' 
# else 
      saldim = 'ppm-0.035' 
# endif

c-----------------------------------------------------------------------
c     open diagnostic file and write meta information 
c-----------------------------------------------------------------------

      if (num_zmbc .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_zmbc = num_zmbc + 1

c-----------------------------------------------------------------------
c       axes (on global domain)
c-----------------------------------------------------------------------

        jbeg = mom_domain%y%global%start_index
        jend = mom_domain%y%global%end_index

        call mpp_write_meta(
     &         io, axis_yt, 'yt_j', 'degrees_N', 'Latitude of T points',
     &         cartesian='Y', domain=mom_domain%y, 
     &         data=yt(jbeg:jend))
        call mpp_write_meta(
     &         io, axis_yu, 'yu_j', 'degrees_N', 'Latitude of U points',
     &         cartesian='Y', domain=mom_domain%y, 
     &         data=yu(jbeg:jend))
        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       attributes of data 
c-----------------------------------------------------------------------

        call mpp_write_meta(
     &      io, field_zmsmf(1), (/axis_yu,axis_t/), 
     &      'taux', 'dyn/cm**2', 'zonal mean taux',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_zmsmf(2), (/axis_yu,axis_t/), 
     &      'tauy', 'dyn/cm**2', 'zonal mean tauy',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_zmsm(1), (/axis_yu,axis_t/), 
     &      'surf_u', 'cm/s', 'zonal mean  zonal velocity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_zmsm(2), (/axis_yu,axis_t/), 
     &      'surf_v', 'cm/s', 'zonal mean meridional velocity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_zmsff, (/axis_yt,axis_t/), 
     &      'sff', 'Sv', 'zonal mean water flux',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

        do n=1,nt
	  if (n .eq. 1) then
	    var   = 'hflx'
	    unit  = 'W/m**2'
            namelong = 'zonal mean heat flux'
	  elseif (n .eq. 2) then
	    var = 'sflx'
	    unit = saldim//'*cm/s'
            namelong = 'zonal mean salt flux'
	  else
            write (var,'(a,i3.3)') 'stf_',n
	    write (unit,'(a)') 'tracer units'
            namelong = 'zonal mean tracer flux'
	  endif
          call mpp_write_meta(
     &      io, field_zmstf(n), (/axis_yt,axis_t/), 
     &      var, unit, namelong,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
	enddo

        do n=1,nt
	  if (n .eq. 1) then
	    var = 'surf_T'
	    unit = 'deg C'
            namelong = 'zonal mean SST'
	  elseif (n .eq. 2) then
	    var = 'surf_S'
	    unit = saldim
            namelong = 'zonal mean SSS'
	  else
            write (var,'(a,i3.3)') 'tracer_',n
	    write (unit,'(a)') 'unknown'
            write (namelong,'(a,i3.3)') 'zonal mean surface tracer_',n
	  endif
          call mpp_write_meta(
     &      io, field_zmst(n), (/axis_yt,axis_t/), 
     &      var, unit, namelong,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
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

c 
c       write axes after all meta data has been written
c      
        call mpp_write (io,axis_yt)
        call mpp_write (io,axis_yu)

c-----------------------------------------------------------------------
c     endif for opening of diagnostic file
c-----------------------------------------------------------------------
      endif

c-----------------------------------------------------------------------
c     point domain to y-domain from mom_domain 
c-----------------------------------------------------------------------

      domainy => mom_domain%y

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

      call mpp_write(io, field_zmsmf(1), domainy, zmsmf(:,1), timrec)
      call mpp_write(io, field_zmsmf(2), domainy, zmsmf(:,2), timrec)
      call mpp_write(io, field_zmsm(1), domainy, zmsm(:,1), timrec)
      call mpp_write(io, field_zmsm(2), domainy, zmsm(:,2), timrec)
      do n=1,nt
        call mpp_write(io, field_zmstf(n), domainy, zmstf(:,n), timrec)
      enddo
      do n=1,nt
        call mpp_write(io, field_zmst(n), domainy, zmst(:,n), timrec)
      enddo
      call mpp_write(io, field_period, period, timrec)

      write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8)')
     &'  ==>  Zonal mean SBCs written to file '
     &,  file_name,' on ts = ', itt,' ', stamp
     &,' between latitudes:',jscomp,' and ', jecomp

      return
      end
#endif


#if defined meridional_tracer_budget
      subroutine mtb_out
c-----------------------------------------------------------------------
c     save meridional tracer budget
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c     updated with mpp_io routines July 1999
c                  s. m. griffies        e-mail=> smg@gfdl.gov  
c-----------------------------------------------------------------------
      use arrays_2d, only : msktmb
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "coord.h"
# include "diag.h"      
# include "switch.h"
# include "taskrows.h"
# include "tmngr.h"
# include "vers.h"
      character*128 file_name, file_stamp, gvtitle
      character*15 tname, units
      character*3 basin_num
      character*20 id
      dimension bufmtb(imt,jmt)
      data num_mtb /0/
      save num_mtb
      save io

c-----------------------------------------------------------------------
c     domain pointer 
c-----------------------------------------------------------------------

      type(domain1D), pointer :: domainy
      type(domain2D), pointer :: domainxy

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_tstor(nt,ntmbb)
      type(fieldtype), save :: field_tdiv(nt,nhreg)
      type(fieldtype), save :: field_tflux(nt,nhreg)
      type(fieldtype), save :: field_tdif(nt,nhreg)
      type(fieldtype), save :: field_tsorc(nt,nhreg)
      type(fieldtype), save :: field_smdvol(nt,nhreg)
      type(fieldtype), save :: field_bufmtb
      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_xt
      type (axistype), save :: axis_yt
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where meridional tracer budget is written
c-----------------------------------------------------------------------

      file_name = file_stamp ('tracer_bud', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, 
c-----------------------------------------------------------------------

      period  = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen - 0.5*period
      else
        timrec = prelyr - 0.5*period/yrlen
      endif
      gvtitle = 'Meridional Tracer Budget' // momver

c-----------------------------------------------------------------------
c     open diagnostic file and write meta information 
c-----------------------------------------------------------------------

      if (num_mtb .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_mtb = num_mtb + 1

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

	do n=1,nt
	  if (n .eq. 1) then
	    tname = "_T"
	    units = "PW"
	  elseif (n .eq. 2) then
	    tname = "_S"
	    units = "giga cm**3/s"
	  else
            write (tname,'(a,i3.3)') "_",n
	    units = "tracer/s"
	  endif
	  do m=1,ntmbb
	    if (m .eq. 1) then
	      basin_num = " "
	    elseif (m .gt. 1 .and. m .le. 9) then
              write (basin_num,'(a,i1)') '_basin_',m
	    else
              write (basin_num,'(a,i2)') '_basin',m
	    endif
	    id = tname//basin_num

            call mpp_write_meta(
     &      io, field_tstor(n,m), (/axis_yt,axis_t/), 
     &      'tstor'//id, units, 'tstor'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

            call mpp_write_meta(
     &      io, field_tdiv(n,m), (/axis_yt,axis_t/), 
     &      'tdiv'//id, units, 'tdiv'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

            call mpp_write_meta(
     &      io, field_tflux(n,m), (/axis_yt,axis_t/), 
     &      'tflux'//id, units, 'tflux'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

            call mpp_write_meta(
     &      io, field_tdif(n,m), (/axis_yt,axis_t/), 
     &      'tdiff'//id, units, 'tdiff'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

            call mpp_write_meta(
     &      io, field_tsorc(n,m), (/axis_yt,axis_t/), 
     &      'tsource'//id, units, 'tsource'//id,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

            if (n .eq. 1) then
              call mpp_write_meta(
     &        io, field_smdvol(n,m), (/axis_yt,axis_t/), 
     &        'volume'//basin_num, 'cm**3', 'basin volume',
     &        -1.e6, 1.e6, -1.0E+34, +1.0E+34)
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
     &      io, field_bufmtb, (/axis_xt,axis_yt,axis_t/), 
     &      'mskhr', 'basin numbers', 'numbers',
     &      0.0, 100.0, -1.0E+34, +1.0E+34)
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
c     point domain to y-domain from mom_domain 
c-----------------------------------------------------------------------

      domainy  => mom_domain%y
      domainxy => mom_domain

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

      do n=1,nt
	do m=1,ntmbb
         call mpp_write(io,field_tstor(n,m),domainy,tstor(:,n,m),timrec)
         call mpp_write(io,field_tdiv(n,m),domainy,tdiv(:,n,m),timrec)
         call mpp_write(io,field_tflux(n,m),domainy,tflux(:,n,m),timrec)
         call mpp_write(io,field_tdif(n,m),domainy,tdif(:,n,m),timrec)
         call mpp_write(io,field_tsorc(n,m),domainy,tsorc(:,n,m),timrec)
	 if (n .eq. 1) then
           call mpp_write(io,field_smdvol(n,m),domainy,smdvol(:,m)
     &          ,timrec)
	 endif
        enddo
      enddo
      call mpp_write(io, field_period, period, timrec)

      if (ntmbb .gt. 1) then
        if (itmb .and. num_mtb .eq. 1) then
          do i=1,imt
            do j=jscomp,jecomp 
              bufmtb(i,j) = float(msktmb(i,j))
            enddo
          enddo
          call mpp_write(io, field_bufmtb, domainxy, bufmtb, timrec)
        endif
      endif

      write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8)')
     &'  ==>  Meridional tracer budget written to file '
     &,  file_name,' on ts = ', itt,' ', stamp
     &,' between latitudes:',jscomp,' and ', jecomp

      return
      end
#endif






