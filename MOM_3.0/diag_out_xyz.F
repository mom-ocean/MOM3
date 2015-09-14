      subroutine cmix_out (joff, js, je, within_window_loop) 
#if defined save_mixing_coeff 

c-----------------------------------------------------------------------
c     save estimated mixing coefficients data
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c     updates to mpp_io made July 1999 by 
c                  s. m. griffies        e-mail=> smg@gfdl.gov 
c-----------------------------------------------------------------------
      use memory_window
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "coord.h"
# include "derived_options.h"
# include "diag.h"
# include "mw.h"
# include "switch.h"
# include "taskrows.h"
# include "tmngr.h"
# include "vers.h"
# ifdef isoneutralmix
#  include "isopyc.h"
      dimension aibuf(imt,km)
# endif
# ifdef kppvmix
#  include "vmixc.h"
# endif
      character*128 file_name, file_stamp, gvtitle
      logical :: within_window_loop
      dimension bufzt(km), bufzw(km)
      data num_cmix /0/
      save num_cmix
      save io
      dimension xue(imt), yun(jmt)

# if defined diagnostics_himem
c
c     domain pointer 
      type(domain2D), pointer :: domain
c
c     himem buffers all the memory windows
c
      real, dimension(:,:,:), save, allocatable :: buf_cet
      real, dimension(:,:,:), save, allocatable :: buf_cnt
      real, dimension(:,:,:), save, allocatable :: buf_cbt
      real, dimension(:,:,:), save, allocatable :: buf_ceu
      real, dimension(:,:,:), save, allocatable :: buf_cnu
      real, dimension(:,:,:), save, allocatable :: buf_cbu
#  if defined gm_skew || defined redi_diffusion
      real, dimension(:,:,:), save, allocatable :: buf_Ai_ez
      real, dimension(:,:,:), save, allocatable :: buf_Ai_nz
      real, dimension(:,:,:), save, allocatable :: buf_Ai_bx
      real, dimension(:,:,:), save, allocatable :: buf_Ai_by
#   if defined redi_diffusion
      real, dimension(:,:,:), save, allocatable :: buf_K11
      real, dimension(:,:,:), save, allocatable :: buf_K22
      real, dimension(:,:,:), save, allocatable :: buf_K33
#   endif
#  endif
#  if defined kppvmix
      real, dimension(:,:,:), save, allocatable :: buf_ghats
#  endif
# endif

# if defined diagnostics_lomem
      type(domain2D) :: domain
c
c     lomem buffers only within a single memory window
c
      real :: buf_cet(1:imt,js+joff:je+joff,km)
      real :: buf_cnt(1:imt,js+joff:je+joff,km)
      real :: buf_cbt(1:imt,js+joff:je+joff,km)
      real :: buf_ceu(1:imt,js+joff:je+joff,km)
      real :: buf_cnu(1:imt,js+joff:je+joff,km)
      real :: buf_cbu(1:imt,js+joff:je+joff,km)
#  if defined gm_skew || defined redi_diffusion
      real :: buf_Ai_ez(1:imt,js+joff:je+joff,km)
      real :: buf_Ai_nz(1:imt,js+joff:je+joff,km)
      real :: buf_Ai_bx(1:imt,js+joff:je+joff,km)
      real :: buf_Ai_by(1:imt,js+joff:je+joff,km)
#   if defined redi_diffusion
      real :: buf_K11(1:imt,js+joff:je+joff,km)
      real :: buf_K22(1:imt,js+joff:je+joff,km)
      real :: buf_K33(1:imt,js+joff:je+joff,km)
#   endif
#  endif
#  if defined kppvmix
      real :: buf_ghats(1:imt,js+joff:je+joff,km)
#  endif
# endif

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_cet
      type(fieldtype), save :: field_cnt
      type(fieldtype), save :: field_cbt
      type(fieldtype), save :: field_ceu
      type(fieldtype), save :: field_cnu
      type(fieldtype), save :: field_cbu
#  if defined gm_skew || defined redi_diffusion
      type(fieldtype), save :: field_Ai_ez
      type(fieldtype), save :: field_Ai_nz
      type(fieldtype), save :: field_Ai_bx
      type(fieldtype), save :: field_Ai_by
#   if defined redi_diffusion
      type(fieldtype), save :: field_K11
      type(fieldtype), save :: field_K22
      type(fieldtype), save :: field_K33
#   endif
#  endif
#  if defined kppvmix
      type(fieldtype), save :: field_ghats
#  endif
      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_xt
      type (axistype), save :: axis_yt
      type (axistype), save :: axis_zt
      type (axistype), save :: axis_xu
      type (axistype), save :: axis_yu
      type (axistype), save :: axis_zw
      type (axistype), save :: axis_xue
      type (axistype), save :: axis_yun
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where mixing coefficients are written
c-----------------------------------------------------------------------

      file_name = file_stamp ('cmix', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, and ocean depth in meters
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Mixing coefficients' // momver
      bufzt(:) = zt(:)*0.01
      bufzw(:) = zw(:)*0.01

c-----------------------------------------------------------------------
c     check whether diagnostic file(s) should be opened
c-----------------------------------------------------------------------

      if (num_cmix .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_cmix = num_cmix + 1

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
     &      io, axis_zt, 'zt_k', 'm', 'Depth of T grid point',
     &      cartesian='Z', sense=-1, data=bufzt)

        call mpp_write_meta(
     &      io, axis_xu, 'xu_i', 'degrees_E', 'Longitude of U points',
     &      cartesian='X', domain=mom_domain%x, 
     &      data=xu(ibeg:iend))

        call mpp_write_meta(
     &      io, axis_yu, 'yu_j', 'degrees_N', 'Latitude of U points',
     &      cartesian='Y', domain=mom_domain%y, 
     &      data=yu(jbeg:jend))

        call mpp_write_meta(
     &      io, axis_zw, 'zw_k', 'm', 'Depth of cell bottom',
     &      cartesian='Z', sense=-1, data=bufzw)
 
        do i=1,imt
	  if (i. lt. imt) then
	    xue(i) = xt(i+1)
	  else
	    xue(i) = xue(i-1) + (xt(imt) - xt(imt-1))
	  endif
        enddo
        call mpp_write_meta(
     &      io, axis_xue, 'xue_i', 'degrees_E', 'East face of U points',
     &      cartesian='X', domain=mom_domain%x, data=xue(ibeg:iend))

        do j=1,jmt
	  if (j. lt. jmt) then
	    yun(j) = yt(j+1)
	  else
	    yun(j) = yun(j-1) + (yt(jmt) - yt(jmt-1))
	  endif
        enddo
        call mpp_write_meta(
     &      io, axis_yun, 'yun_j', 'degrees_N','North face of U points',
     &      cartesian='Y', domain=mom_domain%y, data=yun(jbeg:jend))

        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       attributes of data 
c-----------------------------------------------------------------------
c
c       tracer diffusivities
c
        call mpp_write_meta(
     &      io, field_cet, (/axis_xu,axis_yt,axis_zt,axis_t/), 
     &      'diff_cet', 'cm**2/s', 'diff_cet',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_cnt, (/axis_xt,axis_yu,axis_zt,axis_t/), 
     &      'diff_cnt', 'cm**2/s', 'diff_cnt',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_cbt, (/axis_xt,axis_yt,axis_zw,axis_t/), 
     &      'diff_cbt', 'cm**2/s', 'diff_cbt',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c
c       viscosities 
c
        call mpp_write_meta(
     &      io, field_ceu, (/axis_xue,axis_yu,axis_zt,axis_t/), 
     &      'diff_ceu', 'cm**2/s', 'diff_ceu',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_cnu, (/axis_xu,axis_yun,axis_zt,axis_t/), 
     &      'diff_cnu', 'cm**2/s', 'diff_cnu',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_cbu, (/axis_xu,axis_yu,axis_zw,axis_t/), 
     &      'diff_cbu', 'cm**2/s', 'diff_cbu',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

# if defined gm_skew || defined redi_diffusion
c
c       isoneutral mixing coefficients and tensor components 
c
        call mpp_write_meta(
     &      io, field_Ai_ez, (/axis_xu,axis_yt,axis_zt,axis_t/), 
     &      'Ai_ez', 'cm**2/s', 'Aie*slope',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_Ai_nz, (/axis_xt,axis_yu,axis_zt,axis_t/), 
     &      'Ai_nz', 'cm**2/s', 'Ain*slope',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_Ai_bx, (/axis_xt,axis_yt,axis_zw,axis_t/), 
     &      'Ai_bx', 'cm**2/s', 'Aibx*slope',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_Ai_by, (/axis_xt,axis_yt,axis_zw,axis_t/), 
     &      'Ai_by', 'cm**2/s', 'Aiby*slope',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
#  if defined redi_diffusion
        call mpp_write_meta(
     &      io, field_K11, (/axis_xu,axis_yt,axis_zt,axis_t/), 
     &      'K11', 'cm**2/s', 'K11 from Redi tensor',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_K22, (/axis_xt,axis_yu,axis_zt,axis_t/), 
     &      'K22', 'cm**2/s', 'K22 from Redi tensor',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_K33, (/axis_xt,axis_yt,axis_zw,axis_t/), 
     &      'K33', 'cm**2/s', 'K33 from Redi tensor',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
#  endif
# endif
# if defined kppvmix
c
c       non-local term from KPP
c
        call mpp_write_meta(
     &      io, field_ghats, (/axis_xt,axis_yt,axis_zw,axis_t/), 
     &      'ghats', 's/cm**2', 'nonlocal transfer term in KPP',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
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
      
        call mpp_write (io,axis_xt)
        call mpp_write (io,axis_yt)
        call mpp_write (io,axis_zt)
        call mpp_write (io,axis_xu)
        call mpp_write (io,axis_yu)
        call mpp_write (io,axis_zw)
        call mpp_write (io,axis_xue)
        call mpp_write (io,axis_yun)

c-----------------------------------------------------------------------
c       allocate the himem buffer variables
c-----------------------------------------------------------------------

# ifdef diagnostics_himem
        allocate( buf_cet(1:imt,jstask:jetask,km) )
        allocate( buf_cnt(1:imt,jstask:jetask,km) )
        allocate( buf_cbt(1:imt,jstask:jetask,km) )
        allocate( buf_ceu(1:imt,jstask:jetask,km) )
        allocate( buf_cnu(1:imt,jstask:jetask,km) )
        allocate( buf_cbu(1:imt,jstask:jetask,km) )
#  if defined gm_skew || defined redi_diffusion
        allocate( buf_Ai_ez(1:imt,jstask:jetask,km) )
        allocate( buf_Ai_nz(1:imt,jstask:jetask,km) )
        allocate( buf_Ai_bx(1:imt,jstask:jetask,km) )
        allocate( buf_Ai_by(1:imt,jstask:jetask,km) )
#   if defined redi_diffusion
        allocate( buf_K11(1:imt,jstask:jetask,km) )
        allocate( buf_K22(1:imt,jstask:jetask,km) )
        allocate( buf_K33(1:imt,jstask:jetask,km) )
#   endif
#  endif
#  if defined kppvmix
        allocate( buf_ghats(1:imt,jstask:jetask,km) )
#  endif
# endif

c-----------------------------------------------------------------------
c     endif for opening the diagnostic file
c-----------------------------------------------------------------------
      endif


c-----------------------------------------------------------------------
c     return if outside mw and using lomem 
c-----------------------------------------------------------------------

# if defined diagnostics_lomem
      if( .NOT. within_window_loop )return
# endif

c-----------------------------------------------------------------------
c     point "domain" to the x-y-domain from "mom_domain" 
c-----------------------------------------------------------------------

# ifdef diagnostics_himem
!domain was declared pointer because it needs to inherit linked list
      domain => mom_domain
# else
!domain was not declared pointer
      domain = mom_domain
      domain%y%global = domain%y%compute
      domain%y%compute%start_index = js+joff
      domain%y%compute%end_index   = je+joff
      domain%y%data = domain%y%compute
# endif

c-----------------------------------------------------------------------
c     data indexed by (i,k,j) must be re-arranged to match (i,j,k). 
c-----------------------------------------------------------------------

      if( within_window_loop )then

c-----------------------------------------------------------------------
c     define variables (on data domain)
c     incoming data is already masked to 0.
c     add (1-mask)*missing_val for Ferret. 
c-----------------------------------------------------------------------

        do k=1,km
          do j=js,je
            do i=1,imt
              buf_cet(i,j+joff,k) = ce(i,k,j,2)
     &            + (1-tmask(i,k,j))*field_cet%missing
              buf_cnt(i,j+joff,k) = cn(i,k,j,2)
     &            + (1-tmask(i,k,j))*field_cnt%missing
              buf_cbt(i,j+joff,k) = cb(i,k,j,2)
     &            + (1-tmask(i,k,j))*field_cbt%missing
              buf_ceu(i,j+joff,k) = ce(i,k,j,1)
     &            + (1-umask(i,k,j))*field_ceu%missing
              buf_cnu(i,j+joff,k) = cn(i,k,j,1)
     &            + (1-umask(i,k,j))*field_cnu%missing
              buf_cbu(i,j+joff,k) = cb(i,k,j,1)
     &            + (1-umask(i,k,j))*field_cbu%missing
# if defined gm_skew || defined redi_diffusion
	      buf_Ai_ez(i,j+joff,k) = 0.25*
     &                        (Ai_ez(i,k,j,0,1) + Ai_ez(i,k,j,0,0)
     &                       + Ai_ez(i,k,j,1,1) + Ai_ez(i,k,j,1,0))
     &            + (1-tmask(i,k,j))*field_Ai_ez%missing
	      buf_Ai_nz(i,j+joff,k) = 0.25*
     &                        (Ai_nz(i,k,j,0,1) + Ai_nz(i,k,j,0,0)
     &                       + Ai_nz(i,k,j,1,1) + Ai_nz(i,k,j,1,0))
     &            + (1-tmask(i,k,j))*field_Ai_nz%missing
	      buf_Ai_bx(i,j+joff,k) = 0.25*
     &                         (Ai_bx(i,k,j,0,0) + Ai_bx(i,k,j,1,0)
     &                        + Ai_bx(i,k,j,0,1) + Ai_bx(i,k,j,1,1))
     &            + (1-tmask(i,k,j))*field_Ai_bx%missing
	      buf_Ai_by(i,j+joff,k) = 0.25*
     &                         (Ai_by(i,k,j,0,0) + Ai_by(i,k,j,1,0)
     &                        + Ai_by(i,k,j,0,1) + Ai_by(i,k,j,1,1))
     &            + (1-tmask(i,k,j))*field_Ai_by%missing
#  if defined redi_diffusion
              buf_K11(i,j+joff,k) = K11(i,k,j)
     &            + (1-tmask(i,k,j))*field_K11%missing
              buf_K22(i,j+joff,k) = K22(i,k,j)
     &            + (1-tmask(i,k,j))*field_K22%missing
              buf_K33(i,j+joff,k) = K33(i,k,j)
     &            + (1-tmask(i,k,j))*field_K33%missing
#  endif
# endif
# if defined kppvmix
              buf_ghats(i,j+joff,k) = ghats(i,k,j)
     &            + (1-tmask(i,k,j))*field_ghats%missing
# endif
            enddo
	  enddo
	enddo

      endif  !end of within_window_loop if-test

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

c     For himem, only write data when finished with window loop
c     For lomem, write data within window loop
c
# if defined diagnostics_himem
      if( .NOT. within_window_loop )then
# endif
        call mpp_write(io, field_cet, domain, buf_cet, timrec)
        call mpp_write(io, field_cnt, domain, buf_cnt, timrec)
        call mpp_write(io, field_cbt, domain, buf_cbt, timrec)
        call mpp_write(io, field_ceu, domain, buf_ceu, timrec)
        call mpp_write(io, field_cnu, domain, buf_cnu, timrec)
        call mpp_write(io, field_cbu, domain, buf_cbu, timrec)
# if defined gm_skew || defined redi_diffusion
        call mpp_write(io, field_Ai_ez, domain, buf_Ai_ez, timrec)
        call mpp_write(io, field_Ai_nz, domain, buf_Ai_nz, timrec)
        call mpp_write(io, field_Ai_bx, domain, buf_Ai_bx, timrec)
        call mpp_write(io, field_Ai_by, domain, buf_Ai_by, timrec)
#  if defined redi_diffusion
        call mpp_write(io, field_K11, domain, buf_K11, timrec)
        call mpp_write(io, field_K22, domain, buf_K22, timrec)
        call mpp_write(io, field_K33, domain, buf_K33, timrec)
#  endif
# endif
# if defined kppvmix
        call mpp_write(io, field_ghats, domain, buf_ghats, timrec)
# endif
        call mpp_write(io, field_period, period, timrec)

# if defined diagnostics_himem
      endif
# endif

      if(js+joff .eq. jscomp) then
          write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8,a,i8)')
     &'    ==>  Mixing coefficients written to file '
     &,    file_name, 'on ts = ', itt,' ', stamp
     &,'   between latitudes:',jscomp,' and ', jecomp
     &,'   and from the surface down to depth level ',km
      endif 

      return
      end
#else
      return
      end
#endif



#if defined save_convection
      subroutine conv_out (joff, js, je, within_window_loop)
c-----------------------------------------------------------------------
c     save convection data
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c     updates to mpp_io made July 1999 by 
c                  s. m. griffies        e-mail=> smg@gfdl.gov 
c-----------------------------------------------------------------------
      use memory_window
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod

# include "param.h"
# include "coord.h"
# include "derived_options.h"
# include "diag.h"
# include "mw.h"
# include "switch.h"
# include "taskrows.h"
# include "tmngr.h"
# include "vers.h"
      character*128 file_name, file_stamp, gvtitle
      logical :: within_window_loop
      dimension bufzt(km)
      data num_conv /0/
      save num_conv
      save io

# if defined diagnostics_himem
      type(domain2D), pointer :: domain
c
c     himem buffers all the memory windows
c
      real, dimension(:,:,:), save, allocatable :: buf_conv
# endif

# if defined diagnostics_lomem
      type(domain2D) :: domain
c
c     lomem buffers only within a single memory window
c
      real :: buf_conv(1:imt,js+joff:je+joff,km)
# endif

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_conv
      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_xt
      type (axistype), save :: axis_yt
      type (axistype), save :: axis_zt
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where convection diagnostic is written
c-----------------------------------------------------------------------

      file_name = file_stamp ('cvct', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, and ocean depth in meters
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Rate of Convection' // momver
      bufzt(:) = zt(:)*0.01

c-----------------------------------------------------------------------
c     check whether diagnostic file(s) should be opened
c-----------------------------------------------------------------------

      if (num_conv .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_conv = num_conv + 1

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
     &      io, axis_zt, 'zt_k', 'm', 'Depth of T grid point',
     &      cartesian='Z', sense=-1, data=bufzt)
 
        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       attributes of data 
c-----------------------------------------------------------------------

        call mpp_write_meta(
     &      io, field_conv, (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'convect', 'deg C/s', 'Rate of convection',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

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
      
        call mpp_write (io,axis_xt)
        call mpp_write (io,axis_yt)
        call mpp_write (io,axis_zt)

c-----------------------------------------------------------------------
c       allocate the himem buffer variables on data domain
c-----------------------------------------------------------------------

# ifdef diagnostics_himem
        allocate( buf_conv(1:imt,jstask:jetask,km) )
# endif

c-----------------------------------------------------------------------
c     endif for opening the diagnostic file
c-----------------------------------------------------------------------
      endif

c-----------------------------------------------------------------------
c     return if ouside mw and using lomem
c-----------------------------------------------------------------------

# ifdef diagnostics_lomem
      if( .NOT. within_window_loop ) return
# endif


c-----------------------------------------------------------------------
c     point "domain" to the x-y-domain from "mom_domain" 
c-----------------------------------------------------------------------

# ifdef diagnostics_himem
!domain was declared pointer because it needs to inherit linked list
      domain => mom_domain
# else
!domain was not declared pointer
      domain = mom_domain
      domain%y%global = domain%y%compute
      domain%y%compute%start_index = js+joff
      domain%y%compute%end_index   = je+joff
      domain%y%data = domain%y%compute
# endif

c-----------------------------------------------------------------------
c     data indexed by (i,k,j) must be re-arranged to match (i,j,k). 
c-----------------------------------------------------------------------

      if( within_window_loop )then

c-----------------------------------------------------------------------
c     define variables (on data domain)
c     incoming data is already masked to 0.
c     will add (1-mask)*missing_val for Ferret. 
c-----------------------------------------------------------------------

        do k=1,km
          do j=js,je
            do i=1,imt
              buf_conv(i,j+joff,k) = excnv1(i,k,j)
     &            + (1-tmask(i,k,j))*field_conv%missing
            enddo
          enddo
        enddo

      endif 

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

c     Data is written only outside window loop for himem
c     It is written within window loop for lomem

# if defined diagnostics_himen
      if( .NOT. within_window_loop )then
# endif

        call mpp_write(io, field_conv, domain, buf_conv, timrec)
        call mpp_write(io, field_period, period, timrec)

# if defined diagnostics_himen
      endif
# endif

        if(js+joff .eq. jscomp) then
           write (stdout,'(a,,a,a,i10,a,a/a,i8,a,i8,a,i8)')
     &'    ==>  Writing convection rate to file '
     &,    file_name, 'on ts = ', itt,' ', stamp
     &,'   between latitudes:',jscomp,' and ', jecomp
     &,'   and from the surface down to depth level ',km
         endif

      return
      end
#endif


#if defined save_pressure
      subroutine pressure_out (joff, js, je, within_window_loop)
c-----------------------------------------------------------------------
c     save pressure and pressure gradient 
c
c     author:  s. m. griffies            e-mail=> smg@gfdl.gov 
c     updates to mpp_io made July 1999
c-----------------------------------------------------------------------
      use memory_window
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "coord.h"
# include "derived_options.h"
# include "diag.h"      
# include "mw.h"
# include "switch.h"
# include "taskrows.h"
# include "tmngr.h"
# include "vers.h"
      character*128 file_name, file_stamp, gvtitle
      logical :: within_window_loop
      dimension bufzt(km)
      data num_press /0/
      save num_press
      save io
# if defined diagnostics_himem
c
c     domain pointer 
c
      type(domain2D), pointer :: domain
c
c     himem buffers all the memory windows
c
      real, dimension(:,:,:), save, allocatable :: buf_pressure
      real, dimension(:,:,:,:), save, allocatable :: buf_grad_p
# endif

# if defined diagnostics_lomem
      type(domain2D) :: domain
c
c     lomem buffers only within a single memory window
c
      real :: buf_pressure(1:imt,js+joff:je+joff,km)
      real :: buf_grad_p(1:imt,js+joff:je+joff,km,2)
# endif

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_pressure
      type(fieldtype), save :: field_grad_px
      type(fieldtype), save :: field_grad_py
      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_xt
      type (axistype), save :: axis_yt
      type (axistype), save :: axis_zt
      type (axistype), save :: axis_xu
      type (axistype), save :: axis_yu
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where pressure is written
c-----------------------------------------------------------------------

      file_name = file_stamp ('pressure', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, and ocean depth in meters
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Pressure diagnostic' // momver
      bufzt(:) = zt(:)*0.01

c-----------------------------------------------------------------------
c     check whether diagnostic file(s) should be opened
c-----------------------------------------------------------------------

      if (num_press .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_press = num_press + 1

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
     &      io, axis_xu, 'xu_i', 'degrees_E', 'Longitude of U points',
     &      cartesian='X', domain=mom_domain%x, 
     &      data=xu(ibeg:iend))

        call mpp_write_meta(
     &      io, axis_yu, 'yu_j', 'degrees_N', 'Latitude of U points',
     &      cartesian='Y', domain=mom_domain%y, 
     &      data=yu(jbeg:jend))

        call mpp_write_meta(
     &      io, axis_zt, 'zt_k', 'm', 'Depth of T grid point',
     &      cartesian='Z', sense=-1, data=bufzt)

        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       attributes of data 
c-----------------------------------------------------------------------

        call mpp_write_meta(
     &      io, field_pressure, (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'pressure', 'g/cm/s**2', 'baroclinic pressure',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_grad_px, (/axis_xu,axis_yu,axis_zt,axis_t/), 
     &      'grad_px', 'g/cm**2/s**2', 'x-baroclinic press grad',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_grad_py, (/axis_xu,axis_yu,axis_zt,axis_t/), 
     &      'grad_py', 'g/cm**2/s**2', 'y-baroclinic press grad',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
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
      
        call mpp_write (io,axis_xt)
        call mpp_write (io,axis_yt)
        call mpp_write (io,axis_zt)
        call mpp_write (io,axis_xu)
        call mpp_write (io,axis_yu)

c-----------------------------------------------------------------------
c       allocate the himem buffer variables on data domain
c-----------------------------------------------------------------------

# ifdef diagnostics_himem
        allocate( buf_pressure(1:imt,jstask:jetask,km) )
        allocate( buf_grad_p(1:imt,jstask:jetask,km,2) )
# endif

c-----------------------------------------------------------------------
c     endif for opening the diagnostic file
c-----------------------------------------------------------------------
      endif


c-----------------------------------------------------------------------
c     return if outside mw and using lomem 
c-----------------------------------------------------------------------

# ifdef diagnostics_lomem
      if( .NOT. within_window_loop )return
# endif

c-----------------------------------------------------------------------
c     point "domain" to the x-y-domain from "mom_domain" 
c-----------------------------------------------------------------------

# ifdef diagnostics_himem
!domain was declared pointer because it needs to inherit linked list
      domain => mom_domain
# else
!domain was not declared pointer
      domain = mom_domain
      domain%y%global = domain%y%compute
      domain%y%compute%start_index = js+joff
      domain%y%compute%end_index   = je+joff
      domain%y%data = domain%y%compute
# endif

c-----------------------------------------------------------------------
c     data indexed by (i,k,j) must be re-arranged to match (i,j,k). 
c-----------------------------------------------------------------------

      if( within_window_loop )then

c-----------------------------------------------------------------------
c     define variables (on data domain)
c     incoming data is already masked to 0.
c     add (1-mask)*missing_val for Ferret. 
c-----------------------------------------------------------------------

        do k=1,km
          do j=js,je
            do i=1,imt
              buf_pressure(i,j+joff,k) = pressure(i,k,j)
     &            + (1-tmask(i,k,j))*field_pressure%missing
              buf_grad_p(i,j+joff,k,1) = grad_p(i,k,j,1)
     &            + (1-umask(i,k,j))*field_grad_px%missing
              buf_grad_p(i,j+joff,k,2) = grad_p(i,k,j,2)
     &            + (1-umask(i,k,j))*field_grad_py%missing
            enddo
	  enddo
	enddo

      endif  !end of within_window_loop if-test

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

c     Data is written only outside window loop for himem
c     It is written within window loop for lomem

# if defined diagnostics_himem
      if( .NOT. within_window_loop )then
# endif

        call mpp_write(io, field_pressure, domain, buf_pressure, timrec)
        call mpp_write(io, field_grad_px, domain, buf_grad_p(:,:,:,1)
     &       , timrec)
        call mpp_write(io, field_grad_py, domain, buf_grad_p(:,:,:,2)
     &       , timrec)
        call mpp_write(io, field_period, period, timrec)

# if defined diagnostics_himem
      endif
# endif

       if(js+joff .eq. jscomp) then
          write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8,a,i8)')
     &'    ==>  Writing pressure field to file '
     &,    file_name, 'on ts = ', itt,' ', stamp
     &,'   between latitudes:',jscomp,' and ', jecomp
     &,'   and from the surface down to depth level ',km
       endif

      return
      end

#endif 

#if defined cross_flow && (defined redi_diffusion || defined gent_mcwilliams)
      subroutine cross_out (joff, js, je, within_window_loop)
c-----------------------------------------------------------------------
c     save dianeutral and isoneutral components of velocity field 
c
c     author:      r. c. pacanowski and s.m. griffies
c     e-mail=>     rcp@gfdl.gov     and smg@gfdl.gov
c
c     updates to mpp_io made July 1999 by s. m. griffies    
c-----------------------------------------------------------------------
      use memory_window
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "coord.h"
# include "derived_options.h"
# include "diag.h"
# include "grdvar.h"
# include "isopyc.h"
# include "mw.h"
# include "switch.h"
# include "state.h"
# include "taskrows.h"
# include "tmngr.h"
# include "vers.h"
      logical :: within_window_loop
      character*128 file_name, file_stamp, gvtitle
      dimension bufzt(km)
      dimension sx(imt,km), sy(imt,km), sz(imt,km)
      dimension slmag(imt,km)
      dimension ut(imt,km), vt(imt,km), wt(imt,km)
      data num_cross /0/
      save num_cross
      save io

# if defined diagnostics_himem
c
c     domain pointer 
      type(domain2D), pointer :: domain
c
c     himem buffers all the memory windows
c
      real, dimension(:,:,:,:), save, allocatable :: buf_u_dia
      real, dimension(:,:,:,:), save, allocatable :: buf_u_iso
# endif

# if defined diagnostics_lomem
      type(domain2D) :: domain
c
c     lomem buffers only within a single memory window
c
      real :: buf_u_dia(1:imt,js+joff:je+joff,km,3)
      real :: buf_u_iso(1:imt,js+joff:je+joff,km,3)
# endif

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_u_dia(3)
      type(fieldtype), save :: field_u_iso(3)
      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_xt
      type (axistype), save :: axis_yt
      type (axistype), save :: axis_zt
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where mixing coefficients are written
c-----------------------------------------------------------------------

      file_name = file_stamp ('cross', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, and ocean depth in meters
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'Diapycnal and isopycnal component of flow' // momver
      bufzt(:) = zt(:)*0.01

c-----------------------------------------------------------------------
c     check whether diagnostic file(s) should be opened
c-----------------------------------------------------------------------

      if (num_cross .eq. 0) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_cross = num_cross + 1

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
     &      io, axis_zt, 'zt_k', 'm', 'Depth of T grid point',
     &      cartesian='Z', sense=-1, data=bufzt)

        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       attributes of data 
c-----------------------------------------------------------------------

        call mpp_write_meta(
     &      io, field_u_dia(1), (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'u_dia', 'cm/s', 'zonal dianeutral velocity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_u_dia(2), (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'v_dia', 'cm/s', 'merid dianeutral velocity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_u_dia(3), (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'w_dia', 'cm/s', 'vert dianeutral velocity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_u_iso(1), (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'u_iso', 'cm/s', 'zonal isoneutral velocity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_u_iso(2), (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'v_iso', 'cm/s', 'merid isoneutral velocity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_u_iso(3), (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      'w_iso', 'cm/s', 'vert isoneutral velocity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
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
      
        call mpp_write (io,axis_xt)
        call mpp_write (io,axis_yt)
        call mpp_write (io,axis_zt)

c-----------------------------------------------------------------------
c       allocate the himem buffer variables on data domain 
c-----------------------------------------------------------------------

# ifdef diagnostics_himem
        allocate( buf_u_dia(1:imt,jstask:jetask,km,3) )
        allocate( buf_u_iso(1:imt,jstask:jetask,km,3) )
# endif

c-----------------------------------------------------------------------
c     endif for opening the diagnostic file
c-----------------------------------------------------------------------
      endif

c-----------------------------------------------------------------------
c     return if outside mw and using lomem
c-----------------------------------------------------------------------

# ifdef diagnostics_lomem
      if( .NOT. within_window_loop )return
# endif

c-----------------------------------------------------------------------
c     point "domain" to the x-y-domain from "mom_domain" 
c-----------------------------------------------------------------------

# ifdef diagnostics_himem
!domain was declared pointer because it needs to inherit linked list
      domain => mom_domain
# else
!domain was not declared pointer
      domain = mom_domain
      domain%y%global = domain%y%compute
      domain%y%compute%start_index = js+joff
      domain%y%compute%end_index   = je+joff
      domain%y%data = domain%y%compute
# endif

c-----------------------------------------------------------------------
c     data indexed by (i,k,j) must be re-arranged to match (i,j,k). 
c-----------------------------------------------------------------------

      if( within_window_loop )then

        do j=js,je
          jrow = joff + j
c
c         compute velocity at "T" points
c
          do k=1,km
            do i=2,imt-1
	      ut(i,k) = p25*(u(i,k,j,1,tau)   + u(i-1,k,j,1,tau)
     &                   + u(i,k,j-1,1,tau) + u(i-1,k,j-1,1,tau))
	      vt(i,k) = p25*(u(i,k,j,2,tau)   + u(i-1,k,j,2,tau)
     &                   + u(i,k,j-1,2,tau) + u(i-1,k,j-1,2,tau))
              wt(i,k) = p5*(adv_vbt(i,k,j)    + adv_vbt(i,k-1,j))
            enddo
          enddo
c
c         compute the dianeutral direction
c
          do k=1,km
	    km1 = max(k-1,1)
	    kp1 = min(k+1,km)
            do i=2,imt-1

  	      sx(i,k) = c0
	      sy(i,k) = c0
	      sz(i,k) = c0

              do n=1,2
  	        sx(i,k) =
     &          sx(i,k) +
     &          0.5*drho(i,k,j,n)
     &                *(
     &                 (t(i+1,k,j,n,tau)-t(i,k,j,n,tau))
     &                *cstr(jrow)*dxur(i)*tmask(i+1,k,j)*tmask(i,k,j)
     &                 + (t(i,k,j,n,tau)-t(i-1,k,j,n,tau))
     &                *cstr(jrow)*dxur(i-1)*tmask(i,k,j)*tmask(i-1,k,j)
     &                 )
 
	        sy(i,k) =
     &          sy(i,k) +
     &          0.5*drho(i,k,j,n)
     &                 *(
     &                 (t(i,k,j+1,n,tau)-t(i,k,j,n,tau))
     &                *dyur(jrow)*tmask(i,k,j+1)*tmask(i,k,j)
     &                 + (t(i,k,j,n,tau)-t(i,k,j-1,n,tau))
     &                *dyur(jrow-1)*tmask(i,k,j)*tmask(i,k,j-1)
     &                  ) 

	         sz(i,k) =
     &           sz(i,k) +
     &           0.5*drho(i,k,j,n)
     &                 *(
     &                  (t(i,km1,j+1,n,tau)-t(i,k,j,n,tau))
     &                *dzwr(km1)*tmask(i,km1,j)*tmask(i,k,j)
     &                +  (t(i,k,j+1,n,tau)-t(i,kp1,j,n,tau))
     &                *dzwr(k)*tmask(i,k,j)*tmask(i,kp1,j)
     &                  )

              enddo
            enddo
          enddo

          do k=1,km
            do i=2,imt-1
	      slmag(i,k) = 1.0/(sqrt(sx(i,k)**2+sy(i,k)**2+sz(i,k)**2)
     &                      + epsln)
            enddo
          enddo
c
          do k=1,km
            do i=2,imt-1
  	      sx(i,k) = sx(i,k)*slmag(i,k)
	      sy(i,k) = sy(i,k)*slmag(i,k)
	      sz(i,k) = sz(i,k)*slmag(i,k)
            enddo
          enddo
c
c         compute diapycnal and isopycnal velocity components
c
          do k=1,km
            do i=2,imt-1
	      d = (ut(i,k)*sx(i,k) + vt(i,k)*sy(i,k) + wt(i,k)*sz(i,k))
              buf_u_dia(i,jrow,k,1) = d*sx(i,k)
              buf_u_dia(i,jrow,k,2) = d*sy(i,k)
              buf_u_dia(i,jrow,k,3) = d*sz(i,k)
              buf_u_iso(i,jrow,k,1) = ut(i,k) - buf_u_dia(i,jrow,k,1)
	      buf_u_iso(i,jrow,k,2) = vt(i,k) - buf_u_dia(i,jrow,k,2)
	      buf_u_iso(i,jrow,k,3) = wt(i,k) - buf_u_dia(i,jrow,k,3)
            enddo
          enddo
# if defined cyclic
          buf_u_dia(1,jrow,:,:)   = buf_u_dia(imt-1,jrow,:,:)
          buf_u_dia(imt,jrow,:,:) = buf_u_dia(2,jrow,:,:)
          buf_u_iso(1,jrow,:,:)   = buf_u_iso(imt-1,jrow,:,:)
          buf_u_iso(imt,jrow,:,:) = buf_u_iso(2,jrow,:,:)
# else 
          buf_u_dia(1,jrow,:,:)   = c0
          buf_u_dia(imt,jrow,:,:) = c0
          buf_u_iso(1,jrow,:,:)   = c0
          buf_u_iso(imt,jrow,:,:) = c0
# endif

        enddo !end of j-loop  
      endif  !end of within_window_loop if-test

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

c     Data is written only outside window loop for himem
c     It is written within window loop for lomem
c
# if defined diagnostics_himem
      if( .NOT. within_window_loop )then
# endif
        do n=1,3
          call mpp_write(io, field_u_dia(n), domain, buf_u_dia(:,:,:,n)
     &       , timrec)
          call mpp_write(io, field_u_iso(n), domain, buf_u_iso(:,:,:,n)
     &       , timrec)
        enddo
        call mpp_write(io, field_period, period, timrec)

# if defined diagnostics_himem
      endif
# endif

      if(js+joff .eq. jscomp) then
          write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8,a,i8)')
     &'    ==>  Cross flow components written to file '
     &,    file_name, 'on ts = ', itt,' ', stamp
     &,'   between latitudes:',jscomp,' and ', jecomp
     &,'   and from the surface down to depth level ',km
      endif

      return
      end
#endif



#if defined save_fct
       subroutine fct_diag (joff, js, je, n, within_window_loop)
c-----------------------------------------------------------------------
c     save the change in tracer per time step due to fct 
c     (i.e. fct - centered diff) to see where fct is important
c
c     author:      Michael Eby      e-mail=> eby@uvic.ca
c
c     updates to save flux components 
c                  s. m. griffies        e-mail=> smg@gfdl.gov 
c
c     updates to mpp_io made July 1999
c                  s. m. griffies        e-mail=> smg@gfdl.gov 
c-----------------------------------------------------------------------
      use memory_window
      use mpp_domains_mod
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "coord.h" 
# include "diag.h" 
# include "mw.h" 
# include "switch.h"
# include "taskrows.h"
# include "tmngr.h"
# include "vers.h"
      logical :: within_window_loop
      character*128 file_name, file_stamp, gvtitle
      dimension bufzt(km)
      data num_fct /0/
      save num_fct
      save io

# if defined diagnostics_himem
c
c     domain pointer 
      type(domain2D), pointer :: domain
c
c     himem buffers all the memory windows
c
      real, dimension(:,:,:,:), save, allocatable :: buf_fct_minus_ctr
      real, dimension(:,:,:,:), save, allocatable :: buf_flux_diffX
      real, dimension(:,:,:,:), save, allocatable :: buf_flux_diffY
      real, dimension(:,:,:,:), save, allocatable :: buf_flux_diffZ
# endif

# if defined diagnostics_lomem
      type(domain2D) :: domain
c
c     lomem buffers only within a single memory window
c
      real :: buf_fct_minus_ctr(1:imt,js+joff:je+joff,km,2)
      real :: buf_flux_diffX(1:imt,js+joff:je+joff,km,2)
      real :: buf_flux_diffY(1:imt,js+joff:je+joff,km,2)
      real :: buf_flux_diffZ(1:imt,js+joff:je+joff,km,2)
# endif

      character*15 textfct1
      character*12 textfct2
# if defined salinity_psu
      character*3 saldim
# else 
      character*9 saldim
# endif

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
    
      type(fieldtype), save :: field_fct_minus_ctr(2)
      type(fieldtype), save :: field_flux_diffX(2)
      type(fieldtype), save :: field_flux_diffY(2)     
      type(fieldtype), save :: field_flux_diffZ(2)     
      type(fieldtype), save :: field_period 

c-----------------------------------------------------------------------
c     axis types
c-----------------------------------------------------------------------

      type (axistype), save :: axis_xt
      type (axistype), save :: axis_yt
      type (axistype), save :: axis_zt
      type (axistype), save :: axis_t

c-----------------------------------------------------------------------
c     specify file where FCT diagnostic is written
c-----------------------------------------------------------------------

      file_name = file_stamp ('fct', xstamp, '.dta')

c-----------------------------------------------------------------------
c     set averaging period, time, title, and ocean depth in meters
c-----------------------------------------------------------------------

      period = c0
      if (rununits.eq.'days') then
        timrec = prelyr*yrlen
      else
        timrec = prelyr
      endif
      gvtitle = 'FCT diagnostic' // momver
      bufzt(:) = zt(:)*0.01
# if defined salinity_psu
      saldim = 'psu' 
# else 
      saldim = 'ppm-0.035' 
# endif

c-----------------------------------------------------------------------
c     check whether diagnostic file(s) should be opened
c-----------------------------------------------------------------------

      if (num_fct .eq. 0 .and. n .eq. 1) then

        call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_fct = num_fct + 1

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
     &      io, axis_zt, 'zt_k', 'm', 'Depth of T grid point',
     &      cartesian='Z', sense=-1, data=bufzt)

        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       attributes of data 
c-----------------------------------------------------------------------

        do nn=1,2
          if(nn .eq. 1) then

           write (textfct1,'(a14,i1)') 'fct_minus_ctr_',nn
           call mpp_write_meta(
     &      io, field_fct_minus_ctr(nn)
     &          , (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      textfct1, 'deg C', textfct1,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

           write (textfct2,'(a11,i1)') 'flux_diffX_',nn
           call mpp_write_meta(
     &      io, field_flux_diffX(nn),(/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      textfct2, '(cm/s)*deg C', textfct2,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

           write (textfct2,'(a11,i1)') 'flux_diffY_',nn
           call mpp_write_meta(
     &      io, field_flux_diffY(nn),(/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      textfct2, '(cm/s)*deg C', textfct2,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

           write (textfct2,'(a11,i1)') 'flux_diffZ_',nn
           call mpp_write_meta(
     &      io, field_flux_diffZ(nn),(/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      textfct2, '(cm/s)*deg C', textfct2,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

          endif

          if(nn .eq. 2) then

           write (textfct1,'(a14,i1)') 'fct_minus_ctr_',nn
           call mpp_write_meta(
     &      io, field_fct_minus_ctr(nn)
     &          , (/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      textfct1, saldim, textfct1,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

           write (textfct2,'(a11,i1)') 'flux_diffX_',nn
           call mpp_write_meta(
     &      io, field_flux_diffX(nn),(/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      textfct2, '(cm/s)*'//saldim, textfct2,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

           write (textfct2,'(a11,i1)') 'flux_diffY_',nn
           call mpp_write_meta(
     &      io, field_flux_diffY(nn),(/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      textfct2, '(cm/s)*'//saldim, textfct2,
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

           write (textfct2,'(a11,i1)') 'flux_diffZ_',nn
           call mpp_write_meta(
     &      io, field_flux_diffZ(nn),(/axis_xt,axis_yt,axis_zt,axis_t/), 
     &      textfct2, '(cm/s)*'//saldim, textfct2,
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
c       write axes now that all meta data has been written
c-----------------------------------------------------------------------
      
        call mpp_write (io,axis_xt)
        call mpp_write (io,axis_yt)
        call mpp_write (io,axis_zt)

c-----------------------------------------------------------------------
c       allocate the himem buffer variables on data domain
c-----------------------------------------------------------------------

# ifdef diagnostics_himem
        allocate( buf_fct_minus_ctr(1:imt,jstask:jetask,km,2) )
        allocate( buf_flux_diffX(1:imt,jstask:jetask,km,2) )
        allocate( buf_flux_diffY(1:imt,jstask:jetask,km,2) )
        allocate( buf_flux_diffZ(1:imt,jstask:jetask,km,2) )
# endif

c-----------------------------------------------------------------------
c     endif for opening the diagnostic file
c-----------------------------------------------------------------------
      endif


c-----------------------------------------------------------------------
c     return if outside mw and using lomem 
c-----------------------------------------------------------------------

# ifdef diagnostics_lomem
      if( .NOT. within_window_loop )return
# endif

c-----------------------------------------------------------------------
c     point "domain" to the x-y-domain from "mom_domain" 
c-----------------------------------------------------------------------

# ifdef diagnostics_himem
!domain was declared pointer because it needs to inherit linked list
      domain => mom_domain
# else
!domain was not declared pointer
      domain = mom_domain
      domain%y%global = domain%y%compute
      domain%y%compute%start_index = js+joff
      domain%y%compute%end_index   = je+joff
      domain%y%data = domain%y%compute
# endif

c-----------------------------------------------------------------------
c     data indexed by (i,k,j) must be re-arranged to match (i,j,k). 
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     define variables (on data domain)
c     incoming data is already masked to 0.
c     add (1-mask)*missing_val for Ferret. 
c-----------------------------------------------------------------------

      if( within_window_loop )then

        do k=1,km
          do j=js,je
            jrow=j+joff
            do i=1,imt
              buf_fct_minus_ctr(i,jrow,k,n) = fct_minus_ctr(i,k,j)
     &            + (1-tmask(i,k,j))*field_fct_minus_ctr(n)%missing
              buf_flux_diffX(i,jrow,k,n) = flux_diffX(i,k,j)
     &            + (1-tmask(i,k,j))*field_flux_diffX(n)%missing
              buf_flux_diffY(i,jrow,k,n) = flux_diffY(i,k,j)
     &            + (1-tmask(i,k,j))*field_flux_diffY(n)%missing
              buf_flux_diffZ(i,jrow,k,n) = flux_diffZ(i,k,j)
     &            + (1-tmask(i,k,j))*field_flux_diffZ(n)%missing
            enddo
	  enddo
	enddo

      endif  !end of within_window_loop if-test

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

c     For himem, only write data when finished with window loop
c     For lomem, write data within window loop

# if defined diagnostics_himem
      if( .NOT. within_window_loop )then
# endif
        call mpp_write(io,field_fct_minus_ctr(n),domain
     &       , buf_fct_minus_ctr(:,:,:,n),timrec)
        call mpp_write(io,field_flux_diffX(n),domain
     &       , buf_flux_diffX(:,:,:,n),timrec)
        call mpp_write(io,field_flux_diffY(n),domain
     &       , buf_flux_diffY(:,:,:,n),timrec)
        call mpp_write(io,field_flux_diffZ(n),domain
     &       , buf_flux_diffZ(:,:,:,n),timrec)

        call mpp_write(io, field_period, period, timrec)

# if defined diagnostics_himem
      endif
# endif

        if(js+joff .eq. jscomp) then
          write (stdout,'(a,a,a,i10,a,a/a,i8,a,i8,a,i8)')
     &'  ==>  FCT diagnostics written to file '
     &,  file_name, 'on ts = ', itt,' ', stamp
     &,' between latitudes:',jscomp,' and ', jecomp
     &,' and from the surface down to depth level ',km
        endif

      return
      end

#endif



#if defined save_xbts
      subroutine xbt_out (fname, reltim, period, dstamp)
c
c-----------------------------------------------------------------------
c     save xbt data
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c     updates to mpp_io made August 1999 by 
c                  s. m. griffies        e-mail=> smg@gfdl.gov
c-----------------------------------------------------------------------
c
      use mom3_mpp_mod
      use mpp_io_mod
#  include "param.h"
#  include "coord.h"
#  include "diag.h"      
#  include "tmngr.h"
#  include "diag_xbt.h"
#  include "vers.h"
      dimension buf2(maxxbt*kmxbt), buf(maxxbt)
      dimension bufzt(km), bufzw(km)
      character*128 gvtitle
      character*(*) fname, dstamp
      character*30 tname
      data num_xbt /0/
      save num_xbt
      save io
# if defined salinity_psu
      character*3 saldim
# else 
      character*9 saldim
# endif

      type (axistype), save :: axis_n     
      type (axistype), save :: axis_zt
      type (axistype), save :: axis_zw
      type (axistype), save :: axis_t

      type(fieldtype), save :: field_lons      
      type(fieldtype), save :: field_lats
      type(fieldtype), save :: field_uxbtsf(2)
      type(fieldtype), save :: field_txbtsf(nt)
      type(fieldtype), save :: field_uxbt(nuxbt)
      type(fieldtype), save :: field_vxbt(nuxbt)
      type(fieldtype), save :: field_txbt(nt,ntxbt)
      type(fieldtype), save :: field_xbtw
      type(fieldtype), save :: field_period
      
      timrec   = reltim - 0.5*period/yrlen      
      gvtitle  = 'XBTs' // momver
      bufzt(:) = 0.01*zt(:)
      bufzw(:) = 0.01*zw(:)
      buf(:)   = 0.0
      buf2(:)  = 0.0      
# if defined salinity_psu
      saldim = 'psu' 
# else 
      saldim = 'ppm-0.035' 
# endif

      
      if (num_xbt .eq. 0) then

        call mpp_open(io, fname, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

        num_xbt = num_xbt + 1
        
c-----------------------------------------------------------------------
c       axes
c-----------------------------------------------------------------------

        do n=1,numxbt
          buf(n) = n
        enddo
        call mpp_write_meta(
     &      io, axis_n, 'N', 'station', 'XBT station',
     &      cartesian='X',  data=buf(1:numxbt))

        call mpp_write_meta(
     &      io, axis_zt, 'zt_k', 'm', 'Depth of T points',
     &      cartesian='Z', sense=-1, data=bufzt(1:kmxbt))

        call mpp_write_meta(
     &      io, axis_zw, 'zw_k', 'm', 'Depth of W points',
     &      cartesian='Z', sense=-1, data=bufzw(1:kmxbt))
 
        call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c       attributes of data 
c-----------------------------------------------------------------------

c
c       xbt positions 
c
        call mpp_write_meta(
     &      io, field_lons, (/axis_n,axis_t/), 
     &      'lons', 'degrees_E', 'XBT longitudes',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_lats, (/axis_n,axis_t/), 
     &      'lats', 'degrees_N', 'XBT latitudes',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c
c       surface fluxes 
c
        call mpp_write_meta(
     &      io, field_uxbtsf(1), (/axis_n,axis_t/), 
     &      'taux', 'dyn/cm**2', 'zonal wind stress',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbtsf(2), (/axis_n,axis_t/), 
     &      'tauy', 'dyn/cm**2', 'merid wind stress',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

        call mpp_write_meta(
     &      io, field_txbtsf(1), (/axis_n,axis_t/), 
     &      'hflx', 'cal/cm**2/s', 'heat flux',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbtsf(2), (/axis_n,axis_t/), 
     &      'sflx', 'g/cm**2/s', 'salt flux',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        do n=3,nt
	  write (tname,'(a,i1,a)') 'tracer_',n,'_flux'
          call mpp_write_meta(
     &    io, field_txbtsf(n), (/axis_n,axis_t/), 
     &    tname, 'tracer flux units', tname,
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34) 
	enddo
c
c       zonal velocity terms 
c
        call mpp_write_meta(
     &      io, field_uxbt(1), (/axis_n,axis_zt,axis_t/), 
     &      'dudt', 'cm/s**2', 'du/dt',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(2), (/axis_n,axis_zt,axis_t/), 
     &      'neg_px', 'cm/s**2', '-zonal pressure gradient',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(3), (/axis_n,axis_zt,axis_t/), 
     &      'neg_uux', 'cm/s**2', '-d(uu)/dx',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(4), (/axis_n,axis_zt,axis_t/), 
     &      'neg_vuy', 'cm/s**2', '-d(vu)/dy',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(5), (/axis_n,axis_zt,axis_t/), 
     &      'neg_wuz', 'cm/s**2', '-d(wu)/dz',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(6), (/axis_n,axis_zt,axis_t/), 
     &      'zon_diff_u', 'cm/s**2', 'zonal diffusion of u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(7), (/axis_n,axis_zt,axis_t/), 
     &      'mer_diff_u', 'cm/s**2', 'merid diffusion of u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(8), (/axis_n,axis_zt,axis_t/), 
     &      'ver_diff_u', 'cm/s**2', 'vert diffusion of u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(9), (/axis_n,axis_zt,axis_t/), 
     &      'metric_diff_u', 'cm/s**2', 'metric diffusion of u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(10), (/axis_n,axis_zt,axis_t/), 
     &      'fv', 'cm/s**2', 'u coriolis term',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(11), (/axis_n,axis_zt,axis_t/), 
     &      'source_u', 'cm/s**2', 'u source term',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(12), (/axis_n,axis_zt,axis_t/), 
     &      'neg_surfpx', 'cm/s**2', '-zonal surface pressure grad',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(13), (/axis_n,axis_zt,axis_t/), 
     &      'adv_metric_u', 'cm/s**2', 'advection metric for u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(14), (/axis_n,axis_zt,axis_t/), 
     &      'neg_zonadv_u', 'cm/s**2', '-zonal advection of u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(15), (/axis_n,axis_zt,axis_t/), 
     &      'neg_meradv_u', 'cm/s**2', '-merid advection of u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(16), (/axis_n,axis_zt,axis_t/), 
     &      'neg_veradv_u', 'cm/s**2', '-vert advection of u',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_uxbt(17), (/axis_n,axis_zt,axis_t/), 
     &      'zon_vel', 'cm/s', 'u cell zonal velocity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c
c       meridional velocity terms 
c
        call mpp_write_meta(
     &      io, field_vxbt(1), (/axis_n,axis_zt,axis_t/), 
     &      'dvdt', 'cm/s**2', 'dv/dt',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(2), (/axis_n,axis_zt,axis_t/), 
     &      'neg_py', 'cm/s**2', '-merid pressure gradient',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(3), (/axis_n,axis_zt,axis_t/), 
     &      'neg_uvx', 'cm/s**2', '-d(uv)/dx',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(4), (/axis_n,axis_zt,axis_t/), 
     &      'neg_vvy', 'cm/s**2', '-d(vv)/dy',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(5), (/axis_n,axis_zt,axis_t/), 
     &      'neg_wvz', 'cm/s**2', '-d(wv)/dz',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(6), (/axis_n,axis_zt,axis_t/), 
     &      'zon_diff_v', 'cm/s**2', 'zonal diffusion of v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(7), (/axis_n,axis_zt,axis_t/), 
     &      'mer_diff_v', 'cm/s**2', 'merid diffusion of v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(8), (/axis_n,axis_zt,axis_t/), 
     &      'ver_diff_v', 'cm/s**2', 'vert diffusion of v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(9), (/axis_n,axis_zt,axis_t/), 
     &      'metric_diff_v', 'cm/s**2', 'metric diffusion of v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(10), (/axis_n,axis_zt,axis_t/), 
     &      'neg_fu', 'cm/s**2', '-v coriolis term',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(11), (/axis_n,axis_zt,axis_t/), 
     &      'source_v', 'cm/s**2', 'v source term',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(12), (/axis_n,axis_zt,axis_t/), 
     &      'neg_surfpy', 'cm/s**2', '-merid surface pressure grad',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(13), (/axis_n,axis_zt,axis_t/), 
     &      'adv_metric_v', 'cm/s**2', 'advection metric for v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(14), (/axis_n,axis_zt,axis_t/), 
     &      'neg_zonadv_v', 'cm/s**2', '-zonal advection of v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(15), (/axis_n,axis_zt,axis_t/), 
     &      'neg_meradv_v', 'cm/s**2', '-merid advection of v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(16), (/axis_n,axis_zt,axis_t/), 
     &      'neg_veradv_v', 'cm/s**2', '-vert advection of v',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_vxbt(17), (/axis_n,axis_zt,axis_t/), 
     &      'mer_vel', 'cm/s', 'u cell meridional velocity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c
c       vertical velocity
c
        call mpp_write_meta(
     &      io, field_xbtw, (/axis_n,axis_zw,axis_t/), 
     &      'adv_vbu', 'cm/s', 'u cell vertical velocity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c
c       temperature terms 
c
        call mpp_write_meta(
     &      io, field_txbt(1,1), (/axis_n,axis_zt,axis_t/), 
     &      'dTdt', 'deg C/s', 'd(Temperature)/dt',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(1,2), (/axis_n,axis_zt,axis_t/), 
     &      'neg_utx', 'deg C/s', '-d(u*Temperature)/dx',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(1,3), (/axis_n,axis_zt,axis_t/), 
     &      'neg_vty', 'deg C/s', '-d(v*Temperature)/dy',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(1,4), (/axis_n,axis_zt,axis_t/), 
     &      'neg_wtz', 'deg C/s', '-d(w*Temperature)/dz',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(1,5), (/axis_n,axis_zt,axis_t/), 
     &      'diff_tx', 'deg C/s', 'zonal diffusion of temperature',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(1,6), (/axis_n,axis_zt,axis_t/), 
     &      'diff_ty', 'deg C/s', 'merid diffusion of temperature',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(1,7), (/axis_n,axis_zt,axis_t/), 
     &      'diff_tz', 'deg C/s', 'vert diffusion of temperature',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(1,8), (/axis_n,axis_zt,axis_t/), 
     &      'source_T', 'deg C/s', 'temperature chg due to sources',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(1,9), (/axis_n,axis_zt,axis_t/), 
     &    'convect_T', 'deg C/s', 'temperature chg due to convection',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(1,10), (/axis_n,axis_zt,axis_t/), 
     &      'filter_T', 'deg C/s', 'temperature chg due to filtering',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(1,11), (/axis_n,axis_zt,axis_t/), 
     &      'neg_zonadv_T', 'deg C/s', '-zonal advection of temp',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(1,12), (/axis_n,axis_zt,axis_t/), 
     &      'neg_meradv_T', 'deg C/s', '-merid advection of temp',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)       
        call mpp_write_meta(
     &      io, field_txbt(1,13), (/axis_n,axis_zt,axis_t/), 
     &      'neg_veradv_T', 'deg C/s', '-vert advection of temp',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)  
        call mpp_write_meta(
     &      io, field_txbt(1,14), (/axis_n,axis_zt,axis_t/), 
     &      'var_T', 'deg C**2', 'variance of temperature',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)  
        call mpp_write_meta(
     &      io, field_txbt(1,15), (/axis_n,axis_zt,axis_t/), 
     &      'temp', 'deg C', 'temperature',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c
c       salinity terms 
c
        call mpp_write_meta(
     &      io, field_txbt(2,1), (/axis_n,axis_zt,axis_t/), 
     &      'dSdt', saldim//'/s', 'd(Salinity)/dt',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(2,2), (/axis_n,axis_zt,axis_t/), 
     &      'neg_usx', saldim//'/s', '-d(u*Salinity)/dx',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(2,3), (/axis_n,axis_zt,axis_t/), 
     &      'neg_vsy', saldim//'/s', '-d(v*Salinity)/dy',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(2,4), (/axis_n,axis_zt,axis_t/), 
     &      'neg_wsz', saldim//'/s', '-d(w*Salinity)/dz',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(2,5), (/axis_n,axis_zt,axis_t/), 
     &      'diff_sx', saldim//'/s', 'zonal diffusion of salinity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(2,6), (/axis_n,axis_zt,axis_t/), 
     &      'diff_sy', saldim//'/s', 'merid diffusion of salinity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(2,7), (/axis_n,axis_zt,axis_t/), 
     &      'diff_sz', saldim//'/s', 'vert diffusion of salinity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(2,8), (/axis_n,axis_zt,axis_t/), 
     &      'source_S', saldim//'/s', 'salinity chg due to sources',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(2,9), (/axis_n,axis_zt,axis_t/), 
     &    'convect_S', saldim//'/s', 'salinity chg due to convection',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(2,10), (/axis_n,axis_zt,axis_t/), 
     &      'filter_S', saldim//'/s', 'salinity chg due to filtering',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(2,11), (/axis_n,axis_zt,axis_t/), 
     &   'neg_zonadv_S', saldim//'/s', '-zonal advection of salinity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
        call mpp_write_meta(
     &      io, field_txbt(2,12), (/axis_n,axis_zt,axis_t/), 
     &   'neg_meradv_S', saldim//'/s', '-merid advection of salinity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)       
        call mpp_write_meta(
     &      io, field_txbt(2,13), (/axis_n,axis_zt,axis_t/), 
     &   'neg_veradv_S', saldim//'/s', '-vert advection of salinity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)  
        call mpp_write_meta(
     &      io, field_txbt(2,14), (/axis_n,axis_zt,axis_t/), 
     &      'var_S', saldim//'**2', 'variance of salinity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)  
        call mpp_write_meta(
     &      io, field_txbt(2,15), (/axis_n,axis_zt,axis_t/), 
     &      'salinity', saldim, 'salinity',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)
c
c       passive tracer terms 
c
        do n=3,nt
	  write (tname,'(a,i1)') '_',n

          call mpp_write_meta(
     &    io, field_txbt(n,1), (/axis_n,axis_zt,axis_t/), 
     &    'dPdt'//tname, 'tracer units/s',
     &    'd(passive traver)/dt'//tname,
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &    io, field_txbt(n,2), (/axis_n,axis_zt,axis_t/), 
     &    'neg_uPx'//tname, 'tracer units/s',
     &    '-d(u*P)/dx'//tname,
     &     -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &    io, field_txbt(n,3), (/axis_n,axis_zt,axis_t/), 
     &    'neg_vPy'//tname, 'tracer units/s',
     &    '-d(v*P)/dy'//tname,
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &    io, field_txbt(n,4), (/axis_n,axis_zt,axis_t/), 
     &    'neg_wPz'//tname, 'tracer units/s',
     &    '-d(w*P)/dz'//tname,
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &    io, field_txbt(n,5), (/axis_n,axis_zt,axis_t/), 
     &    'diff_Px'//tname, 'tracer units/s',
     &    'zonal diffusion of tracer'//tname,
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &    io, field_txbt(n,6), (/axis_n,axis_zt,axis_t/), 
     &    'diff_Py'//tname, 'tracer units/s',
     &    'merid diffusion of tracer'//tname,
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &    io, field_txbt(n,7), (/axis_n,axis_zt,axis_t/), 
     &    'diff_Pz'//tname, 'tracer units/s',
     &    'vert diffusion of tracer'//tname,
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &    io, field_txbt(n,8), (/axis_n,axis_zt,axis_t/), 
     &    'source_P'//tname, 'tracer units/s',
     &    'tracer chg due to sources'//tname,
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &    io, field_txbt(n,9), (/axis_n,axis_zt,axis_t/), 
     &    'convect_P'//tname, 'tracer units/s',
     &    'tracer chg due to convection'//tname,
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &    io, field_txbt(n,10), (/axis_n,axis_zt,axis_t/), 
     &    'filter_P'//tname, 'tracer units/s',
     &    'tracer chg due to filtering'//tname,
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &    io, field_txbt(n,11), (/axis_n,axis_zt,axis_t/), 
     &    'neg_zonadv_P'//tname, 'tracer units/s',
     &    '-zonal advection of tracer'//tname,
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &    io, field_txbt(n,12), (/axis_n,axis_zt,axis_t/), 
     &    'neg_meradv_P'//tname, 'tracer units/s',
     &    '-merid advection of tracer'//tname,
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)       
          call mpp_write_meta(
     &    io, field_txbt(n,13), (/axis_n,axis_zt,axis_t/), 
     &    'neg_veradv_P'//tname, 'tracer units/s',
     &    '-vert advection of tracer'//tname,
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)  
          call mpp_write_meta(
     &    io, field_txbt(n,14), (/axis_n,axis_zt,axis_t/), 
     &    'var_P'//tname, 'tracer units**2',
     &    'variance of tracer'//tname,
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)  
          call mpp_write_meta(
     &    io, field_txbt(n,15), (/axis_n,axis_zt,axis_t/), 
     &    'tracer'//tname, 'tracer units', 'tracer'//tname,
     &    -1.e6, 1.e6, -1.0E+34, +1.0E+34)

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
      
        call mpp_write (io,axis_n)
        call mpp_write (io,axis_zt)
        call mpp_write (io,axis_zw)

c-----------------------------------------------------------------------
c     endif for opening the diagnostic file
c-----------------------------------------------------------------------
      endif


c-----------------------------------------------------------------------
c     write out the data 
c-----------------------------------------------------------------------
c
c     station coordinates
c
      do n=1,numxbt
        buf(n) = xt(ixbt(n))
      enddo
      call mpp_write(io, field_lons, buf, timrec)
      do n=1,numxbt
        buf(n) = yt(jxbt(n))
      enddo
      call mpp_write(io, field_lats, buf, timrec)
c
c     surface fluxes
c
      do i=1,2
        do n=1,numxbt
          buf(n) = uxbtsf(i,n)
        enddo
        call mpp_write(io, field_uxbtsf(i), buf, timrec)
      enddo
      do i=1,nt
        do n=1,numxbt
          buf(n) = txbtsf(i,n)
        enddo
        call mpp_write(io, field_txbtsf(i), buf, timrec)
      enddo
c
c     terms for u equation
c
      do m=1,nuxbt      
        index = 0
        do k=1,kmxbt
          do n=1,numxbt
            index = index + 1
	    buf2(index) = uxbt(k,m,1,n)
          enddo
        enddo
        call mpp_write(io, field_uxbt(m), buf2, timrec)
      enddo
c
c     terms for v equation
c
      do m=1,nuxbt      
        index = 0
        do k=1,kmxbt
          do n=1,numxbt
            index = index + 1
	    buf2(index) = uxbt(k,m,2,n)
          enddo
        enddo
        call mpp_write(io, field_vxbt(m), buf2, timrec)
      enddo
c
c     vertical velocity on u-cell
c
      index = 0
      do k=1,kmxbt
        do n=1,numxbt
          index = index + 1
	  buf2(index) = xbtw(k,n)
        enddo
      enddo
      call mpp_write(io, field_xbtw, buf2, timrec)
c
c     terms for tracer equation
c
      do i=1,nt
        do m=1,ntxbt      
          index = 0
          do k=1,kmxbt
            do n=1,numxbt
              index = index + 1
	      buf2(index) = txbt(k,m,i,n)
            enddo
          enddo
          call mpp_write(io, field_txbt(i,m), buf2, timrec)
        enddo
      enddo
c
c     averaging period 
c
      call mpp_write(io, field_period, period, timrec)

      write (stdout,'(a,i6,a,a,a,i10,a,a)')
     &   ' => ',numxbt,' XBT stations written to file '
     &,     fname, ' on ts=',itt, ' ',dstamp

      return
      end
#endif

