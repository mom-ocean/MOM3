#ifdef makeecmwf
# define driver_only
# include "iomngr.F"
# include "grids.F"
# include "size_check.F"
# include "util.F"
      program ecmwf
c
c-----------------------------------------------------------------------
c     ecmwf prepares surface boundary conditions using ...
c     Trenberths ECMWF wind stress monthly climatology
c     for MOM
c
c     author:      b. l. samuels      e-mail=> bls@gfdl.gov
c-----------------------------------------------------------------------
c
      use mpp_io_mod
      use mom3_mpp_mod
# include "size.h"
      character*60 name
# include "stdunits.h"
# include "coord.h"
c
c     intrph  = 1 when model resolution is coarser than ECMWF
c               2 when model resolution is finer than ECMWF
c
      data  name/' MOM boundary conditions '/
      data intrph /1/
      namelist /input/ name, intrph
c
# ifdef cyclic
      write (stdout,*) '  ==> cyclic option set'
# else
      write (stdout,*) '  ==> cyclic option NOT set'
# endif
#if defined cray_t90 && !defined f90
      write (stdout,*) 
     & '=>Error: option "f90" must be enabled with "cray_t90" in "ecmwf"'
      if (imt .gt. 0) stop
#endif
      read  (stdin,input)
      write (stdout,input)
c
      if (intrph .eq. 1) then
        write (stdout,*)
     &   ' Assuming model resolution is coarser than ECMWF'
      elseif (intrph .eq. 2) then
        write (stdout,*)
     &   ' Assuming model resolution is finer than ECMWF'
      else
        write (stdout,*) ' intrph must be 1 or 2'
	stop
      endif
c
c-----------------------------------------------------------------------
c     initialize the units
c-----------------------------------------------------------------------
c
      call ioinit
      call mpp_io_init
c     
c-----------------------------------------------------------------------
c     construct the grid.
c-----------------------------------------------------------------------
c
      write (stdout,'(//1x,a,//)') 'Construct the MOM grid'
c
      maxlen = max(imt,jmt,km)
      call gcoord (maxlen, imt2, jmt2, km2, dxtdeg, dytdeg, dxudeg
     &,            dyudeg, dzt, dzw, xt, xu, yt, yu, zt, zw)
c
      call size_check (imt2, jmt2, km2, 'ecmwf', 'stop')
c
c-----------------------------------------------------------------------
c     interpolate trenberth`s monthly wind stress to MOM "u,v" grid
c-----------------------------------------------------------------------
c
      call stress (intrph, xu, yu)
c
# ifdef save_ECMWF
c
c-----------------------------------------------------------------------
c     save data for visualization 
c-----------------------------------------------------------------------
c
      call ecmwf_data
# endif
      call mpp_io_exit
      
      stop
      end 


      subroutine stress (intrp, xu, yu)
c
c-----------------------------------------------------------------------
c     linearly interpolate Trenberth's ECMWF monthly wind stress
c     climatology to the MOM grid
c-----------------------------------------------------------------------
c
# include "size.h"
      parameter (ib=146, jb=73, lenw=imt*jmt)
# include "stdunits.h"
      character*32 stamps
      character*80 iotext
      dimension xu(imt), yu(jmt)
      dimension stamps(13), dpm(13)
      dimension at(ib), cxt(ib)
      common/tds/ d(imt,jmt), a(ib,jb), b(ib,jb)
      common /tau/ cx(ib), cy(jb), clmtau(ib,jb,13,2)
      common /tau1/ work(lenw)
c
c-----------------------------------------------------------------------
c     read in Trenberth`s ECMWF monthly wind stress climatology
c     (the 13th record is the annual mean. k=(1,2) = (taux,tauy))
c-----------------------------------------------------------------------
c
      call getunit (nutauc, 'ecmwf.tau'
     &,             'unformatted sequential rewind ieee')
      do n=1,13
        read (nutauc)
        read (nutauc) stamps(n), dpm(n), lons, lats, nn
     &, cx, cy, (((clmtau(i,j,n,k),i=1,ib),j=1,jb),k=1,2)
      enddo        
c
c-----------------------------------------------------------------------
c     set up the units
c     lrec is in bytes for the cray (words*8)
c-----------------------------------------------------------------------
c
      call getunit (nutaux, 'ecmwf_taux.mom'
     &,             'unformatted sequential rewind ieee')
      call getunit (nutauy, 'ecmwf_tauy.mom'
     &,             'unformatted sequential rewind ieee')
c
      mmm    = 0
      do n=1,13
        do m=1,2
	  mmm = mmm + 1
c
c-----------------------------------------------------------------------
c         translate data grid so it entirely encloses the MOM "u,v"
c         grid. interpolate stress data to model "u,v" grid.
c         set last row (=jmt) to zero stress
c-----------------------------------------------------------------------
c
# ifdef cyclic
          call tranlon (clmtau(1,1,n,m), ib, ib-2, jb, at
     &,                 cx, xu(2), imt-2, cxt)
c
          if (intrp .eq. 1) then
            call ftc (clmtau(1,1,n,m), ib, jb, cxt, cy, d, imt, jmt, 2
     &,         imt-1, 1, jmt-1, xu, yu, mmm, work, lenw)
          else if (intrp .eq. 2) then
            call ctf (clmtau(1,1,n,m), ib, jb, cxt, cy, d, imt, jmt, 2
     &,         imt-1, 1, jmt-1, xu, yu, mmm, work, lenw)
          else
            write (stdout,*) ' Error: intrp must = 1 or 2'
	    stop
          endif
          do j=1,jmt
            d(1,j)   = d(imt-1,j)
            d(imt,j) = d(2,j)
          enddo
# else
          call tranlon (clmtau(1,1,n,m), ib, ib-2, jb, at
     &,                cx, xu, imt, cxt)
c
          if (intrp .eq. 1) then
            call ftc (clmtau(1,1,n,m), ib, jb, cxt, cy, d, imt, jmt, 1
     &,         imt, 1, jmt-1, xu, yu, mmm, work, lenw)
          else if (intrp .eq. 2) then
            call ctf (clmtau(1,1,n,m), ib, jb, cxt, cy, d, imt, jmt, 1
     &,         imt, 1, jmt-1, xu, yu, mmm, work, lenw)
          else
            write (stdout,*) ' Error: intrp must = 1 or 2'
	    stop
          endif
# endif
c
          if (n .eq. 1 .and. m .eq. 1) then
            write(stdout,'(///,1x, a38,/)') 
     &       'ECMWF wind stress grid definition:'
            write (stdout,'(//,1x,a30,/)') '            latitudes:'
            write (stdout,'(1x,10f10.5)') cy
            write (stdout,'(//,1x,a30,/)') '  original longitudes:'
            write (stdout,'(1x,10f10.5)') cx
            write (stdout,'(//,1x,a30,/)') 'translated longitudes:'
            write (stdout,'(1x,10f10.5)') cxt
          endif
c
          do i=1,imt
            d(i,jmt) = 0.0
          enddo
c
          if (m .eq. 1) then
         iotext='dim x(i),y(j),tx(i,j);read(nu)stamp,avg,i,j,mon,x,y,tx'
            write (stdout,*) iotext
            write (nutaux) iotext
            write (nutaux) stamps(n), dpm(n), imt, jmt, n, xu, yu, d
            write (stdout,*) ' ==> taux done ',stamps(n),'  ', dpm(n)
     &,     ' day average'
            call scope (d, imt, imt, jmt, 'taux')
          else
         iotext='dim x(i),y(j),ty(i,j);read(nu)stamp,avg,i,j,mon,x,y,ty'
            write (stdout,*) iotext
            write (nutauy) iotext
            write (nutauy) stamps(n), dpm(n), imt, jmt, n, xu, yu, d
            write (stdout,*) ' ==> tauy done ',stamps(n),'  ', dpm(n)
     &,     ' day average'
            call scope (d, imt, imt, jmt, 'tauy')
          endif
        enddo
      enddo
      call relunit (nutaux)
      call relunit (nutauy)
      return
      end


# ifdef save_ECMWF
      subroutine ecmwf_data
c
c-----------------------------------------------------------------------
c     Saves an output file with the ECMWF winds
c
c     author: SM Griffies  e-mail==>smg@gfdl.gov
c-----------------------------------------------------------------------
c
      use mom3_mpp_mod
      use mpp_io_mod
#  include "size.h"
#  include "coord.h"
#  include "diag.h"
#  include "stdunits.h"
#  include "tmngr.h"
#  include "vers.h"
      character*120 file_name, gvtitle
      dimension buf(imt,jmt), xb(imt), yb(jmt)
      character iotext*80
      character*6 taux_name,tauy_name
      dimension taux_name(13),tauy_name(13)
      character*30 taux_title,tauy_title
      dimension taux_title(13),tauy_title(13)

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_sbc(13,2)
      type (axistype), save :: axis_xu
      type (axistype), save :: axis_yu

      data taux_name/
     &               'taux1', 'taux2', 'taux3'
     &,              'taux4', 'taux5', 'taux6'
     &,              'taux7', 'taux8', 'taux9'
     &,              'taux10', 'taux11', 'taux12','taux13'/      
      data tauy_name/
     &               'tauy1', 'tauy2', 'tauy3'
     &,              'tauy4', 'tauy5', 'tauy6'
     &,              'tauy7', 'tauy8', 'tauy9'
     &,              'tauy10', 'tauy11', 'tauy12','tauy13'/      

      data taux_title/
     &               'ECMWF jan zonal windstress'
     &,              'ECMWF feb zonal windstress'
     &,              'ECMWF mar zonal windstress'
     &,              'ECMWF apr zonal windstress'
     &,              'ECMWF may zonal windstress'
     &,              'ECMWF jun zonal windstress'
     &,              'ECMWF jul zonal windstress'
     &,              'ECMWF aug zonal windstress'
     &,              'ECMWF sep zonal windstress'
     &,              'ECMWF oct zonal windstress'
     &,              'ECMWF nov zonal windstress'
     &,              'ECMWF dec zonal windstress'
     &,              'ECMWF ann zonal windstress'/

      data tauy_title/
     &               'ECMWF jan merid windstress'
     &,              'ECMWF feb merid windstress'
     &,              'ECMWF mar merid windstress'
     &,              'ECMWF apr merid windstress'
     &,              'ECMWF may merid windstress'
     &,              'ECMWF jun merid windstress'
     &,              'ECMWF jul merid windstress'
     &,              'ECMWF aug merid windstress'
     &,              'ECMWF sep merid windstress'
     &,              'ECMWF oct merid windstress'
     &,              'ECMWF nov merid windstress'
     &,              'ECMWF dec merid windstress'
     &,              'ECMWF ann merid windstress'/
      
      
      file_name = 'ecmwf.dta'
      gvtitle   = 'ECMWF wind stress ' // momver

      call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

c-----------------------------------------------------------------------
c     axes 
c-----------------------------------------------------------------------

      ibeg = 1
      iend = imt
      jbeg = 1
      jend = jmt

      call mpp_write_meta(
     &       io, axis_xu, 'xu_i', 'degrees_E', 'Longitude of U points',
     &       cartesian='X',  data=xu(ibeg:iend))

      call mpp_write_meta(
     &         io, axis_yu, 'yu_j', 'degrees_N', 'Latitude of U points',
     &         cartesian='Y', data=yu(jbeg:jend))

      do n=1,13
        call mpp_write_meta( 
     &         io, field_sbc(n,1), (/axis_xu,axis_yu/), 
     &        taux_name(n), 'dynes/cm**2', taux_title(n),  
     &        -1.0e6, 1.e6, -1.0E+34, +1.0E+34)
      enddo
      
      do n=1,13
        call mpp_write_meta( 
     &         io, field_sbc(n,2), (/axis_xu,axis_yu/), 
     &        tauy_name(n), 'dynes/cm**2', taux_title(n),  
     &        -1.0e6, 1.e6, -1.0E+34, +1.0E+34)
      enddo
c
c     figure title 
c
      call mpp_write_meta(io, 'title', cval=gvtitle )

c-----------------------------------------------------------------------
c     write axes now that all meta information has been written
c-----------------------------------------------------------------------
      
      call mpp_write (io,axis_xu)
      call mpp_write (io,axis_yu)

c-----------------------------------------------------------------------
c     get data from .mom files and then mpp_write
c-----------------------------------------------------------------------
      
      call getunit (nutaux, 'ecmwf_taux.mom'
     &,             'unformatted sequential rewind ieee')
      call getunit (nutauy, 'ecmwf_tauy.mom'
     &,             'unformatted sequential rewind ieee')
      
      rewind nutaux
      rewind nutauy

      do n=1,13
        buf(:,:) = field_sbc(n,1)%missing
        read (nutaux) iotext
        read (nutaux) stamp, dpm, iimt, jjmt, nn, xb, yb, buf
        call mpp_write(io, field_sbc(n,1), buf)
      enddo 
        
      do n=1,13
        buf(:,:) = field_sbc(n,2)%missing        
        read (nutauy) iotext
        read (nutauy) stamp, dpm, iimt, jjmt, nn, xb, yb, buf
        call mpp_write(io, field_sbc(n,2), buf)
      enddo

      write (stdout,'(a,a)')
     &   ' => ECMWF data written in output file ', file_name
      return
      end
# endif

#else
      subroutine ecmwf
      return
      end
#endif
