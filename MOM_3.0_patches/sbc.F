#ifdef makesbc
# define driver_only
# include "iomngr.F"
# include "grids.F"
# include "size_check.F"
# include "util.F"
      program sbc
c
c-----------------------------------------------------------------------
c     sbc prepares surface boundary conditions using ...
c     Hellerman wind stress monthly climatology and 
c     Oort air temperature monthly climatology
c     for MOM
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      use mpp_io_mod
      use mom3_mpp_mod
# include "size.h"
      character*60 name
# include "stdunits.h"
# include "coord.h"
c
c     intrph  = 1 when model resolution is coarser than Hellerman
c               2 when model resolution is finer than Hellerman
c     intrpo  = 1 when model resolution is coarser than Oort
c               2 when model resolution is finer than Oort
c
      data  name/' MOM boundary conditions '/
      data intrph /1/, intrpo /1/
      namelist /input/ name, intrph, intrpo
c
      stdout = 6
# ifdef cyclic
      write (stdout,*) '  ==> cyclic option set'
# else
      write (stdout,*) '  ==> cyclic option NOT set'
# endif
      read  (stdin,input)
      write (stdout,input)
c
      if (intrph .eq. 1) then
        write (stdout,*)
     &   ' Assuming model resolution is coarser than Hellerman'
      elseif (intrph .eq. 2) then
        write (stdout,*)
     &   ' Assuming model resolution is finer than Hellerman'
      else
        write (stdout,*) ' intrph must be 1 or 2'
	stop
      endif
c
      if (intrpo .eq. 1) then
        write (stdout,*)
     &  ' Assuming model resolution is coarser than Oort'
      elseif (intrpo .eq. 2) then
        write (stdout,*)
     &  ' Assuming model resolution is finer than Oort'
      else
        write (stdout,*) ' intrpo must be 1 or 2'
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
      call size_check (imt2, jmt2, km2, 'sbc', 'stop')
c
c-----------------------------------------------------------------------
c     interpolate sol`s monthly wind stress to MOM "u,v" grid
c-----------------------------------------------------------------------
c
      call stress (intrph, xu, yu)
c
c-----------------------------------------------------------------------
c     interpolate oorts monthly air temp to MOM "t" grid
c-----------------------------------------------------------------------
c
      call airtem (intrpo, xt, yt)
c
# if defined save_SBC
c-----------------------------------------------------------------------
c     save data into a file for visualization
c-----------------------------------------------------------------------
c
      call sbc_data
# endif
      call mpp_io_exit
c
      stop
      end 


      subroutine stress (intrp, xu, yu)
c
c-----------------------------------------------------------------------
c     linearly interpolate Sol Hellerman`s monthly wind stress
c     climatology to the MOM grid
c-----------------------------------------------------------------------
c
# include "size.h"
      parameter (ib=182, jb=92, lenw=imt*jmt)
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
c     read in Sol Hellerman`s monthly wind stress climatology
c     (the 13th record is the annual mean. k=(1,2) = (taux,tauy))
c-----------------------------------------------------------------------
c
      call getunit (nutauc, 'hellerman.tau'
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
      call getunit (nutaux, 'taux.mom'
     &,             'unformatted sequential rewind ieee')
      call getunit (nutauy, 'tauy.mom'
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
     &       'Hellerman wind stress grid definition:'
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


      subroutine airtem (intrp, xt, yt)
c
c-----------------------------------------------------------------------
c     interp Oort`s climatological air temperature to MOM "t" grid.
c-----------------------------------------------------------------------
c
# include "size.h"
      parameter (iw=73, jw=75, lenw=imt*jmt)
# include "stdunits.h"
      character*32 stamps
      character*80 iotext
      dimension xt(imt), yt(jmt)
      dimension stamps(13), dpm(13)
      common /airs/ cx(iw), cy(jw), clmair(iw,jw,13)
      common /air1/ work(lenw)
      dimension d(imt,jmt), airt(iw,jw)
      dimension airtp(iw), cxt(iw)
c
c-----------------------------------------------------------------------
c     read in Oort`s monthly air temperature climatology
c     (the 13th record is the annual mean)
c-----------------------------------------------------------------------
c
      call getunit (nuairc, 'oorts.air'
     &,             'unformatted sequential rewind ieee')
      do n=1,13
        read (nuairc)
        read (nuairc) stamps(n), dpm(n), lons, lats, nn
     &, cx, cy, ((clmair(i,j,n),i=1,iw),j=1,jw)
      enddo        
c
c-----------------------------------------------------------------------
c     set up the units
c     lrec is in bytes for the cray (words*8)
c-----------------------------------------------------------------------
c
      call getunit (nuair, 'airtemp.mom'
     &,             'unformatted sequential rewind ieee')
c
      do n=1,13
c
c-----------------------------------------------------------------------
c       translate data then interp to MOM "t" grid
c-----------------------------------------------------------------------
c
# ifdef cyclic
        call tranlon (clmair(1,1,n), iw, iw-1, jw, airtp
     &,               cx, xt(2), imt-2, cxt)
c
        if (intrp .eq. 1) then
          call ftc (clmair(1,1,n), iw, jw, cxt, cy, d, imt, jmt, 2
     &,         imt-1, 1, jmt, xt, yt, n, work, lenw)
        else if (intrp .eq. 2) then
          call ctf (clmair(1,1,n), iw, jw, cxt, cy, d, imt, jmt, 2
     &,         imt-1, 1, jmt, xt, yt, n, work, lenw)
        else
          write (stdout,*) ' Error: intrp must = 1 or 2'
	  stop
        endif
c
        do j=1,jmt
          d(1,j)   = d(imt-1,j)
          d(imt,j) = d(2,j)
        enddo
# else
        call tranlon (clmair(1,1,n), iw, iw-1, jw, airtp
     &,               cx, xt, imt, cxt)
        if (intrp .eq. 1) then
          call ftc (clmair(1,1,n), iw, jw, cxt, cy, d, imt, jmt, 1
     &,         imt, 1, jmt, xt, yt, n, work, lenw)
        else if (intrp .eq. 2) then
          call ctf (clmair(1,1,n), iw, jw, cxt, cy, d, imt, jmt, 1
     &,         imt, 1, jmt, xt, yt, n, work, lenw)
        else
          write (stdout,*) ' Error: intrp must = 1 or 2'
	  stop
        endif
# endif
        if (n .eq. 1) then
          write(stdout,'(///,1x, a30,/)') 
     &     'Oort air temp grid definition:'
          write (stdout,'(//,1x,a30,/)') '            latitudes:'
          write (stdout,'(1x,10f10.5)') cy
          write (stdout,'(//,1x,a30,/)') '  original longitudes:'
          write (stdout,'(1x,10f10.5)') cx
          write (stdout,'(//,1x,a30,/)') 'translated longitudes:'
          write (stdout,'(1x,10f10.5)') cxt
        endif
c
      iotext='dim x(i),y(j),at(i,j);read(nu)stamp,avg,i,j,mon,x,y,at'
        write (stdout,*) iotext
        write (nuair) iotext
        write (nuair) stamps(n), dpm(n), imt, jmt, n, xt, yt, d
        write (stdout,*) ' ==> air temp done ',stamps(n), '  ',dpm(n)
     &, ' day average'
        call scope (d, imt, imt, jmt, 'air temp')
      enddo
      call relunit (nuair)
      return
      end

# if defined save_SBC
      subroutine sbc_data
c
c-----------------------------------------------------------------------
c     Saves an output file with the SBC data 
c
c     author: SM Griffies  e-mail==>smg@gfdl.gov
c-----------------------------------------------------------------------
c
      use mom3_mpp_mod
      use mpp_io_mod
#  include "param.h"
#  include "coord.h"
#  include "diag.h"      
#  include "tmngr.h"
#  include "vers.h"
      character*120 file_name, gvtitle
      dimension buf(imt,jmt), xb(imt), yb(jmt)
      character iotext*80
      character*6 taux_name,tauy_name,airt_name
      dimension taux_name(13),tauy_name(13),airt_name(13)      
      character*30 taux_title,tauy_title
      dimension taux_title(13),tauy_title(13)
      character*17 airt_title
      dimension airt_title(13)

c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_sbc(13,3)
      type (axistype), save :: axis_xt
      type (axistype), save :: axis_yt
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
      data airt_name/
     &               'airt1', 'airt2', 'airt3'
     &,              'airt4', 'airt5', 'airt6'
     &,              'airt7', 'airt8', 'airt9'
     &,              'airt10', 'airt11', 'airt12','airt13'/      

      data taux_title/
     &               'Hellerman jan zonal windstress'
     &,              'Hellerman feb zonal windstress'
     &,              'Hellerman mar zonal windstress'
     &,              'Hellerman apr zonal windstress'
     &,              'Hellerman may zonal windstress'
     &,              'Hellerman jun zonal windstress'
     &,              'Hellerman jul zonal windstress'
     &,              'Hellerman aug zonal windstress'
     &,              'Hellerman sep zonal windstress'
     &,              'Hellerman oct zonal windstress'
     &,              'Hellerman nov zonal windstress'
     &,              'Hellerman dec zonal windstress'
     &,              'Hellerman ann zonal windstress'/

      data tauy_title/
     &               'Hellerman jan merid windstress'
     &,              'Hellerman feb merid windstress'
     &,              'Hellerman mar merid windstress'
     &,              'Hellerman apr merid windstress'
     &,              'Hellerman may merid windstress'
     &,              'Hellerman jun merid windstress'
     &,              'Hellerman jul merid windstress'
     &,              'Hellerman aug merid windstress'
     &,              'Hellerman sep merid windstress'
     &,              'Hellerman oct merid windstress'
     &,              'Hellerman nov merid windstress'
     &,              'Hellerman dec merid windstress'
     &,              'Hellerman ann merid windstress'/
      
      data airt_title/
     &               'Oort jan air temp'
     &,              'Oort feb air temp'
     &,              'Oort mar air temp'
     &,              'Oort apr air temp'
     &,              'Oort may air temp'
     &,              'Oort jun air temp'
     &,              'Oort jul air temp'
     &,              'Oort aug air temp'
     &,              'Oort sep air temp'
     &,              'Oort oct air temp'
     &,              'Oort nov air temp'
     &,              'Oort dec air temp' 
     &,              'Oort ann air temp' /
      
      file_name = 'sbc.dta'
      gvtitle   = 'Surface Boundary Conditions ' // momver

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
     &       io, axis_xt, 'xt_i', 'degrees_E', 'Longitude of T points',
     &       cartesian='X', data=xt(ibeg:iend))

      call mpp_write_meta(
     &         io, axis_yt, 'yt_j', 'degrees_N', 'Latitude of T points',
     &         cartesian='Y', data=yt(jbeg:jend))

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
      
      do n=1,13
        call mpp_write_meta( 
     &         io, field_sbc(n,3), (/axis_xt,axis_yt/), 
     &        airt_name(n), 'deg C', airt_title(n),  
     &        -1.0e6, 1.e6, -1.0E+34, +1.0E+34)
      enddo
c
c     figure title 
c
      call mpp_write_meta(io, 'title', cval=gvtitle )

c-----------------------------------------------------------------------
c     write axes now that all meta information has been written
c-----------------------------------------------------------------------
      
      call mpp_write (io,axis_xt)
      call mpp_write (io,axis_yt)
      call mpp_write (io,axis_xu)
      call mpp_write (io,axis_yu)

c-----------------------------------------------------------------------
c     get data from .mom files and then mpp_write
c-----------------------------------------------------------------------
      
      call getunit (nutaux, 'taux.mom'
     &,             'unformatted sequential rewind ieee')
      call getunit (nutauy, 'tauy.mom'
     &,             'unformatted sequential rewind ieee')
      call getunit (nuair, 'airtemp.mom'
     &,             'unformatted sequential rewind ieee')
      
      rewind nutaux
      rewind nutauy
      rewind nuair      

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
      
      do n=1,13        
        buf(:,:) = field_sbc(n,3)%missing
        read (nuair) iotext
        read (nuair) stamp, dpm, iimt, jjmt, nn, xb, yb, buf
        call mpp_write(io, field_sbc(n,3), buf)
      enddo

      write (stdout,'(a,a)')
     &   ' => SBC data written in to file ', file_name
      return
      end
# endif
      
#else
      subroutine sbc
      return
      end
#endif
