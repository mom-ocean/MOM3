#if defined obc_south || defined obc_north
# define obc1
#endif
#if defined obc_west || defined obc_east
# define obc2
#endif
#ifdef makeobc
      program mkobc
c
c-----------------------------------------------------------------------
c     Prepare sponge layers for damping T and S back to data at
c     northern and southern open boundaries.
c
c     If option "readferret" is set, open boundary data will be read
c     from an external file, written by FERRET. Subroutine "getferdata"
c     has to be adapted.
c     If the above option is not set, open boundary data will be
c     constructed from the initial data. Therefore this routine must
c     be run after run_ic.
c
c     author:   Arne Biastoch   e-mail: abiastoch@ifm.uni-kiel.de
c     based on sponge.F in PREP_DATA
c-----------------------------------------------------------------------
c
#include "size.h"
      integer stdout, stdin
      parameter (stdout=6, stdin=5)
      parameter (spval= -9.999998e+33)
      parameter (itm=13)
      character*32 stamp
      character*80 iotext
c
#ifdef readferret
      dimension stamp(itm), dpm(itm)
#endif
c
      dimension dxtdeg(imt), dytdeg(jmt), dzt(km)
      dimension dxudeg(imt), dyudeg(jmt), dzw(0:km)
      dimension xt(imt),xu(imt),yt(jmt),yu(jmt),zt(km),zw(km)
      dimension kmt(imt,jmt), map(imt,jmt)
      dimension ht(imt,jmt)
c
#ifndef readferret
      character*12 tempin, saltin
      dimension tempin(itm), saltin(itm)
      dimension iotin(itm), iosin(itm)
c
      data tempin/
     &               'jan.temp.mom', 'feb.temp.mom', 'mar.temp.mom'
     &,              'apr.temp.mom', 'may.temp.mom', 'jun.temp.mom'
     &,              'jul.temp.mom', 'aug.temp.mom', 'sep.temp.mom'
     &,              'oct.temp.mom', 'nov.temp.mom', 'dec.temp.mom'
     &,              'ann.temp.mom' /
c
      data saltin/
     &               'jan.salt.mom', 'feb.salt.mom', 'mar.salt.mom'
     &,              'apr.salt.mom', 'may.salt.mom', 'jun.salt.mom'
     &,              'jul.salt.mom', 'aug.salt.mom', 'sep.salt.mom'
     &,              'oct.salt.mom', 'nov.salt.mom', 'dec.salt.mom'
     &,              'ann.salt.mom' /
#endif
      dimension temps(imt,km), salts(imt,km)
      dimension tempn(imt,km), saltn(imt,km) 
      dimension tempw(jmt,km), saltw(jmt,km) 
      dimension tempe(jmt,km), salte(jmt,km) 
      dimension    tempsouth(imt,km,itm), saltsouth(imt,km,itm)
     &            ,tempnorth(imt,km,itm), saltnorth(imt,km,itm) 
     &            ,tempwest(jmt,km,itm), saltwest(jmt,km,itm) 
     &            ,tempeast(jmt,km,itm), salteast(jmt,km,itm) 
c
#ifdef readferret
      data stamp /
     &               'm/d/y= 1/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 2/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 3/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 4/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 5/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 6/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 7/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 8/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 9/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y=10/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y=11/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y=12/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 1/ 1/0000, h:m:s= 0: 0: 0'
     &            /
c
      data dpm /31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0
     &, 30.0, 31.0, 30.0, 31.0, 365.0/
#endif
c-----------------------------------------------------------------------
c     initialize arrays         
c-----------------------------------------------------------------------
c
      do it=1,itm
        do k=1,km 
          do i=1,imt
            tempsouth(i,k,it) = spval
            tempnorth(i,k,it) = spval
            saltsouth(i,k,it) = spval
            saltnorth(i,k,it) = spval
          enddo
          do j=1,jmt
            tempwest(j,k,it) = spval
            tempeast(j,k,it) = spval
            saltwest(j,k,it) = spval
            salteast(j,k,it) = spval
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     initialize the i/o units
c-----------------------------------------------------------------------
c
      call ioinit
#ifdef readferret
# ifdef obc_south
        call getunit (ioget1, 'obc_tracer_south.dat'
     &,               'formatted sequential rewind')
# endif
# ifdef obc_south
        call getunit (ioget2, 'obc_tracer_north.dat'
     &,               'formatted sequential rewind')
# endif
# ifdef obc_west
        call getunit (ioget3, 'obc_tracer_west.dat'
     &,               'formatted sequential rewind')
# endif
# ifdef obc_east
        call getunit (ioget4, 'obc_tracer_east.dat'
     &,               'formatted sequential rewind')
# endif
#else
      do n=1,itm
        call getunit (iotin(n), tempin(n)
     &,               'unformatted sequential rewind ieee')
        call getunit (iosin(n), saltin(n)
     &,               'unformatted sequential rewind ieee')
      enddo
#endif
#ifdef obc1
      call getunit (iosp1, 'obc1.mom'
     &,             'unformatted sequential rewind ieee')
#endif
#ifdef obc2
      call getunit (iosp2, 'obc2.mom'
     &,             'unformatted sequential rewind ieee')
#endif
c
c-----------------------------------------------------------------------
c     construct the grid.
c-----------------------------------------------------------------------
c
      write (stdout,'(//1x,a,//)') 'Construct the MOM grid'
c
      maxlen = max(imt,jmt,km)
      call gcoord (maxlen, iimt, ijmt, ikm, dxtdeg, dytdeg, dxudeg
     &,            dyudeg, dzt, dzw, xt, xu, yt, yu, zt, zw)
c
      if (imt .ne. iimt .or. jmt .gt. ijmt .or. km .gt. ikm) then
        write (stdout,*) ' => parameter mismatch in mkobc'
	write (stdout,*) 'imt=',imt,' iimt=',iimt
	write (stdout,*) 'jmt=',jmt,' ijmt=',ijmt
	write (stdout,*) 'km=',km,' ikm=',ikm
	stop '=> mkobc'
      endif
c
c-----------------------------------------------------------------------
c     generate the kmt field on the grid defined by grids.F
c-----------------------------------------------------------------------
c
      write (stdout,'(/,20x,a/)') ' Driving the topog module...'
      call topog (kmt, map, xt, yt, zt, xu, yu, zw, dzt
     &,           iimt, ijmt, ikm, ht)
c
c-----------------------------------------------------------------------
c     apply "sponge" regions adjacent to southern & northern walls by
c     setting damping factor of "1/spngn" days at jrow=jmt and 
c     "1/spngs" days at jrow=1.
c     
c     damping will be of newtonian form: -spng(j)*(T-sponge)
c     where spng(j) = 1/amp, and amp = some factor * spngts
c     and "sponge" is the Levitus T & S from jrows = 1 and jmt
c-----------------------------------------------------------------------
c
c   USER INPUT
c
      spngts = 25.0*86400.0
c
      spngn  = 1.0/spngts
      spngs  = 1.0/spngts
      spngw  = 1.0/spngts
      spnge  = 1.0/spngts
c
c   END USER INPUT
c
c     show sponge damping factors
c
      write (stdout,'(/a/)') ' Obc damping factors (1/sec):'
        write (stdout,'(1x,a,f10.5,a,e14.7,a, a, e14.7, a)')
     &  'jrow=   2, lat=',yt(2),' spngs=',spngs,' (1/sec)'
     &,' . 1/spngs = ',1.0/(86400.0*spngs), 'days'
        write (stdout,'(1x,a,f10.5,a,e14.7,a, a, e14.7, a)')
     &  'jrow= jmt-1, lat=',yt(jmt-1),' spngn=',spngn,' (1/sec)'
     &,' . 1/spngn = ',1.0/(86400.0*spngn), 'days'
        write (stdout,'(1x,a,f10.5,a,e14.7,a, a, e14.7, a)')
     &  '   i=   2, lon=',xt(2),' spngw=',spngw,' (1/sec)'
     &,' . 1/spngw = ',1.0/(86400.0*spngw), 'days'
        write (stdout,'(1x,a,f10.5,a,e14.7,a, a, e14.7, a)')
     &  '   i= imt-1, lon=',xt(imt-1),' spnge=',spnge,' (1/sec)'
     &,' . 1/spnge = ',1.0/(86400.0*spnge), 'days'
c
#ifdef readferret
# ifdef obc_south
      call getferdata(ioget1, imt, km, itm, tempsouth, saltsouth) 
# endif
# ifdef obc_north
      call getferdata(ioget2, imt, km, itm, tempnorth, saltnorth) 
# endif
# ifdef obc_west
      call getferdata(ioget3, jmt, km, itm, tempwest, saltwest) 
# endif
# ifdef obc_east
      call getferdata(ioget4, jmt, km, itm, tempeast, salteast) 
# endif
c
c       save the sponge
c
        iotext =
     &'        data written by FERRET                             '
c
      do month=1,13
c
        do k=1,km
          do i=1,imt
             if (tempsouth(i,k,month) .eq. 0.) then
                temps(i,k) = 0.
                salts(i,k) = 0.
             else
                temps(i,k) = tempsouth(i,k,month)
                salts(i,k) = saltsouth(i,k,month)*0.001 - 0.035
             endif
c
             if (tempnorth(i,k,month) .eq. 0.) then
                tempn(i,k) = 0.
                saltn(i,k) = 0.
             else
                tempn(i,k) = tempnorth(i,k,month)
                saltn(i,k) = saltnorth(i,k,month)*0.001 - 0.035
             endif
c
           enddo
          do j=1,jmt
             if (tempwest(j,k,month) .eq. 0.) then
                tempw(j,k) = 0.
                saltw(j,k) = 0.
             else
                tempw(j,k) = tempwest(j,k,month)
                saltw(j,k) = saltwest(j,k,month)*0.001 - 0.035
             endif
c
             if (tempeast(j,k,month) .eq. 0.) then
                tempe(j,k) = 0.
                salte(j,k) = 0.
             else
                tempe(j,k) = tempeast(j,k,month)
                salte(j,k) = salteast(j,k,month)*0.001 - 0.035
             endif
c
           enddo
        enddo
c
        js=0
        jn=0
        write (iosp1) iotext
        write (iosp1) stamp(month), dpm(month), imt, km, jmt
     &,               js, jn, month
#ifdef obc_south     
     & ,spngs,((temps(i,k),i=1,imt),k=1,km),((salts(i,k),i=1,imt),k=1,km)
#endif
#ifdef obc_north
     & ,spngn,((tempn(i,k),i=1,imt),k=1,km),((saltn(i,k),i=1,imt),k=1,km)
#endif
c
        js=0
        jn=0
        write (iosp2) iotext
        write (iosp2) stamp(month), dpm(month), imt, km, jmt
     &,               js, jn, month
#ifdef obc_west     
     &               ,spngw,((tempw(j,k),j=1,jmt),k=1,km)
     &                     ,((saltw(j,k),j=1,jmt),k=1,km)
#endif
#ifdef obc_east     
     &               ,spnge,((tempe(j,k),j=1,jmt),k=1,km)
     &                     ,((salte(j,k),j=1,jmt),k=1,km)
#endif
c	 
         write (stdout,'(1x,a,a,a,g14.7,a)')
     &   'Writing OBC for ',stamp(month)
     &,  ' averaged over previous ',dpm(month), ' days'
c
      enddo
c
#else
c
#ifdef obc1
c
c     get the sponge rows
c
      do month=1,itm
c
c       skip row 1
c
        read (iotin(month))
        read (iotin(month))
        read (iosin(month))
	read (iosin(month))        
c
c       read jrow = 2
c
        read (iotin(month)) 
        read (iotin(month))stamp, dpm, im, kk, js, alt, m, xt, zt, temps
c
        read (iosin(month))
        read (iosin(month))stamp, dpm, im, kk, js, alt, m, xt, zt, salts
c
        do i=1,imt
          do k=1,km
             if (k .gt. kmt(i,1)) then
                tempsouth(i,k,month) = spval
                saltsouth(i,k,month) = spval
             else
                tempsouth(i,k,month) = temps(i,k)
                saltsouth(i,k,month) = (salts(i,k)+0.035)*1000.
             endif
          enddo
        enddo
c
c       skip rows 3 .. jmt-2
c
        do jrow=3,jmt-2
          read (iotin(month))
          read (iotin(month))
          read (iosin(month))
	  read (iosin(month))        
	enddo
c
c       read jrow = jmt-1
c
        read (iotin(month)) 
        read (iotin(month))stamp, dpm, im, kk, jn, alt, m, xt, zt, tempn
c
        read (iosin(month))
        read (iosin(month))stamp, dpm, im, kk, jn, alt, m, xt, zt, saltn
c
        do i=1,imt
          do k=1,km
             if (k .gt. kmt(i,jmt)) then
                tempnorth(i,k,month) = spval
                saltnorth(i,k,month) = spval
             else
                tempnorth(i,k,month) = tempn(i,k)
                saltnorth(i,k,month) = (saltn(i,k)+0.035)*1000.
             endif
          enddo
        enddo
c
c       skip row jmt
c
        read (iotin(month))
        read (iotin(month))
        read (iosin(month))
	read (iosin(month))        
c
c       save the sponge
c
        iotext =
     &'dim ss(j),sn(j),d(i,k,4);read(n)stmp,avg,i,k,j,js,jn,mon,ss,sn,d'
        write (iosp1) iotext
        write (iosp1) stamp,dpm,imt,km,jmt,js,jn,month
#ifdef obc_south     
     & ,spngs,((temps(i,k),i=1,imt),k=1,km),((salts(i,k),i=1,imt),k=1,km)
#endif
#ifdef obc_north
     & ,spngn,((tempn(i,k),i=1,imt),k=1,km),((saltn(i,k),i=1,imt),k=1,km)
#endif
        write (stdout,'(1x,a,a,a,g14.7,a)')
     &   'Writing OBC for ',stamp
     &,  ' averaged over previous ',dpm, ' days'
      enddo
#endif
c
c
#ifdef obc2
c
      do month=1,13
        rewind iotin(month)
        rewind iosin(month)
      enddo 
c
c     get the sponge rows for the meridional sponges
c
      do month=1,13
c
c       read jrow = 1 to most southern row of segment 1
c
         do jrow=1,jmt
c	 
          read (iotin(month)) 
          read (iotin(month))stamp, dpm, im, kk, jxs, alt, m, xt, zt
     &  , ((temps(i,k), i=1,imt),k=1,km)
          do k=1,km
            tempw(jrow,k)= temps(2,    k)
            tempe(jrow,k)= temps(imt-1,k)
             if (k .gt. kmt(1,jrow)) then
                tempwest(jrow,k,month) = spval
             else
                tempwest(jrow,k,month) = tempw(jrow,k)
             endif
             if (k .gt. kmt(imt,jrow)) then
                tempeast(jrow,k,month) = spval
             else
                tempeast(jrow,k,month) = tempe(jrow,k)
             endif
          enddo
c
          read (iosin(month))
          read (iosin(month))stamp, dpm, im, kk, jxs, alt, m, xt, zt
     &  , ((salts(i,k), i=1,imt),k=1,km)
c
          do k=1,km
            saltw(jrow,k)= salts(2,    k)
            salte(jrow,k)= salts(imt-1,k)
             if (k .gt. kmt(1,jrow)) then
                saltwest(jrow,k,month) = spval
             else
                saltwest(jrow,k,month) = (saltw(jrow,k)+0.035)*1000.
             endif
             if (k .gt. kmt(imt,jrow)) then
                salteast(jrow,k,month) = spval
             else
                salteast(jrow,k,month) = (salte(jrow,k)+0.035)*1000.
             endif
          enddo
c
        enddo
c
c       save the sponge
c
        js=0
	jn=0
        iotext =
     &'dim ss(j),sn(j),d(i,k,4);read(n)stmp,avg,i,k,j,js,jn,mon,ss,sn,d'
        write (iosp2) iotext
        write (iosp2) stamp,dpm,imt,km,jmt,js,jn,month
#ifdef obc_west     
     &               ,spngw,((tempw(j,k),j=1,jmt),k=1,km)
     &                     ,((saltw(j,k),j=1,jmt),k=1,km)
#endif
#ifdef obc_east
     &               ,spnge,((tempe(j,k),j=1,jmt),k=1,km)
     &                     ,((salte(j,k),j=1,jmt),k=1,km)
#endif
c      do j=1,jmt
c      do k=1,km
c      if (k .le. kmt(2,j)) then
c      write(stdout,*)'j,k,tempw(j,k),saltw(j,k): ',
c     &           j,k,tempw(j,k),(saltw(j,k)+0.035)*1000.
c      endif
c      if (k .le. kmt(imt-1,j)) then
c      write(stdout,*)'j,k,tempe(j,k),salte(j,k): ',
c     &           j,k,tempe(j,k),(salte(j,k)+0.035)*1000.
c      endif
c      enddo
c      enddo
c
        write (stdout,'(1x,a,a,a,g14.7,a)')
     &   'Writing meridional sponges for ',stamp
     &,  ' averaged over previous ',dpm, ' days'
      enddo
      call relunit (iosp2)
#endif

#endif
c
#ifdef write_netcdf
c
c      write open boundaries in NETCDF file
c
      call  write1_netcdf(itm,xt,xu,yt,yu,zt,zw
     &           , tempsouth,saltsouth,tempnorth,saltnorth
     &           , tempwest,saltwest,tempeast,salteast)
#endif
c
      stop
      end

#define driver_only
#include "iomngr.F"
#include "util.F"
#include "topog.F"
#include "grids.F"
#include "size_check.F"

#else
      subroutine sponge
      return
      end
#endif

c
#ifdef write_netcdf
      subroutine write1_netcdf(itm,xt,xu,yt,yu,zt,zw
     &           , tempsouth,saltsouth,tempnorth,saltnorth
     &           , tempwest,saltwest,tempeast,salteast)
c
      include 'netcdf.inc'
c
# include "size.h"
# include "switch.h"
# include "tmngr.h"
c
      integer stdin, stdout, stderr
      parameter (stdin = 5, stdout = 6, stderr = 6)
      parameter (spval= -9.999998e+33)
c
      character iotext*80,expnam*60
c
c error return code and netCDF id
      integer  iret, ncid
c dimension ids
      integer  Longitude_tdim, Longitude_udim, Depth_tdim, Depth_wdim,
     &         Latitude_tdim,  Latitude_udim,  Timedim
c variable ids
      integer  Longitude_tid,  Longitude_uid, Depth_tid, Depth_wid,
     &         Latitude_tid,   Latitude_uid,  
     &         tempnid,tempsid,saltnid,saltsid, 
     &         tempwid,tempeid,saltwid,salteid, Timeid
c
c variable shapes, corner and edge lengths
      integer dims(4), corner(4), edges(4)
c data variables
      real Longitude_t(imt)
      real Longitude_u(imt)
      real Depth_t(km)
      real Depth_w(km+1)
      real Latitude_t(jmt)
      real Latitude_u(jmt)
      real         Time(itm)
c attribute vectors
      integer  longval(1)
      real  floatval(2),Monat(13)
c time information
      character*20 timestamp
      character*24 obcname1,obcname2
c
      dimension xt(imt),xu(imt),yt(jmt),yu(jmt),zt(km),zw(km)
      dimension    tempsouth(imt,km,itm), saltsouth(imt,km,itm)
     &           , tempnorth(imt,km,itm), saltnorth(imt,km,itm) 
     &           , tempwest(jmt,km,itm), saltwest(jmt,km,itm) 
     &           , tempeast(jmt,km,itm), salteast(jmt,km,itm) 
      data Monat / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13/
c
      if (itm .eq. 1) then
          Time(itm) = 13
      else 
        do m=1,itm
          Time(m) = Monat(m)
        enddo
      endif
c
c-----------------------------------------------------------------------
c     define netCDF variables
c-----------------------------------------------------------------------
c
c enter define mode: NCCLOB=overwrite, NCNOCLOB=don't overwrite
c
      ncid = nccre ('obc.cdf', NCCLOB, iret)
c
c define dimensions
c
      Longitude_tdim = ncddef(ncid, 'Longitude_t', imt, iret)
      Latitude_tdim  = ncddef(ncid, 'Latitude_t', jmt, iret)
      Depth_tdim     = ncddef(ncid, 'Depth_t',      km, iret)
      Depth_wdim     = ncddef(ncid, 'Depth_w',    km+1, iret)
      Timedim        = ncddef(ncid, 'Time',     NCUNLIM, iret)
c
c define variables
c
      dims(1) = Longitude_tdim
      Longitude_tid =  ncvdef (ncid,'Longitude_t',NCFLOAT,1,dims,iret)
      dims(1) = Latitude_tdim
      Latitude_tid =  ncvdef (ncid,'Latitude_t',NCFLOAT,1,dims,iret)
      dims(1) = Depth_tdim
      Depth_tid =      ncvdef (ncid,'Depth_t',    NCFLOAT,1,dims,iret)
      dims(1) = Depth_wdim
      Depth_wid =      ncvdef (ncid,'Depth_w',    NCFLOAT,1,dims,iret)
      dims(1) = Timedim
      Timeid  =        ncvdef (ncid,'Time',       NCFLOAT,1,dims,iret)
      dims(1) = Longitude_tdim
      dims(2) = Depth_tdim
      dims(3) = Timedim
#ifdef obc_north
      tempnid=     ncvdef (ncid,'tempn',  NCFLOAT,3,dims,iret)
      saltnid=     ncvdef (ncid,'saltn',  NCFLOAT,3,dims,iret)
#endif
#ifdef obc_south
      tempsid=     ncvdef (ncid,'temps',  NCFLOAT,3,dims,iret)
      saltsid=     ncvdef (ncid,'salts',  NCFLOAT,3,dims,iret)
#endif
      dims(1) = Latitude_tdim
      dims(2) = Depth_tdim
      dims(3) = Timedim
#ifdef obc_west
      tempwid=     ncvdef (ncid,'tempw',  NCFLOAT,3,dims,iret)
      saltwid=     ncvdef (ncid,'saltw',  NCFLOAT,3,dims,iret)
#endif
#ifdef obc_east
      tempeid=     ncvdef (ncid,'tempe',  NCFLOAT,3,dims,iret)
      salteid=     ncvdef (ncid,'salte',  NCFLOAT,3,dims,iret)
#endif
c
c-----------------------------------------------------------------------
c     assign attributes
c-----------------------------------------------------------------------
c
      call ncaptc(ncid, Latitude_tid, 'long_name',  NCCHAR,  24,
     &            'Latitude on T grid      ', iret)
      call ncaptc(ncid, Latitude_tid, 'units',      NCCHAR,  16, 
     &            'degrees_N       ', iret)  
      call ncaptc(ncid, Longitude_tid, 'long_name', NCCHAR,  24,
     &            'Longitude on T grid     ', iret)
      call ncaptc(ncid, Longitude_tid, 'units',     NCCHAR,  16,
     &            'degrees_W       ', iret)
      call ncaptc(ncid, Depth_tid, 'long_name',     NCCHAR,  24,
     &            'Depth of T grid points  ', iret)
      call ncaptc(ncid, Depth_tid, 'units',         NCCHAR,  16,
     &            'meters          ', iret)
      call ncaptc(ncid, Depth_tid, 'positive',      NCCHAR,   4,
     &            'down'            , iret)
      call ncaptc(ncid, Depth_tid, 'point_spacing', NCCHAR,   6,
     &            'uneven'          , iret)
      call ncaptc(ncid, Depth_tid, 'edges'        , NCCHAR,   7,
     &            'Depth_w'         , iret)
      call ncaptc(ncid, Depth_wid, 'long_name',     NCCHAR,  24,
     &            'Depth at bottom of T box', iret)
      call ncaptc(ncid, Depth_wid, 'units',         NCCHAR,  16,
     &            'meters          ', iret)
      call ncaptc(ncid, Depth_wid, 'positive',      NCCHAR,   4,
     &            'down'            , iret)
      call ncaptc(ncid,       Timeid, 'long_name',  NCCHAR,  24,
     &            'Time                    ', iret)
      call ncaptc(ncid,       Timeid, 'units',      NCCHAR,  16,
     &            'months          ', iret)
      call ncaptc(ncid, Timeid,'time_origin' ,      NCCHAR,  20,
     &  timestamp)
      call ncaptc(ncid, Timeid,'time_origin' ,      NCCHAR,  20,
     &  '01-JAN-0000 00:00:00', iret)
c
#ifdef obc_north
      if (yt(jmt).ge.0.) then
        write (obcname1,'("TEMPERATURE AT ",F6.2," N ")') yt(jmt)
        write (obcname2,'("SALINITY AT ",F6.2," N    ")') yt(jmt)
      else
        write (obcname1,'("TEMPERATURE AT ",F6.2," S ")') yt(jmt)
        write (obcname2,'("SALINITY AT ",F6.2," S    ")') yt(jmt)
      endif
      call ncaptc(ncid, tempnid,'long_name',NCCHAR,24,obcname1,iret) 
      call ncaptc(ncid, tempnid, 'units', NCCHAR,  16,
     &            'deg C           ', iret) 
      floatval(1) = spval
      call ncapt(ncid, tempnid, 'missing_value', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = spval
      call ncapt(ncid, tempnid, '_FillValue', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = -5
      floatval(2) = 50
      call ncapt(ncid, tempnid, 'valid_range', NCFLOAT, 2, 
     &           floatval, iret)
      call ncaptc(ncid, saltnid,'long_name',NCCHAR,24,obcname2,iret) 
      call ncaptc(ncid, saltnid, 'units', NCCHAR,  16,
     &            'PSU             ', iret) 
      floatval(1) = spval
      call ncapt(ncid, saltnid, 'missing_value', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = spval
      call ncapt(ncid, saltnid, '_FillValue', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = -100
      floatval(2) = 100
      call ncapt(ncid, saltnid, 'valid_range', NCFLOAT, 2, 
     &           floatval, iret)
#endif
c
#ifdef obc_south
      if (yt(1).ge.0.) then
        write (obcname1,'("TEMPERATURE AT ",F6.2," N ")') yt(1)
        write (obcname2,'("SALINITY AT ",F6.2," N    ")') yt(1)
      else
        write (obcname1,'("TEMPERATURE AT ",F6.2," S ")') yt(1)
        write (obcname2,'("SALINITY AT ",F6.2," S    ")') yt(1)
      endif
      call ncaptc(ncid, tempsid,'long_name',NCCHAR,24,obcname1,iret) 
      call ncaptc(ncid, tempsid, 'units', NCCHAR,  16,
     &            'deg C           ', iret) 
      floatval(1) = spval
      call ncapt(ncid, tempsid, 'missing_value', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = spval
      call ncapt(ncid, tempsid, '_FillValue', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = -5
      floatval(2) = 50
      call ncapt(ncid, tempsid, 'valid_range', NCFLOAT, 2, 
     &           floatval, iret)
      call ncaptc(ncid, saltsid,'long_name',NCCHAR,24,obcname2,iret) 
      call ncaptc(ncid, saltsid, 'units', NCCHAR,  16,
     &            'PSU             ', iret) 
      floatval(1) = spval
      call ncapt(ncid, saltsid, 'missing_value', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = spval
      call ncapt(ncid, saltsid, '_FillValue', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = -100
      floatval(2) = 100
      call ncapt(ncid, saltsid, 'valid_range', NCFLOAT, 2, 
     &           floatval, iret)
#endif
#ifdef obc_east
      if (xt(imt).ge.0.) then
        write (obcname1,'("TEMPERATURE AT ",F6.2," W ")') xt(imt)
        write (obcname2,'("SALINITY AT ",F6.2," W    ")') xt(imt)
      else
        write (obcname1,'("TEMPERATURE AT ",F6.2," E ")') xt(imt)
        write (obcname2,'("SALINITY AT ",F6.2," E    ")') xt(imt)
      endif
      call ncaptc(ncid, tempeid,'long_name',NCCHAR,24,obcname1,iret) 
      call ncaptc(ncid, tempeid, 'units', NCCHAR,  16,
     &            'deg C           ', iret) 
      floatval(1) = spval
      call ncapt(ncid, tempeid, 'missing_value', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = spval
      call ncapt(ncid, tempeid, '_FillValue', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = -5
      floatval(2) = 50
      call ncapt(ncid, tempeid, 'valid_range', NCFLOAT, 2, 
     &           floatval, iret)
      call ncaptc(ncid, salteid,'long_name',NCCHAR,24,obcname2,iret) 
      call ncaptc(ncid, salteid, 'units', NCCHAR,  16,
     &            'PSU             ', iret) 
      floatval(1) = spval
      call ncapt(ncid, salteid, 'missing_value', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = spval
      call ncapt(ncid, salteid, '_FillValue', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = -100
      floatval(2) = 100
      call ncapt(ncid, salteid, 'valid_range', NCFLOAT, 2, 
     &           floatval, iret)
#endif
c
#ifdef obc_west
      if (xt(1).ge.0.) then
        write (obcname1,'("TEMPERATURE AT ",F6.2," E ")') xt(1)
        write (obcname2,'("SALINITY AT ",F6.2," E    ")') xt(1)
      else
        write (obcname1,'("TEMPERATURE AT ",F6.2," W ")') xt(1)
        write (obcname2,'("SALINITY AT ",F6.2," W    ")') xt(1)
      endif
      call ncaptc(ncid, tempwid,'long_name',NCCHAR,24,obcname1,iret) 
      call ncaptc(ncid, tempwid, 'units', NCCHAR,  16,
     &            'deg C           ', iret) 
      floatval(1) = spval
      call ncapt(ncid, tempwid, 'missing_value', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = spval
      call ncapt(ncid, tempwid, '_FillValue', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = -5
      floatval(2) = 50
      call ncapt(ncid, tempwid, 'valid_range', NCFLOAT, 2, 
     &           floatval, iret)
      call ncaptc(ncid, saltwid,'long_name',NCCHAR,24,obcname2,iret) 
      call ncaptc(ncid, saltwid, 'units', NCCHAR,  16,
     &            'PSU             ', iret) 
      floatval(1) = spval
      call ncapt(ncid, saltwid, 'missing_value', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = spval
      call ncapt(ncid, saltwid, '_FillValue', NCFLOAT, 1,
     &            floatval, iret)
      floatval(1) = -100
      floatval(2) = 100
      call ncapt(ncid, saltwid, 'valid_range', NCFLOAT, 2, 
     &           floatval, iret)
#endif
c
c
      expnam='OPEN BOUNDARY VALUES'
      call ncaptc(ncid, NCGLOBAL, 'title', NCCHAR, 60, expnam, iret)
c      call ncaptc(ncid, NCGLOBAL, 'history', NCCHAR, 32,
c     &     stamp(itm), iret)
c
c-----------------------------------------------------------------------
c     leave define mode
c-----------------------------------------------------------------------
c
      call ncendf(ncid, iret)
c
c-----------------------------------------------------------------------
c     convert variables
c     change depth to m
c     put land mask
c-----------------------------------------------------------------------
c
      if (xt(1).gt.360 .or. xt(imt).gt.360 .or.
     &    xu(1).gt.360 .or. xu(imt).gt.360) then
      do i=1,imt
         Longitude_t(i)=xt(i)-360.
         Longitude_u(i)=xu(i)-360.
      enddo
      else
      do i=1,imt
         Longitude_t(i)=xt(i)
         Longitude_u(i)=xu(i)
      enddo
      endif
c
      do j=1,jmt
         Latitude_t(j)=yt(j)
         Latitude_u(j)=yu(j)
      enddo
c
         Depth_w(1)  =0.
      do k=1,km
         Depth_t(k)=zt(k)/100.
         Depth_w(k+1)=zw(k)/100.
      enddo
c
c-----------------------------------------------------------------------
c     put variables in netCDF file
c-----------------------------------------------------------------------
c
c store Longitude_t
      corner(1) = 1
      edges(1) = imt
      call ncvpt(ncid, Longitude_tid, corner, edges, Longitude_t, iret)
c store Latitude_t
      corner(1) = 1
      edges(1) = jmt
      call ncvpt(ncid, Latitude_tid, corner, edges, Latitude_t, iret)
c store Depth_t
      corner(1) = 1
      edges(1) = km
      call ncvpt(ncid, Depth_tid, corner, edges, Depth_t, iret)
c store Depth_w
      corner(1) = 1
      edges(1) = km+1
      call ncvpt(ncid, Depth_wid, corner, edges, Depth_w, iret)
c store Time in hours    
      corner(1) = 1
      edges(1) = itm
      call ncvpt(ncid, Timeid, corner, edges, Time, iret)
c store data
      corner(1) = 1
      corner(2) = 1
      corner(3) = 1
      edges(1) = imt
      edges(2) = km
      edges(3) = itm
#ifdef obc_north
      call ncvpt(ncid, tempnid, corner, edges, tempnorth, iret)
      call ncvpt(ncid, saltnid, corner, edges, saltnorth, iret)
#endif
#ifdef obc_south
      call ncvpt(ncid, tempsid, corner, edges, tempsouth, iret)
      call ncvpt(ncid, saltsid, corner, edges, saltsouth, iret)
#endif
      corner(1) = 1
      corner(2) = 1
      corner(3) = 1
      edges(1) = jmt
      edges(2) = km
      edges(3) = itm
#ifdef obc_west
      call ncvpt(ncid, tempwid, corner, edges, tempwest, iret)
      call ncvpt(ncid, saltwid, corner, edges, saltwest, iret)
#endif
#ifdef obc_east
      call ncvpt(ncid, tempeid, corner, edges, tempeast, iret)
      call ncvpt(ncid, salteid, corner, edges, salteast, iret)
#endif
c
      call ncclos (ncid, iret)
c
      return
      end
#endif
c
#ifdef readferret
      subroutine getferdata(ioget, imt, km, itm, data1, data2)
c
c     reads data written by FERRET
c
      parameter (spval= -9.999998e+33)
c
      character*70 dummy(8),dummy2
      dimension data1(imt,km,itm),data2(imt,km,itm)
c
      do 1 m=1,itm
      do 1 k=1,km
      do 1 i=1,imt
         data1(i,k,m)=spval
         data2(i,k,m)=spval
1     continue
c
      do 10 m=1,itm
      do l=1,9
         read(ioget,100)  dummy(l)
100      format(A70)
      enddo
c
      do k=1,km
          read(ioget,100)  dummy2
          read(ioget,100)  dummy2
        do i=2,imt-1
          read(ioget,200)  data1(i,k,m),data2(i,k,m)
200       format(11X,2F7.2)
210       format(2F7.3)
        enddo
          read(ioget,100)  dummy2
      enddo
          rewind ioget
10    continue
c
      close(1)
      return
      end
#endif

