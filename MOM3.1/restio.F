      subroutine read_restart (js, je)
c
c=======================================================================
c     read the restart file
c     
c     author:        a. rosati             e-mail=> ar@gfdl.gov
c                    r. c. pacanowski      e-mail=> rcp@gfdl.gov
c    (mpp items)     v. balaji             e-mail=> vb@gfdl.gov
c=======================================================================
c
      use mom3_mpp_mod
      use arrays_2d
      use memory_window
      character*128 file_name, file_stamp
#include "param.h"
#include "coord.h"
#ifdef trajectories
# include "diag_ptraj.h"
#endif
#include "emode.h"
#include "grdvar.h"
#include "iounit.h"
#ifdef nonconst_diffusivity
# include "isopyc.h"
#endif
# include "mw.h"
#include "taskrows.h"
#include "tmngr.h"
      integer(kind=8) :: akmt
      common /restbuf/ apsi(imt,jmt)
      dimension akmt(imt,jmt)
#if defined rigid_lid_surface_pressure || defined free_surf_on
# if defined time_averaged_surface
      common /restbuf/ aps(imt,jmt,3)
# else
      common /restbuf/ aps(imt,jmt,2)
# endif
#endif
c
#if !defined partial_cell
      dimension hup(1,1), htp(1,1)
#endif

      
      character*32 ostamp, seqrecl
#if !defined sgi_mipspro || defined use_libMPI
      common / restbuf / bufsl(imt,km,nvar)
# if defined pressure_depth
      common / restbuf / bufr(imt,km)
# endif
#else
      dimension bufsl(imt,km,nvar)
# if defined pressure_depth
      dimension bufr(imt,km)
# endif
# ifdef parallel_1d
      pointer( ptrsl, bufsl )
      integer, save :: lensl=0
      save :: ptrsl
#  if defined pressure_depth
      pointer( ptrr, bufr )
      integer, save :: lenr=0
      save :: ptrr
#  endif
# endif

#endif
c
c mjh: explicitly define 8 byte integers for restart I/O
c      this ensures restart portability between CRAY and 
c      non-CRAY machines
      
      integer (kind=8) :: itt8, irstdy8, msrsdy8, kkm8
     &     ,                   iimt8,ijmt8,ilensl8,npart8

c     
c-----------------------------------------------------------------------
c     read restart file (unformatted)
c
c     record 1 on the restart file is the control block (containing
c              the time step and relative year)
c     record`s 2 & 3 are for stream function time levels
c     record`s 4 & 5 are for the d(stream function)/dt fields used
c              to construct a guess for the elliptic solver
c     record 6 is for "kmt" indicating the number of ocean T-cells
c              from the surface to the ocean floor.
#if defined partial_cell
c     record 7 is for "htp" depth to bottom of partial T-cells
#endif
#if defined nonconst_diffusivity
c     record 8 is for "growth2a" squared Eady growth rate sec^-2.
c     record 9 is for "agm" nonconstant Redi/GM diffusivity.
#endif
c     the next "jmt" records contain prognostic variables at "tau"
c     the last "jmt" records contain prognostic variables at "tau+1"
c-----------------------------------------------------------------------
c
#ifdef timing
      call tic ('ocean', 'reading restart')
#endif

#ifdef parallel_1d
# if defined sgi_mipspro && !defined use_libMPI 
      call mpp_malloc( ptrsl, size(bufsl), lensl )
#  ifdef pressure_depth
      call mpp_malloc( ptrr, size(bufr), lenr )
#  endif
# endif
#endif
c
      if (pn .eq. 1) then
c      
        write (seqrecl,'(a,i12)') ' words=',nvar*imt*km
c
c       use generic form of restart file because time stamp is unknown
c
        file_name = 'restart.dta'
        call getunit (iorest, trim(file_name)
     &,             'unformatted sequential rewind'//seqrecl)
c
c-----------------------------------------------------------------------
c       read the time step and stamp from the restart file
c       timestep "itt" at this point corresponds to "tau+1" which was
c       the condition at the end of the previous timestep
c-----------------------------------------------------------------------
c
        read (iorest) stamp, iotext
        read (iorest) itt8, irstdy8, msrsdy8, iimt8, ijmt8, kkm8
        itt=itt8
        irstdy=irstdy8
        msrsdy=msrsdy8
        im=iimt8
        jm=ijmt8
        kkm=kkm8

        write (stdout,'(a,i10,1x,a)')
     &  '===> Reading MOM restart from restart.dta on ts=',itt, stamp
        if (im .ne. imt .or. jm .ne. jmt .or. kkm .ne. km) then
          write (stdout,*) '=>Error restio. im=',im,' jm=',jm,' km=',kkm
        endif
c
c-----------------------------------------------------------------------
c       verify that this restart file is the right one. The time step
c       and stamp must match those in the "archive.time" file which the
c       run script should restore upon restarting. (security measure)
c-----------------------------------------------------------------------
c
        ittarc = 0
        ostamp = 'no archive.temp file'
        call getunit (ioarch, 'archive.time'
     &,             'formatted sequential rewind')
        n = 0
        io_status = 0
        do while (n .le. 100000 .and. io_status .eq. 0)
          n = n + 1 
          read (ioarch,8888, IOSTAT=io_status) ittarc, ostamp
          if (io_status .eq. 0) then
	    write (stdout,'(1x,a,1x,i10,1x,a32)')
     &      ' read archive.time record for:', ittarc, ostamp
          endif
        enddo
        if (io_status .gt. 0) then
          write(*,*)
     &     'ERROR in reading archive.time - IOSTAT=',io_status
          STOP
        endif

        call relunit (ioarch)
8888    format (1x,i10,1x,a32)
        if (itt .ne. ittarc .or. stamp .ne. ostamp) then
          if (ittarc .eq. 0) then
            write (stdout,*) '=>Error: could not find "archive.time"'
     &,   ' to verify this is the correct restart.'
            write (stdout,*) ' restart itt=',itt,' restart stamp=',stamp
	  else
	    write (stdout,*) '=>Error: archive.time itt=',ittarc
     &,                  ' but restart has itt = ',itt	  
	    write (stdout,*) '         archive.time stamp=',ostamp
     &,                    ' but restart stamp = ',stamp
	  endif
          write (stdout,*) '=>Error in restio.F'
          call abort()
        endif
      endif
c
#ifdef parallel_1d
      itt8    = itt
      irstdy8 = irstdy
      msrsdy8 = msrsdy
      call mpp_transmit (itt8, 1, ALL_PES, itt8, 1, 0)
      call mpp_transmit (irstdy8, 1, ALL_PES, irstdy8, 1, 0)
      call mpp_transmit (msrsdy8, 1, ALL_PES, msrsdy8, 1, 0)
      call mpp_sync()
      itt    = itt8
      irstdy = irstdy8
      msrsdy = msrsdy8
#endif
c
c-----------------------------------------------------------------------
c     update pointers to tau-1, tau, & tau+1 data on disk based on itt
c     for latitude rows they point to latdisk(1) or latdisk(2)
c     for 2D fields they point to records on kflds
c-----------------------------------------------------------------------
c
      taum1disk = mod(itt+1,2) + 1
      taudisk   = mod(itt  ,2) + 1
      taup1disk = taum1disk
c
c-----------------------------------------------------------------------
c     update pointers to tau-1, tau, & tau+1 data in the MW based on itt
c-----------------------------------------------------------------------
c
#if defined max_window
c
c     rotate time levels instead of moving data
c
      taum1 = mod(itt+0,3) - 1
      tau   = mod(itt+1,3) - 1
      taup1 = mod(itt+2,3) - 1
#endif
c
#if defined rigid_lid_surface_pressure || defined free_surf_on
# ifdef restnosp
c
c-----------------------------------------------------------------------
c     when restarting a prognostic surface pressure configuration from a
c     stream function restart...
c     initialize ubar and vbar from the stream function case
c     first read psi at "tau" then psi at "tau+1"
c     warning: ubar and vbar will be multiplied by "hr" below when
c     "hr" is available
c-----------------------------------------------------------------------
c
      do m=2,1,-1
        if (pn .eq. 1) then
          read (iorest) stamp, iotext
          read (iorest) iimt8,ijmt8,apsi
          iimt=iimt8
          ijmt=ijmt8
          if (iimt .ne. imt .or. ijmt .ne. jmt) then
            write (stdout,*) ' => Error restio #0,imt=',iimt,'
     &,     jmt=',ijmt, ' m=',m
          endif
        endif
# ifdef parallel_1d
        call mpp_transmit
     &     (apsi(1,1), imt*jmt, ALL_PES, apsi(1,1), imt*jmt, 0)
# endif
        do jrow=js,je
	  jp1 = min(jrow+1,jmt)
          do i=2,imtm1
	    ptd(i,jrow) = apsi(i,jrow)
            diag1       = apsi(i+1,jp1) - apsi(i  ,jrow)
            diag0       = apsi(i  ,jp1) - apsi(i+1,jrow)
	    if (m .eq. 1) then
              ubar(i,jrow,1) = -(diag1+diag0)*dyu2r(jrow)
              ubar(i,jrow,2) =  (diag1-diag0)*dxu2r(i)*csur(jrow)
	    else
              ubarm1(i,jrow,1) = -(diag1+diag0)*dyu2r(jrow)
              ubarm1(i,jrow,2) =  (diag1-diag0)*dxu2r(i)*csur(jrow)
	    endif
          enddo
        enddo
      enddo
      call setbcx (ptd(1,js), imt, je-js+1)
      call setbcx (ubar(1,js,1), imt, je-js+1)
      call setbcx (ubar(1,js,2), imt, je-js+1)
      call setbcx (ubarm1(1,js,1), imt, je-js+1)
      call setbcx (ubarm1(1,js,2), imt, je-js+1)
c
c     skip guess fields. then set the surface pressure to zero
c
      if (pn .eq. 1) then
        read (iorest) stamp, iotext
        read (iorest)
        read (iorest) stamp, iotext
        read (iorest)
      endif
#  if defined explicit_free_surface && defined time_averaged_surface
      etat(:,max(1,js-1):min(jmt,je+1),:) = c0
#  endif
      ps(:,js:je,:)                       = c0
# else
c
c     read surface pressure or surface height
c
      if (pn .eq. 1) then
        read (iorest) stamp, iotext
        read (iorest) iimt8,ijmt8,aps
        iimt=iimt8
        ijmt=ijmt8
        if (iimt .ne. imt .or. ijmt .ne. jmt) then
          write (stdout,*) ' => Error restio #A,imt=',iimt,',jmt=',ijmt
        endif
      endif
#   ifdef parallel_1d
      call mpp_transmit
     &     (aps(1,1,1), size(aps), ALL_PES, aps(1,1,1), size(aps), 0)
      call mpp_sync()
#   endif
#   if defined debug_restio
        write (stdout,*) 'restio: tau jrow='
     &,                 'ps cksum=',checksum(aps(1,1,1),imt,jmt)
        write (stdout,*) 'restio: taup1 jrow='
     &,                 'ps cksum=',checksum(aps(1,1,2),imt,jmt)
#    if defined time_averaged_surface        
        write (stdout,*) 'restio: taup1 jrow='
     &,                 'ps cksum=',checksum(aps(1,1,3),imt,jmt)
#    endif
#   endif

#   if defined explicit_free_surface && defined time_averaged_surface
      jsstrt = max(1,js-1)
      jeend  = min(jmt,je+1)
      call extract (etat, jsstrt, jeend, 3, aps, imt, jmt)
#   else
      call extract (ps, js, je, 2, aps, imt, jmt)
#   endif
c
c     read "tau" (ubarm1)
c
      if (pn .eq. 1) then
        read (iorest) stamp, iotext
        read (iorest) iimt8,ijmt8,aps
        iimt=iimt8
        ijmt=ijmt8
        if (iimt .ne. imt .or. ijmt .ne. jmt) then
          write (stdout,*) ' => Error restio #B,imt=',iimt,',jmt=',ijmt
        endif
      endif
#   ifdef parallel_1d
      call mpp_transmit
     &     (aps(1,1,1), size(aps), ALL_PES, aps(1,1,1), size(aps), 0)
      call mpp_sync()
#   endif
#   if defined debug_restio
        write (stdout,*) 'restio: tau jrow='
     &,                 'ubarm1 cksum=',checksum(aps(1,1,1),imt,jmt)
        write (stdout,*) 'restio: tau jrow='
     &,                 'vbarm1 cksum=',checksum(aps(1,1,2),imt,jmt)
#   endif
      call extract (ubarm1, js, je, 2, aps, imt, jmt)
c
c     read  "tau+1" (ubar)
c
      if (pn .eq. 1) then
        read (iorest) stamp, iotext
        read (iorest) iimt8,ijmt8,aps
        iimt=iimt8
        ijmt=ijmt8
        if (iimt .ne. imt .or. ijmt .ne. jmt) then
          write (stdout,*) ' => Error restio #C,imt=',iimt,',jmt=',ijmt
        endif
      endif
#   ifdef parallel_1d
      call mpp_transmit
     &     (aps(1,1,1), size(aps), ALL_PES, aps(1,1,1), size(aps), 0)
      call mpp_sync()
#   endif
#   if defined debug_restio
        write (stdout,*) 'restio: tau jrow='
     &,                 'ubarp1 cksum=',checksum(aps(1,1,1),imt,jmt)
        write (stdout,*) 'restio: tau jrow='
     &,                 'vbarp1 cksum=',checksum(aps(1,1,2),imt,jmt)
#   endif
      call extract (ubar, js, je, 2, aps, imt, jmt)
c
c
#  if !defined explicit_free_surface
c
c     read  guess field
c
      if (pn .eq. 1) then
        read (iorest) stamp, iotext
        read (iorest) iimt8,ijmt8,apsi
        iimt=iimt8
        ijmt=ijmt8
        if (iimt .ne. imt .or. ijmt .ne. jmt) then
          write (stdout,*) ' => Error restio #D,imt=',iimt,',jmt=',ijmt
        endif
      endif
#   ifdef parallel_1d
      call mpp_transmit
     &     (apsi(1,1), imt*jmt, ALL_PES, apsi(1,1), imt*jmt, 0)
      call mpp_sync()
#   endif

      call extract (pguess, js, je, 1, apsi, imt, jmt)
#  endif
c
# endif restnosp
c
c     initialize the kflds disk which is used for surface pressure. 
c
# if !defined explicit_free_surface
      do jrow=js,je
        do i=1,imt
          ptd(i,jrow)    = ps(i,jrow,1)
        enddo
      enddo
# endif
      num_wrds = (je-js+1)*imt
      call oput (kflds, num_wrds, nkflds-1, ps(1,js,1))
      call oput (kflds, num_wrds, nkflds, ps(1,js,1))
c
#endif
c
#ifdef stream_function
c
c
c     read "tau" then "tau+1" stream function
c
      read (iorest) stamp, iotext
      read (iorest) iimt8,ijmt8,apsi
      iimt=iimt8
      ijmt=ijmt8
      if (iimt .ne. imt .or. ijmt .ne. jmt) then
        write (stdout,*)'=> Error restart #1, imt=',iimt,', jmt=',ijmt
      endif
      call extract (psi(1,js,2), js, je, 1, apsi, imt, jmt)
c
      read (iorest) stamp, iotext
      read (iorest) iimt8,ijmt8,apsi
      iimt=iimt8
      ijmt=ijmt8
      if (iimt .ne. imt .or. ijmt .ne. jmt) then
        write (stdout,*)'=> Error restart #2, imt=',iimt,', jmt=',ijmt
      endif
      call extract (psi(1,js,1), js, je, 1, apsi, imt, jmt)
c
c     read the two guess fields
c
      read (iorest) stamp, iotext
      read (iorest) iimt8,ijmt8,apsi
      iimt=iimt8
      ijmt=ijmt8
      if (iimt .ne. imt .or. ijmt .ne. jmt) then
        write (stdout,*)'=> Error restart #3, imt=',iimt,', jmt=',ijmt
      endif
      call extract (ptd(1,js), js, je, 1, apsi, imt, jmt)
c
      num_wrds = (je-js+1)*imt
      call oput (kflds, num_wrds, nkflds-1, ptd(1,js))
c
      read (iorest) stamp, iotext
      read (iorest) iimt8,ijmt8,apsi
      iimt=iimt8
      ijmt=ijmt8

      if (iimt .ne. imt .or. ijmt .ne. jmt) then
        write (stdout,*)'=> Error restart #4, imt=',iimt,', jmt=',ijmt
      endif
      call extract (ptd(1,js), js, je, 1, apsi, imt, jmt)
c
      num_wrds = (je-js+1)*imt
      call oput (kflds, num_wrds, nkflds, ptd(1,js))
#endif
c
c     read "kmt" (the number of ocean T-cells to the ocean floor)
c
      if (pn .eq. 1) then
        read (iorest) stamp, iotext
        read (iorest) iimt8,ijmt8,akmt
        iimt=iimt8
        ijmt=ijmt8
        if (iimt .ne. imt .or. ijmt .ne. jmt) then
          write (stdout,*)'=> Error restart #5, imt=',iimt,', jmt=',ijmt
        endif
      endif
#ifdef parallel_1d
      call mpp_transmit
     &     (akmt(1,1), imt*jmt, ALL_PES, akmt(1,1), imt*jmt, 0)
      call mpp_sync()
#endif
      do jrow=js,je
        if( jrow.EQ.1 .OR. jrow.EQ.jmt )cycle
	do i=1,imt
	  if (kmt(i,jrow) .ne. akmt(i,jrow)) then
	    write (stdout,*) '=>Error: read_restart i=',i,' jrow=',jrow
     &,     'kmt=',kmt(i,jrow),' akmt=',akmt(i,jrow)
            call abort()
	  endif
	enddo
      enddo
#ifdef partial_cell
c
c     read "htp" which is the depth to bottom of partial T cells and
c
      if (pn .eq. 1) then
        read (iorest) stamp, iotext
        read (iorest) iimt8,ijmt8,apsi
        iimt=iimt8
        ijmt=ijmt8
        if (iimt .ne. imt .or. ijmt .ne. jmt) then
          write (stdout,*)'=>Error restart #5a, imt=',iimt,', jmt=',ijmt
        endif
      endif
# ifdef parallel_1d
      call mpp_transmit
     &     (apsi(1,1), imt*jmt, ALL_PES, apsi(1,1), imt*jmt, 0)
      call mpp_sync()
# endif
      if(je+1 .GT. jmt) then
        jeend  = je
      else
        jeend  = je+1
      endif
      call extract (htp(1,js), js, jeend, 1, apsi, imt, jmt)
#endif
c
#if defined nonconst_diffusivity && !defined nonconst_diffusivity_initial
c
c     read "growth2a" which is the squared Eady growth rate
c
      if (pn .eq. 1) then
        read (iorest) stamp, iotext
        read (iorest) iimt8,ijmt8,apsi
        iimt=iimt8
        ijmt=ijmt8
        if (iimt .ne. imt .or. ijmt .ne. jmt) then
          write (stdout,*)'=>Error restart #5b, imt=',iimt,', jmt=',ijmt
        endif
      endif
# ifdef parallel_1d
      call mpp_transmit
     &     (apsi(1,1), imt*jmt, ALL_PES, apsi(1,1), imt*jmt, 0)
      call mpp_sync()
# endif
      call extract (growth2a(1,js), js, je, 1, apsi, imt, jmt)
c
c     read "agm" which is the nonconstant Redi/GM diffusivity
c      
      if (pn .eq. 1) then
        read (iorest) stamp, iotext
        read (iorest) iimt8,ijmt8,apsi
        iimt=iimt8
        ijmt=ijmt8
        if (iimt .ne. imt .or. ijmt .ne. jmt) then
          write (stdout,*)'=>Error restart #5c, imt=',iimt,', jmt=',ijmt
        endif
      endif
# ifdef parallel_1d
      call mpp_transmit
     &     (apsi(1,1), imt*jmt, ALL_PES, apsi(1,1), imt*jmt, 0)
      call mpp_sync()
# endif
      call extract (agm(1,js), js, je, 1, apsi, imt, jmt)
c
#endif
#ifdef parallel_1d
      call mpp_sync()
#endif
c
#if defined rigid_lid_surface_pressure || defined free_surf_on
# ifdef restnosp
c
c     complete external mode calculation using "hr"
c
      do m=1,2
        do jrow=js,je
	  do i=1,imt
	    ubar(i,jrow,m)   = ubar(i,jrow,m)*hr(i,jrow)
	    ubarm1(i,jrow,m) = ubarm1(i,jrow,m)*hr(i,jrow)
	  enddo
	enddo
	call setbcx (ubar(1,js,m), imt, je-js+1)
	call setbcx (ubarm1(1,js,m), imt, je-js+1)
      enddo        
# endif
#endif
c
c-----------------------------------------------------------------------
c     read the "tau" latitude rows
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
#if defined max_window
	j = jrow + jbuf - 1
#else
	j = jmw
#endif
c
        if (pn .eq. 1) then
          read (iorest) stamp, iotext
          read (iorest,IOSTAT=io_status) ilensl8,bufsl
          ilensl=ilensl8
          if (ilensl .ne. nslab) then
            write (stdout,*) ' => Error rest #6, lensl=',ilensl
     &,                         ',jrow=',jrow
            write (stdout,*)
     &     'possibly "nt" on restart differs from parameter nt =',nt
            call abort()
          endif
#ifdef tcvmix
c         How are the tcvmix variables read in???????
          uu = 1.0/undef
	  uu = uu**2
#endif
        endif
#ifdef parallel_1d
        if( pe .eq. 0 ) then
            if (pe_of_lat(jrow) .ne. 0) call mpp_transmit (bufsl(1,1,1)
     &,         size(bufsl), pe_of_lat(jrow), bufsl(1,1,1), size(bufsl)
     &,          NULL_PE )
        elseif (pe .eq. pe_of_lat(jrow)) then
            call mpp_transmit (bufsl(1,1,1), size(bufsl), NULL_PE
     &,          bufsl(1,1,1),size(bufsl), 0)
        endif
        call mpp_sync()
#endif

#if defined debug_restio
        write (stdout,*) 'restio: tau jrow=',jrow
     &,                 'u cksum=',checksum(bufsl(1,1,nvar-1),imt,km)
        write (stdout,*) 'restio: tau jrow=',jrow
     &,                 'v cksum=',checksum(bufsl(1,1,nvar),imt,km)
        write (stdout,*) 'restio: tau jrow=',jrow
     &,                 't cksum=',checksum(bufsl(1,1,1),imt,km)
#endif

c
        if (jrow .ge. js .and. jrow .le. je) then
c
c         velocity is internal mode component only
c         note that tracers are stored first, then velocities on restart
c         for compatibility with previous versions
c
	  do k=1,km
	    do i=1,imt
	      u(i,k,j,1,tau) = bufsl(i,k,nvar-1)
	      u(i,k,j,2,tau) = bufsl(i,k,nvar)
	      do n=1,nvar-2
	        t(i,k,j,n,tau) = bufsl(i,k,n)
	      enddo
	    enddo
	  enddo
c
c         initialize every latitude
c
#if defined max_window
          if (jrow .eq. 1) then
c
c           when jrow=jmt, initialize buffer rows below jrow=1
c
            do jj=1,jbuf
	      do k=1,km
	        do i=1,imt
	          do n=1,2
	            u(i,k,jj,n,taum1) = c0
	            u(i,k,jj,n,tau)   = c0
	            u(i,k,jj,n,taup1) = c0
	          enddo
	          do n=1,nvar-2
	            t(i,k,jj,n,taum1) = c0
	            t(i,k,jj,n,tau)   = c0
	            t(i,k,jj,n,taup1) = c0
	          enddo
	          tmask(i,k,jj) = c0
	          umask(i,k,jj) = c0
	        enddo
	      enddo
	    enddo
	  elseif (jrow .eq. jmt) then
c
c           when jrow=jmt, initialize buffer rows above jrow=jmt
c
            do jj=jmw-jbuf+1,jmw
	      do k=1,km
	        do i=1,imt
	          do n=1,2
	            u(i,k,jj,n,taum1) = c0
	            u(i,k,jj,n,tau)   = c0
	            u(i,k,jj,n,taup1) = c0
	          enddo
	          do n=1,nvar-2
	            t(i,k,jj,n,taum1) = c0
	            t(i,k,jj,n,tau)   = c0
	            t(i,k,jj,n,taup1) = c0
	          enddo
	          tmask(i,k,jj) = c0
	          umask(i,k,jj) = c0
	        enddo
	      enddo
	    enddo
	  else
c
c           do nothing since "tau" data is in place in the MW
c
	  endif
#else
          call putrow (latdisk(taudisk), nslab, jrow, u(1,1,j,1,tau)
     &,                                               t(1,1,j,1,tau))
#endif
        endif
      enddo
c
c
c-----------------------------------------------------------------------
c     read the "tau+1" latitude rows
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
c
#if defined max_window
        j = jrow + jbuf - 1
#else
	j = jmw
#endif
c
        if (pn .eq. 1) then
          read (iorest) stamp, iotext
          read (iorest) ilensl8,bufsl
          ilensl=ilensl8
          if (ilensl .ne. nslab) then
            write (stdout,*) ' => Error restart #7, ilensl=',ilensl
     &,                    ',jrow=',jrow
          endif
#ifdef tcvmix
c         How are the tcvmix variables read in???????
          uu = 1.0/undef
	  uu = uu**2
#endif
        endif
#ifdef parallel_1d
        if ( pe .eq. 0 )then
          if (pe_of_lat(jrow) .ne. 0) call mpp_transmit (bufsl(1,1,1)
     &,      size(bufsl), pe_of_lat(jrow), bufsl(1,1,1), size(bufsl)
     &,      NULL_PE)
        else if (pe .eq. pe_of_lat(jrow)) then
            call mpp_transmit (bufsl(1,1,1), size(bufsl), NULL_PE
     &,          bufsl(1,1,1),size(bufsl), 0)
        endif
        call mpp_sync()
#endif

#if defined debug_restio
        write (stdout,*) 'restio: taup1 jrow=',jrow
     &,                 'u cksum=',checksum(bufsl(1,1,nvar-1),imt,km)
        write (stdout,*) 'restio: taup1 jrow=',jrow
     &,                 'v cksum=',checksum(bufsl(1,1,nvar),imt,km)
        write (stdout,*) 'restio: taup1 jrow=',jrow
     &,                 't cksum=',checksum(bufsl(1,1,1),imt,km)
#endif
c
        if (jrow .ge. js .and. jrow .le. je) then
c
c         velocity is internal mode component only
c         note that tracers are stored first, then velocities on restart
c         for compatibility with previous versions
c
	  do k=1,km
	    do i=1,imt
	      u(i,k,j,1,taup1) = bufsl(i,k,nvar-1)
	      u(i,k,j,2,taup1) = bufsl(i,k,nvar)
	      do n=1,nvar-2
	        t(i,k,j,n,taup1) = bufsl(i,k,n)
	      enddo
	    enddo
	  enddo
c
c         initialize every latitude
c
#if defined max_window
c
c         do nothing since "tau+1" data is in place in the MW
c
#else
          call putrow (latdisk(taup1disk), nslab, jrow
     &,                u(1,1,j,1,taup1), t(1,1,j,1,taup1))
#endif
        endif
      enddo
c
#ifdef pressure_depth
      do jrow=1,jmt
c
c       read rho
c
        if (pn .eq. 1) then
          read (iorest) stamp, iotext
          read (iorest, IOSTAT=io_status) bufr
        endif
# ifdef parallel_1d
        if( pe .eq. 0 ) then
            if (pe_of_lat(jrow) .ne. 0) call mpp_transmit (bufr(1,1)
     &,         size(bufr), pe_of_lat(jrow), bufr(1,1), size(bufr)
     &,          NULL_PE )
        elseif (pe .eq. pe_of_lat(jrow)) then
            call mpp_transmit (bufr(1,1), size(bufr), NULL_PE
     &,          bufr(1,1),size(bufr), 0)
        endif
        call mpp_sync()
# endif
        if (jrow .ge. js .and. jrow .le. je) then
	  rho(:,:,jrow-js+1+jbuf) = bufr(:,:)
	endif
      enddo
#endif
c
#ifdef trajectories
      read (iorest, err=210) stamp, iotext
      read (iorest)
      read (iorest) stamp, iotext
# ifdef lyapunov
      read (iorest,err=205) npart8,pxyz,pijk,em
      npart=npart8
# else
      read (iorest) npart8,pxyz,pijk
      npart=npart8
# endif
      write (stdout, '(a,i8,a)') 
     &'===> also read trajectory restart for ',npart,' particles'
      if (npart .ne. nptraj) then
        write (stdout,*) ' => Particle trajectory restart error:'
        write (stdout,*) '    read ',npart,' particles but expected '
     &,                     nptraj 
        call abort()
      endif
      go to 220
# ifdef lyapunov
205   continue
      write (stdout,'(/a/a/)')
     & '==>Note: Did not find lyapunov part of trajectory restart...'
     &,'         Setting all particles to Initial Positions'   
# endif
210   continue
c
c     if no restart data is available. initialize the particles
c
      call ptraji
220   continue
#endif
      if (pn .eq. 1) then
        call relunit (iorest)
      endif
#ifdef timing
      call toc ('ocean', 'reading restart')
#endif
#ifdef parallel_1d
      call mpp_sync()
#endif
c
      return
      end



      subroutine write_restart (js, je)
c
c=======================================================================
c     write the restart file
c     
c     author:        a. rosati             e-mail=> ar@gfdl.gov
c                    r. c. pacanowski      e-mail=> rcp@gfdl.gov
c    (mpp items)     v. balaji             e-mail=> vb@gfdl.gov
c=======================================================================
c
      use mom3_mpp_mod
      use arrays_2d
      use time_manager_mod
      use memory_window
      character*128 file_name, file_stamp
#include "param.h"
#include "coord.h"
#ifdef trajectories
# include "diag_ptraj.h"
#endif
#include "emode.h"
#include "grdvar.h"
#include "iounit.h"
#ifdef nonconst_diffusivity
# include "isopyc.h"
#endif
#include "mw.h"
#include "taskrows.h"
#include "tmngr.h"
      integer(kind=8) :: akmt
      common /restbuf/ apsi(imt,jmt)
      dimension akmt(imt,jmt)
#if defined rigid_lid_surface_pressure || defined free_surf_on
# if defined time_averaged_surface
      common /restbuf/ aps(imt,jmt,3)
# else
      common /restbuf/ aps(imt,jmt,2)
# endif
#endif
c
#if !defined partial_cell
      dimension hup(1,1), htp(1,1)
#endif
      character*32 ostamp, seqrecl
#if !defined sgi_mipspro || defined use_libMPI
      common / restbuf / bufsl(imt,km,nvar)
# if defined pressure_depth
      common / restbuf / bufr(imt,km)
# endif
#else
      dimension bufsl(imt,km,nvar)
# if defined pressure_depth
      dimension bufsr(imt,km)
# endif
# if defined parallel_1d
      pointer( ptrsl, bufsl )
      integer, save :: lensl=0
      save :: ptrsl
#  if defined pressure_depth
      pointer( ptrr, bufr )
      integer, save :: lenr=0
      save :: ptrr
#  endif
# endif
#endif
      
#ifdef stream_function
      dimension ext(imt,2)
#endif

      integer(kind=8) :: itt8, irstdy8, msrsdy8, imt8, jmt8, km8
     &,                  nslab8,nptraj8
c
c-----------------------------------------------------------------------
c     write restart file using unformatted fortran i/o
c
c     record 1 on the restart file is the control block (containing
c                the time step and relative year)
c     record`s 2 & 3 are for stream function time levels
c     record`s 4 & 5 are for the d(stream function)/dt fields used
c                to construct a guess for the elliptic solver
c     record 6 is for "kmt" indicating the number of ocean T cells
c                from the surface to the ocean floor.
#if defined partial_cell
c     record 7 is for "htp" depth to bottom of partial T-cells
#endif
#if defined nonconst_diffusivity
c     record 8 is for "growth2a" squared Eady growth rate sec^-2.
c     record 9 is for "agm" nonconstant Redi/GM diffusivity.
#endif
c     the next "jmt" records contain prognostic variables at "tau"
c     the last "jmt" records contain prognostic variables at "tau+1"
c-----------------------------------------------------------------------
c
#ifdef timing
      call tic ('ocean', 'writing restart')
#endif
#ifdef parallel_1d
# if defined sgi_mipspro && !defined use_libMPI 
      call mpp_malloc( ptrsl, size(bufsl), lensl )
#  ifdef pressure_depth
      call mpp_malloc( ptrr, size(bufr), lenr )
#  endif
# endif
#endif
      apsi(:,:) = 0.
#if defined rigid_lid_surface_pressure || defined free_surf_on
      aps(:,:,:) = 0.
#endif
c
      if (pn .eq. 1) then
c
c       reduce size of old restart.dta to save space
c
        file_name = 'restart.dta'
        call getunit (iorest, trim(file_name)
     &,             'unformatted sequential rewind')
        write (iorest) 'restart.dta removed'
        call relunit (iorest)
c
c       open new restart file with date stamp on filename
c	
        write (seqrecl,'(a,i12)') ' words=',nvar*imt*km
        file_name = file_stamp ('restart', xstamp, '.dta')
        call getunit (iorest, trim(file_name)
     &,             'unformatted sequential rewind'//seqrecl)
c
        write (stdout,'(a,i8,1x,a)')
     & '==> Writing MOM restart to file restart.dta at ts=', itt, stamp
c
c       note that timestep "itt" corresponds to "tau+1"
c
        iotext = ' read (iorest) itt, irstdy, msrsdy, imt, jmt, km'
        write (iorest) stamp, iotext, expnam
        call get_time (model_time, isec, irstdy)
	msrsdy = 1000*isec
        itt8=itt
        irstdy8=irstdy
        msrsdy8=msrsdy
        imt8=imt
        jmt8=jmt
        km8=km
        write (iorest) itt8, irstdy8, msrsdy8, imt8, jmt8, km8
      endif
c
#if defined rigid_lid_surface_pressure || defined free_surf_on
c
# ifdef parallel_1d
      call mpp_sync()
# endif
c
c     ps or etat are at (tau,tau+1) here
c
      aps(:,:,:) = 0.
# if defined explicit_free_surface && defined time_averaged_surface
      aps(:,js:je,:) = etat(:,js:je,:)
# else
      aps(:,js:je,:) = ps(:,js:je,:)
# endif
# ifdef parallel_1d
      call mpp_update_domains (aps, mom_global_domain)
# endif
# if defined debug_restio
        write (stdout,*) 'restio: tau jrow='
     &,                 'ps cksum=',checksum(aps(1,1,1),imt,jmt)
        write (stdout,*) 'restio: taup1 jrow='
     &,                 'ps cksum=',checksum(aps(1,1,2),imt,jmt)
#  if defined time_averaged_surface
        write (stdout,*) 'restio: taup1 jrow='
     &,                 'ps cksum=',checksum(aps(1,1,3),imt,jmt)
#  endif
# endif 
      if (pn .eq. 1) then
# if defined explicit_free_surface && defined time_averaged_surface
        iotext = 
     &' read (iorest) imt,jmt,((etat(i,j,1),i=1,imt),j=1,jmt),etat(,,2)
     &       ,etat(,,3)'
# else
        iotext = 
     &' read (iorest) imt,jmt, ((ps(i,j,1),i=1,imt),j=1,jmt), ps(,,2)'
# endif
        write (iorest) stamp, iotext, expnam
        imt8=imt
        jmt8=jmt
        write (iorest) imt8, jmt8, aps
      endif
c
c     ubarm1 is at tau here
c
       aps(:,:,:) = 0.
       aps(:,js:je,1:2) = ubarm1(:,js:je,1:2)
# ifdef parallel_1d
      call mpp_update_domains (aps, mom_global_domain)
# endif
# if defined debug_restio
        write (stdout,*) 'restio: tau jrow='
     &,                 'ubarm1 cksum=',checksum(aps(1,1,1),imt,jmt)
        write (stdout,*) 'restio: tau jrow='
     &,                 'vbarm1 cksum=',checksum(aps(1,1,2),imt,jmt)
# endif
      if (pn .eq. 1) then
        iotext = 
     &' read (iorest) imt,jmt,((ubarm1(i,j),i=1,imt),j=1,jmt), vbarm1()'
        write (iorest) stamp, iotext, expnam
        imt8=imt
        jmt8=jmt
        write (iorest) imt8, jmt8, aps
      endif
c
c     ubar is at tau+1 here
c
      aps(:,:,:) = 0.
      aps(:,js:je,1:2) = ubar(:,js:je,1:2)
# ifdef parallel_1d
      call mpp_update_domains (aps, mom_global_domain)
# endif
# if defined debug_restio
        write (stdout,*) 'restio: taup1 jrow='
     &,                 'ubar cksum=',checksum(aps(1,1,1),imt,jmt)
        write (stdout,*) 'restio: taup1 jrow='
     &,                 'vbar cksum=',checksum(aps(1,1,2),imt,jmt)
# endif
      if (pn .eq. 1) then
        iotext = 
     &' read (iorest) imt, jmt, ((ubar(i,j),i=1,imt),j=1,jmt), vbar()'
        write (iorest) stamp, iotext, expnam
        imt8=imt
        jmt8=jmt
        write (iorest) imt8, jmt8, aps
      endif
c
# if !defined explicit_free_surface
c
c     guess field
c

      aps(:,:,:) = 0.
      apsi(:,js:je) = pguess(:,js:je)
#  ifdef parallel_1d
      call mpp_update_domains (apsi, mom_global_domain)
#  endif
      if (pn .eq. 1) then
        iotext = 
     &' read (iorest) imt, jmt, ((pguess(i,j),i=1,imt),j=1,jmt)'
        write (iorest) stamp, iotext, expnam
        imt8=imt
        jmt8=jmt
        write (iorest) imt8, jmt8, apsi
      endif
# endif
#endif
c
#ifdef stream_function
c
c     psi at "tau"
c

      apsi(:,:) = 0.
      apsi(:,js:je) = psi(:,js:je,2)
# ifdef parallel_1d
      call mpp_update_domains (apsi, mom_global_domain)
# endif
# if defined debug_restio
      write (stdout,*) 'restio: tau jrow='
     &,                 'psi cksum=',checksum(apsi(1,1),imt,jmt)
# endif
      iotext = ' read (iorest) imt,jmt,((psi(i,j,2),i=1,imt),j=1,jmt)'
      write (iorest) stamp, iotext, expnam
      imt8=imt
      jmt8=jmt
      write (iorest) imt8, jmt8, apsi
c
c     psi at "tau+1"
c

      apsi(:,:) = 0.
      apsi(:,js:je) = psi(:,js:je,2)
# ifdef parallel_1d
      call mpp_update_domains (apsi, mom_global_domain)
# endif
# if defined debug_restio
      write (stdout,*) 'restio: taup1 jrow='
     &,                 'psi cksum=',checksum(apsi(1,1),imt,jmt)
# endif
      iotext = ' read (iorest) imt,jmt,((psi(i,j,1),i=1,imt),j=1,jmt)'
      write (iorest) stamp, iotext, expnam
      imt8=imt
      jmt8=jmt
      write (iorest) imt8, jmt8, apsi
c
c     guess fields
c
      call oget (kflds, nwds, nkflds-1, ptd)
      iotext = ' read (iorest) imt, jmt, ((g1(i,j),i=1,imt),j=1,jmt)'
      write (iorest) stamp, iotext, expnam
      write (iorest) imt8, jmt8, ptd

c
      call oget (kflds, nwds, nkflds, ptd)
      iotext = ' read (iorest) imt, jmt, ((g2(i,j),i=1,imt),j=1,jmt)'
      write (iorest) stamp, iotext, expnam
      imt8=imt
      jmt8=jmt
      write (iorest) imt8, jmt8, ptd
#endif
c
c     "kmt"
c
      akmt(:,:) = 0.0
      do jrow=js,je
	do i=1,imt
	  akmt(i,jrow) = kmt(i,jrow)
	enddo
      enddo
c
#ifdef parallel_1d
!     use apsi instead of akmt, since mpp_update_domains expects reals
      call mpp_update_domains (akmt, mom_global_domain)
      call mpp_sync()
#endif
      if (pn .eq. 1) then
        iotext = ' read (iorest) imt, jmt, ((kmt(i,j),i=1,imt),j=1,jmt)'
        write (iorest) stamp, iotext, expnam
        imt8=imt
        jmt8=jmt
        write (iorest) imt8, jmt8, akmt
      endif
c
#ifdef partial_cell
c
c     "htp"
c
      apsi(:,:) = 0.
      apsi(:,js:je) = htp(:,js:je)
#ifdef parallel_1d
      call mpp_update_domains (apsi, mom_global_domain)
#endif
      if (pn .eq. 1) then
        iotext = ' read (iorest) imt, jmt, ((htp(i,j),i=1,imt),j=1,jmt)'
        write (iorest) stamp, iotext, expnam
        imt8=imt
        jmt8=jmt
        write (iorest) imt8, jmt8, apsi
      endif
#endif
c
#if defined nonconst_diffusivity && !defined nonconst_diffusivity_final
c
c     "growth2a"
c

      apsi(:,:) = 0.
      apsi(:,js:je) = growth2a(:,js:je)
#ifdef parallel_1d
      call mpp_update_domains (apsi, mom_global_domain)
#endif
      if (pn .eq. 1) then
      iotext=' read (iorest) imt, jmt,((growth2a(i,j),i=1,imt),j=1,jmt)'
        write (iorest) stamp, iotext, expnam
        imt8=imt
        jmt8=jmt
        write (iorest) imt8, jmt8, apsi
      endif
c
c     "agm"
c
      apsi(:,:) = 0.
      apsi(:,js:je) = agm(:,js:je)
#ifdef parallel_1d
      call mpp_update_domains( apsi, mom_global_domain )
#endif
      if (pn .eq. 1) then
      iotext=' read (iorest) imt, jmt,((agm(i,j),i=1,imt),j=1,jmt)'
        write (iorest) stamp, iotext, expnam
        imt8=imt
        jmt8=jmt
        write (iorest) imt8, jmt8, apsi
      endif
#endif
c
c-----------------------------------------------------------------------
c     save the "tau" latitude rows
c     (disk pointers need re-setting to handle euler backward)
c-----------------------------------------------------------------------
c
      taum1disk = mod(itt+1,2) + 1
      taudisk   = mod(itt  ,2) + 1
      taup1disk = taum1disk
c
      do jrow=1,jmt
        bufsl(:,:,:) = c0
c
#if defined max_window
	j = jrow + jbuf - 1
c
c       remove external mode from "tau". since psi has been updated
c       psi(,,2) is at "tau" 
c
        if (jrow .ge. js .and. jrow .le. je) then
# ifdef stream_function
#  if !defined robert_time_filter
	  do i=2,imt-1
            diag1    = psi(i+1,jrow+1,2) - psi(i  ,jrow,2)
            diag0    = psi(i  ,jrow+1,2) - psi(i+1,jrow,2)
	    ext(i,1) = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
	    ext(i,2) =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)
     &                       *csur(jrow)
          enddo
#  endif
	  do k=1,km
            do i=2,imt-1
              if (k .le. kmu(i,jrow)) then
#  if !defined robert_time_filter
                bufsl(i,k,nvar-1) = (u(i,k,j,1,tau) - ext(i,1))
                bufsl(i,k,nvar)   = (u(i,k,j,2,tau) - ext(i,2))
#  else
                bufsl(i,k,nvar-1) = u(i,k,j,1,tau)
                bufsl(i,k,nvar)   = u(i,k,j,2,tau)
#  endif
              else
                bufsl(i,k,nvar-1) = c0
                bufsl(i,k,nvar)   = c0
              endif
            enddo
	  enddo
# endif
# if defined rigid_lid_surface_pressure || defined free_surf_on
c
c         since external mode has been updated, ubarm1 is at "tau"
c
	  do k=1,km
            do i=2,imt-1
              if (k .le. kmu(i,jrow)) then
#  if !defined robert_time_filter
                bufsl(i,k,nvar-1) = (u(i,k,j,1,tau)-ubarm1(i,jrow,1))
                bufsl(i,k,nvar)   = (u(i,k,j,2,tau)-ubarm1(i,jrow,2))
#  else
                bufsl(i,k,nvar-1) = u(i,k,j,1,tau)
                bufsl(i,k,nvar)   = u(i,k,j,2,tau)
#  endif
              else
                bufsl(i,k,nvar-1) = c0
                bufsl(i,k,nvar)   = c0
              endif
            enddo
          enddo
# endif
	  call setbcx (bufsl(1,1,nvar-1), imt, km)
	  call setbcx (bufsl(1,1,nvar), imt, km)
        endif
c
#else
c
	j = jmw
        if (jrow .ge. js .and. jrow .le. je) then
          call getrow (latdisk(taudisk), nslab, jrow, u(1,1,j,1,tau)
     &,                                               t(1,1,j,1,tau))
	  do k=1,km
	    do i=1,imt
	      bufsl(i,k,nvar-1) = u(i,k,j,1,tau)
	      bufsl(i,k,nvar)   = u(i,k,j,2,tau)
	    enddo
	  enddo
	endif
#endif
c
        if (jrow .ge. js .and. jrow .le. je) then
	  do k=1,km
	    do i=1,imt
	      do n=1,nvar-2
	        bufsl(i,k,n) = t(i,k,j,n,tau)
	      enddo
	    enddo
	  enddo
	endif

#ifdef parallel_1d
        call mpp_sync()
        if (pe .eq. 0) then
            if (pe_of_lat(jrow) .ne. 0) call mpp_transmit (bufsl(1,1,1),
     &           size(bufsl), NULL_PE, bufsl(1,1,1),
     &           size(bufsl), pe_of_lat(jrow))
        else if (pe .eq. pe_of_lat(jrow)) then
            call mpp_transmit (bufsl(1,1,1), size(bufsl), 0,
     &           bufsl(1,1,1), size(bufsl), NULL_PE)
        endif
        call mpp_sync()
#endif

        if (jrow.eq.1.or.jrow.eq.jmt) bufsl(:,:,:)=0.0

#if defined debug_restio
        write (stdout,*) 'restio: tau jrow=',jrow
     &,                 'u cksum=',checksum(bufsl(1,1,nvar-1),imt,km)
        write (stdout,*) 'restio: tau jrow=',jrow
     &,                 'v cksum=',checksum(bufsl(1,1,nvar),imt,km)
        write (stdout,*) 'restio: tau jrow=',jrow
     &,                 't cksum=',checksum(bufsl(1,1,1),imt,km)
#endif

c
        if (pn .eq. 1) then
          write (iotext,'(a15,i4)') ' taudisk jrow =',jrow
          iotext(20:) = ' read(iorest) lensl, (bufsl(i),i=1,lensl)'
          write (iorest) stamp, iotext, expnam
          nslab8=nslab
          write (iorest) nslab8, bufsl
        endif
      enddo
c
c-----------------------------------------------------------------------
c     save the "tau+1" latitude rows
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        bufsl(:,:,:) = c0
c
#if defined max_window
	j = jrow + jbuf - 1
c
c       do nothing since data is already in "tau+1" position in MW
c       and velocity is already internal mode only
c
#else
	j = jmw
        if (jrow .ge. js .and. jrow .le. je) then
          call getrow (latdisk(taup1disk), nslab, jrow
     &,                  u(1,1,j,1,taup1), t(1,1,j,1,taup1))
	endif
#endif
c
        if (jrow .ge. js .and. jrow .le. je) then
	  do k=1,km
	    do i=1,imt
	      bufsl(i,k,nvar-1) = u(i,k,j,1,taup1)
	      bufsl(i,k,nvar) = u(i,k,j,2,taup1)
	      do n=1,nvar-2
	        bufsl(i,k,n) = t(i,k,j,n,taup1)
	      enddo
	    enddo
	  enddo
	endif

#ifdef parallel_1d
        call mpp_sync()
        if (pe .eq. 0) then
            if (pe_of_lat(jrow) .ne. 0)call mpp_transmit( bufsl(1,1,1),
     &           size(bufsl), NULL_PE, bufsl(1,1,1),
     &           size(bufsl), pe_of_lat(jrow))
        else if (pe .eq. pe_of_lat(jrow)) then
            call mpp_transmit( bufsl(1,1,1), size(bufsl), 0,
     &           bufsl(1,1,1), size(bufsl), NULL_PE)
        endif
        call mpp_sync()
#endif

#if defined debug_restio
        write (stdout,*) 'restio: taup1 jrow=',jrow
     &,                 'u cksum=',checksum(bufsl(1,1,nvar-1),imt,km)
        write (stdout,*) 'restio: taup1 jrow=',jrow
     &,                 'v cksum=',checksum(bufsl(1,1,nvar),imt,km)
        write (stdout,*) 'restio: taup1 jrow=',jrow
     &,                 't cksum=',checksum(bufsl(1,1,1),imt,km)
#endif
c
        if (pn .eq. 1) then
          write (iotext,'(a17,i4)') ' taup1disk jrow =',jrow
          iotext(22:) = ' read(iorest) lensl, (bufsl(i),i=1,lensl)'
          write (iorest) stamp, iotext, expnam
          nslab8=nslab
          write (iorest) nslab8, bufsl
	endif
      enddo
c
#ifdef pressure_depth
      do jrow=1,jmt
        bufr(:,:) = c0
        if (jrow .ge. js .and. jrow .le. je) then
	
          bufr(:,:) = rho(:,:,jrow-js+1+jbuf)
# ifdef parallel_1d
          call mpp_sync()
          if (pe .eq. 0) then
            if (pe_of_lat(jrow) .ne. 0) call mpp_transmit (bufr(1,1),
     &           size(bufr), NULL_PE, bufr(1,1),
     &           size(bufr), pe_of_lat(jrow))
          else if (pe .eq. pe_of_lat(jrow)) then
            call mpp_transmit (bufr(1,1), size(bufr), 0,
     &           bufr(1,1), size(bufr), NULL_PE)
          endif
          call mpp_sync()
# endif
c
          if (pn .eq. 1) then
            write (iotext,'(a15,i4)') ' rho jrow =',jrow
            iotext(20:) = ' read(iorest) (bufr(i),i=1,imt*km)'
            write (iorest) stamp, iotext, expnam
            write (iorest) bufr
          endif
        endif
      enddo
#endif
c
#ifdef trajectories
      reltim = relyr
      iotext = 'read (iorest) reltim'
      write (iorest) stamp, iotext, expnam
      write (iorest) reltim
# ifdef lyapunov
      iotext ='read (iorest) nptraj, pxyz, pijk, em'
      write (iorest) stamp, iotext, expnam
      nptraj8=nptraj
      write (iorest) nptraj8, pxyz, pijk, em
# else
      iotext ='read (iorest) nptraj, pxyz, pijk'
      write (iorest) stamp, iotext, expnam
      nptraj8=nptraj
      write (iorest) nptraj8, pxyz, pijk
# endif
      write (stdout, '(a,i8,a)')
     &'===> also writing trajectory restart for ',nptraj, 'particles'
#endif
c
c-----------------------------------------------------------------------
c     write time step and stamp to file archive.time for purpose of
c     verifying that the restart file "restart.dta" is the correct one
c     when beginning the next job. The run script should save this file
c     and restore it when restarting. (security measure)
c-----------------------------------------------------------------------
c     
      if (pn .eq. 1) then
        call getunit (ioarch, 'archive.time'
     &,             'formatted sequential append')
        write (ioarch,8888) itt, stamp
        call relunit (ioarch)
c
        call relunit (iorest)
      endif
#ifdef parallel_1d
      call mpp_sync()
#endif
c
#ifdef timing
      call toc ('ocean', 'writing restart')
#endif
c
      return
8888  format (1x,i10,1x,a32)
      end


      subroutine extract (a, js, je, itm, b, im, jm )
c
c     extract "a" from "b" on each processor
c
      dimension a(im,js:je,itm), b(im,jm,itm)
      do it=1,itm
        do jrow=js,je
	  do i=1,im
	    a(i,jrow,it) = b(i,jrow,it)
	  enddo
        enddo
      enddo
      return
      end
