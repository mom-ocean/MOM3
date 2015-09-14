#if defined ramdrive
      subroutine ostart (lu, fname, ntot, nwrec, nbuf)
c
c=======================================================================
c
c     ramdrive section uses memory as a ramdisk
c
c     check for validity of unitnumber and show memory requirement
c
c     lu    = i/o unit number associated with "fname"
c     fname = filename 
c     ntot  = total length of unit, in words
c     nwrec = number of words in the record on the unit
c     nbuf  = number of buffers for the unit
c=======================================================================
c
      character*(*) fname
# include "param.h"
# include "iounit.h"
      if (lu .ne. kflds .and. lu .ne. latdisk(1) .and. 
     &                        lu .ne. latdisk(2)) then
        write (stdout,'(/1x,a,/,11x,a,/11x,a,i10)')
     &  '=>Error: ostart... when using the "ramdrive" option, only'
     &,' units "kflds","latdisk(1)","latdisk(2)" are              '
     &,' allowed and you were trying to open unit #',lu 
	call abort()
      endif
      write (stdout,'(/a,a,a,1pg10.3,a)')
     & ' Memory resource for file ',trim(fname),' = ',ntot*1.e-6,' (MW)'
c
c     attach "lu" to "fname" to remove "lu" from the list of possible
c     unit numbers. This is just a dummy file.
c
      call getunit (lu, trim(fname), 'sequential, unformatted, rewind')
      return
      end
c
      subroutine ofind (lu, nwrs, nrec)
      return
      end
c
      subroutine oclose (lu)
      return
      end


      subroutine getrow (lu, nwrs, jrow, u, t)
c
c=======================================================================
c     get velocity and tracer data from ramdisk latitude row "jrow"
c
c         lu    = i/o unit number associated with "fname"
c         nwrs  = number of words to read
c         jrow  = latitude row
c         u     = velocities
c         t     = tracers
c         ntrac = number of tracers
c=======================================================================
c
      use datamod
# include "param.h"
# include "iounit.h"
# include "taskrows.h"
c
      parameter (ntrac=nvar-2)
      dimension u(imt,km,jmw,2), t(imt,km,jmw,ntrac)
c
      if (lu .eq. latdisk(1) .or. lu .eq. latdisk(2)) then
	if (lu .eq. latdisk(1)) then
          ntl = 1
	elseif (lu .eq. latdisk(2)) then
          ntl = 2
	endif
        call getlat (ramdisk3d(1,jrow,ntl), u, t)
      else
        write (stdout,*) '=>Error: unit ',lu,' not allowed. jrow=', jrow
	call abort()
      endif
      return
      end

      subroutine putrow (lu, nwrs, jrow, u, t)
c
c=======================================================================
c     put velocity and tracer data to ramdisk latitude row "jrow"
c
c         lu    = i/o unit number associated with "fname"
c         nwrs  = number of words to read
c         jrow  = latitude row
c         u     = velocities
c         t     = tracers
c         ntrac = number of tracers
c=======================================================================
c
      use datamod
# include "param.h"
# include "iounit.h"
# include "taskrows.h"
c
      parameter (ntrac=nvar-2)
      dimension u(imt,km,jmw,2), t(imt,km,jmw,ntrac)
c
      if (lu .eq. latdisk(1) .or. lu .eq. latdisk(2)) then
	if (lu .eq. latdisk(1)) then
          ntl = 1
	elseif (lu .eq. latdisk(2)) then
          ntl = 2
	endif
        call putlat (ramdisk3d(1,jrow,ntl), u, t)
      else
        write (stdout,*) '=>Error: unit ',lu,' not allowed. jrow=', jrow
	call abort()
      endif
      return
      end

      subroutine getlat (disk, u, t)
c
c-----------------------------------------------------------------------
c     copy 3-D data from ramdisk to 3D array in memory.
c-----------------------------------------------------------------------
c
#include "param.h"
      parameter (ntrac=nvar-2)
      dimension u(imt,km,jmw,2), t(imt,km,jmw,ntrac)
      dimension disk(imt,km,nvar)
      do n=1,2
        do k=1,km
	  do i=1,imt
	    u(i,k,1,n) = disk(i,k,n)
	  enddo
	enddo
      enddo
      do n=1,ntrac
        do k=1,km
	  do i=1,imt
	    t(i,k,1,n) = disk(i,k,n+2)
	  enddo
	enddo
      enddo
      return
      end
      
      subroutine putlat (disk, u, t)
c
c-----------------------------------------------------------------------
c     copy 3D array in memory to 3-D data on ramdisk. 
c-----------------------------------------------------------------------
c
#include "param.h"
      parameter (ntrac=nvar-2)
      dimension u(imt,km,jmw,2), t(imt,km,jmw,ntrac), disk(imt,km,nvar)
      do n=1,2
        do k=1,km
	  do i=1,imt
	    disk(i,k,n) = u(i,k,1,n)
	  enddo
	enddo
      enddo
      do n=1,ntrac
        do k=1,km
	  do i=1,imt
	    disk(i,k,n+2) = t(i,k,1,n)
	  enddo
	enddo
      enddo
      return
      end


      subroutine oget (lu, nwrs, nrec, a)
c
c=======================================================================
c     get 2-D data from ramdisk record "nrec"
c
c         lu    = i/o unit number associated with "fname"
c         nwrs  = number of words to read
c         nrec  = record
c         a     = two dimensional field
c=======================================================================
c
      use datamod
# include "param.h"
# include "iounit.h"
c
      dimension a(nwds)
c
      if (lu .ne. kflds) then
        write (stdout,*) '=>Error: unit # ',lu,' not allowed. rec=',nrec
	call abort()
      else
        do n=1,nwrs
          a(n) = ramdisk2d(n,nrec)
        enddo
      endif
      return
      end

      subroutine oput (lu, nwrs, nrec, a)
c
c=======================================================================
c     put 2D field data to ramdisk record "nrec"
c
c         lu    = i/o unit number associated with "fname"
c         nwrs  = number of words to read
c         nrec  = record
c         a     = two dimensional field
c=======================================================================
c
      use datamod
# include "param.h"
# include "iounit.h"
c
      dimension a(nwds)
c
      if (lu .ne. kflds) then
        write (stdout,*) '=>Error: unit # ',lu,' not allowed. rec=',nrec
	call abort()
      else
        do n=1,nwrs
          ramdisk2d(n,nrec) = a(n)
        enddo
      endif
      return
      end
#endif




#ifdef max_window
      subroutine copy_all_rows (ifrom, ito)
c
c-----------------------------------------------------------------------
c     copy all latitude rows from time level "ifrom" to time level "ito"
c     for prognostic variables in the MW
c     
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      use mom3_mpp_mod
      use memory_window
# include "param.h"
# include "mw.h"
# if defined parallel_1d
      real,dimension(imt*km,jmw,nt) :: tmp_var_t
      real,dimension(imt*km,jmw,2) :: tmp_var_u
      pointer (ptr_t,tmp_var_t)
      pointer (ptr_u,tmp_var_u)
# endif
# ifdef timing
      call tic ('copy_all_rows', 'within the MW')
# endif
      do n=1,2
        do j=1,jmw
	  do k=1,km
	    do i=1,imt
	      u(i,k,j,n,ito) = u(i,k,j,n,ifrom)
	    enddo
	  enddo
	enddo
      enddo
      do n=1,nvar-2
        do j=1,jmw
	  do k=1,km
	    do i=1,imt
	      t(i,k,j,n,ito) = t(i,k,j,n,ifrom)
	    enddo
	  enddo
	enddo
      enddo
# ifdef timing
      call toc ('copy_all_rows', 'within the MW')
# endif
      return
      end
#endif
c
      subroutine euler_shuffle
c
c-----------------------------------------------------------------------
c     after the second pass of an euler backward step, exchange "tau"
c     and "tau+1" data, after shuffling, data will be in proper position
c     for the next time step.
c     
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      use memory_window
#include "param.h"
#include "mw.h"
#ifdef timing
      call tic ('euler_shuffle', 'within the MW')
#endif
      do n=1,2
        do j=1,jmw
	  do k=1,km
	    do i=1,imt
	      temp             = u(i,k,j,n,tau)
	      u(i,k,j,n,tau)   = u(i,k,j,n,taup1)
	      u(i,k,j,n,taup1) = temp
	    enddo
	  enddo
	enddo
      enddo
      do n=1,nt
        do j=1,jmw
	  do k=1,km
	    do i=1,imt
	      temp             = t(i,k,j,n,tau)
	      t(i,k,j,n,tau)   = t(i,k,j,n,taup1)
	      t(i,k,j,n,taup1) = temp
	    enddo
	  enddo
	enddo
      enddo
#ifdef timing
      call toc ('euler_shuffle', 'within the MW')
#endif
      return
      end


#if defined crayio || defined ssread_sswrite
      subroutine  ostart (lu, fname, ntot, nwrec, nbuf)
c
c=======================================================================
c
c     cray i/o section using getwa/putwa
c
c     assign a unit number ot "fname" and show resource requirements
c     
c     fname = filename 
c     lu    = i/o unit number associated with "fname"
c     ntot  = total length of unit, in words
c     nwrec = number of words in the record on the unit
c     nbuf  = number of buffers supplied to the unit
c=======================================================================
c
      character*(*) fname
      character*80 optlst
# include "stdunits.h"
# ifdef ssread_sswrite
      nrows   = ntot/nwrec
      nwdprow = (nwrec/512+1)*512
      ntw     = nrows * nwdprow
      write (optlst,'(a,i10,a)')
     &  'ssdrw buffer:words=',ntw,' sds scratch novfl'
      call getunit (lu, trim(fname), optlst) 
      write (stdout,'(/a,a,a,1pg10.3,a)')
     & ' Cray SSREAD/SSWRITE i/o disk resource for file ',trim(fname)
     &,' = ',ntw*1.e-6,' (MW)'
# else 
      write (optlst,'(a,i10,a)')
     &  'word buffer:words=',nwrec,' sds scratch novfl'
      call getunit (lu, trim(fname), optlst) 
      write (stdout,'(/a,a,a,1pg10.3,a)')
     & ' Cray Wordio disk resource for file ',trim(fname)
     &,' = ',ntot*1.e-6,' (MW)'
# endif
      return
c
      entry oclose (lu)
# ifndef ssread_sswrite
      call wclose (lu)
# endif
      return
c
      entry ofind (lu, nwrs, nrec)
# ifndef ssread_sswrite
      call seek (lu, (nrec-1)*nwrs+1, nwrs)
# endif
      return
      end
      
      subroutine getrow (lu, nwrs, nrec, u, t)
c
c=======================================================================
c     get velocity and tracer data from latitude row "nrec"
c
c         lu    = i/o unit number associated with "fname"
c         nwrs  = number of words to read
c         nrec  = latitude row
c         u     = velocities
c         t     = tracers
c         ntrac = number of tracers
c=======================================================================
c
# include "param.h"
      parameter (ntrac=nvar-2)
# include "iounit.h"
c
      dimension u(imt,km,jmw,2), t(imt,km,jmw,ntrac)
# ifdef ssread_sswrite
#  include "iomngr.h"
c
      dimension ssbuf(imt*km*nvar+512)
c
      nblocks = (nwrs/512+1)
      nwprow = nblocks*512
      sdsaddr = sdsbase(lu) + (nrec-1)*nblocks + 1
      iret = ssread (ssbuf, sdsaddr, nwprow)
      if (iret .ne. 0) then
        write (stdout,*) 'getrow: ssread error is ', iret,' file=',lu
        call abort()
      endif
      m = 0
      lenik = imt*km
      do k=1,km
        do i=1,imt
          m = m + 1
          u(i,k,1,1) = ssbuf(m)
          u(i,k,1,2) = ssbuf(m + lenik)
        enddo
      enddo
      m = 2*lenik
      do n=1,ntrac            
        do k=1,km
          do i=1,imt
            m = m + 1
            t(i,k,1,n) = ssbuf(m)
          enddo
        enddo            
      enddo            
# else
c
      nfrst = (nrec-1)*nwrs + 1
      if (lu .eq. latdisk(1) .or. lu .eq. latdisk(2)) then
        call getlat (lu, u, t, nfrst, nwrs)
      else
        write (stdout,*) '=>Error: unit #',lu,' not allowed. rec=',nrec
	call abort()
      endif
# endif
      return
      end

      subroutine getlat (lu, u, t, nfrst, nwrs)
c
c-----------------------------------------------------------------------
c     transfer contiguous portions of disk to non-contiguous
c     portions of memory.
c-----------------------------------------------------------------------
c
# include "param.h"
      parameter (ntrac=nvar-2)
      dimension u(imt,km,jmw,2), t(imt,km,jmw,ntrac)
      nwords = nwrs / nvar
      do n=1,2
        nfirst = nfrst + (n-1)*imtkm
        call getwa (lu, u(1,1,1,n), nfirst, nwords)
      enddo
      do n=1,ntrac
        nfirst = nfrst + (2+n-1)*imtkm
        call getwa (lu, t(1,1,1,n), nfirst, nwords)
      enddo
      return
      end

      subroutine oget (lu, nwrs, nrec, a)
c
c=======================================================================
c     get data from record "nrec"
c
c         lu    = i/o unit number associated with "fname"
c         nwrs  = number of words to read
c         nrec  = record
c         a     = array where data will go
c=======================================================================
c
# include "param.h"
# include "iounit.h"
c
      dimension a(nwrs)
# ifdef ssread_sswrite
#  include "iomngr.h"
      dimension ssbuf(nwrs+512)
c
      nblocks = (nwrs/512+1)
      nwprow = nblocks*512
      sdsaddr = sdsbase(lu) + (nrec-1)*nblocks + 1
      iret = ssread (ssbuf, sdsaddr, nwprow)
      if (iret .ne. 0) then
        write (stdout,*) 'oget: ssread error is ', iret,' file=',lu
        call abort()
      endif
      do n=1,nwrs
        a(n) = ssbuf(n)
      enddo
# else
c
      nfrst = (nrec-1)*nwrs + 1
      if (lu .eq. latdisk(1) .or. lu .eq. latdisk(2)) then
	write (stdout,*) '=>Error: unit #',lu, ' not allowed. rec=',nrec
        call abort()
      else
        call getwa (lu, a, nfrst, nwrs)
      endif
# endif
      return
      end


      subroutine putrow (lu, nwrs, nrec, u, t)
c
c=======================================================================
c     put velocity and tracer data to latitude row "nrec"
c
c         lu    = i/o unit number associated with "fname"
c         nwrs  = number of words to read
c         nrec  = latitude row
c         u     = velocities
c         t     = tracers
c         ntrac = number of tracers
c=======================================================================
c
# include "param.h"
      parameter (ntrac=nvar-2)
# include "iounit.h"
c
      dimension u(imt,km,jmw,2), t(imt,km,jmw,ntrac)
# ifdef ssread_sswrite
#  include "iomngr.h"
c
      dimension ssbuf(imt*km*nvar+512)
c
      nblocks = (nwrs/512+1)
      nwprow = nblocks*512
      sdsaddr = sdsbase(lu) + (nrec-1)*nblocks + 1
      m = 0
      lenik = imt*km
      do k=1,km
        do i=1,imt
          m = m + 1
          ssbuf(m) = u(i,k,1,1)
          ssbuf(m + lenik) = u(i,k,1,2) 
        enddo
      enddo
      m = 2*lenik
      do n=1,ntrac            
        do k=1,km
          do i=1,imt
            m = m + 1
            ssbuf(m) = t(i,k,1,n)
          enddo
        enddo            
      enddo            
      iret = sswrite (ssbuf, sdsaddr, nwprow)
      if (iret .ne. 0) then
        write (stdout,*) 'putrow: sswrite error is ', iret,' file=',lu
        call abort()
      endif
# else
c
      nfrst = (nrec-1)*nwrs + 1
      if (lu .eq. latdisk(1) .or. lu .eq. latdisk(2)) then
        call putlat (lu, u, t, nfrst, nwrs)
      else
        write (stdout,*) '=>Error: unit #',lu,' not allowed. rec=',nrec
	call abort()
      endif
# endif
      return
      end

      subroutine putlat (lu, u, t, nfrst, nwrs)
c
c-----------------------------------------------------------------------
c     transfer non-contiguous portions of memory to contiguous portions
c     of disk. 
c-----------------------------------------------------------------------
c
# include "param.h"
      parameter (ntrac=nvar-2)
      dimension u(imt,km,jmw,2), t(imt,km,jmw,ntrac)
      nwords = nwrs / nvar
      do n=1,2
        nfirst = nfrst + (n-1)*imtkm
        call putwa (lu, u(1,1,1,n), nfirst, nwords)
      enddo
      do n=1,ntrac
        nfirst = nfrst + (2+n-1)*imtkm
        call putwa (lu, t(1,1,1,n), nfirst, nwords)
      enddo
      return
      end

      subroutine oput (lu, nwrs, nrec, a)
c
c=======================================================================
c     put data to record "nrec"
c
c         lu    = i/o unit number associated with "fname"
c         nwrs  = number of words to read
c         nrec  = record
c         a     = array where data will go
c=======================================================================
c
# include "param.h"
# include "iounit.h"
c
      dimension a(nwrs)
# ifdef ssread_sswrite
#  include "iomngr.h"
c
      dimension ssbuf(nwrs+512)
c
      nblocks = (nwrs/512+1)
      nwprow = nblocks*512
      sdsaddr = sdsbase(lu) + (nrec-1)*nblocks + 1
      do n=1,nwrs
        ssbuf(n) = a(n)
      enddo
      iret = sswrite (ssbuf, sdsaddr, nwprow)
      if (iret .ne. 0) then
        write (stdout,*) 'oput: sswrite error is ', iret,' file=',lu
        call abort()
      endif
# else
c
      nfrst = (nrec-1)*nwrs + 1
      if (lu .eq. latdisk(1) .or. lu .eq. latdisk(2)) then
	write (stdout,*) '=>Error: unit #',lu, ' not allowed. rec=',nrec
        call abort()
      else
        call putwa (lu, a, nfrst, nwrs)
      endif
# endif
      return
      end
#endif




#if defined fio
      subroutine  ostart (lu, fname, ntot, nwrec, nbuf)
c
c=======================================================================
c
c     fortran i/o section
c
c     assign a unit number ot "fname" and show resource requirements
c     
c     fname = filename 
c     lu    = i/o unit number associated with "fname"
c     ntot  = total length of unit, in words
c     nwrec = number of words in the record on the unit
c     nbuf  = number of buffers supplied to the unit
c=======================================================================
c
      character*(*) fname
      character*80 optlst
# include "param.h"
# include "iounit.h"
      if (fname .eq. 'latdisk1' .or. fname .eq. 'latdisk2') then
        lenblk = nwrec / nvar
      else
        lenblk = nwrec
      endif
      write (optlst,'(a,i10)') 'direct recl:words = ',lenblk
      call getunit (lu, fname, optlst) 
      write (stdout,'(/a,a,a,1pg10.3,a)')
     & ' Direct access disk resource for file ',trim(fname),' = '
     &, ntot*1.e-6,' (MW)'
      return
c
      entry oclose (lu)
      return
c
      entry ofind (lu, nwrs, nrec)
      return
      end

      
      subroutine getrow (lu, nwrs, nrec, u, t)
c
c=======================================================================
c     get velocity and tracer data from latitude row "nrec" on disk
c
c         lu    = i/o unit number associated with "fname"
c         nwrs  = number of words to read
c         nrec  = latitude row
c         u     = velocities
c         t     = tracers
c         ntrac = number of tracers
c=======================================================================
c
# include "param.h"
      parameter (ntrac=nvar-2)
# include "iounit.h"
c
      dimension u(imt,km,jmw,2), t(imt,km,jmw,ntrac)
c
      nfrst = (nrec-1)*nwrs + 1
      if (lu .eq. latdisk(1) .or. lu .eq. latdisk(2)) then
        call getlat (lu, u, t, nfrst, nwrs)
      else
        write (stdout,*) '=>Error: unit #',lu,' not allowed. rec=',nrec
	call abort()
      endif
      return
      end


      subroutine getlat (lu, u, t, nfrst, nwrs)
c
c-----------------------------------------------------------------------
c     transfer contiguous portions of disk to non-contiguous
c     portions of memory.
c-----------------------------------------------------------------------
c
#include "param.h"
      parameter (ntrac=nvar-2)
      dimension u(imt,km,jmw,2), t(imt,km,jmw,ntrac)
      dimension buf(imt,km)
      nwords = nwrs / nvar
      nfirst = (nfrst-1)/nwords+1
      do n=1,2
        nrec = nfirst + n - 1
	read (lu, rec=nrec) buf
        do k=1,km
          do i=1,imt
            u(i,k,1,n) = buf(i,k)
          enddo
        enddo
      enddo
      nfirst = nfirst + 2  - 1
      do n=1,ntrac
        nrec = nfirst + n
	read (lu, rec=nrec) buf
        do k=1,km
          do i=1,imt
            t(i,k,1,n) = buf(i,k)
          enddo
        enddo
      enddo
      return
      end

      subroutine oget (lu, nwrs, nrec, a)
c
c=======================================================================
c     get data from record "nrec" on disk
c
c         lu    = i/o unit number associated with "fname"
c         nwrs  = number of words to read
c         nrec  = record
c         a     = array where data will go
c=======================================================================
c
# include "param.h"
# include "iounit.h"
c
      dimension a(nwrs)
c
      if (lu .eq. latdisk(1) .or. lu .eq. latdisk(2)) then
	write (stdout,*) '=>Error: unit #',lu, ' not allowed. rec=',nrec
        call abort()
      else
        read (lu, rec=nrec) a
      endif
      return
      end

      
      subroutine putrow (lu, nwrs, nrec, u, t)
c
c=======================================================================
c     put velocity and tracer data to latitude row "nrec" on disk
c
c         lu    = i/o unit number associated with "fname"
c         nwrs  = number of words to read
c         nrec  = latitude row
c         u     = velocities
c         t     = tracers
c         ntrac = number of tracers
c=======================================================================
c
# include "param.h"
      parameter (ntrac=nvar-2)
# include "iounit.h"
c
      dimension u(imt,km,jmw,2), t(imt,km,jmw,ntrac)
c
      nfrst = (nrec-1)*nwrs + 1
      if (lu .eq. latdisk(1) .or. lu .eq. latdisk(2)) then
        call putlat (lu, u, t, nfrst, nwrs)
      else
        write (stdout,*) '=>Error: unit #',lu,' not allowed. rec=',nrec
	call abort()
      endif
      return
      end


      subroutine putlat (lu, u, t, nfrst, nwrs)
c
c-----------------------------------------------------------------------
c     transfer non-contiguous portions of memory to
c     contiguous portions of disk. 
c-----------------------------------------------------------------------
c
#include "param.h"
      parameter (ntrac=nvar-2)
      dimension u(imt,km,jmw,2), t(imt,km,jmw,ntrac)
      dimension buf(imt,km)
      nwords = nwrs / nvar
      nfirst = (nfrst-1)/nwords+1
      do n=1,2
        nrec = nfirst + n - 1
        do k=1,km
          do i=1,imt
            buf(i,k) = u(i,k,1,n)
          enddo
        enddo
	write (lu, rec=nrec) buf
      enddo
      nfirst = nfirst + 2  - 1
      do n=1,ntrac
        nrec = nfirst + n
        do k=1,km
          do i=1,imt
            buf(i,k) = t(i,k,1,n)
          enddo
        enddo
	write (lu, rec=nrec) buf
      enddo
      return
      end

      subroutine oput (lu, nwrs, nrec, a)
c
c=======================================================================
c     put data to record "nrec" on disk
c
c         lu    = i/o unit number associated with "fname"
c         nwrs  = number of words to read
c         nrec  = record
c         a     = array where data will go
c=======================================================================
c
# include "param.h"
# include "iounit.h"
c
      dimension a(nwrs)
c
      if (lu .eq. latdisk(1) .or. lu .eq. latdisk(2)) then
	write (stdout,*) '=>Error: unit #',lu, ' not allowed. rec=',nrec
        call abort()
      else
        write (lu, rec=nrec) a
      endif
      return
      end

#endif
      

