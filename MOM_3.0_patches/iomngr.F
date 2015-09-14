
!-----------------------------------------------------------------------
!     Modified for parallel execution
!
!     author: V. Balaji      e-mail => vb@gfdl.gov
!-----------------------------------------------------------------------

#include "derived_options.h"

      subroutine ioinit
c
c-----------------------------------------------------------------------
c     setup allowable unit numbers. Files on solid state disk cannot
c     be interrogated by the inquire statement so explicit handling
c     of unit numbers is needed.
c-----------------------------------------------------------------------
c
      logical allowable_unit
      character*128 string
      character*32 unit_name
      common /iounums/ allowable_unit(100), unit_name(100), ioinitx
c
      do n=1,10
        allowable_unit(n) = .false.
        unit_name(n)      = ' '
      enddo
      do n=11,99
        allowable_unit(n) = .true.
        unit_name(n)      = ' '
      enddo
      ioinitx = 1      
c
#if defined _CRAY
# if defined CRAYPVP
        string = 'assign -P thread p:%'
# else
        string = 'assign -P private p:%'
# endif
        call do_assign (string)
#endif
      return
      end

      
      subroutine getunit (nu, filename, attributes)
c
c-----------------------------------------------------------------------
c     Associate unit number "nu" with "filename" and open the file with
c     the requested "attributes"
c
c     input:  "filename" and "attributes"
c             
c     output: nu = unit number
c
c     author: R.C. Pacanowski      e-mail => rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      character*(*) filename, attributes
      character*20 attrib, form, position, access
      character*2 unit
      character*128 string, text1
      logical find_substring
      logical keep_going
#include "stdunits.h"
      logical allowable_unit
      character*32 unit_name
      common /iounums/ allowable_unit(100), unit_name(100), ioinitx
#if defined SGICRAY
# include "iomngr.h"
      dimension info(10)
      integer oldbase, base
      data sdsoff /0/
      save oldbase, base
#endif
c
c     check to see if the filename is already attached to a unit
c     if it is, use the old unit number and return
c
      keep_going = .true. 
      n = 0
      do while (n .le. 99 .and. keep_going)
        n = n + 1 
        if (.not. allowable_unit(n)) then
          if (unit_name(n) .eq. trim(filename)) then
            nu = n
            allowable_unit(nu) = .false.
#if defined test_iomngr || defined debug_iomngr
            write (stdout,*) '-->Re-using unit ',nu
     &,     ' already attached to filename=',trim(filename)
#endif
            keep_going = .false.
          endif
        endif
      enddo
c
c     get a unit number and and remember the filename
c
      if (keep_going) then

      call getunitnumber (nu)
      write (unit,'(i2)') nu
      unit_name(nu) = trim(filename)
c
c     determine the attributes of the file and open the file using
c     the previously found unit number
c
#if defined test_iomngr || defined debug_iomngr
      write (stdout,*) 'Binding unit ',nu,' to filename=',trim(filename)
#endif
c
c              S E Q U E N T I A L      F I L E S
c
      if (find_substring(attributes, 'sequential')) then
        access = 'sequential'
c
        if (find_substring(attributes, 'unformatted')) then
          form = 'unformatted'
        else
          form = 'formatted'
        endif
c
        if (find_substring(attributes, 'rewind')) then
          position = 'rewind'
        else
          position = 'append'
        endif
c
#if defined SGICRAY
c
        if (find_substring(attributes, 'ieee')) then
c
c         assign for 32 bit IEEE sequential files
c
          if (find_substring(attributes, 'words')) then
            call find_number (len_words, attributes, 'words')
            nblock = (len_words/512) + 1
            call block_2_string (nblock, text1)
!            string = 'assign -F f77,cachea::::1 -N ieee f:'//trim(filename)
!          else
!            string = 'assign -F f77,cachea::::1 -N ieee f:'//trim(filename)
          endif
# ifdef _CRAY
         string = 'assign -F f77,cachea::::1 -N ieee f:'//trim(filename)
# endif
# ifdef sgi_mipspro
!         assign -N not yet implemented in Irix, Balaji Sep 1998.
          string = 'assign -F f77 f:'//trim(filename)
# endif
          call do_assign (string)
        elseif (find_substring(attributes, 'cray_float')) then
c
c         assign for reading older C90 non 64 bit IEEE sequential files
c
          string = 'assign -N cray f:'//trim(filename)
          call do_assign (string)
        elseif (find_substring(attributes, 'words')) then
c
c         assign for plain vanilla sequential files
c
          call find_number (len_words, attributes, 'words')
          nblock = (len_words/512) + 1
          call block_2_string (nblock, text1)
          string = 'assign  -F cos,cachea::::1 f:'//trim(filename)
          call do_assign (string)
        elseif (form .eq. 'unformatted') then
          string = 'assign  -F cos,cachea::::1 f:'//trim(filename)
          call do_assign (string)
        endif
#endif
c
#if defined test_iomngr || defined debug_iomngr
        write (stdout,*) trim(filename), ' attributes: ',form
     &,                  ' ', access,' ', position 
#endif
c
        open (nu, FILE=trim(filename), FORM=form, ACCESS=access
     &,       POSITION=position)
c
c
c              D I R E C T    A C C E S S     F I L E S
c
      elseif (find_substring(attributes, 'direct')) then
c
        call find_number (len_words, attributes, 'words')
c
c       convert to bytes assuming 8 bytes per word
c
        len_bytes = 8*len_words
#if defined test_iomngr || defined debug_iomngr
        write (stdout,*) ' file attributes: direct, record length='
     &,         len_words,' words or '
     &,         len_bytes,' bytes assuming 8 bytes per word'
#endif
c
#if defined CRAYPVP
        if (find_substring(attributes, 'sds')) then
c          recall that on ssd, "scr" implies once closed, data is gone!
          nblock = (len_words/512) + 1
          call block_2_string (nblock, text1)
          string = 'assign -F cachea::::1,sds.scr.novfl u:'//unit
          call do_assign (string)
        endif
#endif
        open (nu, FILE=trim(filename), FORM='UNFORMATTED'
     &,       ACCESS='DIRECT',  RECL=len_bytes)
c
#if defined CRAYPVP
c
c              C R A Y    S S R E A D / S S W R I T E     F I L E S
c
      elseif (find_substring(attributes, 'ssdrw ')) then
c
        ierr = 0
        call find_number (len_words, attributes, 'words')
c
c       SDS allocation is zero based for ssreads and writes.
c
        nblkpf = len_words/512 + 1
        nblock = sdsoff + nblkpf - 1
        if (sdsoff .eq. 0) then
          base = sdsalloc (nblock, ierr)
        else
          base = sdsrealc (oldbase, nblock, ierr)
        endif
# if defined test_iomngr || defined debug_iomngr
        write (stdout,'(/a,a,a,i4,3(a,i10))')
     &'=> filename: ',trim(filename),', unit: ',nu,' allocated '
     &, nblkpf,' blocks from ',sdsoff,' to ',sdsoff + nblkpf - 1
c
        call sdsinfo (istat, info, 9)
        write (stdout,*) ' SDS space available=',info(1)
        write (stdout,*) ' SDS arena size=',info(2)
        write (stdout,*) ' SDS tot allocated space in arena=',info(7)
        write (stdout,*) ' SDS # allocations in arena=',info(9)
        if (istat .ne. 0) then
          write (stdout,*) ' Error: sdsinfo istat =',istat
          call abort()
        endif
# endif
        if (ierr .ne. 0) then
          write (stdout,*) ' Error: sdsalloc ierr =',ierr
          call abort()
        endif
c
        sdsbase(nu) = sdsoff
        sdsoff      = nblock + 1
        oldbase     = base
# if defined test_iomngr || defined debug_iomngr
        write (stdout,*)
     &' file attributes: CRAY ssread/sswrite, file length ='
     &,         len_words,'words. Tot blocks on SSD = ', nblock
     &,         ', sds base address = ',sdsbase(nu)
# endif
c
c              C R A Y    W O R D I O     F I L E S
c
      elseif (find_substring(attributes, 'word ')) then
c
        if (find_substring(attributes, 'sds')) then
          string = 'assign -F sds.scr.novfl u:'//unit
          call do_assign (string)
        endif
c
        call find_number (len_words, attributes, 'words')
        nblock = (len_words/512) + 1
        call wopen (nu, nblock, 1, ierr)
        if (ierr .ne. 0) then
          write (stdout,*) ' Error: wopen ierr =',ierr
          call abort()
        endif
# if defined test_iomngr || defined debug_iomngr
        write (stdout,*)' file attributes: CRAY wordio, record length ='
     &,         len_words,'words or '
     &,         nblock,' 512 word blocks'
# endif
#endif
        
      else
        write (stdout,*) 'Error: filename = ',trim(filename)
     &,' has illegal attribute.  requested attributes = ',attributes
        call abort()
      endif
      endif
      return
      end


      subroutine getunitnumber (nu)
c
c-----------------------------------------------------------------------
c     find a unit number that is not currently attached to a file
c
c     input:  nothing
c     output: nu = the unit number
c
c     author: R.C. Pacanowski      e-mail => rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      logical in_use
      character*32 string
      character*64 name
      logical allowable_unit
      character*32 unit_name
      common /iounums/ allowable_unit(100), unit_name(100), ioinitx
#include "stdunits.h"
c
      if (ioinitx .ne. 1) then
        write (stdout,*) '=>Error: must make a call to "ioinit" before'
        write (stdout,*) '         using the iomngr.'
        call abort()
      endif
c
c     find an available unit number. "inquire" is used to catch files
c     written without using the iomngr (i.e. suppose an atmos model
c     were to use: write (30) xx. "inquire" would not use unit 30)
c
      do n=1,99
        if (allowable_unit(n)) then
          inquire (n, OPENED=in_use, NAME=name)
#if defined test_iomngr
          write (stdout,*) ' inquire: unit=',n,' in_use=',in_use
     &,                    ' name=',name
#endif
          if (.not. in_use) then
            nu = n
            allowable_unit(n) = .false.
            unit_name(n)      = ' '
#if defined test_iomngr || defined debug_iomngr
            write (stdout,*) 'Acquiring unit n=',n
#endif
            exit
          else
            allowable_unit(n) = .false.
            unit_name(n)      = name
          endif
        endif
#if defined test_iomngr || defined debug_iomngr
        if (unit_name(n) .ne. ' ') then
          write (stdout,*) 'Searching units: unit n=',n
     &,                    ' is in use. name=', unit_name(n)
        endif
#endif
        if (n .eq. 99) then
          write (stdout,*) 'Error: exhausted available unit numbers'
          write (stdout,*) '       Are you forgetting to close units?'
          call abort()
        endif
      enddo
      return
      end

      subroutine relunit (nu)
c
c-----------------------------------------------------------------------
c     release unit "nu" and close associated file
c-----------------------------------------------------------------------
c
      logical allowable_unit
      character*32 unit_name
      common /iounums/ allowable_unit(100), unit_name(100), ioinitx
#include "stdunits.h"
      character*2 unit
#if defined test_iomngr || defined debug_iomngr
      write (stdout,*) 'relunit: Releasing unit=',nu
     &,                ' and closing file=',trim(unit_name(nu))
#endif
      
      close (nu)
      allowable_unit(nu) = .true.
      unit_name(nu)      = ' '
#if defined SGICRAY
      write (unit,'(i2)') nu
      call assign ('assign -R u:'//unit, ierr)
      if (ierr .ne. 0) then
        write (stdout,*)' assign -R failed in relunit: err=',ierr
        call assign ('assign -V', ierr)
        call abort()
      endif
#endif
      return
      end



#if defined SGICRAY
      subroutine do_assign (string)
# include "stdunits.h"
      character*(*) string
# ifdef parallel_1d
#  if defined test_iomngr || defined debug_iomngr
      write (stdout,*) ' string in do_assign = ', trim(string)
#  endif
# endif
      call assign (trim(string), ierr)
      if (ierr .ne. 0) then
        write (stdout,*)' assign failed:err=',ierr,' string=',string
        call assign ('assign -V', ierr)
        call abort()
      endif
# if defined test_iomngr || defined debug_iomngr
      call assign ('assign -V', ierrr)
# endif
      return
      end
#endif

      function find_substring (string, substring)
      character*(*) substring, string
      logical find_substring
      if (index (string, substring) .eq. 0) then
        find_substring = .false.
      else
        find_substring = .true.
      endif
      return
      end
          
      subroutine block_2_string (nblock, text1)
      character*(*) text1
      character*10 text2
      write (text2,'(i10)') nblock
      text1 = trim(adjustl(text2))
      return
      end


      subroutine find_number (number, string, target1)
#include "stdunits.h"
      character*(*) string, target1
      character*1 c
c
      lens = len_trim(string)
c
      itarget  = index (string, target1)
      ilocal = index (string(itarget:lens), '=')
      if (ilocal .ne. 0) then
        ist = ilocal + itarget
        do i=ist,lens
          c = string(i:i)
          if (c .ge. '0' .and. c .le. '9') then
            inum = i
            exit
          endif
        enddo
      else
        write (stdout,*) 
     &  ' Error: need "=" sign when specifying word length'
     &, ' attributes = ',string
        call abort()
      endif
c
# if defined test_iomngr || defined debug_iomngr
      write (stdout,*)
     &' searching for word: "',target1,'" in string: "',string,'"'
      write (stdout,*)
     &' itarget=',itarget, ' inum=',inum 
# endif
      is = inum
      ie = lens
      do i=is,lens
        c = string(i:i)
        if (c .lt. '0' .or. c .gt. '9') then
          ie = i-1
          exit
        endif
      enddo
      if (ie .ne. 0) then
        read (string(is:ie), *) number 
      else
        write (stdout,*) 
     &  ' Error: cannot find word length in string=',string
        call abort()
      endif
      return
      end      



#ifdef test_iomngr
      program driver
c
c----------------------------------------------------------------------
c     test the I/O manager
c     author:  R. C. Pacanowski     e-mail => rcp@gfdl.gov
c----------------------------------------------------------------------
c
c      parameter (imt=92,jmt=61,km=15,nt=2,num2df=3)
c      parameter (imt=182,jmt=86,km=22,nt=4,num2df=3)
c      parameter (imt=73,jmt=58,km=15,nt=2,num2df=3)
c      parameter (imt=500,jmt=100,km=15,nt=2,num2df=3)
c      parameter (imt=512,jmt=100,km=10,nt=2,num2df=3)
      parameter (imt=192,jmt=62,km=50,nt=2,num2df=3)
#include "stdunits.h"
#if defined SGICRAY
# include "iomngr.h"
#endif
c
      parameter (len2d =((imt*jmt)/512+1)*512)
      parameter (len3d = ((imt*km*(2+nt))/512+1)*512) 
      dimension ssbuf2d(len2d), ssbuf3d(len3d)
c
      dimension a(imt,jmt), b(imt,jmt), buf(imt)
      dimension var(imt,km,nt+2), tdf(imt,jmt,num2df)
      character*64 st1, st2, reclen, opts, filename, totlen, string
      namelist /name1/ a1, i1, b1
c
      stdout = 6
c
c     set allowable unit numbers for files
c
      call ioinit
c      
      a1 = 1.0
      i1 = 5
      b1 = 3.0
      a2 = a1
      i2 = i1
      b2 = b1
c
      write(stdout,'(/,10x,a,/)') 'Testing the I/O manager'
# ifdef sgi_mipspro
      write (stdout,*) 'Platform = sgi mipspro'
# endif
# ifdef YMP
      write (stdout,*) 'Platform = cray_ymp'
# endif
# ifdef _CRAYC90
      write (stdout,*) 'Platform = cray_c90'
# endif
# ifdef _CRAYT90
      write (stdout,*) 'Platform = cray_t90'
# endif
# ifdef _CRAYT3E
      write (stdout,*) 'Platform = cray_t3e'
# endif
c----------------------------------------------------------------------
      write (stdout,'(/,a,/)') 'Testing: formatted sequential files'
      filename = 'iotest.doc'
      call getunit (iou, filename, 'formatted sequential rewind')
      write (iou,name1)
      call relunit (iou)
      call getunit (iou, filename, 'formatted sequential append')
      write (iou,'(1x,a,i10)') 'imt=',imt
      call relunit (iou)
      call getunit (iou, filename, 'formatted sequential rewind')
      read (iou,name1)
      read (iou,'(1x,a,i10)') st1(1:4), ii
      call relunit (iou)
      if ((a1 .ne. a2) .or. (i1 .ne. i2) .or. (b1 .ne. b2) .or. 
     &    (st1(1:4) .ne. 'imt=') .or. (ii .ne. imt)) then
       write (stdout,*) ' Error on formatted sequential file test'
       write (stdout,*) ' a1, i1, b1, st1, ii=',a1, i1, b1, st1(1:4), ii
       call abort()
      endif
      write (stdout,*) ' ok'
c----------------------------------------------------------------------
      write (stdout,'(/,a,/)') 'Testing: unformatted files'
      do j=1,jmt
        do i=1,imt
          a(i,j) = i*j
          b(i,j) = i+j
        enddo
      enddo
      filename = 'ioregular.doc'
      write (string,'(a,i10)')
     & 'unformatted sequential rewind words=',imt*jmt
      call getunit (iou, filename, string)
      write (iou)  a
      call relunit (iou)
      write (string,'(a,i10)')
     & 'unformatted sequential append words=',imt*jmt
      call getunit (iou, filename, string)
      write (iou)  b
      call relunit (iou)
      write (string,'(a,i10)')
     & 'unformatted sequential rewind words=',imt*jmt
      call getunit (iou, filename, string)
      read (iou)   a      
      read (iou)   a      
      call relunit (iou)
      do j=1,jmt
        do i=1,imt
          if (a(i,j) .ne. b(i,j)) then
            write (stdout,*) ' Error reading unformatted file'
            call abort()
          endif
        enddo
      enddo
      write (stdout,*) ' ok'
c----------------------------------------------------------------------
      write (stdout,'(/,a,/)') 'Testing: unformatted 32bit ieee files'
      do j=1,jmt
        do i=1,imt
          a(i,j) = i*j
          b(i,j) = i+j
        enddo
      enddo
      filename = 'ioieee.doc'
      call getunit (iou, filename, 'unformatted sequential rewind ieee')
      write (iou) st1, imt, a1, i1, b1, a
      call relunit (iou)
      call getunit (iou, filename, 'unformatted sequential append ieee')
      write (iou) st1, imt, a1, i1, b1, b
      call relunit (iou)
      call getunit (iou, filename, 'unformatted sequential rewind ieee')
      read (iou)  st2, ii, a2, i2, b2, a      
      read (iou)  st2, ii, a2, i2, b2, a      
      call relunit (iou)
      do j=1,jmt
        do i=1,imt
          if (a(i,j) .ne. b(i,j)) then
            write (stdout,*) ' Error reading ieee file'
            call abort()
          endif
        enddo
      enddo
      write (stdout,*) ' ok'
c----------------------------------------------------------------------
      write (stdout,'(/,a,/)') 'Testing: direct access files'
c      allow a conservative 4 bytes per word for character variables
      st1 = 'abcdef'
      write(reclen,'(i10)') imt + len(st1)/4 + 1
      opts = 'direct words='// reclen
      filename = 'iod.dta'
      call getunit (iou, filename, opts)
      do j=1,jmt
        do i=1,imt
          buf(i) = a(i,j)
        enddo
        write (iou,rec=j) st1, buf
      enddo
c      call relunit (iou)
      call getunit (iou, filename, opts)
      do j=1,jmt
        read (iou,rec=j) st2, buf
        do i=1,imt
          if (buf(i) .ne. a(i,j)) then
            write (stdout,*) ' Error reading direct access file'
            write (stdout,*) ' Character string written =',st1
            write (stdout,*) ' Character string read =',st2
            write (stdout,*) ' i=',i,' j=',j,' wrote a(i,j)=',a(i,j)
     &,             ' read a(i,j)=',buf(i)
            call abort()
          endif
        enddo
      enddo
      call relunit (iou)
      write (stdout,*) ' ok'

# if defined CRAYPVP
c----------------------------------------------------------------------
      write (stdout,'(/,a,/)') 
     &'Testing: CRAY "wordio" on rotating disk (no sds)'
c      allow a conservative 4 bytes per word for character variables
      st1 = 'abcdef'
      filename = 'iowordr.dta'
      write(reclen,'(i10)') imt
      opts = 'word  words='// reclen
      call getunit (iou, filename, opts)
      do j=1,jmt
        do i=1,imt
          buf(i) = a(i,j)
        enddo
        call putwa (iou, buf, (j-1)*imt+1, imt)
      enddo
      do j=1,jmt
        call getwa (iou, buf, (j-1)*imt+1, imt)
        do i=1,imt
          if (buf(i) .ne. a(i,j)) then
            write (stdout,*) ' Error reading crayio file'
            write (stdout,*) ' i=',i,' j=',j,' wrote a(i,j)=',a(i,j)
     &,             ' read a(i,j)=',buf(i)
            call abort()
          endif
        enddo
      enddo
      call wclose (iou)
      call getunit (iou, filename, opts)
      do j=1,jmt
        call getwa (iou, buf, (j-1)*imt+1, imt)
        do i=1,imt
          if (buf(i) .ne. a(i,j)) then
            write (stdout,*) ' Error reading re-opened crayio file'
            write (stdout,*) ' i=',i,' j=',j,' wrote a(i,j)=',a(i,j)
     &,             ' read a(i,j)=',buf(i)
            call abort()
          endif
        enddo
      enddo
      call wclose (iou)
      write (stdout,*) ' ok'
c
c----------------------------------------------------------------------
c     TESTING SOLID STATE DISK (can't release then reopen units because
c     units are assigned as "scratch files" and data is gone once
c     file is released
c----------------------------------------------------------------------
c
      write (stdout,'(//,a,//)')
     & '  =========> NOTE: must have CRAY SSD for the following to work'
      write (stdout,*) ' '
      write (stdout,'(/,a,/)') 'Testing: direct access on CRAY SSD'
c      allow a conservative 4 bytes per word for character variables
      st1 = 'abcdef'
      filename = 'iodssd.dta'
      write(reclen,'(i10)') imt + len(st1)/4 + 1
      opts = 'direct sds words='// reclen
      call getunit (iou, filename, opts)
      do j=1,jmt
        do i=1,imt
          buf(i) = a(i,j)
        enddo
        write (iou,rec=j) st1, buf
      enddo
      do j=1,jmt
        read (iou,rec=j) st2, buf
        do i=1,imt
          if (buf(i) .ne. a(i,j)) then
            write (stdout,*) ' Error reading direct access file'
            write (stdout,*) ' Character string written =',st1
            write (stdout,*) ' Character string read =',st2
            write (stdout,*) ' i=',i,' j=',j,' wrote a(i,j)=',a(i,j)
     &,             ' read a(i,j)=',buf(i)
            call abort()
          endif
        enddo
      enddo
      call relunit (iou)
      write (stdout,*) ' ok'
c----------------------------------------------------------------------
      write (stdout,'(/,a,/)') 'Testing: CRAY "wordio" on CRAY SSD'
      filename = 'ioword.dta'
      write(reclen,'(i10)') imt
      opts = 'word sds words='// reclen
      call getunit (ioun, filename, opts)
      do j=1,jmt
        do i=1,imt
          buf(i) = a(i,j)
        enddo
        call putwa (ioun, buf, (j-1)*imt+1, imt)
      enddo
      do j=1,jmt
        call getwa (ioun, buf, (j-1)*imt+1, imt)
        do i=1,imt
          if (buf(i) .ne. a(i,j)) then
            write (stdout,*) ' Error reading wordio file'
            write (stdout,*) ' i=',i,' j=',j,' wrote a(i,j)=',a(i,j)
     &,             ' read a(i,j)=',buf(i)
            call abort()
          endif
        enddo
      enddo
      call wclose (ioun)
      write (stdout,*) ' ok'
c
c----------------------------------------------------------------------
c     SSREAD/SSWRITE Section (the ultimate in speed for CRAY SSD)
c----------------------------------------------------------------------
c
      write (stdout,'(/,a,/)')
     & 'Testing: CRAY "ssread/sswrite" on CRAY SSD'
c
      do n=1,num2df
        do j=1,jmt
          do i=1,imt
            tdf(i,j,n) = (i+0.5)*j*n
          enddo
        enddo
      enddo
      do n=1,nt+2
        do k=1,km
          do i=1,imt
            var(i,k,n) = i*k*n
          enddo
        enddo
      enddo
c
      nblocks2 = ((imt*jmt)/512+1)
      nwprow2 = nblocks2*512
      ntw2    = num2df * nwprow2
      write(totlen,'(i10)') ntw2
      opts = 'ssdrw sds words='// totlen
      filename = 'ioss2d.dta'
      call getunit (ioun2, filename, opts)
      do i=1,nwprow2
        ssbuf2d(i) = 0.0
      enddo
c
      nblocks3 = ((imt*km*(nt+2))/512+1)
      nwprow3 = nblocks3*512
      ntw3    = jmt * nwprow3
      write(totlen,'(i10)') ntw3
      opts = 'ssdrw sds words='// totlen
      filename = 'ioss3d.dta'
      call getunit (ioun3, filename, opts)
      filename = 'ioss3d2.dta'
      call getunit (ioun4, filename, opts)
      do i=1,nwprow3
        ssbuf3d(i) = 0.0
      enddo
c
c-------------------------- write 2D fields-------------
c
      do nfield=num2df,1,-1
        n = 0
        do j=1,jmt
          do i=1,imt
            n = n + 1
            ssbuf2d(n) = tdf(i,j,nfield)
          enddo
        enddo
        sdsaddr = sdsbase(ioun2) + (nfield-1)*nblocks2 + 1
        iret = sswrite (ssbuf2d, sdsaddr, nwprow2)
        if (iret .ne. 0) write (stdout,*) 'err1. sswrite iret is ',iret
      enddo
c
c-------------------------- write 3D fields-------------
c
      do j=1,jmt
        n = 0
        do m=1,nt+2
          do k=1,km
            do i=1,imt
              n = n + 1
              ssbuf3d(n) = var(i,k,m)*j
            enddo
          enddo
        enddo
        sdsaddr =  sdsbase(ioun3) + (j-1)*nblocks3 + 1
        iret = sswrite (ssbuf3d, sdsaddr, nwprow3)
        if (iret .ne. 0) write (stdout,*)'err2. sswrite iret is ',iret
        if (iret .ne. 0) call abort()
	sdsaddr =  sdsbase(ioun4) + (j-1)*nblocks3 + 1
        iret = sswrite (ssbuf3d, sdsaddr, nwprow3)
        if (iret .ne. 0) write (stdout,*)'err3. sswrite iret is ',iret
        if (iret .ne. 0) call abort()
      enddo
c
c-------------------------- read 2D fields-------------
c
      do nfield=1,num2df
        sdsaddr = sdsbase(ioun2) + (nfield-1)*nblocks2 + 1
        iret = ssread (ssbuf2d, sdsaddr, nwprow2)
        if (iret .ne. 0) write (stdout,*) 'err1a. ssread iret is ',iret
        n = 0
        do j=1,jmt
          do i=1,imt
            n = n + 1
            if (ssbuf2d(n) .ne. tdf(i,j,nfield)) then
              write (stdout,*) ' Error reading ssread/write file',nfield
              write (stdout,*) ' i=',i,' j=',j,' wrote b(i,j)='
     &,             tdf(i,j,nfield), ' read b(i,j)=',ssbuf2d(n)
              call abort()
            endif
          enddo
        enddo
      enddo
      write (stdout,*) ' ok: reading 2d fields is correct'
c
c-------------------------- read 3D fields-------------
      do io=1,2
      if (io .eq. 1) ioun = ioun3
      if (io .eq. 2) ioun = ioun4
      do j=1,jmt
        sdsaddr = sdsbase(ioun) + (j-1)*nblocks3 + 1
        iret = ssread (ssbuf3d, sdsaddr, nwprow3)
        if (iret .ne. 0) write (stdout,*) 'err2. ssread return is ',iret
        n = 0
        do m=1,nt+2
          do k=1,km
            do i=1,imt
              n = n + 1
              if (var(i,k,m)*j .ne. ssbuf3d(n)) then
                write (stdout,*) ' Error reading ssread/write file 2'
                write (stdout,*) 'i=',i,' j=',j, 'k=',k,' wrote b(i,j)='
     &,             var(i,k,m)*j, ' read b(i,k,j)=',ssbuf3d(n)
                call abort()
              endif 
            enddo
          enddo
        enddo
      enddo
      write (stdout,*) ' ok: reading 3d fields is correct', io
      enddo
# endif      
      end
#endif







