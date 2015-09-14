      subroutine nor1st (nunit, normsk)
#ifdef meridional_overturning_basin
c
c=======================================================================
c
c     Subroutine nor1st is i/o routine for the masks used in northward 
c     meridional overturning calculation.
c
c     author:      k. dixon      e-mail=> kd@gfdl.gov
c=======================================================================
c
      use arrays_2d
      use memory_window
c
#include "param.h"
c
c     nunit = unit to be written to or read from
c
      integer normsk(imt,jmt)
c
      character *1 chrow(imt), chnum(0:9)
      data (chnum(n),n=0,9) /'0','1','2','3','4','5','6','7','8','9'/
c
      dimension ncol(imt)
      logical calculate_mask_internally
c
c-----------------------------------------------------------------------
c     read in (write out) horizontal masks
c     set linel to length of desired formatted printout line
c-----------------------------------------------------------------------
c
      linemx = 100
      linel  = 105
      line   = linel - 5
      if (line .gt. linemx) line = linemx
      nwr = (imt/line) + 1
c
      do i=1,imt
        ncol(i) = mod(i,10)
      enddo
c
      kntnor = 0

c     if calculate_mask_internally = .false., the mask is calculated 
c     internally only if the mask file does not exist;
c     if calculate_mask_internally = .true., the mask is calculated 
c     internally in any case:
c
      calculate_mask_internally = .false.
c
      read (nunit,9099,err=100,end=100)
      go to 101
 100  calculate_mask_internally = .true.
 101  continue
c
      if (calculate_mask_internally) then
        call setglobalmasks(normsk,mskvr,'northmask')
      end if
      if (.not. calculate_mask_internally) read (nunit,9099)
      do n=1,nwr
        ia = 1 + (line*(n-1))
        ib = ia + line - 1
        if (ib .gt. imt) ib = imt
        write(stdout,9021) (ncol(i),i=ia,ib)
        if (.not. calculate_mask_internally) read (nunit,9099)
        do jj=1,jmt
          jjj = jmt - jj + 1
          if (.not. calculate_mask_internally) then
            read (nunit,9022) jr , (normsk(i,jjj),i=ia,ib)
          else
            jr = jjj
          end if
          if (jr .ne. jjj) then
            write (stdout,999) nunit, jjj, jr
            write (stderr,999) nunit, jjj, jr
            call abort()
          endif
          write(stdout,'(" nor=",1x,i3,1x,100(i1))')
     &         jjj,(normsk(i,jjj),i=ia,ib)
          do i=ia,ib
            if (i.ne.1 .and. i.ne.imt) then
              if (kmt(i,jjj).GT.0 .and. normsk(i,jjj).EQ.0) then
                write(stdout,9116) i,jjj,kmt(i,jjj),normsk(i,jjj)
              elseif (kmt(i,jjj).EQ.0 .and. normsk(i,jjj).ne.0) then
                write(stdout,9116) i,jjj,kmt(i,jjj),normsk(i,jjj)
              else
                kntnor = kntnor + 1
              endif
            endif
	  enddo
	enddo
      enddo
c
      write(stdout,9123) kntnor
c
#endif
      return
c
999   format(/' error => bad j-row when reading normsk from unit ',
     &       i3,/'   expected',i4,'    read in',i4)
9021  format('  i=>',100(i1))
9022  format(1x,i3,1x,100(i1))
9116  format(' WARNING: CHECK normsk at i,j=',2i4
     &     ,'  kmt=',i3,'  normsk=',i3)
9123    format(' A total of ',i6,' normsk points checked out OK')     
9099  format(1x)
      end
