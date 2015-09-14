      subroutine reg1st (nunit, wrform, wrvol, wrk, wrmask, rmask)

#if defined tracer_averages || defined term_balances
# define traceravgbal
#endif
#if defined traceravgbal || defined gyre_components
c
c=======================================================================
c
c     Subroutine reg1st is i/o routine for the user defined horizontal
c     and vertical regional masks used in MOM.
c     It is also the i/o routine that can be used to write volume
c     information for the horizontal regions used in calculating
c     volume weighted averages of tracers and surface tracer fluxes.
c     Both formatted and unformatted i/o is supported.
c     (see "diag_region.h" for more details on variables).
c
c     when "readrmsk" is defined, "mskhr","mskvr", "hregnm" &"vregnm"
c     will be read from file "regionmasks" which is assigned to
c     "iormsk" in the same manner as this routine writes that
c     information. ("iormsk" is set in "blkdta.F")
c     alternatively, the horizontal region masks can be specified in
c     "setocn.F", and although a simple, latitude dependent default
c     setting is provided, the user should consider specifying a more
c     appropriate set of masks.
c     
c     author:      k. dixon      e-mail=> kd@gfdl.gov
c
c=======================================================================
c
      use arrays_2d
#include "param.h"
#include "diag_region.h"
#include "iounit.h"
#include "taskrows.h"
c
c  nunit = unit to be written to or read from
c  wrform= true(false) switch for formatted(unformatted) writes
c  wrvol = switch to write volume & area information
c  wrk   = switch to write k-level volume information
c  wrmask= switch to write horizontal region masks field
c  rmask = switch to read horizontal region masks field from specified
c          unit
c  settop, setbot are used in defining level limits for vertical regions
c
      logical  wrform, wrvol, wrk, wrmask, rmask
      logical  settop, setbot, calculate_mask_internally
c
      dimension ncol(imt)
c
c-----------------------------------------------------------------------
c     write out regional volume information
c-----------------------------------------------------------------------
c
      if (pn .ne. 1) return
      if (wrvol) then
        if(wrform) then
          write(nunit,9000)
          if(wrk) then
            write(nunit,9001) (mask,mask=1,nhreg)
            do 700 k=1,km
              write(nunit,9002) k, volgk(k),
     &             (volbk(mask,k),mask=1,nhreg)
 700        continue
          endif
          write(nunit,9003) volgt, (volbt(mask),mask=1,nhreg)
          write(nunit,9004) areag, (areab(mask),mask=1,nhreg)
        else
c         
          iotext ='read (iotavg) km, nhreg, volgt, areag'
          write (nunit) 'no stamp available in reg1st    ', iotext
     &         ,                   expnam       
          write (nunit) km, nhreg, volgt, areag
c         
          iotext ='read (iotavg) (volgk(k),k=1,km)'
          write (nunit) 'no stamp available in reg1st    ', iotext
     &         ,                   expnam       
          call wrufio (nunit, volgk, km)
c         
          iotext ='read (iotavg) ((volbk(l,k),l=1,nhreg),k=1,km)'
          write (nunit) 'no stamp available in reg1st    ', iotext
     &         ,                   expnam       
          call wrufio (nunit, volbk, nhreg*km)
c         
          iotext ='read (iotavg) (volbt(l),l=1,nhreg)'
          write (nunit) 'no stamp available in reg1st    ', iotext
     &         ,                   expnam       
          call wrufio (nunit, volbt, nhreg)
c         
          iotext ='read (iotavg) (areab(l),l=1,nhreg)'
          write (nunit) 'no stamp available in reg1st    ', iotext
     &         ,                   expnam       
          call wrufio (nunit, areab, nhreg)
        endif
      endif
c     
c-----------------------------------------------------------------------
c     read in (write out) horizontal & vertical region masks
c     set linel to length of desired formatted printout line
c-----------------------------------------------------------------------
c
c     if calculate_mask_internally = .false., the mask is calculated 
c     internally only if the mask file does not exist;
c     if calculate_mask_internally = .true., the mask is calculated 
c     internally in any case:
      calculate_mask_internally = .false.

      if (wrmask .or. rmask) then
c
        if(wrform) then
          if (wrmask) write(nunit,9011)
          if (wrmask) write(nunit,9012)
     &         (' domain for hor mask #',i,'=',hregnm(i),i=1,nhreg)
          if ( rmask) then
            read (nunit,9099,err=100,end=100)
            go to 101
 100        calculate_mask_internally = .true.
 101        continue
            if (calculate_mask_internally) then
              call setglobalmasks(mskhr,mskvr,'regionmask')
            end if
            if (.not. calculate_mask_internally) read (nunit,9099)
          endif
          if (.not. calculate_mask_internally) then
            if ( rmask) read (nunit,9013) (hregnm(i),i=1,nhreg)
          end if
          if (nunit .eq. stdout) then
            call iplot (mskhr, imt, imt, jmt)
          else
            linemx = 100
            linel  = 105
            line   = linel - 5
            if (line .gt. linemx) line = linemx
            nwr = (imt/line) + 1
c           
            do 900 i=1,imt
              ncol(i) = mod(i,10)
 900        continue
c           
            do 1000 n=1,nwr
              ia = 1 + (line*(n-1))
              ib = ia + line - 1
              if (ib .gt. imt) ib = imt
              if (wrmask) write(nunit,9021) (ncol(i),i=ia,ib)
              if ((rmask) .and. (.not. calculate_mask_internally))
     &             read (nunit,9099)
              do 990 jj=1,jmt
                jjj = jmt - jj + 1
                if (wrmask)write(nunit,9022)jjj,(mskhr(i,jjj),i=ia,ib)
                if ( rmask) then
                  if (.not. calculate_mask_internally) then
                    read (nunit,9022) jr , (mskhr(i,jjj),i=ia,ib)
                  else
                    jr = jjj
                  end if
                  if (jr .ne. jjj) then
                    write (stdout,999) nunit, jjj, jr
                    write (stderr,999) nunit, jjj, jr
                    call abort()
                  endif
                endif
 990          continue
 1000       continue
          endif
c         
          if (wrmask) write(nunit,9031)
          if (wrmask) write(nunit,9032)
     &         (' domain for ver mask #',i,'=',vregnm(i),i=1,nvreg)
          if (wrmask) write(nunit,9034) mskvr
c         
          if ( rmask) then
            if (.not. calculate_mask_internally) read (nunit,9099)
            if (.not. calculate_mask_internally) read (nunit,9099)
          endif
          if (.not. calculate_mask_internally) then
            if ( rmask) read (nunit,9033) (vregnm(i),i=1,nvreg)
            if ( rmask) read (nunit,9034) mskvr
          endif
c         
        else
          if (wrmask) then
            iotext = ' read (nunit) mskhr, mskvr'
            write (nunit) 'no stamp available in reg1st    ', iotext
     &           ,                   expnam       
            write(nunit) mskhr, mskvr
c           
            iotext = ' read (nunit) hregnm, vregnm'
            write (nunit) 'no stamp available in reg1st    ', iotext
     &           ,                   expnam       
            write(nunit) hregnm, vregnm
          endif
          if (rmask) then
            if (.not. calculate_mask_internally) then
              read (nunit)
              read (nunit) mskhr, mskvr
              read (nunit)
              read (nunit) hregnm, vregnm
            endif
          endif
        endif
      endif
c     
c     if vertical masks were read in, set level limits for defining
c     vertical regions in term balance calculations
c
      if (rmask) then
        if (.not. calculate_mask_internally) then
          do 1100 i=1,nvreg
            settop = .false.
            setbot = .false.
            do 1090 k=1,km
              kk = km-k+1
              if (mskvr(k)  .eq. i .and. .not. settop) then
                llvreg(i,1) = k
                settop = .true.
              endif
              if (mskvr(kk) .eq. i .and. .not. setbot) then
                llvreg(i,2) = kk
                setbot = .true.
              endif
 1090       continue
 1100     continue
        endif
      endif
      return
c
 999  format(/' error => bad j-row when reading regionmasks from unit ',
     &     i3,/'   expected',i4,'    read in',i4)
9000  format(/' Horizontal regional volumes [cubic m] and areas [sq m]')
9001  format('    k','  All Regions ',9(1x,i7,5x))
9002  format(1x,i4,10(1x,e12.6))
9003  format('  SUM',10(1x,e12.6))
9004  format(/' AREA',10(1x,e12.6))
9011  format(/' Horizontal regional mask names & domains:')
9012  format(a22,i4,a1,a40)
9013  format(27x,a40)
9021  format('  i=>',100(i1))
9022  format(1x,i3,1x,100(i1))
9031  format(/' Vertical regional mask names & domains:')
9032  format(a22,i4,a1,a20)
9033  format(27x,a20)
9034  format (1x, 42i3)
9099  format(1x)
      end
#else
      return
      end
#endif
