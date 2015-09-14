  module timeinterp_mod
!  The module timeinterp_mod contains a revised version of subroutine
!  timeinterp. The module provides a variable maxinterp  which counts the
!  defined data sets. This allows to call timeinterp independently from 
!  different model components (ocean, river, S.B.C., open lateral boundary 
!  conditions).
!  
!  Changes:
!  A call of timeinterp with n = -1 opens a new data set for interpolation.
!  timeinterp returns the data set index. To open and interpolate a new data 
!  set do something like this
!
!      do m=mstart,mend  
!        n = -1
!        call timeinterp(time, n, ....
!  
!  The return value of n has to be stored for later identification of the 
!  data set.
! 
!        index(m) = n   
!
!  To interpolate a previously opened data set 
!
!      call timeinterp(time, index(i), ....
!
!  index(i) is the index of the data set provided by timeinterp in the 
!  first call. It is valid also if data sets of other model components 
!  have been processed.
!
!  Please note: In FORTRAN90 the type of the procedure parameters is 
!  checked in more detail than in FORTRAN77. tdrec and aprec must be 
!  fields of type real, statements like
!      call timeinterp(time, index, tdrec(1,n) .. 
!  are not accepted by some compilers since tdrec(1,n) is a simple variable
!  of type real. Use 
!      call timeinterp(time, index, tdrec(:,n) .. 
!       

  private
  integer, public :: maxinterp
  data maxinterp /0/      

  public position_within_data, timeinterp
      
  contains
      
      function position_within_data (mom_time, data_start_time&
      ,                              data_end_time, periodic)
      use time_manager_mod
!
!=======================================================================
!     compute where the model time is in relation to starting and ending
!     time for the dataset
!
!     author: r. c. pacanowski      e-mail=> rcp@gfdl.gov
!=======================================================================
!
      type(time_type) ::  mom_time, data_start_time, data_end_time
      type(time_type) ::  temp_time
      logical periodic
!
      if (periodic) then
!
!       perodic assumes a climatology of 12 months
!       use month and day from model time but year from dataset
!
        call get_date (mom_time, iy, imo, id, ih, im, is) 
        call get_date (data_start_time, iy2, imo2, id2, ih2, im2, is2)
	if (imo .eq. 2 .and. id .gt. 28) id = 28
        temp_time = set_date (iy2, imo, id, ih, im, is) 
        call get_time (temp_time - data_start_time, isec, iday)
      else
        if (mom_time .ge. data_start_time .and. &
           mom_time .le. data_end_time) then
          call get_time (mom_time - data_start_time, isec, iday)
        else
	  if (mom_time .lt. data_start_time) then
            call get_time (data_start_time - data_start_time, isec,iday)
	  endif
	  if (mom_time .gt. data_end_time) then
            call get_time (data_end_time - data_start_time, isec, iday)
	  endif
	endif
      endif
      position_within_data = iday + isec/86400.0
      return
      end function position_within_data 
       
      
      subroutine timeinterp (tm, n, tdrec, aprec, ndr, periodic, method&
      ,                     ia, ib, wb, change, inext, iprev)
!
!=======================================================================
!
!     time interpolator ... constructs indices & weight needed for
!     linearly interpolating data defined at arbitrary time intervals
!     (midpoints of years, months, days or  random intervals) to
!     the time of the current model time step.
!
!     inputs:
!
!     tm     = the time at which the data is desired (units of "tdrec")
!
!     tdrec  = the times at which the data records in the dataset are
!              defined. times must be monotonically increasing and are
!              assumed to be at the centers of the averaging periods.
!              (eg: the centers of the months if using monthly averaged
!               climatology. units are arbitrary)
!
!     aprec  = array of averaging periods for the data records
!              (eg: the number of days per month)
!
!     ndr    = number of data records in the dataset. (eg: 12 if using
!              monthly climatology)
!
!     periodic = (true,false) if the dataset is to be treated as
!              (periodic, not periodic). if periodic, then the model
!               time is always mapped into the dataset. It is assumed
!               the dataset is a one year climatology. if not periodic,
!               then all model times < the midpoint of the 1st record 
!               in the dataset use the midpoint of the 1st record and
!               all model times > midpoint of the last record "ndr" use
!               the midpoint of the last record in the dataset.
!               
!     method = interpolation scheme desired.  (0..3)
!                0 = no interpolation; the average value is used
!                    for all times in the entire averaging period.
!                    (preserves the integral over averaging periods,
!                    but is discontinuous at period boundaries.)
!                1 = linear interpolation between the middles of
!                    two adjacent averaging periods.
!                    (continuous but does not preserve integral for
!                    unequal periods.)
!                2 = equal linear interpolation.  Assumes that the
!                    value on the boundary between two adjacent
!                    averaging periods is the unweighted average of
!                    the two average values.  Linearly interpolates
!                    between the midperiod and period boundary.
!                    (continuous but does not preserve integral for
!                    unequal periods.)
!                3 = equal area (midperiod to midperiod) interpolation
!                    chooses a value for the boundary between two
!                    adjacent periods such that linear interpolation 
!                    between the two midperiods and this value will
!                    preserve the integral midperiod to midperiod.
!                Note that methods 1,2, and 3 are equivalent if
!                all periods lengths are equal.
!
!     maxinterp= number of previosly defined datasets
!     n      = a number denoting which dataset is being interpolated
!              (each dataset must be referenced by a unique number
!               starting with 1 for the 1st, 2 for the 2nd, ...etc)
!
!              If timeinterp is called with n=-1, a new dataset is 
!              assumed. maxinterp is incremented and n=maxinterp is the
!              returned value.
!                
!              If timeinterp is called with n=1..maxinterp a previously 
!              defined dataset is assumed. n is not changed and the
!              interpolation is performed.
!
!
!     outputs:
!
!     n      = a number denoting which dataset is being interpolated
!              See input variable n
!     ia     = index for pointing to the next data record which will be
!              reached by the model. (eg: ahead of the model. "ia" would
!              be 3 if "tm" was beyond the  middle of {but still within}
!              february)
!     ib     = index for pointing to the data record which was just
!              passed by the model. (eg: behind the model. "ib" would
!              be 2 if "tm" was beyond the middle of {but still within}
!              february)
!     inext  = index to memory buffer containing data from "ia"
!     iprev  = index to memory buffer containing data from "ib"
!     wb     = interpolation weight for defining data at "tm"
!              schematically the interpolation is defined by:
!
!              data(iprev) <== disk data "ib"              
!              data(inext) <== disk data "ia"
!              data(tm) = wb*data(iprev) + (1-wb)*data(inext)
!
!     change = logical for sensing when "ia" and "ib" change.
!              when change = T then it is time to read the disk
!              and update "inext" and "iprev"
!     
!     author:        r. c. pacanowski      e-mail=> rcp@gfdl.gov
!     revised by:    d. j. goldberg        e-mail=> chg@gfdl.gov
!     revised by:    m. schmidt e-mail=>martin.schmidt@io-warnemuende.de
!=======================================================================
!
      parameter (iflag=-99999)
      logical change, periodic
# include "stdunits.h"
      integer, allocatable :: imethod(:),itemp(:), iaold(:)
      integer              :: method, ndr, ia, ib, inext, iprev
      real, intent(in)     :: tdrec(:), aprec(:)
      real, intent(out)    :: wb
      save imethod, iaold
!
!-----------------------------------------------------------------------
!     statement function
!-----------------------------------------------------------------------
!
      frac(r) = r - int(r)
!
!!      if (n .gt. maxsets) then
!!        write (stdout,'(a,i10,a,i10)') 'Error: n=', n&
!!      ,          ' maxsets=',maxsets
!!        call abort()
!!      end if
!
!     call for a dataset which is not defined
!     The first call

      if (maxinterp.eq.0) then
        if(n.le.-1) then            ! Initialize timeinterp
          write (stdout,*) ' => Initializing timeinterp'
	  allocate(imethod(1),iaold(1)) 
        endif
      endif
      if (n.eq.-1) then ! Define a new dataset
        maxinterp = maxinterp + 1
	n = maxinterp
	if (maxinterp.gt.1) then    ! Save previously defined datasets
	  allocate (itemp(maxinterp))
	  do nd=1, maxinterp-1
	    itemp(nd) = imethod(nd) 
	  enddo
          write (stdout,'(/1x,a,i2,a,i3/)')&
           'Assigning interpolation method ',method, ' to dataset # ',n
	  deallocate(imethod)
	  allocate(imethod(maxinterp))
	  imethod(:) = itemp(:)
!
	  do nd=1, maxinterp-1
	    itemp(nd) = iaold(nd) 
	  enddo
	  deallocate(iaold)
	  allocate(iaold(maxinterp))
	  iaold(:) = itemp(:)
	  deallocate (itemp)
	endif
	imethod(maxinterp) = method
	iaold  (maxinterp) = iflag
      else
        if(n.gt.maxinterp) then
	  write (stdout,'(/a,i2,a/,a,i2,a/)')&
          '=>Error: trying to define or interpolate dataset ',n,'.'&
      ,   ' However, only ',maxinterp,' datasets have been defined.'
          call abort()
	endif      
      endif
!
      if (method .ne. imethod(n)) then
        write (stdout,'(/a,i2,a,i3/a,i2,a/)')&
        'Error: trying to use method ',method, ' on dataset # ',n&
      ,  'originally, method ',imethod(n),' was used in timeinterp'
        call abort()
      endif 
!
      if (periodic) then
!
!       define the position of the dataset in time
!      
        dstart = tdrec(1) - 0.5*aprec(1)
        dend   = tdrec(ndr) + 0.5*aprec(ndr)
        dlen   = dend - dstart
!
!       map the model time into the dataset assuming dataset periodicity
!
        if (tm .lt. dstart) then
          d = dstart - tm
          f = frac(d/dlen)
          time = dend - f*dlen
        elseif (tm .gt. dend) then
          d = tm - dend
          f = frac(d/dlen)
          time = dstart + f*dlen
        else
          time = tm
        end if
      else
!
!       define the position of the dataset in time. no periodicity
!      
        dstart = tdrec(1)
        dend   = tdrec(ndr)
        dlen   = dend - dstart
!
!       map the model time into the dataset. assume data is constant
!       before the beginning and after the end of the dataset
!
        if (tm .lt. dstart) then
          time = dstart
        elseif (tm .gt. dend) then
          time = dend 
        else
          time = tm
        end if
      end if
!
!     calculate record pointers and weighting for interpolation of
!     dataset records to the model time step.
!
      ib = indp (time, tdrec, ndr)
      if (tdrec(ib) .gt. time) ib = ib - 1
      if (periodic) then
        ia = mod(ib, ndr) + 1
        if (ib .lt. 1) ib = ndr
      else
        ia = ib + 1
        if (ia .gt. ndr) ia = ib
        if (ib .lt. 1)   ib = ia
      end if
!
!     find whether "time" is closer to midpoint of record "ia" or "ib"
!     ic is the index of the closest midpoint
!     io is the index of the other midpoint
!
      startaft = tdrec(ia) - 0.5*aprec(ia)
      if (time .ge. startaft .and. time .le. tdrec(ia)) then
        ic = ia
        io = ib
      else
        ic = ib
        io = ia
      end if
!
!     dtmid = distance from "time" to midpoint of closer record
!     dtbnd = distance from "time" to boundary of closer record
!     dtomid = distance from "time" to midpoint of other record
!
      dtmid  = abs(time - tdrec(ic))
      dtbnd  = 0.5*aprec(ic) - dtmid
      dtomid = 0.5*aprec(io) + dtbnd
!
!-----------------------------------------------------------------------
!     3) equal area (midperiod to midperiod) interpolation formula
!-----------------------------------------------------------------------
!
      if (method .eq. 3) then
        wc = 2.0*dtbnd/aprec(ic) + 2.0*dtmid/(aprec(ic) + aprec(io))
!
!-----------------------------------------------------------------------
!     2) equal linear interpolation
!             value on period boundary assumed to be average of values
!             on the two adjacent periods.
!-----------------------------------------------------------------------
!
      elseif (method .eq. 2) then
        wc = (2.0*dtbnd + dtmid)/aprec(ic)
!
!-----------------------------------------------------------------------
!     1) linear interpolation
!-----------------------------------------------------------------------
!
      elseif (method .eq. 1) then
        wc = dtomid/(dtmid + dtomid)
!
!-----------------------------------------------------------------------
!     0) no interpolation
!-----------------------------------------------------------------------
!
      elseif (method .eq. 0) then
        wc = 1.0
      else
!
!-----------------------------------------------------------------------
!     anything else is not allowed for (unless you want to add one!)
!-----------------------------------------------------------------------
!
        write (stdout,*)&
         '=>Error: method = ',method,' not allowed in timeinterp'
        call abort()
      end if
!
      if (ib .eq. ic) then
        wb = wc
      else
        wb = 1.0 - wc
      end if
      if (wc .lt. 0.0 .or. wc .gt. 1.0) then
        write (stdout,*) ' ic=',ic,' io=',io, ' dtmid=',dtmid,' dtbnd='&
      , dtbnd,' dtomid=',dtomid, ' time=',time, ' ia=',ia,' ib=',ib&
      ,' wc=',wc
        write (stdout,*)& 
      ' =>Error: bad interpolation wieght in timeinterp'
        call abort()
      endif
!
!     refresh pointers to memory buffers when reading disk data
!
      if (iaold(n) .ne. ia) then
        change = .true.
        ite   = iprev
        iprev = inext
        inext = ite
      else
        change = .false.
      end if
      iaold(n) = ia
!
#ifdef debug_timeinterp
      write (stdout,*) 'timeinterp: ic=',ic,' io=',io, ' dtmid=',dtmid&
      ,' dtbnd=', dtbnd,' dtomid=',dtomid, ' time=',time, ' ia=',ia&
      ,' ib=',ib, ' wc=',wc, ' change=',change, ' iprev=',iprev&
      , ' inext=',inext
#endif
      return
      end subroutine timeinterp
      
      end module timeinterp_mod
      
