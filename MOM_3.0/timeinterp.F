#ifdef test_timeinterp
      program driver
c
c=======================================================================
c     test the time interpolator module
c
c     author: R. C. Pacanowski        e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use time_manager_mod
#include "stdunits.h"
#include "switch.h"
#include "tmngr.h"
      logical periodic, readit
      character*12 recname
      parameter (nrec=12)
      character*32 dstamp(0:nrec)
      dimension tdrec(nrec), aprec(nrec), recname(nrec)
      type(time_type), dimension(nrec+1) :: a_time
      type(time_type) :: init_time, data_start_time, data_end_time
      type(time_type) :: temp_time
c
      stdout = 6
c
c-----------------------------------------------------------------------
c     set type of calendar
c
c     if (eqyear and eqmon) then thirty_day_months
c     if (eqyear and !eqmon) then no_leap (feb always has 28 days)
c     if (!eqyear and !eqmon) then julian
c-----------------------------------------------------------------------
c
      eqyear   = .true.
      eqmon    = .false.
      call set_calendar (eqyear, eqmon)
c
c-----------------------------------------------------------------------
c     specify an arbitrary time for initial conditions
c-----------------------------------------------------------------------
c
      year0  = 1
      month0 = 1
      day0   = 1
      hour0  = 0
      min0   = 0
      sec0   = 0
      init_time = set_date (year0, month0, day0, hour0, min0, sec0)
      call set_stamp (ic_stamp, year0, month0, day0, hour0, min0, sec0)
c
c-----------------------------------------------------------------------
c     decide on reference date for calculating logical switches
c
c     refrun  = T ==> reference to the start of each run
c     refinit = T ==> reference to initial condition time given by:
c                     year0, month0, day0, hour0, min0, sec0
c     refuser = T ==> reference to user specified time given by:
c                     ryear, rmonth, rday, rhour, rmin, rsec
c
c     choose by setting one of the above to TRUE
c-----------------------------------------------------------------------
c
      refrun  = .false.
      refinit = .true.
      refuser = .false.
c
c     if refuser = T then set date for referencing periodic switches
c
      if (refuser) then
        ryear  = 1900
        rmonth = 1
        rday   = 1
        rhour  = 0
        rmin   = 0
        rsec   = 0
      end if
c
c-----------------------------------------------------------------------
c     set the time step "dt" (sec) and how long to integrate "runlen".
c     "rununits" must be 'days'
c-----------------------------------------------------------------------
c
      dt       = 3600.0
      runlen   = 1.2*365.0
      rununits = 'days'
c
c-----------------------------------------------------------------------
c     for use in MOM, the integration time "runlen" is broken into a
c     number of segments each of length "segtim" days. updated surface
c     boundary conditions are applied to MOM every "segtim" days. this
c     is useful when coupling to atmospheric models in which case
c     both models exchange surface boundary conditions every "segtim"
c     days where "segtim" is 1/(coupling frequency). when MOM is without
c     an atmospheric model, "segtim" is set to the time step (in days)
c     as done here.     
c-----------------------------------------------------------------------
c
      segtim = dt/86400.0
c
c     set number of time steps between mixing time steps
c
      nmix      = 17
c
      write (stdout,'(///,tr5,a,g14.7,a)')
     &'Driver will test the time interpolator by simulating integrating'
     &, runlen, rununits
      write (stdout,*) '     using timesteps of dt = ',dt,' (sec) '
c
c-----------------------------------------------------------------------
c     set an arbitrary starting time for the model
c-----------------------------------------------------------------------
c
      model_year  = 1992
      model_month = 1
      model_day   = 1
      model_hour  = 0
      model_min   = 0
      model_sec   = 0
      model_time = set_date (model_year, model_month, model_day
     &,                      model_hour, model_min, model_sec)
      call get_time (model_time, iseconds, idays)
c
c     set corresponding days and milliseconds for calendar
c     (typically iday and isec would be read from a restart)
c
      iday = idays
      isec = iseconds*1000
      itt  = 0
c
c-----------------------------------------------------------------------
c     initialize the time manager with specified initial condition time,
c     user reference time, and model time (from iday and isec).
c-----------------------------------------------------------------------
c
      call tmngri (year0, month0, day0, hour0, min0, sec0
     &,            ryear, rmonth, rday, rhour, rmin, rsec
     &,            iday, isec
     &,            runlen, rununits, rundays, dt)
c
c-----------------------------------------------------------------------
c     test interpolating data from a climatological monthly mean 
c     dataset to current model time (at time step "itt")
c
c     define the dataset as periodic and choose an interpolation
c     method.  see subroutine "timeinterp" included within this module
c     to explain the various methods of interpolation.
c-----------------------------------------------------------------------
c
      periodic = .true.
      method   = 3
c
c     Build a 12 month climatological dataset.
c     Note that the time stamp marks the end of each month. (1st = jan)
c     dstamp(m) = time at end of month m
c     a_time(m) = time structure at beginning of month m
c     aprec(m)  = number of days in month m
c     tdrec(m)  = center of month m in days counting from a_time(1)
c
      iyr  = 1
      imon = 1
      iday = 1
      ihr  = 0
      imin = 0
      isec = 0
      a_time(1) = set_date (iyr, imon, iday, ihr, imin, isec)
      call set_stamp (dstamp(0), iyr, imon, iday, ihr, imin, isec) 
      aprec(1)  = days_in_month(a_time(1))
      tdrec(1)  = 0.5*aprec(1)
c
      do m=2,nrec+1
	recname(m-1) = monname(imon)
        imon = imon + 1
	if (imon .gt. 12) then
	  imon = 1
	  iyr  = iyr + 1
	endif
        a_time(m) = set_date (iyr, imon, iday, ihr, imin, isec)
        call set_stamp (dstamp(m-1), iyr, imon, iday, ihr, imin, isec)
	if (m .le. nrec) then
          aprec(m)  = days_in_month(a_time(m))
	  tdrec(m)  = tdrec(m-1) + 0.5*(aprec(m-1) + aprec(m))
	endif
      enddo
c
c     construct time structures to mark bounds of dataset
c
      data_start_time = a_time(1)
      data_end_time   = a_time(nrec+1)
c
      do m=1,nrec
	write (stdout,'(a,i3,3a,f10.4,4a)') ' data record # ',m
     &, ' is ',recname(m), ' centered at ',tdrec(m), 'days. from '
     &, dstamp(m-1),' to ',dstamp(m)
      enddo
c
c     express model time as number of days "daysid" into dataset
c
      daysid = position_within_data (model_time, data_start_time
     &,                              data_end_time, periodic)
c
c     initialize pointers to memory buffers
c      iprevm will point to memory holding data before model timestep
c      inextm will point to memory holding data after model timestep
c     also, set the number of datasets being interpolated to 1
c     
      iprevm = 1
      inextm = 2
      ntdrec = 1
c
      call timeinterp (daysid, ntdrec, tdrec, aprec, nrec, periodic
     &,                method, inextd, iprevd, wprev, readit, inextm
     &,                iprevm)
c
c     read data from disk records "inextd" and "iprevd" into
c     memory buffers "inextm" and "iprevm" to initialize them
c
      write (stdout,'(/a,a,a,i2/a,a,a,i2/)')  
     &  ' read disk data record for ', recname(iprevd)
     &, ' into memory buffer ',iprevm 
     &, ' and disk data record for ', recname(inextd)
     &, ' into memory buffer ', inextm
          write (stdout,'(a,i3,a,i3,a,f10.8,a,i2,a,i2)')
     &    ' inextd=',inextd,' iprevd=',iprevd, ' wprev=',wprev
     &,   ' inextm=',inextm,' iprevm=',iprevm
c
c-----------------------------------------------------------------------
c     integrate the time manager for the required time
c     and print when things will happen
c-----------------------------------------------------------------------
c
      do n=1,1000000
c
        itt = itt + 1
        call tmngr (dt)
	call set_time_switches
c
c-----------------------------------------------------------------------
c       test for when to read monthly S.B.C. (surface
c       boundary conditions) for use in interpolating them to the   
c       current model time.
c-----------------------------------------------------------------------
c
        daysid = position_within_data (model_time, data_start_time
     &,                                data_end_time, periodic)
c
        call timeinterp (daysid, ntdrec, tdrec, aprec, nrec, periodic
     &,                  method, inextd, iprevd, wprev, readit, inextm
     &,                  iprevm)
c
c-----------------------------------------------------------------------
c       print when things happen
c-----------------------------------------------------------------------
c
        if (itt .le. 10) then
          write (stdout,9000) itt, stamp, relday, daysid
        endif
c
c       show when to update the memory buffers holding monthly S.B.C.
c       by reading data from disk record "inextd"  into
c       memory buffers "inextm"
c
        if (readit) then
          write (stdout,'(/a,i7,a,g14.7,a,a,a,i2,a,i1,a)')  
     &     ' ==> time to read data on itt=',itt, ', day= ',relday
     &,    ' from ', recname(inextd), ' into buffer ', inextm
     &,   ' (previous data is now in ',iprevm,')'
c            write (stdout,'(1x,a,i5,a,i2,a,i2,a,e14.7,a,i2,a,i2)')
c     &    ' itt=',itt,' inextd=',inextd,' iprevd=',iprevd
c     &,   ' wprev=',wprev, ' inextm=',inextm,' iprevm=',iprevm
        endif
	if (wprev .lt. 0.005 .or. wprev .gt. 0.995) then
	  write (stdout,'(1x,a,i7,a,a,a,f12.8,a,f12.8)')
     &      ' itt=',itt, ' ',stamp, ' weight for (prev,next) data ='
     &,      wprev,',',1.0-wprev
	endif	
c
        call get_date (model_time, iy, im, id, ih, imin, is)
        if (eomon) write (stdout,9001) ' ...end of '//trim(monname(im))
     &//' on itt=', itt, ' ', stamp
        if (midmon) write (stdout,9001) ' ... mid '//trim(monname(im))
     &//' on itt  =',itt, ' ', stamp
        if (eoyear) write (stdout,*) ' ==> end of year ',iy,' on itt='
     &,itt,' ', stamp
c
        if (eorun)  then
	  write (stdout,9001) ' ==> end of run on itt=',itt,' ', stamp
          call abort()
	endif
c
      enddo
      write (stdout,*) ' Increase loop limit on "n" in "driver.F"'
      call abort()
9000  format (1x, ' itt=',i7,1x,a,' model days since I.C.=',f12.4
     &,      '. days since start of S.B.C.=', f12.4)
9001  format (a, i8, a, a)
      end
#endif

      function position_within_data (mom_time, data_start_time
     &,                              data_end_time, periodic)
      use time_manager_mod
c
c=======================================================================
c     compute where the model time is in relation to starting and ending
c     time for the dataset
c
c     author: r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      type(time_type) ::  mom_time, data_start_time, data_end_time
      type(time_type) ::  temp_time
      logical periodic
c
      if (periodic) then
c
c       perodic assumes a climatology of 12 months
c       use month and day from model time but year from dataset
c
        call get_date (mom_time, iy, imo, id, ih, im, is) 
        call get_date (data_start_time, iy2, imo2, id2, ih2, im2, is2)
	if (imo .eq. 2 .and. id .gt. 28) id = 28
        temp_time = set_date (iy2, imo, id, ih, im, is) 
        call get_time (temp_time - data_start_time, isec, iday)
      else
        if (mom_time .ge. data_start_time .and. 
     &      mom_time .le. data_end_time) then
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
      end
       


      subroutine timeinterp (tm, n, tdrec, aprec, ndr, periodic, method
     &,                     ia, ib, wb, change, inext, iprev)
c
c=======================================================================
c
c     time interpolator ... constructs indices & weight needed for
c     linearly interpolating data defined at arbitrary time intervals
c     (midpoints of years, months, days or  random intervals) to
c     the time of the current model time step.
c
c     inputs:
c
c     tm     = the time at which the data is desired (units of "tdrec")
c
c     tdrec  = the times at which the data records in the dataset are
c              defined. times must be monotonically increasing and are
c              assumed to be at the centers of the averaging periods.
c              (eg: the centers of the months if using monthly averaged
c               climatology. units are arbitrary)
c
c     aprec  = array of averaging periods for the data records
c              (eg: the number of days per month)
c
c     ndr    = number of data records in the dataset. (eg: 12 if using
c              monthly climatology)
c
c     periodic = (true,false) if the dataset is to be treated as
c              (periodic, not periodic). if periodic, then the model
c               time is always mapped into the dataset. It is assumed
c               the dataset is a one year climatology. if not periodic,
c               then all model times < the midpoint of the 1st record 
c               in the dataset use the midpoint of the 1st record and
c               all model times > midpoint of the last record "ndr" use
c               the midpoint of the last record in the dataset.
c               
c     method = interpolation scheme desired.  (0..3)
c                0 = no interpolation; the average value is used
c                    for all times in the entire averaging period.
c                    (preserves the integral over averaging periods,
c                    but is discontinuous at period boundaries.)
c                1 = linear interpolation between the middles of
c                    two adjacent averaging periods.
c                    (continuous but does not preserve integral for
c                    unequal periods.)
c                2 = equal linear interpolation.  Assumes that the
c                    value on the boundary between two adjacent
c                    averaging periods is the unweighted average of
c                    the two average values.  Linearly interpolates
c                    between the midperiod and period boundary.
c                    (continuous but does not preserve integral for
c                    unequal periods.)
c                3 = equal area (midperiod to midperiod) interpolation
c                    chooses a value for the boundary between two
c                    adjacent periods such that linear interpolation 
c                    between the two midperiods and this value will
c                    preserve the integral midperiod to midperiod.
c                Note that methods 1,2, and 3 are equivalent if
c                all periods lengths are equal.
c
c     n      = a number denoting which dataset is being interpolated
c              (each dataset should be referenced by a unique number
c               starting with 1 for the 1st, 2 for the 2nd, ...etc)
c
c
c     outputs:
c
c     ia     = index for pointing to the next data record which will be
c              reached by the model. (eg: ahead of the model. "ia" would
c              be 3 if "tm" was beyond the  middle of {but still within}
c              february)
c     ib     = index for pointing to the data record which was just
c              passed by the model. (eg: behind the model. "ib" would
c              be 2 if "tm" was beyond the middle of {but still within}
c              february)
c     inext  = index to memory buffer containing data from "ia"
c     iprev  = index to memory buffer containing data from "ib"
c     wb     = interpolation weight for defining data at "tm"
c              schematically the interpolation is defined by:
c
c              data(iprev) <== disk data "ib"              
c              data(inext) <== disk data "ia"
c              data(tm) = wb*data(iprev) + (1-wb)*data(inext)
c
c     change = logical for sensing when "ia" and "ib" change.
c              when change = T then it is time to read the disk
c              and update "inext" and "iprev"
c     
c     author:        r. c. pacanowski      e-mail=> rcp@gfdl.gov
c     revised by:    d. j. goldberg        e-mail=> chg@gfdl.gov
c=======================================================================
c
      logical change, periodic
# include "stdunits.h"
      parameter (maxsets=15, iflag=-99999)
      dimension iaold(maxsets), tdrec(ndr), aprec(ndr)
      dimension imethod(maxsets)
      data iaold /maxsets*iflag/
      save iaold, imethod
c
c-----------------------------------------------------------------------
c     statement function
c-----------------------------------------------------------------------
c
      frac(r) = r - int(r)
c
      if (n .gt. maxsets) then
        write (stdout,'(a,i10,a,i10)') 'Error: n=', n
     &,          ' maxsets=',maxsets
        call abort()
      end if
c
      if (iaold(n) .eq. iflag) then
        write (stdout,'(/1x,a,i2,a,i3/)')
     &      'Assigning interpolation method ',method, ' to dataset # ',n
        imethod(n) = method
      endif 
c
      if (method .ne. imethod(n)) then
        write (stdout,'(/a,i2,a,i3/a,i2,a/)')
     &   'Error: trying to use method ',method, ' on dataset # ',n
     &,  'originally, method ',imethod(n),' was used in timeinterp'
        call abort()
      endif 
c
      if (periodic) then
c
c       define the position of the dataset in time
c      
        dstart = tdrec(1) - 0.5*aprec(1)
        dend   = tdrec(ndr) + 0.5*aprec(ndr)
        dlen   = dend - dstart
c
c       map the model time into the dataset assuming dataset periodicity
c
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
c
c       define the position of the dataset in time. no periodicity
c      
        dstart = tdrec(1)
        dend   = tdrec(ndr)
        dlen   = dend - dstart
c
c       map the model time into the dataset. assume data is constant
c       before the beginning and after the end of the dataset
c
        if (tm .lt. dstart) then
          time = dstart
        elseif (tm .gt. dend) then
          time = dend 
        else
          time = tm
        end if
      end if
c
c     calculate record pointers and weighting for interpolation of
c     dataset records to the model time step.
c
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
c
c     find whether "time" is closer to midpoint of record "ia" or "ib"
c     ic is the index of the closest midpoint
c     io is the index of the other midpoint
c
      startaft = tdrec(ia) - 0.5*aprec(ia)
      if (time .ge. startaft .and. time .le. tdrec(ia)) then
        ic = ia
        io = ib
      else
        ic = ib
        io = ia
      end if
c
c     dtmid = distance from "time" to midpoint of closer record
c     dtbnd = distance from "time" to boundary of closer record
c     dtomid = distance from "time" to midpoint of other record
c
      dtmid  = abs(time - tdrec(ic))
      dtbnd  = 0.5*aprec(ic) - dtmid
      dtomid = 0.5*aprec(io) + dtbnd
c
c-----------------------------------------------------------------------
c     3) equal area (midperiod to midperiod) interpolation formula
c-----------------------------------------------------------------------
c
      if (method .eq. 3) then
        wc = 2.0*dtbnd/aprec(ic) + 2.0*dtmid/(aprec(ic) + aprec(io))
c
c-----------------------------------------------------------------------
c     2) equal linear interpolation
c             value on period boundary assumed to be average of values
c             on the two adjacent periods.
c-----------------------------------------------------------------------
c
      elseif (method .eq. 2) then
        wc = (2.0*dtbnd + dtmid)/aprec(ic)
c
c-----------------------------------------------------------------------
c     1) linear interpolation
c-----------------------------------------------------------------------
c
      elseif (method .eq. 1) then
        wc = dtomid/(dtmid + dtomid)
c
c-----------------------------------------------------------------------
c     0) no interpolation
c-----------------------------------------------------------------------
c
      elseif (method .eq. 0) then
        wc = 1.0
      else
c
c-----------------------------------------------------------------------
c     anything else is not allowed for (unless you want to add one!)
c-----------------------------------------------------------------------
c
        write (stdout,*)
     &   '=>Error: method = ',method,' not allowed in timeinterp'
        call abort()
      end if
c
      if (ib .eq. ic) then
        wb = wc
      else
        wb = 1.0 - wc
      end if
      if (wc .lt. 0.0 .or. wc .gt. 1.0) then
        write (stdout,*) ' ic=',ic,' io=',io, ' dtmid=',dtmid,' dtbnd='
     &, dtbnd,' dtomid=',dtomid, ' time=',time, ' ia=',ia,' ib=',ib
     &, ' wc=',wc
        write (stdout,*) 
     & ' =>Error: bad interpolation wieght in timeinterp'
        call abort()
      endif
c
c     refresh pointers to memory buffers when reading disk data
c
      if (iaold(n) .ne. ia) then
        change = .true.
        itemp = iprev
        iprev = inext
        inext = itemp
      else
        change = .false.
      end if
      iaold(n) = ia
c
#ifdef debug_timeinterp
      write (stdout,*) 'timeinterp: ic=',ic,' io=',io, ' dtmid=',dtmid
     &,' dtbnd=', dtbnd,' dtomid=',dtomid, ' time=',time, ' ia=',ia
     &,' ib=',ib, ' wc=',wc, ' change=',change, ' iprev=',iprev
     &, ' inext=',inext
#endif
      return
      end


