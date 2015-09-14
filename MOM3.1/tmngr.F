      subroutine tmngr (dt)
c
c=======================================================================
c     increment model time by "dt" seconds and calculate time elapsed 
c     since the specified base time. 
c     time step "dt" should be an integral number of seconds
c
c     input:
c
c       dt = integral number of seconds in a time step
c
c     output:
c
c       model_time = a time structure corresponding to the model time
c       time_since_base = a time structure for time since a base time 
c
c     author: R.C.Pacanowski        e-mail => rcp@gfdl.gov
c=======================================================================
c
      use time_manager_mod
c
      model_time = model_time + set_time (nint(dt), 0)
      time_since_base = model_time - base_time
c
      return
      end



      subroutine tmngri (icyear, icmonth, icday, ichour, icmin, icsec
     &,                  rfyear, rfmonth, rfday, rfhour, rfmin, rfsec
     &,                  idayrestart, msrestart
     &,                  runlen0, rununits0, rundays0, timestep)
c
c=======================================================================
c     initialize calendar and internal time variables for initial 
c     conditions, user reference, and model time at start of run.
c
c     author: R.C.Pacanowski        e-mail => rcp@gfdl.gov
c=======================================================================
c
      use time_manager_mod
      type(time_type) ::  time_since_ic, user_time
#include "stdunits.h"
#include "switch.h"
#include "tmngr.h"
      character*(*) rununits0
c
      integer sday, smon, syr, ssec, smin, shrs
      integer years, months, days, hours, mins, secs
      integer icyear, icmonth, icday, ichour, icmin, icsec
      integer rfyear, rfmonth, rfday, rfhour, rfmin, rfsec
      integer idayrestart, msrestart
      character*32 user_stamp
#ifdef use_netCDF
      character*20 type_of_calendar
#endif
c
      write (stdout,'(//,10x,a,/)') 'Time manager initialization'
c
c-----------------------------------------------------------------------
c     set logical to denote first time step of run
c-----------------------------------------------------------------------
c
      first = .true.
c
c-----------------------------------------------------------------------
c     set the calendar type.
c-----------------------------------------------------------------------
c
      call set_calendar (eqyear, eqmon)
      yrlen = 365.0
c
c-----------------------------------------------------------------------
c     set the names for days and months
c-----------------------------------------------------------------------
c
      dayname(1) = 'Sunday'
      dayname(2) = 'Monday'
      dayname(3) = 'Tuesday'
      dayname(4) = 'Wednesday'
      dayname(5) = 'Thursday'
      dayname(6) = 'Friday'
      dayname(7) = 'Saturday'
c
      monname(1) = 'January'
      monname(2) = 'Febuary'
      monname(3) = 'March'
      monname(4) = 'April'
      monname(5) = 'May'
      monname(6) = 'June'
      monname(7) = 'July'
      monname(8) = 'August'
      monname(9) = 'September'
      monname(10) = 'October'
      monname(11) = 'November'
      monname(12) = 'December'
c
c-----------------------------------------------------------------------
c     set the initial condition time as a time structure
c-----------------------------------------------------------------------
c
      if (icyear .lt. 1) then
        write (stdout,*) '=>Error: icyear must be >= 1 but is ',icyear
        call abort()
      endif
      ic_time = set_date (icyear, icmonth, icday, ichour, icmin, icsec)
      call set_stamp (ic_stamp, icyear, icmonth, icday, ichour, icmin
     &,               icsec)
      write (stdout,'(a,a)') 'Initial Condition time stamp = ',ic_stamp
c
c-----------------------------------------------------------------------
c     set the model time as a time structure
c-----------------------------------------------------------------------
c
      if (init) then
        model_time = ic_time
      else
        model_time = set_time (msrestart/1000, idayrestart)
      endif
      call get_date (model_time, years, months, days, hours, mins, secs)
      call set_stamp (stamp, years, months, days, hours, mins, secs)
      pstamp = stamp
c
c-----------------------------------------------------------------------
c     compute the time at the end of the run as a time structure
c-----------------------------------------------------------------------
c
      if (rununits0 .eq. 'days') then
c
c       when units of "runlen0" are days, allow for fractional days
c
        rundays0 = runlen0
        days = int(runlen0)
        secs = (runlen0-days)*86400.
        stop_time = set_time(secs,days) + model_time
        call get_date (stop_time, syr, smon, sday, shrs, smin, ssec)
        call set_stamp (stop_stamp, syr, smon, sday, shrs, smin, ssec)
      elseif (rununits0 .eq. 'months') then
c
c       when units of "runlen0" are months, stop on a month boundary
c
        call get_date (model_time, iy, imon, id, ih, imin, isec)
	do n=1,nint(runlen0)
	  imon = imon + 1
	  if (imon .gt. 12) then
	    imon = 1
	    iy = iy + 1
	  endif
	enddo
	stop_time = set_date (iy, imon, 1, 0, 0, 0)
	call set_stamp (stop_stamp, iy, imon, 1, 0, 0, 0)
	call get_time (stop_time - model_time, isec, idays)
	rundays0 = idays + isec/86400.0 
      elseif (rununits0 .eq. 'years') then
c
c       when units of "runlen0" are years, stop on a month boundary
c
        call get_date (model_time, iy, imon, id, ih, imin, isec)
	iy = iy + nint(runlen0)
	stop_time = set_date (iy, imon, 1, 0, 0, 0)
	call set_stamp (stop_stamp, iy, imon, 1, 0, 0, 0) 
	call get_time (stop_time - model_time, isec, idays)
	rundays0 = idays + isec/86400.0 
      else
        write (stdout,*)
     & '=>Error: specified integration time units = ',trim(rununits0)
     &,'         units must be days, months, or years'
        call abort()
      endif
c
c-----------------------------------------------------------------------
c     calculate time since initial conditions and other goodies 
c-----------------------------------------------------------------------
c
      time_since_ic = model_time - ic_time
      call get_time (time_since_ic, secs, days)
      relday = days + secs/86400.
      relyr  = relday/yrlen
      prelyr = relyr
      write (stdout,'(/a,i10,a,a,a,1pe14.7,a,e14.7/a,a,/)')
     &' Time step (itt) =',itt,', time stamp=',stamp
     &,', days since IC=',relday,', years since IC=',relyr
     &,' model will stop on ',stop_stamp
      call get_time (model_time - set_date(years,1,1), secs, days)
      dayoyr = days + secs/86400.
c
c-----------------------------------------------------------------------
c     make sure time step is an integral number of seconds 
c-----------------------------------------------------------------------
c
      frac = timestep - nint(timestep)
      if (frac .ne. 0.0) then
      write (stdout,'(/a/a,e14.7,a,i10/)')
     & '-->Warning: the time step must be an integral number of seconds'
     &,'            resetting from ',timestep, ' to ', nint(timestep)
      endif 
c
c----------------------------------------------------------------------
c     select reference time for computing logical switches
c----------------------------------------------------------------------
c
      irefs = 0
      if (refrun) then
        base_time = model_time
        irefs = irefs + 1
        write (stdout,'(/a,a,a/)') ' "refrun = .true." was specified:'
     &,  ' =>logical switch calculations will be referenced to '
     &,  'the beginning of each run (each job submission).'
      endif
      if (refinit) then
        base_time = ic_time
        irefs = irefs + 1
        write (stdout,'(/a,a,a/)') ' "refinit = .true." was specified:'
     &,  ' =>logical switch calculations will be referenced to '
     &,  'the Initial Condition time.'
      endif
      if (refuser) then
        user_time = set_date (rfyear, rfmonth, rfday, rfhour, rfmin
     &,                       rfsec)
        call set_stamp (user_stamp, rfyear, rfmonth, rfday, rfhour,rfmin
     &,                       rfsec)
        write (stdout,'(a,a)') 'User specified time stamp = ',user_stamp
        base_time = user_time
        irefs = irefs + 1
        write (stdout,'(/a,a,a/)') ' "refuser = .true." was specified:'
     &,  ' =>logical switch calculations will be referenced to '
     &, 'the user specified time.'
      endif
c
      if (irefs .ne. 1) then
        write (stdout, *) '=>Error: Only one of the following can be'
     &  // ' set to "true": refrun, refinit, or refuser.'
        call abort()
      endif
c
c----------------------------------------------------------------------
c     convert time step to a time structure
c----------------------------------------------------------------------
c
      secs       = int(timestep)
      days       = 0
      dt_time    = set_time (secs, days)
      start_time = model_time + dt_time
#ifdef use_netCDF
c
c----------------------------------------------------------------------
c     set units of time for netcdf output 
c----------------------------------------------------------------------
c
c     common_year is 365.0 days. If working on planets with
c     different length years, there does not seem to be
c     any provision for this in the NETCDF udunits.dat file.
c     at least as of 7/11/96 
c
      if (rununits.eq.'days') then
        write(timeunits,'(a,i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a)') 
     &        'days since ',icyear,'-', icmonth,'-',icday,' '
     &       ,ichour,':',icmin,':',icsec,'.0'
      
      else
        if (eqyear) then
          if (eqmon) then
            type_of_calendar = 'common_year'
	  else
            type_of_calendar = 'Gregorian_year'
	  endif
        else
          type_of_calendar = 'Julian_year'
        endif 
c
        write(timeunits,'(a,i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a)') 
     &        'Gregorian_year since ',icyear,'-', icmonth,'-',icday,' '
     &       ,ichour,':',icmin,':',icsec,'.0'
      endif
#endif
c
      return
      end


      subroutine set_calendar (equal_years, equal_months)
c
c-----------------------------------------------------------------------
c     choose a calendar.
c-----------------------------------------------------------------------
c
      use time_manager_mod
#include "stdunits.h"
#include "tmngr.h"
      logical equal_years, equal_months
c
      if (equal_years .and. equal_months) then
        call set_calendar_type (thirty_day_months)
      endif
      if (equal_years .and. .not. equal_months) then
        call set_calendar_type (no_leap)
      endif
      if (.not. equal_years .and. .not. equal_months) then
        call set_calendar_type (julian)
      endif
      return
      end 


      subroutine set_stamp (stamp, year, month, day, hour, min, sec)
c
c=======================================================================
c     make a 32 character time stamp from day,month,year,sec,min,hour
c     author:       R.C. Pacanowski     e-mail=>  rcp@gfdl.gov
c=======================================================================
c
      integer year, month, day, hour, min, sec
      character stamp*32
c
      if (year .le. 9999) then
        write (stamp,'(a6,i2,a1,i2,a1,i4,a8,i2,a1,i2,a1,i2)') 
     &                'm/d/y=',month,'/',day,'/'
     &,                year,', h:m:s=', hour,':', min,':', sec
      elseif (year .ge. 10000 .and. year .le. 99999) then
        write (stamp,'(a6,i2,a1,i2,a1,i5,a7,i2,a1,i2,a1,i2)')
     &            'm/d/y=',month,'/',day,'/'
     &,           year,',h:m:s=', hour,':', min,':', sec
      elseif (year .ge. 100000 .and. year .le. 999999) then
        write (stamp,'(a6,i2,a1,i2,a1,i6,a6,i2,a1,i2,a1,i2)')
     &            'm/d/y=',month,'/',day,'/'
     &,           year,'h:m:s=', hour,':',  min,':', sec
      else
        write (*,*) '=>Error: year=',year,' is too large in set_stamp.'
        call abort()
      endif
      return
      end

      subroutine get_stamp (stamp, year, month, day, hour, min, sec)
c
c=======================================================================
c     convert 32 character time stamp into day,month,year,sec,min,hour
c     author:       R.C. Pacanowski     e-mail=>  rcp@gfdl.gov
c=======================================================================
c
      character stamp*32
      integer   year, month, day, hour, min, sec
      character skip6*6, skip7*7, skip8*8, skip1*1
c
c     for years: 0..9999
      if (stamp(17:17) .eq. ',') then 
        read  (stamp, '(a6,i2,a1,i2,a1,i4,a8,i2,a1,i2,a1,i2)')
     &      skip6, month, skip1, day, skip1, year, skip8, hour
     &,     skip1, min, skip1, sec
c     for years:  10000..99999
      elseif (stamp(18:18) .eq. ',') then
        read  (stamp, '(a6,i2,a1,i2,a1,i5,a7,i2,a1,i2,a1,i2)')
     &      skip6, month, skip1, day, skip1, year, skip7, hour
     &,     skip1, min, skip1, sec
c     for years:  100000..999999
      elseif (stamp(19:19) .eq. 'h') then
        read  (stamp, '(a6,i2,a1,i2,a1,i6,a6,i2,a1,i2,a1,i2)')
     &      skip6, month, skip1, day, skip1, year, skip6, hour
     &,     skip1, min, skip1, sec
      else
        write (*,*) '=>Error: year=',year,' is too large in get_stamp.'
        call abort()
      endif
      return
      end


      function file_stamp (fname, stamp, suffix) 
c
c=======================================================================
c     add day, month and year information from "stamp" to "fname"
c     author:       R.C. Pacanowski     e-mail=>  rcp@gfdl.gov
c=======================================================================
c
      character*(*) fname, stamp, suffix, file_stamp
      integer   year, month, day, hour, min, sec
      character skip6*6, skip7*7, skip8*8, skip1*1
c      
      character (len=2), dimension (1:31) :: dd
      character (len=2),dimension (1:12) :: mm
      save dd, mm
      character (len=6) yy
      data (dd(i),i=1,9) /'01','02','03','04','05','06','07','08','09'/
      data (dd(i),i=10,18)/'10','11','12','13','14','15','16','17','18'/
      data (dd(i),i=19,27)/'19','20','21','22','23','24','25','26','27'/
      data (dd(i),i=28,31) /'28','29','30','31'/
      data (mm(i),i=1,9) /'01','02','03','04','05','06','07','08','09'/
      data (mm(i),i=10,12)/'10','11','12'/
c
c     for years: 0..9999
      if (stamp(17:17) .eq. ',') then 
        read  (stamp, '(a6,i2,a1,i2,a1,i4,a8,i2,a1,i2,a1,i2)')
     &      skip6, month, skip1, day, skip1, year, skip8, hour
     &,     skip1, min, skip1, sec
        yy = '000000'
        if (year .le. 9) then
          write (yy(6:6),'(i1)') year
        elseif (year .le. 99) then
          write (yy(5:6),'(i2)') year
        elseif (year .le. 999) then
          write (yy(4:6),'(i3)') year
        else
          write (yy(3:6),'(i4)') year
        endif
c     for years:  10000..99999
      elseif (stamp(18:18) .eq. ',') then
        read  (stamp, '(a6,i2,a1,i2,a1,i5,a7,i2,a1,i2,a1,i2)')
     &      skip6, month, skip1, day, skip1, year, skip7, hour
     &,     skip1, min, skip1, sec
        yy(1:1)='0'
        write (yy(2:),'(i5)') year
c     for years:  100000..999999
      elseif (stamp(19:19) .eq. 'h') then
        read  (stamp, '(a6,i2,a1,i2,a1,i6,a6,i2,a1,i2,a1,i2)')
     &      skip6, month, skip1, day, skip1, year, skip6, hour
     &,     skip1, min, skip1, sec
        write (yy,'(i6)') year
      else
        write (*,*) '=>Error: bad stamp in file_stamp.'
	write (*,*) '         stamp=',stamp
        call abort()
      endif
      file_stamp = trim(fname)//'.'//yy//'.'//mm(month)//'.'//dd(day)
     &//trim(suffix)
      return
      end
