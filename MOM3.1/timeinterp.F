#ifdef test_timeinterp
      program driver
c
c=======================================================================
c     test the time interpolator module
c
c     author :   R. C. Pacanowski        e-mail=> rcp@gfdl.gov
!     revised:   m. schmidt e-mail=>martin.schmidt@io-warnemuende.de
c=======================================================================
c
      use time_manager_mod
      use timeinterp_mod
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
!     Initialize a new data set. timeinterp should return ntdrec=1
      ntdrec = -1
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
#else
      subroutine tdummy
      return
      end
#endif


