!====================== include file "tmngr.h" =========================
!
!                       time manager variables
!
!     time manager inputs:
!
!     refrun  = T ==> logical switches referenced to the start of run
!     refinit = T ==> logical switches referenced to initial condition
!                     time: year0, month0, day0, hour0, min0, sec0
!     refuser = T ==> logical switches referenced to user specified 
!                     time: ryear, rmonth, rday, rhour, rmin, rsec
!
!     dayoyr  = relative day number referenced to the beginning
!               of the current year.  (real)
!     relyr   = number of years (and fractional years) of model
!               integration (for time tau+1 {itt}) relative to 
!               initial condition
!     prelyr  = relyr for previous time step
!
!     stamp   = 32 character date and time for current model timestep
!     pstamp  = 32 character date and time for previous model timestep
!     ic_stamp= 32 character date and time for initial conditions
!     stop_stamp = 32 character date and time for last time step
!
!     itt     = current time step counter (from initial cond.)
!     itt0    = time step at start of current run
!
!               variables used for initialization
!
!     irstdy  = integer number of days at start of run
!     msrsdy  = fractional day in millisec at start of run
!
!     year0   = year of initial conditions
!     month0  = month of initial conditions
!     day0    = day of initial conditions
!     hour0   = hour of initial conditions
!     min0    = minute of initial conditions
!     sec0    = second of initial conditions  
!     
!     ryear   = year of user specified reference time
!     rmonth  = month of user specified reference time
!     rday    = day of user specified reference time
!     rhour   = hour of user specified reference time
!     rmin    = minute of user specified reference time
!     rsec    = second of user specified reference time  
!
!-----------------------------------------------------------------------
!
!
      logical refrun, refinit, refuser, eqmon, eqyear
!
      integer year0, month0, day0, hour0, min0, sec0
      integer ryear, rmonth, rday, rhour, rmin, rsec
      character *32 pstamp, stamp, ic_stamp, stop_stamp, xstamp
!
      common /tmngrl/ refrun, refinit, refuser, eqmon, eqyear
!
      common /tmngrii/ itt0, itt, irstdy, msrsdy
!
      common /tmngrr/ dayoyr, relyr, prelyr, relday
!
      common /tmngrc/ stamp, pstamp, ic_stamp, stop_stamp, xstamp
!
      common /tmngrii/ year0, month0, day0, hour0, min0, sec0
      common /tmngrii/ ryear, rmonth, rday, rhour, rmin, rsec
!
      character*10 dayname
      character*12 monname
      common /tmngrc/ dayname(7), monname(12)
!
      real yrlen
      common /tmngrr/ yrlen
!
      character(len=80) :: timeunits
      common /tmngrc/ timeunits
      
