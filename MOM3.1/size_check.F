      subroutine size_check(imt2, jmt2, km2, sub, option)
c
c-----------------------------------------------------------------------
c     check that array bounds (imt2, jmt2, km2) = (imt, jmt, km)
c
c     inputs:
c       imt2 = input value for imt
c       jmt2 = input value for jmt
c       km2  = input value for km
c       sub  = name of the subroutine requesting the size check
c       option = what to do if size check fails:
c              'stop' is the only option
c
c     author:  c. h. goldberg             e-mail chg@gfdl.gov
c-----------------------------------------------------------------------
c
      character *(*) sub, option
      character *60 msg
#include "size.h"
#include "stdunits.h"
c
      if (imt .ne. imt2 .or. jmt .ne. jmt2 .or. km .ne. km2) then
        write (stdout,'(a/2(a,i4,a,i4,a,i4,a,a,/))')
     & '==>Error:  size_check '
     &,       'imt = ', imt2, '  jmt = ',jmt2,'  km = ',km2
     &,       ' in ', sub
     &,       'imt = ', imt, '  jmt = ',jmt,'  km = ',km
     &,       ' in "size.h"'
        write (stdout,'(/,a,a,a)')
     &        'Sizes in ',sub,' are incompatible with "size.h"'
	call abort()
      end if
      return
      end

