#ifdef time_mean_sbc_data
      subroutine setatm (dtatm)
c
c-----------------------------------------------------------------------
c     set up all quantities that need to be defined only once for this
c     simple atmosphere.
c     
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
# include "stdunits.h"
c
c     The atmospheric time step "dtatm" is not needed since this is
c     not an atmospheric model. it will be defined in mom.F 
c
      write (stdout,'(/,a,/)')
     &'==> S.B.C. for MOM are comming from annual mean S.B.C. data'
c
# ifdef minimize_sbc_memory
      write (stdout,'(/,a,a,/)')
     &'==> Error: "minimize_sbc_memory" has not been extended to work '
     &,'with "time_mean_sbc_data"'
      call abort()
# endif
# ifdef simple_sbc
      write (stdout,'(/,a,a,/)')
     &'==> Error: "time_mean_sbc_data" and '
     &,' "simple_sbc" are incompatable'
      call abort()
# endif
# ifdef coupled
c
c     verify that "coupled" is not enabled.
c
      write (stdout,*) '==> Error: "time_mean_sbc_data" and '
     &, '"coupled" are incompatable'
      call abort()
# endif
       return
       end

#else
      subroutine nodatasetatm
      return
      end
#endif
