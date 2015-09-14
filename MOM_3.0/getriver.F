      subroutine getriver
#ifdef river_inflow
c
c=======================================================================
c      This subroutine gets the mass- and tracer inflow due to rivers.
c      The position of the river is defined in river.h.  
c      rff(nriv)    = river water discharge as volume flux in cm**3/sec
c      rtf(nriv,nt) = river tracer input in data units
c=======================================================================
c
# include "param.h"
# include "river.h"
c
# ifdef simple_sbc
      rff(1) = 3200.e6 / areariv(1)
      rtf(1,1) = 10.e0
      rtf(1,2) = 0.1
# else
c     Here the interpolation of the river data from a database has to be
c     added. Such a tool is available from
c     martin.schmidt@io-warnemuende.de      
# endif
#endif
      return
      end
