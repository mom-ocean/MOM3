!====================== include file "param.h" =========================
!
#include "size.h"
!
!     nvar   = number of prognostic variables
!     lseg   = maximum number of longitudinal stream function segments
!     nlatpr = maximum number of latitudes for matrix printouts
!              on diagnostic time steps
!     nhreg  = number of regions in the horizontal used for averaging
!              tracers.
!     nvreg  = number of regions in the vertical used for term balance
!              calculations. note "nvreg" isn`t used for tracer averages
!     numreg = total number of regions ( = product of nhreg & nvreg)
!              used for term balance calculations
!
      parameter (lseg=5, nlatpr=10)
      parameter (nhreg=5, nvreg=2, numreg=nhreg*nvreg)
      parameter (nvar=nt+2)
!
      parameter (imtm1=imt-1, imtm2=imt-2)
      parameter (kmm1=km-1, kmp1=km+1)                                
      parameter (jmtm1=jmt-1, jmtm2=jmt-2)
      parameter (imtkm=imt*km, nwds=imt*jmt, nkflds=2)
      parameter (nslab=imt*nvar*km, ntmin2=nt+1/nt)
!
!     parameter constants
!
#include "pconst.h"
!
!     standard I/O units
!
#include "stdunits.h"
