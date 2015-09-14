      block data blkdta
      use mpp_io_mod
#include "param.h"
#ifdef neptune
# include "cnep.h"
#endif
#include "cprnts.h"
#include "diag.h"
#include "iounit.h"
#include "switch.h"
#include "tmngr.h"
c
c     data initializations
c
c     initialize diagnostic intervals and averaging periods
c     intervals > 0.0 are in days and set the time between output
c     intervals < 0.0 give no output
c     intervals = 0.0 imply output every time step
c     averaging periods are in days and indicate the time period 
c     for averaging diagnostics. If the averaging period = 0.0,
c     the period is equal to the time of one time step. 
c
      data glenint /-365.25/
      data prxzint /-365.25/
      data extint  /-365.25/
      data tyzint  /-365.25/
      data fctint  /-365.25/
      data vmsfint /-365.25/
      data gyreint /-365.25/
      data igyre   /.true./
      data dspint  /-365.25/
      data dspper  /20.0/
      data stabint /-365.25/
      data zmbcint /-365.25/
      data trajint /-365.25/
c
      data tsiint      /1.0/
      data cmixint     /-365.25/
      data crossint    /-365.25/
      data pressint    /-365.25/
      data exconvint   /-365.25/
      data tavgint     /-365.25/
      data itavg       /.true./
      data tmbint      /-365.25/
      data tmbper      /20.0/
      data itmb        /.true./
      data trmbint     /-365.25/
      data itrmb       /.true./
      data snapint     /-365.25/
      data timavgint   /-365.25/
      data timavgper   /20.0/
      data xbtint      /-365.25/
      data xbtper      /20.0/
      data rhoint      /-365.25/
      data diffsnapint /-365.25/
c
      data cmixts     /.false./
      data crossts    /.false./
      data pressts    /.false./
      data fctts      /.false./
      data exconvts   /.false./
      data snapts     /.false./
      data timavgts   /.false./
      data prxzts     /.false./
      data diffsnapts /.false./  
c
c     default starting & ending coordinates for diagnostic matrix prints
c     (lats and lons in degrees, depths in cm)
c     see "cprnts.h"
c
      data prlat /-60.0, 0.0, 27.0, 55.0, 6*100.0/
      data prslon /180.0, 180.0, 180.0, 180.0, 6*0.0/
      data prelon /230.0, 230.0, 230.0, 230.0, 6*0.0/
      data prsdpt /nlatpr*0.0/
      data predpt /nlatpr*6000.e2/
      data slatxy, elatxy, slonxy, elonxy /-90.0, 90.0, 1.5, 358.5/
c
c     default some control #`s for how diagnostics are to be written
c      (see "iounit.h" for explanation) 
c
      data restrt /.false./
      data iotavg  /-1/
      data iotmb   /-1/
      data ioglen  /-1/
      data iotrmb  /-1/
      data iogyre  /-1/
      data iovmsf  /-1/
      data ioprxz  /-1/
      data ioext   /-1/
      data iodsp   /-1/
      data iotsi   /-1/
      data iozmbc  /-1/
      data ioxbt   /-1/
c
#ifdef neptune
c
c     set length scales for topographic stress
c
      data spnep /3.0e5/
      data senep /12.0e5/
#endif
c
c
      end
