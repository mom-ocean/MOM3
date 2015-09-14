c====================== include file "iomngr.h"=========================
c
c
c     sdsbase  = base block number for ssread/sswrite file on SSD 
c     sdsaddr  = block number within ssread/sswrite file on SSD
c     sdsoff   = offset to handle multiple files on SSD
c
      integer sdsbase, sdsaddr, ssread, sswrite, sdsoff
      character*16 sds_name
      common /ssdrw/ sdsbase(100), sds_name(100), sdsaddr, sdsoff
c


