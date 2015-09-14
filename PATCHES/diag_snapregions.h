!==================== include file "diag_snapregions.h" =====================
!     nregions must be defined in diag_region_out.h
      character regionname(nregions) *16
      integer, dimension(nregions):: isreg, iereg, jsreg, jereg
      integer, dimension(nregions):: ksreg, kereg
      integer, dimension(nregions):: regionindex
      data (regionindex(n),isreg(n), iereg(n), jsreg(n), jereg(n)
     &       ,ksreg(n), kereg(n)       
     &       ,regionname(n), n=1,nregions)	       
     &  / 1,  1, 10, 20, 50, 10, 20,'obc'  
     &  , 2, 40, 70, 30, 60,  1,  6,'gotland'  
     &  /
