!==================== include file "diag_snapregions.h" =====================
      character regionname(nregions) *16
      integer, dimension(nregions):: isreg, iereg, jsreg, jereg
      integer, dimension(nregions):: regionindex
      data (regionindex(n),isreg(n), iereg(n), jsreg(n), jereg(n)
     &       ,regionname(n), n=1,nregions)	       
     &  / 1,  2,  4, 20, 60, 'obc'    
     &  , 2, 40, 70, 30, 60, 'gotland'  
     &  /
