c====================== include file "dens.h" ==========================
c
c-----------------------------------------------------------------------
c     statement functions for 
c     1. Density
c     2. First derivatives of density.
c     3. Second derivatives of density.
c-----------------------------------------------------------------------
c
#ifdef linearized_density
c
      dens(tq,sq,k)       = c(k,1)*tq
      drhodt (tq, sq, k)  = c(k,1) 
      drhods (tq, sq, k)  = 0.0
      ddensdtdt (tq,sq,k) = 0.0
      ddensdtds (tq,sq,k) = 0.0
      ddensdsds (tq,sq,k) = 0.0
c
#else
      dens (tq, sq, k) = (c(k,1) + (c(k,4) + c(k,7)*sq)*sq +
     &                   (c(k,3) + c(k,8)*sq + c(k,6)*tq)*tq)*tq +
     &                   (c(k,2) + (c(k,5) + c(k,9)*sq)*sq)*sq
c
       drhodt (tq, sq, k) = c(k,1) + (c(k,4) + c(k,7)*sq)*sq 
     &                + (2.0*c(k,3) + 2.0*c(k,8)*sq + 3.0*c(k,6)*tq)*tq

       drhods (tq, sq, k) = (c(k,4) + 2.0*c(k,7)*sq + c(k,8)*tq)*tq
     &                   + c(k,2) + (2.0*c(k,5) + 3.0*c(k,9)*sq)*sq
c
       ddensdtdt (tq,sq,k) = 2.0*c(k,3) + 6.0*c(k,6)*tq + 2.0*c(k,8)*sq 
       ddensdtds (tq,sq,k) = c(k,4) + 2.0*c(k,7)*sq + 2.0*c(k,8)*tq
       ddensdsds (tq,sq,k) = 2.0*c(k,5) + 2.0*c(k,7)*tq + 6.0*c(k,9)*sq
#endif
