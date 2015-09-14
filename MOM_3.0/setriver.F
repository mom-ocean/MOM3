      subroutine setriver (js, je)
#ifdef river_inflow
      use arrays_2d 
# include "param.h"
# include "grdvar.h"
# include "river.h"
c----------------------------------------------------------------------
c     calculate the area of the boxes with river inflow
c----------------------------------------------------------------------
      write(stdout,'(a,i2)') "Number of defined rivers : ",nriv
      do n=1,nriv 
        areariv(n)=0.0
        do m=1,nboxriv(n)
           ir = iriv(m,n)
           jr = jriv(m,n)
           if (kmt(ir,jr).lt.2) then
             write(stdout,'(a,i2,a,i2,a)') 
     &         'error => River box',m,' of river',n,' is on land!' 
             write(stdout,'(a,i4,a,i4,a)') 
     &         'ir = ',ir,'jr = ',jr
             call abort()
           endif
           areariv(n) = areariv(n)+dxt(ir)*cst(jr)*dyt(jr)
        enddo
        write(stdout,'(a,i2,a,i2,a)')
     &  'river ',n,' is simulated in ',nboxriv(n),' boxes at'
        do m=1,nboxriv(n)
          write(stdout,'(i4,i4)') iriv(m,n),jriv(m,n)
        enddo
        write(stdout,'(a,1pe12.5,a)') "Surface area of the river boxes"
     &  , areariv(n),' cm**2'
      enddo
#endif      
      return
      end
