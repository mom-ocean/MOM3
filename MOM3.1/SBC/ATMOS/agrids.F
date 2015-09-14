      subroutine agrids
c
c-----------------------------------------------------------------------
c     the atmospheric b.c. grid is assumed to be global from the south
c     pole (j=1) to the north pole (j=jma) and from Greenwich (i=2)
c     eastward to 1 grid point west of Greenwich (i=ima+1). i=1 and
c     i = ima+2 are cyclic boundary points
c     
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      use sbc_info
      use sbc_atmos
#include "param.h"
c
      write (stdout,9000)
c
c-----------------------------------------------------------------------
c      based on this grid layout, define the b.c. grids "abcgx"
c      and "abcgy".
c-----------------------------------------------------------------------
c
c      the atmos model has "ima" longitudes and "jma" latitudes. the
c      S.B.C. grid has 2 extra points in longitude as defined below. 
c      
       dx     = 360.0/ima
c       abcgx0 = 0.0
       abcgx0 = 0.5*dx
       do i=1,ima
         abcgx(i+1) = abcgx0 + (i-1)*dx
         abcgdx(i+1)  = dx
       enddo
c
c      add 2 boundary points necessary to contain the MOM S.B.C. grid 
c      within the atmosphere grid.
c
       abcgx(1)     = abcgx(2)    - dx
       abcgx(imap2) = abcgx(imap2-1) + dx
       abcgdx(1)      = dx
       abcgdx(imap2)  = dx
       write (stdout,9101) ' longitudes for atmos b.c. grid:'
       write (stdout,9001) (abcgx(i),i=1,imap2)
c
c      in the latitude direction, it may be necessary to define "jma"
c      so that the atmosphere grid goes from j=2..jma-1 leaving j=1
c      and j=jma as extra b.c. rows to insure the ocean b.c. grids are
c      contained within the atmosphere b.c. grids.
c
       dy  = 180.0/(jma-1)
       abcgy0 = -90.0
       degrad = 4.0*atan(1.0)/180.0
       do j=1,jma
         abcgy(j)   = abcgy0 + (j-1)*dy
         abcgdy(j)  = dy
	 abcgcs(j)  = cos(abcgy(j)*degrad)
       enddo
       write (stdout,9101) ' latitudes for atmos b.c. grid:'
       write (stdout,9001) (abcgy(j),j=1,jma)
c
       return
9000  format (/,10x,'==> Constructing atmos b.c. grids'/)
9001  format (1x,10g13.5)
9101  format (/,1x,a34,i2)
       end
