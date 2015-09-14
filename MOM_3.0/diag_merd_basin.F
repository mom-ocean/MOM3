      subroutine merid1st 
#ifdef meridional_overturning_basin
c
c-----------------------------------------------------------------------
c     read the basin mask to divide the ocean into various basin
c     numbers for the meridional transport calculation
c
c     author:      k. dixon      e-mail=> kd@gfdl.gov
c---------------------------------------------------------------------
c
      use arrays_2d
      use mpp_io_mod
#include "param.h"
#include "diag.h"
      integer normsk(imt,jmt)
c
c     for test case, number of basins in diag.h is set 
c     to 1 and no mask is read in. In this case, setup the mask for 
c     a single basin which is the entire ocean:
c
      if (nbasin .eq. 1) then
        do j=1,jmt
          do i=1,imt
            if (kmt(i,j) .eq. 0) then
              normsk(i,j) = 0
            else
              normsk(i,j) = 1
            end if
          enddo
        enddo
      else
        call getunit (ionrth, 'northmasks'
     &       ,             'formatted sequential rewind')
        call nor1st (ionrth, normsk)
        call relunit (ionrth)
      end if
c
c     check if all ocean points belong to a basin
c
      icheck = 0
      do j=1,jmt
        do i=1,imt
          if (normsk(i,j) .gt. nbasin) then
            icheck = icheck + 1
          endif
        enddo
      enddo
      if (icheck .gt. 0) then
        write(stdout,900) icheck
900     format(' Warning:  not all ocean points belong to a basin.'
     &,        ' Number of missing points=',i7) 
      endif
c
c     make a mask for each basin (on t points, following kd) 
c
      do j=1,jmt
        do i=1,imt 
c
c         global overturning.
c
          basinmask(i,j,0) = 1.
c
c         by basin 
c
          do l=1,nbasin
            if (normsk(i,j) .eq. l) then            
              basinmask(i,j,l) = 1. 
            else
              basinmask(i,j,l) = 0.
            endif
          enddo
c
        enddo
      enddo
c
#endif
      return
      end
