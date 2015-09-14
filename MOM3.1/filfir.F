      subroutine filfir (t, im, kind, num)
#ifdef firfil
c
c=======================================================================
c     simple finite impulse response filter with [.25, .5, .25] weights
c     modified for assymmetric and symmetric boundary conditions
c
c     input:
c             t    = array of quantity to be filtered
c             im   = length of "t"
c             kind = (1,2,3) = (symmetric, asymmetric, cyclic) b.c.
c             num = number of filter passes
c
c     output:
c             t    = array of filtered quantity
c
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#include "tmngr.h"
c
      dimension t(im), s(imt)
c
      call tic ('filtering', 'filfir (finite impulse)')

      if (kind .eq. 1) then
c
c-----------------------------------------------------------------------
c       apply the filter "num" times using a symmetric (no flux)
c       boundary condition
c-----------------------------------------------------------------------
c
	do n=1,num
          s(1)  = 0.75*t(1)  + p25*t(2)
          s(im) = 0.75*t(im) + p25*t(im-1)
          do i=2,im-1
            s(i) = p25*t(i-1) + p5*t(i) + p25*t(i+1)
          enddo
          do i=1,im
            t(i) = s(i)
          enddo
        enddo
      elseif (kind .eq. 2) then
c
c----------------------------------------------------------------------
c       apply the filter "num" times using an asymmetric (flux)
c       boundary condition
c----------------------------------------------------------------------
c
c       calculate mean
c
        sum = 0.0
        do i=1,im
          sum = sum + t(i)
        enddo
        avg1 = sum/im
c
	do n=1,num
          s(1)  = 0.5*t(1)  + p25*t(2)
          s(im) = 0.5*t(im) + p25*t(im-1)
          do i=2,im-1
            s(i) = p25*t(i-1) + p5*t(i) + p25*t(i+1)
          enddo
          do i=1,im
            t(i) = s(i)
          enddo
        enddo
c
        sum = 0.0
        do i=1,im
          sum = sum + t(i)
        enddo
        avg2 = sum/im
c
c       replace mean
c
        do i = 1,im
          t(i) = t(i) - avg2 + avg1
        enddo
c
      elseif (kind .eq. 3) then
c
c----------------------------------------------------------------------
c       apply the filter "num" times using a cyclic
c       boundary condition
c----------------------------------------------------------------------
c
	do n=1,num
          s(1)  = p25*t(im) + 0.5*t(1)  + p25*t(2)
          s(im) = p25*t(im-1) + 0.5*t(im) + p25*t(1)
          do i=2,im-1
            s(i) = p25*t(i-1) + p5*t(i) + p25*t(i+1)
          enddo
          do i=1,im
            t(i) = s(i)
          enddo
        enddo
      else
        write (stdout,'(/a,i10,a)') ' error=> kind =', kind,' in filfir'
        call abort()
      endif
c
      call toc ('filtering', 'filfir (finite impulse)')
#endif
      return
      end
