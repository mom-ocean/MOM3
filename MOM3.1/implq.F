      subroutine implq (q2, q1, tdt, q2bct, j)
#if defined tcvmix && defined implicitvmix
c
c======================================================================
c      solve the vertical diffusion equation of tke implicitly
c      using the method of inverting a tridiagonal matrix
c      as described in richtmyer and morton.
c      dissipation is included implicitly
c
c      inputs:
c
c      q1      = right hand side terms
c      tdt     = 2 * timestep
c      q2bct   = top boundary condition
c      j       = row j
c
c      outputs:
c
c      q2      = returned solution
c
      use arrays_2d
#include "param.h"
#include "ctcmix.h"
#include "vmixc.h"
      common /invt/ a, b, c, d, e, f, g
      dimension a(imt,km), b(imt,km), c(imt,km), d(imt,km)
      dimension e(imt,0:kmp1), f(imt,0:kmp1), g(imt)
      dimension q1(imt,km), q2(imt,km)
      dimension q2bct(imt)
c
c======================================================================
c
      do 100 k=1,kmm1
        do 90 i=1,imt
          a(i,k) = eeq(i,k)*tdt
          c(i,k) = ffq(i,k)*tdt
          e(i,k) = tdt*dissp(i,k)
90      continue
100    continue
c
      do 200 k=1,kmm1
        do 190 i=1,imt
          b(i,k) = c1 + a(i,k) + c(i,k) + e(i,k)
          d(i,k) = q1(i,k)
190     continue
200   continue
c
      do 300 i = 1,imt
        e(i,kmm1) = c0
        f(i,kmm1) = c0
300   continue
c
c     b.c. at bottom (q2=q2l=0)
c
      do 400 i=1,imt
        kz = kmt(i,j)
        if (kz .eq. 0) goto 400
        a(i,kz) = c0
        d(i,kz) = c0
400   continue
c
      do 500 kk=2,km
        k = km + 1 - kk
          do 490 i=1,imt
            g(i) = c1/(b(i,k)-c(i,k)*e(i,k))
            f(i,k-1) = (d(i,k)+c(i,k)*f(i,k))*g(i)
            e(i,k-1) = a(i,k)*g(i)
490       continue
500   continue
c
      do 600 i=1,imt
        q2(i,1) = e(i,0)*q2bct(i) + f(i,0)
600   continue
c
      do 700 k=2,km
        do 690 i=1,imt
          q2(i,k) = e(i,k-1)*q2(i,k-1) + f(i,k-1)
690     continue
700   continue
c
#endif
      return
      end
