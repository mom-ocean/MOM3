      subroutine gosbc (ncall)
c
c=======================================================================
c     interpolate the ocean S.B.C. (surface boundary conditions 
c     which were prepared by the atmosphere) to the ocean grid
c
c     inputs:
c
c     ncall = number of times this routine was called
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#ifdef coupled
      use arrays_2d
      use sbc_info
      use sbc_ocean
      use sbc_atmos
# include "param.h"
# include "coord.h"
# include "grdvar.h"
c
      parameter (lenw=10*imt)
      common /gosbci/ nc1, nc2
      common /gosbcr/ work1(lenw), work2(lenw)
      dimension sor(imap2,jma), resa(imap2,jma)
      dimension average(maxsbc)
c
      call tic ('mom', 'get S.B.C. for ocean (gosbc)')
c
      write (stdout,8900)
c
c-----------------------------------------------------------------------
c     initialize counters to store interpolation weights
c-----------------------------------------------------------------------
c
      if (ncall .eq. 1) then
        nc1 = 0
	nc2 = 0
c
c-----------------------------------------------------------------------
c       compute initial checksums of ocean S.B.C.'s
c       "n" is the S.B.C. and "m" refers to its ordering within arrays
c-----------------------------------------------------------------------
c
        do m=1,numosbc
	  cksum = checksum (sbcatm(1,1,m), imap2, jma)
          write (stdout,*) sbcname(m),' S.B.C. checksum =',cksum
        enddo
      endif
c
c-----------------------------------------------------------------------
c     prepare each ocean S.B.C. one at a time
c     "n" is the S.B.C. and "m" refers to its ordering within arrays
c     This assumes that the S.B.C. have been initialized, accumulated,
c     and averaged in the Atmosphere model
c-----------------------------------------------------------------------
c
      do m=1,numosbc
# ifdef cyclic
c
c       set cyclic conditions on the atmosphere grid 
c
        do j=1,jma
          sbcatm(1,j,m)     = sbcatm(imap2-1,j,m)
          sbcatm(imap2,j,m) = sbcatm(2,j,m)
        enddo
# endif
c
c       compute the global mean from the heat flux at bottom of
c       atmosphere before interpolating to ocean
c
        if (sbcname(m) .eq. ' heatflux') then
          avg  = c0
          anum = c0
          do j=2,jma-1
            cosdy = abcgcs(j)*abcgdy(j)
            do i=2,imap2-1
              weight = aland(i,j)*abcgdx(i)*cosdy
              anum = anum + weight
              avg  = avg + weight*sbcatm(i,j,m)
	    enddo
	  enddo
          if (anum .ne. c0) avg = avg/anum
	  fluxcor = avg*4.186e4
	  write (stdout,8400) fluxcor
	endif
c
# ifdef trace_coupled_fluxes
        write (stdout,*) ' ===> values on entering gosbc.F:'
        call scope (sbcatm(1,1,m), imap2, imap2, jma, sbcname(m))
# endif
c
c-----------------------------------------------------------------------
c       extrapolate values into land areas on the atmospheric grid
c       to accommodate mismatches in ocean and atmospheric land masks
c       when interpolating to ocean grid
c-----------------------------------------------------------------------
c
        call extrap (sbcatm(1,1,m), aland,  sor, resa,imap2, jma, numpas
     &,              crits(m), sbcname(m), 2, 0)
# ifdef trace_coupled_fluxes
        write (stdout,*) ' ===> after extrapolating into land:'
        call scope (sbcatm(1,1,m), imap2, imap2, jma, sbcname(m))
# endif
c
c-----------------------------------------------------------------------
c       linearly interpolate to the ocean grid assuming the
c       atmos grid is coarse relative to the ocean grid
c-----------------------------------------------------------------------
c
        if (sbcname(m) .eq. ' taux     ' .or. 
     &      sbcname(m) .eq. ' tauy     ') then
          nc1 = nc1 + 1
          call ctf (sbcatm(1,1,m), imap2, jma, abcgx, abcgy 
     &,       sbcocn(1,1,m), imt, jmt, 2, imtm1, 2, jmt-1, xu, yu, nc1
     &,       work1, lenw)
        else
	  nc2 = nc2 + 1
          call ctf (sbcatm(1,1,m), imap2, jma, abcgx, abcgy 
     &,       sbcocn(1,1,m), imt, jmt, 2, imtm1, 2, jmt-1, xt, yt, nc2
     &,       work2, lenw)
	endif
c
        do j=1,jmt
# ifdef cyclic
          sbcocn(1,j,m)   = sbcocn(imtm1,j,m)
          sbcocn(imt,j,m) = sbcocn(2,j,m)
# else
          sbcocn(1,j,m)   = sbcocn(2,j,m)
          sbcocn(imt,j,m) = sbcocn(imtm1,j,m)
# endif
	enddo
c
        do i=1,imt
          sbcocn(i,jmt,m) = sbcocn(i,jmtm1,m)
          sbcocn(i,1,m)   = sbcocn(i,2,m)
        enddo
c
# ifdef trace_coupled_fluxes
        write (stdout,*) ' ===> after interpolating to "sbcocn":'
        call scope (sbcocn(1,1,m), imt, imt, jmt, sbcname(m))
# endif
c
c-----------------------------------------------------------------------
c       convert to units expected by ocean
c
c       note: a westerly wind means taux is positive (toward the east)
c             a southerly wind means tauy is positive (toward the north)
c             a positive heat flux means the ocean is warming
c             a positive fresh water flux means the ocean is getting
c               fresher (less dense)
c-----------------------------------------------------------------------
c
        do j=1,jmt
          do i=1,imt
            sbcocn(i,j,m) = coabc(m)*sbcocn(i,j,m)
	  enddo
	enddo
# ifdef trace_coupled_fluxes
        write (stdout,*) ' ===> after converting units:'
        call scope (sbcocn(1,1,m), imt, imt, jmt, sbcname(m))
# endif
c
c-----------------------------------------------------------------------
c      calculate averages of the ocean S.B.C.
c      (the S.B.C. are assumed to be defined on the same grid as "kmt")
c-----------------------------------------------------------------------
c
        average(m)  = c0
        anum        = c0
        do j=1,jmt
          cosdy = cst(j)*dyt(j)
          do i=2,imt-1
            weight = min(kmt(i,j),1)*dxt(i)*cosdy
            anum = anum + weight
            average(m)  = average(m) + weight*sbcocn(i,j,m)
          enddo
        enddo
        if (anum .ne. c0) average(m) = average(m)/anum
c
      enddo
c
c-----------------------------------------------------------------------
c      show averages of the ocean boundary conditions
c-----------------------------------------------------------------------
c
      write (stdout,9100)
      do n=1,numosbc
        write (stdout,9200) n, sbcname(n), average(n), dunits(n)       
        if (sbcname(n) .eq. ' heatflux') then
c
c         also show in units of watts/m**2
c
	  write (stdout,8500) average(n)*4.186e4, ' watts/m**2'
        endif
        if (sbcname(n) .eq. ' saltflux') then
c
c         also so equivalent fresh water flux
c
	  write (stdout,8500) average(n)*86400.0/0.035, ' cm/day    '
        endif
      enddo
c
      write (stdout,9400)
      call toc ('mom', 'get S.B.C. for ocean (gosbc)')
8400  format (1x,'Note: not removing a net imbalance of ',1pe14.7
     &,       'watts/m**2',/,'        from the atmospheric heat flux.') 
8500  format (57x,1pe14.7,18x,a15)
8900  format (/,10x, ' ==> Getting ocean S.B.C.')
9100  format (/10x,' ==> S.B.C. averages for the ocean follow:'/)
9200  format (17x,'for S.B.C. #',i2,', the average ',a10,' is ',1pe14.7
     &,       1x, ' after converting from ',a15)
9400  format (/10x,' ==> S.B.C. prepared for this ocean segment.'/)
#endif
      return
      end
