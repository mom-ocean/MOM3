      subroutine gasbc (ncall)
c
c=======================================================================
c     interpolate the atmospheric S.B.C. (surface boundary conditions 
c     which were prepared by the ocean) to the atmosphere grid
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
c
      logical errors
      parameter (lenw=10*imt)
      common /gasbcr/ work1(lenw)
      dimension sor(imt,jmt), reso(imt,jmt), average(maxsbc)
      save errors
      data errors /.false./
c
      call tic ('mom', 'get S.B.C. for atmos (gasbc)')
      write (stdout,8900)
c
c-----------------------------------------------------------------------
c     set up the necessary things on the first entry
c     "n" is the atmos S.B.C. and "m" refers to its ordering in arrays
c-----------------------------------------------------------------------
c
      if (ncall .eq. 1) then
c
c       find ocean grid domain in terms of atmosphere grid indices
c
        isocn = indp (xt(1), abcgx, imap2)
	if (abcgx(isocn) .lt. xt(1)) isocn = isocn + 1
	ieocn = indp (xt(imt), abcgx, imap2)
	if (abcgx(ieocn) .gt. xt(imt)) ieocn = ieocn - 1
	jsocn = indp (yt(1), abcgy, jma)
	if (abcgy(jsocn) .lt. yt(1)) jsocn = jsocn + 1
	jeocn = indp (yt(jmt), abcgy, jma)
	if (abcgy(jeocn) .gt. yt(jmt)) jeocn = jeocn - 1
c
        do m=numosbc+1,numsbc
c
c-----------------------------------------------------------------------
c         if ocean and atmospheric domains are not coincident, define
c         a blending zone to provide a smooth transition between SST
c         from within the ocean model and prescribed SST outside of 
c         the ocean`s domain. This is for the case of a global
c         atmosphere with a regional ocean. bzone is controlled by
c         setting the blending width "bwidth" and is defined on the
c         atmospheric grid as follows:
c
c         bzone=0 ==> region of width "bwidth" where SST will be
c                       blended
c         bzone=1 ==> domain of the ocean model
c         bzone=2 ==> domain outside the ocean where SST is prescribed
c
c         note: if "bwidth" = 0 then there is no "bzone"
c-----------------------------------------------------------------------
c
          if (sbcname(m) .eq. ' sst') then
            if (bwidth .ne. c0) then
              do j=1,jma
                do i=2,imap2-1
                  if (abcgx(i) .le. xt(2) - bwidth .or.
     &              abcgx(i) .ge. xt(imt-1) + bwidth .or.
     &              abcgy(j) .le. yt(2)     - bwidth .or.
     &              abcgy(j) .ge. yt(jmt-1) + bwidth) then
                    bzone(i,j) = 2
                  elseif (i .le. isocn .or. i .ge. ieocn .or.
     &              j .le. jsocn .or. j .ge. jeocn) then
                    bzone(i,j) = 0
		  else
                    bzone(i,j) = 1
                  endif
                enddo
                bzone(1,j)  = bzone(imap2-1,j)
                bzone(imap2,j) = bzone(2,j)
              enddo
c
              write (stdout,9100)
              call iplot (bzone, imap2, imap2, jma)
            endif
          endif
        enddo
c
c---------------------------------------------------------------------
c       compute initial checksums of the atmos S.B.C.
c---------------------------------------------------------------------
c
        do m=numosbc+1,numsbc
	  cksum = checksum (sbcocn(1,1,m), imt, jmt)
          write (stdout,*) sbcname(m),' S.B.C. checksum =',cksum
        enddo
      endif
c
c-----------------------------------------------------------------------
c     prepare each atmosphere S.B.C. one at a time.
c-----------------------------------------------------------------------
c 
      n = 0
      do m=numosbc+1,numsbc
        n = n + 1      
c
c-----------------------------------------------------------------------
c       apply prescribed "SST" outside of the ocean domain.
c       this is needed only when the ocean domain is not global. in the
c       case where the ocean domain is a cyclic strip between two
c       latitudes, "sstpre" can be set to -2 deg C at the poles and the
c       blending zone (bzone=0) can be defined poleward of the ocean
c       domain. in the case where the ocean domain is a basin, "sstpre"
c       should really be changed to be a function of latitude and
c       longitude.
c-----------------------------------------------------------------------
c
        if (bwidth .ne. c0 .and. sbcname(m) .eq. ' sst') then
          do j=1,jma
            do i=1,imap2
              if (bzone(i,j) .eq. 2) sbcatm(i,j,m) = sstpre
            enddo
          enddo
        endif
c
c-----------------------------------------------------------------------
c       set conditions on ocean grid. (cyclic or no flux)
c-----------------------------------------------------------------------
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
        do i=1,imt
          sbcocn(i,jmt,m) = sbcocn(i,jmtm1,m)
          sbcocn(i,1,m)   = sbcocn(i,2,m)
        enddo
# ifdef trace_coupled_fluxes
        write (stdout,*) ' ===> values on entering gasbc.F:'
        call scope (sbcocn(1,1,m), imt, imt, jmt, sbcname(m))
# endif
c
c-----------------------------------------------------------------------
c       extrapolate values into land areas on the ocean grid
c       to accommodate mismatches in ocean and atmospheric land masks
c       when interpolating to atmosphere grid
c-----------------------------------------------------------------------
c
        call extrap (sbcocn(1,1,m), kmt, sor, reso, imt, jmt, numpas
     &,              crits(m), sbcname(m), 1, 0)
c
# ifdef trace_coupled_fluxes
        write (stdout,*) ' ===> after extrapolation into land:'
        call scope (sbcocn(1,1,m), imt, imt, jmt, sbcname(m))
# endif
c
c-----------------------------------------------------------------------
c       interpolate (by averaging) to the atmosphere grid assuming the
c       atmos grid is coarse relative to the ocean grid. The "sbcocn"
c       will be zeroed in "tracer" and "baroclinic" before accumulating
c       at the start of the segment based on switch "osegs"
c-----------------------------------------------------------------------
c
        call ftc (sbcocn(1,1,m), imt, jmt, xt, yt, sbcatm(1,1,m)
     &,           imap2, jma, isocn, ieocn, jsocn, jeocn, abcgx, abcgy
     &,           ncall+n-1, work1, lenw)
c
c       set cyclic conditions on the atmosphere grid
c
# ifdef cyclic
	do j=1,jma
	  sbcatm(1,j,m)     = sbcatm(imap2-1,j,m)
	  sbcatm(imap2,j,m) = sbcatm(2,j,m)
        enddo
# endif
c
c-----------------------------------------------------------------------
c       blend SST between ocean domain and domain of prescribed SST
c-----------------------------------------------------------------------
c
        if (bwidth .ne. c0 .and. sbcname(m) .eq. ' sst') then
          call extrap (sbcatm(1,1,m), bzone, sor, res, imap2, jma
     &,                numpas, crits(m), 'blending SST', 2, 0)
        endif
# ifdef trace_coupled_fluxes
        write (stdout,*) ' ===> after spatially averaging:'
        call scope (sbcatm(1,1,m), imap2, imap2, jma, sbcname(m))
# endif
c
c-----------------------------------------------------------------------
c       convert to units expected by atmospheric model
c       also knock out values over atmospheric land grid.
c-----------------------------------------------------------------------
c
        if (sbcname(m) .eq. ' sst') then
          do j=1,jma
            do i=1,imap2
              sbcatm(i,j,m) = (sbcatm(i,j,m) + coabc(m))*aland(i,j)
            enddo
          enddo
        else
          do j=1,jma
            do i=1,imap2
              sbcatm(i,j,m) = (coabc(m)*sbcatm(i,j,m))*aland(i,j) 
            enddo
          enddo
        endif
# ifdef trace_coupled_fluxes
        write (stdout,*) ' ===> after converting units:'
        call scope (sbcatm(1,1,m), imap2, imap2, jma, sbcname(m))
# endif
c
c-----------------------------------------------------------------------
c      calculate averages of the atmosphere S.B.C.
c      (the S.B.C. are assumed to be defined on the same grid as "aland"
c-----------------------------------------------------------------------
c
        average(m)  = c0
        anum        = c0
        do j=1,jma
          cosdy = abcgcs(j)*abcgdy(j)
          do i=1,ima
            weight = aland(i,j)*abcgdx(i)*cosdy
            anum = anum + weight
            average(m)  = average(m) + weight*sbcatm(i,j,m)
          enddo
        enddo
        if (anum .ne. c0) average(m) = average(m)/anum
c
       enddo
c
c-----------------------------------------------------------------------
c      show averages of the atmosphere S.B.C. (they are assumed to be
c      defined on the same grid as "aland")
c-----------------------------------------------------------------------
c
      write (stdout,9700)
      do m=numosbc+1,numsbc
       write (stdout,9800) m, sbcname(m), average(m), dunits(m)       
      enddo
c
      if (errors) then
        write (stdout,9500)
        call abort()
      else
        write (stdout,9600)
      endif
c
      call toc ('mom', 'get S.B.C. for atmos (gasbc)')
8900  format (/,10x, ' ==> Getting atmosphere S.B.C.')
9100  format (/' Summary of the blending region "bzone" on atmos S.B.C.'
     &,      ' grid:'/,' 1 => indicates domain of ocean model'
     &,       ' 2 => indicates domain of prescribed SST'
     &,       /,' 0 => indicates region where SST is blended')
9500  format (/1x,' ==> Error in gasbc.F')
9600  format (/10x,' ==> S.B.C. prepared for this atmosphere segment.'/)
9700  format (/10x,' ==> S.B.C. averages for the atmosphere follow:'/)
9800  format (17x,'for S.B.C. #',i2,', the average ',a10,' is ',1pe14.7
     &,       1x, ' after converting from ',a15)
#endif
      return
      end
