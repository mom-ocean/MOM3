      subroutine atmos
c
c=======================================================================
c     this "atmos" model is an EXAMPLE of how to supply surface
c     boundary conditions which are a function of latitude, longitude
c     and time to MOM. The "atmos" model grid is assumed to be of
c     different resolution from MOM so interpolations of S.B.C. are
c     needed. This is accomplished by specifying the "-Dcoupled" option
c     when compiling MOM.
c
c     This model is for illustrative purposes only but could be replaced
c     by an atmospheric GCM, intermediate, or anomoly model. 
c
c      "atmos" simply prescribes annual mean zonally averaged wind 
c      stress along with a newtonian type heat flux which drives the
c      ocean back towards annual mean zonally averaged sst. the  
c      salt flux into the ocean is arbitrarily prescribed zero.
c      if it wasn`t, then the precip - evap + runoff of the "atmos"
c      would have to be converted to ...
c      salt flux  = -(P-E+R) * rho * Sref
c      where...
c      P-E+R represents a precipitation minus evaporation
c            plus runoff rate in cm of water/sec
c      rho is the density if water taken to be 1 g/cm**3, and
c      Sref is a reference salinity in units of grams of salt per
c            gram of water (in units of "parts per part" such
c            as 0.035 ... not "parts per thousand")
c      Depending upon the application of interest, users may wish
c      to set Sref to either be a constant over the entire model
c      domain or for it to be the locally predicted salinity of
c      the uppermost model level (S(1)).  If one desires to ensure
c      that a global average P-E+R flux of zero translates to zero
c      trend in the salt content of the ocean, then a constant Sref
c      is needed.
c      An alternative is to prescribe a fresh water flux. This is
c      under development. Refer to the manual.      
c
c      the atmospheric model must do the following:
c      
c
c       1) predict variables for one atmospheric time step per call
c
c       2) accumulate the ocean surface boundary conditions every call
c
c     
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use sbc_info
      use sbc_atmos
#include "param.h"
c
      dimension atmvar(4,jma)
c
c-----------------------------------------------------------------------
c     predict the state of the atmosphere for one atmospheric time step.
c     typically, the atmosphere boundary condition SST would influence
c     the atmosphere. However, in the test case atmosphere nothing
c     happens.
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c     after the time step is complete, calculate the ocean S.B.C.
c     use the "abcgy" array to fetch some zonal meansfor this test case.
c     assume atmvar(,1) = taux
c            atmvar(,2) = tauy
c            atmvar(,3) = tstar (some atmospheric surface temperature)  
c-----------------------------------------------------------------------
c
      do j=1,jma
        y1 = abcgy(j)
        call bcest (y1, y1, atmvar(1,j))
      enddo
c
c
c-----------------------------------------------------------------------
c     accumulate atmospheric output as a S.B.C. for driving the
c     ocean during the next ocean segment.
c
c     ordering of S.B.C. was set in driver.F 
c
c     note: only the interior longitudes (2 through imap2-1) need to be
c     set.  recall that "sbcatm" is offset by one because of adding two
c     points for the east-west cyclic condition.
c
c     note: accumulate the b.c. only over non land areas to preserve the
c     guess for extrapolation in "gasbc.F". assume all s.b.c. for the 
c     ocean are defined on the land/sea mask grid.
c
c     note: if the atmosphere grid goes from North Pole (j=1) to South
c           Pole (j=jma), then all references to "sbcatm"
c           should use inverted indices. For example:
c           jj = jma+1-j
c           sbcatm(i,jj,.) = sbcatm(i,jj,.) ... for accumulating
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     initialize all S.B.C. for the ocean at the beginning of each
c     segment. ordering of S.B.C. was set in driver.F
c-----------------------------------------------------------------------
c
      if (asegs) then
        do n=1,numosbc
          do j=1,jma
            do i=1,imap2
	      if (aland(i,j) .ne. 0) sbcatm(i,j,n) = c0
            enddo
	  enddo
	enddo
      endif
c
c-----------------------------------------------------------------------
c     accumulate all S.B.C. on every time step within the segment
c
c         note: "i" in "sbcatm"  corresponds to "i-1" in the
c              atmosphere model grid because of the cyclic points which
c              were added to "bctatm". it does not matter here because
c              zonal means are used.
c-----------------------------------------------------------------------
c
      itaux = index_of_sbc ('taux')
      itauy = index_of_sbc ('tauy')
      ihflx = index_of_sbc ('heatflux')
      ipme  = index_of_sbc ('saltflux')
      isst  = index_of_sbc ('sst')
c
      if (itaux .ne. 0) then
        do j=1,jma
          do i=2,imap2-1
            if (aland(i,j) .ne. 0) then
              sbcatm(i,j,itaux) = sbcatm(i,j,itaux) + atmvar(1,j)
            endif
          enddo
	enddo
      endif
c
      if (itauy .ne. 0) then
        do j=1,jma
          do i=2,imap2-1
            if (aland(i,j) .ne. 0) then
              sbcatm(i,j,itauy) = sbcatm(i,j,itauy) + atmvar(2,j)
            endif
          enddo
	enddo
      endif
c
c     for heatflux, choose an arbitrary mixed layer depth and e-folding
c     time. no fresh water flux
c
      if (ihflx .ne. 0) then
        daysec = 86400.0
        thick  = 30.e2
        efold  = 50.0*daysec
        fx     = thick/efold
        do j=1,jma
          do i=2,imap2-1
            if (aland(i,j) .ne. 0) then
              sbcatm(i,j,ihflx) = sbcatm(i,j,ihflx) +fx*( atmvar(3,j)
     &                                         -sbcatm(i,j,isst))
            endif
          enddo
	enddo
      endif
c
      if (ipme .ne. 0) then
        do j=1,jma
          do i=2,imap2-1
            if (aland(i,j) .ne. 0) then
              sbcatm(i,j,ipme) = c0
            endif
          enddo
	enddo
      endif
c
c-----------------------------------------------------------------------
c     average all S.B.C. at the end of each atmospheric segment
c     "n" is the S.B.C. and "m" refers to its ordering within arrays
c-----------------------------------------------------------------------
c
      if (asege) then
        rts = c1/ntspas
        do n=1,numosbc
          do j=1,jma
            do i=2,imap2-1
              if (aland(i,j) .ne. 0) sbcatm(i,j,n) = rts*sbcatm(i,j,n)
            enddo
          enddo
	enddo
      endif
c
      return
      end
