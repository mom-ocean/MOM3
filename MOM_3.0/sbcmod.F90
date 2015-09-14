      module sbc_info
!=======================================================================
!
!                    S.B.C. (surface boundary conditions) info
!
!     numosbc = actual number of S.B.C. used for the ocean.
!
!     numasbc =  actual number of S.B.C. used for the atmos.
!
!     numsbc  = total number of S.B.C. used by the ocean and atmos
!
!     maxsbc  = maximum number of allowable S.B.C. fields used.
!
!
!     sbcname = S.B.C. character string names (eg: "sst", "taux", etc)
!               note that these character strings should be changed with
!               care since they are used to identify fields for special
!               treatment
!
!     dunits  = dimensional units for S.B.C. in the model
!               where the S.B.C. are applied. (eg: "deg C" for "sst" 
!               which is used as a surface boundary condition for
!               the atmos model)
!
!     ntspos  = number of time steps per ocean segment. The S.B.C are
!               held constant in time over this period. If not
!               "coupled" then ntspos = ntspas = 1
!
!     ntspas  = number of time steps per atmos segment. The S.B.C are
!               held constant in time over this period. If not
!               "coupled" then ntspas = ntspos = 1
!
#ifdef restorst
!     damp_factor = damping factor for restoring surface tracers
!                   to data. damp_factor = dzt(1)/(dampts(n)*86400)
!                   for tracer n. dampts(n) is input through namelist.
#endif
#ifdef coupled
!
!     coabc   = conversion factors for converting S.B.C.
!               from their dimensional units (in the model where they
!               were constructed) to the other model`s units
!               note: if not "coupled" then data is assumed to be 
!               in units required by MOM so... no conversion is done.
!
!     crits   = convergence criteria for forcing S.B.C. into land on the
!               model grid where they were constructed.
!
!     numpas  = number of passes used to extrapolate data into land
!               regions on the model grids where it was constructed.
!
#endif
#if !defined simple_sbc
!
!     asegs   = true on the 1st time step of an atmosphere segment.
!               otherwise false.
!
!     asege  =  true on the last time step of an atmosphere segment.
!               otherwise false.
!
!     afirst =  true on the 1st atmospheric time step of the run
!
!     alast  =  true on the last atmospheric time step of the run
!
!     iombc   = unit for reading/writing MOM boundary condition restart
!               data
!
#endif
!     
!     author:   r. c. pacanowski      e-mail=> rcp@gfdl.gov
!               m. schmidt    e-mail=> martin.schmidt@io-warnemuende.de         
!=======================================================================
!
      public setup_sbc, index_of_sbc 
!
      integer, parameter :: maxsbc = 12
!
      character (len=20), dimension(maxsbc) :: sbcname, dunits
      character (len=20), dimension(maxsbc) :: ocean_sbc, atmos_sbc
      real, dimension(:), allocatable :: damp_factor
      real, dimension(maxsbc) :: coabc, crits
      integer :: numsbc, numosbc, numasbc, ntspos, ntspas
#if !defined simple_sbc
      logical :: asegs, asege, afirst, alast
#endif
#if defined coupled
      integer :: numpas, iombc
#endif
!
      contains
!
      subroutine setup_sbc (name, m)
!
!=======================================================================
!     Ordering of Surface Boundary Conditions (S.B.C.). The order is
!     given by the order of declaration.
!
!     input:
!     name = name of S.B.C. appearing in list of possible S.B.C. defined
!            below     
!     m    = index given to previous "name" (or zero to initialize)
!
!     output:
!     m          = index given to "name".
!     sbcname(m) = name of S.B.C. #m. Other properties of "name"
!     dunits(m)  = units of "name" in the model where "name" is used. 
!                  eg: heatflux is used as a S.B.C. in the ocean with
!                      expected units of cal/cm**2/sec.
!     coabc(m)   = conversion factor which converts units of S.B.C in
!                  model where S.B.C is defined to units of S.B.C
!                  expected by model where S.B.C is used. eq: if the
!                  atmos calculates heat flux in watts/m**2 then the
!                  conversion factor is 0.2389e-4 to go from watts/m**2
!                  to cal/cm**2/sec. Conversion factors are
!                  multiplicative except for name="sst" in which case it
!                  is additive.
!     crit(m)    = convergence criteria for extrapolating S.B.C. into
!                  land on model grid where S.B.C. was constructed.
!                  As a guide, the criteria should only be large enough
!                  to give reasonable values near coastlines. This is
!                  to account for interpolating between atmos and ocean
!                  grids which are not coincident and therefore may have
!                  mismatched coastlines.
!
!     author:   r. c. pacanowski      e-mail=> rcp@gfdl.gov
!=======================================================================
!
      character (len=*) :: name
      integer :: m
      m = m + 1
!
!     list of allowable S.B.C. and associated properties
!
      if (index('taux',name) .ne. 0) then
        sbcname(m) = ' taux'
        dunits(m)  = ' dynes/cm**2'
        coabc(m)   = 1.0
        crits(m)   = 1.e-3 
      elseif (index('tauy',name) .ne. 0) then
        sbcname(m) = ' tauy'
        dunits(m)  = ' dynes/cm**2'
        coabc(m)   = 1.0
        crits(m)   = 1.e-3
      elseif (index('air press',name) .ne. 0) then
	sbcname(m) = ' air press'
	dunits(m)  = ' hPa (mbar)'
	coabc(m)   = 1.e3
	crits(m)   = 1.e-5
      elseif (index('fresh wtr',name) .ne. 0) then
	sbcname(m) = ' fresh wtr'
	dunits(m)  = ' cm/sec'
	coabc(m)   = 1.0
	crits(m)   = 1.e-8
      elseif (index('short wave',name) .ne. 0) then
        sbcname(m) = ' short wave'
        dunits(m)  = ' cal/cm**2/sec '
        coabc(m)   = 0.2389e-4
        crits(m)   = 1.e-5
      elseif (index('sun angle',name) .ne. 0) then
	sbcname(m) = ' sun angle'
	dunits(m)  = ' radian'
	coabc(m)   = 1.0
	crits(m)   = 1.e-5
      elseif (index('heatflux',name) .ne. 0) then
        sbcname(m) = ' heatflux'
        dunits(m)  = ' cal/cm**2/sec'
        coabc(m)   = 0.2389e-4
        crits(m)   = 1.e-5
      elseif (index('saltflux',name) .ne. 0) then
        sbcname(m) = ' saltflux'
        dunits(m)  = ' gm/cm**2/sec'
        coabc(m)   = 1.0
        crits(m)   = 1.e-6
      elseif (index('sst',name) .ne. 0) then
        sbcname(m) = ' sst'
        dunits(m)  = ' deg C'
        coabc(m)   = 0.0
        crits(m)   = 0.05
      elseif (index('sss',name) .ne. 0) then
        sbcname(m) = ' sss'
        dunits(m)  = ' (S-35.0)/1000'
        coabc(m)   = 1.0
        crits(m)   = 1.e-5
      elseif (index('surface u',name) .ne. 0) then
        sbcname(m) = ' surface u'
        dunits(m)  = ' cm/sec'
        coabc(m)   = 1.0
        crits(m)   = 0.1
      elseif (index('surface v',name) .ne. 0) then
        sbcname(m) = ' surface v'
        dunits(m)  = ' cm/sec'
        coabc(m)   = 1.0
        crits(m)   = 0.1
      else
        write (*,*) '=>Error: ',name,' is not in the list.'
	stop 'setup_sbc'
      endif
      end subroutine setup_sbc


      function index_of_sbc (name)
!
!=======================================================================
!     input:
!     name    = name of S.B.C. being sought
!
!     output:
!     index_of_sbc = index of "name" within "sbcname"
!
!     author:   r. c. pacanowski      e-mail=> rcp@gfdl.gov
!=======================================================================
!
      integer :: index_of_sbc
      character (len=*) name
      index_of_sbc = 0
      do n=1,numsbc
        if (index(sbcname(n),name) .ne. 0) then
	  index_of_sbc = n
	  exit
	endif
      enddo
      if (index_of_sbc .eq. 0) then
        write (*,*) '=>Error: S.B.C. ',name,' has not been setup'
        call abort()
      endif
      end function index_of_sbc 

      end module sbc_info





      module sbc_ocean
!
!=======================================================================
!     sbcocn  = S.B.C. fields (n=1..numsbc) on the ocean grid.
!                
!     author:   r. c. pacanowski      e-mail=> rcp@gfdl.gov
!               m. schmidt    e-mail=> martin.schmidt@io-warnemuende.de         
!=======================================================================
!
      real, dimension(:,:,:), allocatable :: sbcocn     
      end module



      module sbc_atmos
#ifdef coupled
!=======================================================================
!
!     sbcatm = S.B.C. fields (n=1..numsbc) on the atmospheric grid.
!
!
!     abcgx   = longitudes (degs) for the boundary condition fields on
!               the atmosphere grid. (similar to "xt" in the ocean)
!
!     abcgy   = latitudes (degs) for the boundary condition fields on
!               the atmosphere grid. (similar to "yt" in the ocean)
!
!     abcgcs  = cosine of "abcgy"
!
!     abcgdx  = width (degs) of the boundary condition grid boxes
!               in the longitude direction on the atmosphere grid.
!
!     abcgdy  = height (degs) of the boundary condition grid boxes
!               in the latitude direction on the atmosphere grid.
!
!     aland   = mask of (0,1) for (land,ocean) on the atmospheric grid
!
!     isocn, ieocn, jsocn, jeocn locate the ocean grid domain in terms
!     of atmosphere grid indices
!
!     isocn   = index of the first atmos grid point (in longitude) that
!               falls within the ocean domain
!     ieocn   = index of the last atmos grid point (in longitude) that
!               falls within the ocean domain
!     jsocn   = index of the first atmos grid point (in latitude) that
!               falls within the ocean domain
!     jeocn   = index of the last atmos grid point (in latitude) that
!               falls within the ocean domain
!
!-----------------------------------------------------------------------
!     this section is only for embedding regional ocean domains within
!     global atmosheric domains
!
!     sstpre  = prescribed value of "sst" outside of the ocean model
!               domain. this is used only when the ocean model is NOT
!               global in extent (the atmospheric model needs "sst"
!               globally). see "bwidth" below for further explanation.
!
!     bwidth  = when the ocean domain is NOT global, "bwidth" specifies
!               the width of a blending zone (in degrees) surrounding
!               the ocean domain. the purpose of this is to provide a
!               smooth transition between "sst" inside the ocean and
!               the prescribed sst "sstpre" outside of the ocean (since
!               the atmosphere always requires global "sst" as a b.c.).
!               for example:
!               when the ocean domain is a cyclic strip between two
!               latitudes, "sstpre" can be set to a polar temperature
!               (outside the ocean domain) and "bwidth" can be set to
!               extend to within one atmospheric grid point of the pole.
!               when the ocean domain is a basin, "sstpre" should
!               really be changed to a function of latitude and
!               and longitude.
!               setting "bwidth" = 0 implies a global ocean with no
!               blending zone.
!
!     bzone   = is defined by setting a non zero "bwidth". this assumes
!               the ocean is NOT global (see "bwidth"). "bzone"
!               defines buffer zones exterior to the ocean domain when
!               the ocean domain is less than global. "bzone" is
!               defined on the atmosphere b.c. grid. it defines
!               three regions as shown below:
!
!                 bzone=0 is the blending zone of width "bwidth"
!                 bzone=1 is the non global ocean domain
!                 bzone=2 is everything external to the blending zone
!
!               Note: these zones are automatically defined by the
!               grid definitions of atmos & ocean models. the only
!               freedom is to choose "bwidth" (described above)
!
!                x--------------------------------------------x
!                |      global atmos domain                   |
!                |        (bzone=2)                           |
!                |    x----------------------------------x    |
!                |    |          blending zone           |    |
!                |    |          (bzone=0)               |    |
!                |    |    x------------------------x    |    |
!                |    |    |                        |    |    |
!                |    |    |       ocean domain     |    |    |
!                |    |    |        (bzone=1)       |    |    |
!                |    |    |                        |    |    |
!                |    |    x------------------------x    |    |
!                |    |                                  |    |
!                |    |                                  |    |
!                |    x----------------------------------x    |
!                |                                            |
!                |                                            |
!                x--------------------------------------------x
!
!     author:   r. c. pacanowski      e-mail=> rcp@gfdl.gov
!               m. schmidt    e-mail=> martin.schmidt@io-warnemuende.de         
!=======================================================================
!
# include "aparam.h"
      integer, dimension (imap2,jma) :: aland, bzone
      integer :: isocn, ieocn, jsocn, jeocn
      real, dimension (imap2) :: abcgx,  abcgdx
      real, dimension (jma) ::abcgy, abcgdy, abcgcs
      real :: sstpre, bwidth
      real, dimension(:,:,:), allocatable :: sbcatm     
#endif
      end module


