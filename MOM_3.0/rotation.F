#ifdef rot_grid
# ifdef drive_rotation
#  define driver_only
#  include "util.F"
#  include "grids.F"
#  include "size_check.F"
      program driver
c
c=======================================================================
c
c                    SPHERICAL GRID ROTATION MODULE
c
c     To generate latitude data needed when using a rotated grid:
c
c        1) specify a geographic latitude and longitude for the position
c           of the rotated north pole and a prime meridian point.     
c
c        2) setup the rotated grid. (see grids.F)
c    
c        3) compile and run this module using the "run_rotation" script
c
c     running this module generates the latitude data file needed to
c     correct the Coriolis term when running MOM with the "rot_grid"
c     option. other output shows details of the rotated model grid and 
c     examples of how to interpolate scalar and vector quantities from 
c     a geographic grid to the rotated model grid (required when 
c     providing forcing data). double precision is strongly 
c     recommended when using these routines.
c
c     author:   M. Eby            e-mail eby@uvic.ca
c=======================================================================
c
#  include "size.h"
#  include "coord.h"
#  include "stdunits.h"
c
c     dimension local data
c
      dimension glat(imt,jmt), tr(imt,jmt), ur(imt,jmt,2)
c
c-----------------------------------------------------------------------
c     USER INPUT ===>  read this section and specify the rotation.
c-----------------------------------------------------------------------
c
c     any spherical grid rotation can be achieved by specifying three, 
c     Euler angle, solid body rotations. the first rotation angle phir
c     is defined as a rotation about the original z axis. the second 
c     rotation angle thetar is defined as a rotation about the new x 
c     axis. the final rotation angle psir is defined as a rotation about
c     the new z axis. one method of deciding on the values of these 
c     angles is to define the geographic (unrotated) latitude and 
c     longitude of the north pole of the rotated grid. this defines 
c     the angles phir and thetar. The latitude of the new grid is now
c     specified but the longitude is still arbitrary. to define the new 
c     grid longitude, specify the geographic location of any point on 
c     the rotated grid`s prime meridian. the angle psir is then 
c     calculated as the rotation angle necessary to place this point on 
c     the prime meridian of the rotated grid.
c
c     for example, to set up a model with an equatorial grid over the 
c     Arctic and North Atlantic, one could specify the rotated grid pole
c     to be at 0 N, 110 W and a prime meridian point to be at 0 N O E. 
c     This specifies a grid rotation in which the new grid equator is 
c     defined along the 20 W and 160 E meridians. The rotated grid 
c     longitude would be east, north of the geographic equator and west 
c     to the south. The geographic pole would be at 0 N and 90 E.
c
c     latitudes and longitudes should be specified such that N and E are
c     positive and S and W are negative.
c 
c     gltnp = geographic latitude of the rotated grid north pole 
c     glnnp = geographic longitude of the rotated grid north pole
c     gltpp = geographic latitude of a prime point on the rotated grid
c     glnpp = geographic longitude of a prime point on the rotated grid
c
      gltnp =    0.
      glnnp = -120.
      gltpp =    0.
      glnpp =  -30.
c
c-----------------------------------------------------------------------
c     ==> end of USER INPUT
c-----------------------------------------------------------------------
c
      stdout = 6
c
      write(stdout,'(/a/)') '=>Warning: 64bit precision is recommended'
c
      rad = acos(-1.)/180.
c
c-----------------------------------------------------------------------
c     calculate the Euler angles required for the rotation
c-----------------------------------------------------------------------
c
      phir   = (glnnp - 90.)*rad
      thetar = (gltnp - 90.)*rad
c
c     set rotation angle psir to zero and rotate prime meridian point to
c     find the required angle psir to give the point a longitude of zero
c
      psir = 0.
      call rotate (gltpp, glnpp, phir, thetar, psir, rltpp, rlnpp)
      psir = rlnpp*rad
c
c-----------------------------------------------------------------------
c     write out the grid definitions and rotation angles
c-----------------------------------------------------------------------
c
      call rotate (gltnp, glnnp, phir, thetar, psir, rltnp, rlnnp)
      call rotate (gltpp, glnpp, phir, thetar, psir, rltpp, rlnpp)
c
      write(stdout,'(/,a,/)') 'Rotated Model Grid Definitions:'
      write(stdout,'(16x,a)') '    Geographic              Rotated     '
      write(stdout,'(16x,a)')'Latitude  Longitude   Latitude  Longitude'
      write(stdout,'(a,4f11.3)') 'North Pole: ', gltnp, glnnp
     &,                           rltnp, rlnnp
      write(stdout,'(a,4f11.3)') 'Prime point:', gltpp, glnpp
     &,                           rltpp, rlnpp
c
      write(stdout,'(//,a,/)') 'Rotation Angles:'
      write(stdout,'(a,f11.3,a)') 'phir:       ', phir/rad, ' deg'
      write(stdout,'(a,f11.3,a)') 'thetar:     ', thetar/rad, ' deg'
      write(stdout,'(a,f11.3,a)') 'psir:       ', psir/rad, ' deg'
c
c-----------------------------------------------------------------------
c     set up model grids in x (longitude), y (latitude), and z (depth)
c     corresponding to Arakawa "b" gird system
c-----------------------------------------------------------------------
c
      maxlen = max(imt, jmt, km)
      call gcoord (maxlen, imt2, jmt2, km2
     &,            dxtdeg, dytdeg, dxudeg, dyudeg
     &,            dzt, dzw, xt, xu, yt, yu, zt, zw)
c
c     check that returned grid sizes match those in file "size.h"
c
      call size_check (imt2, jmt2, km2, 'rotation', 'stop')
c
c-----------------------------------------------------------------------
c     calculate and write the geographic latitude of the model grid 
c-----------------------------------------------------------------------
c
      open (unit=10, file='lat.dta', status='UNKNOWN',
     &      form='UNFORMATTED')
c
c     latitude of velocity grid is needed to calculate Coriolis
c
      do j=1,jmt
        do i=1,imt
          call rotate (yu(j), xu(i), -psir, -thetar, -phir, glt, gln)
          glat(i,j) = glt
        enddo
      enddo
c
      write (10) imt,jmt
      write (10) ((glat(i,j),i=1,imt),j=1,jmt)
#  ifdef neptune
c
c     latitude of tracer grid is needed to calculate neptune velocities
c
      do j=1,jmt
        do i=1,imt
          call rotate (yt(j), xt(i), -psi, -theta, -phi, glt, gln)
          glat(i,j) = glt
        enddo
      enddo
c
      write (10) ((glat(i,j),i=1,imt),j=1,jmt)
      
#  endif
      close (10)
c
c-----------------------------------------------------------------------
c     write the locations of the corners of the rotated model grid
c-----------------------------------------------------------------------
c
      write(stdout,'(/,/,a,/)') 'Rotated Model Grid Corners:'
      write(stdout,'(16x,a)') '    Geographic              Rotated  '
      write(stdout,'(16x,a)') 
     &                'Latitude  Longitude   Latitude  Longitude'
      call rotate (yt(jmt), xt(1), -psir, -thetar, -phir, glt, gln)
      write(stdout,'(a,4f11.3)') 'NW corner:  ', glt, gln, yt(jmt)
     &,                             xt(1)
      call rotate (yt(jmt), xt(imt), -psir, -thetar, -phir, glt, gln)
      write(stdout,'(a,4f11.3)') 'NE corner:  ', glt, gln, yt(jmt)
     &,                             xt(imt)
      call rotate (yt(1), xt(1), -psir, -thetar, -phir, glt, gln)
      write(stdout,'(a,4f11.3)') 'SW corner:  ', glt, gln, yt(1), xt(1)
      call rotate (yt(1), xt(imt), -psir, -thetar, -phir, glt, gln)
      write(stdout,'(a,4f11.3,//)') 'SE corner:  ', glt, gln, yt(1)
     &,                             xt(imt)
c
c-----------------------------------------------------------------------
c     show how to do interpolation of scalar and vector quantities from
c     a geographic grid to the rotated model grid
c-----------------------------------------------------------------------
c
      call eg (phir, thetar, psir, xt, yt, xu, yu, imt, jmt, tr, ur)
c
      end

      subroutine eg (phir, thetar, psir, xtr, ytr, xur, yur, imtr, jmtr
     &,              tr, ur)
c
c=======================================================================
c     e.g. of how to interpolate from a geographic global data set
c    (with resolution del) to a rotated model grid.
c=======================================================================
c
      parameter (del=5, imtg=360/del, jmtg=180/del)
# include "stdunits.h"
c
c     geographic data grid variables
c
      dimension tg(imtg,jmtg), xtg(imtg), ytg(jmtg)
      dimension ug(imtg,jmtg,2), xug(imtg), yug(jmtg)
c
c     rotated model grid variables
c
      dimension tr(imtr,jmtr), xtr(imtr), ytr(jmtr)
      dimension ur(imtr,jmtr,2), xur(imtr), yur(jmtr)
c
      pi = 4.*atan(1.)
c
c     define data grid latitude and longitude offsets
c
      ytgoff = -90. + del/2.
      xtgoff = del/2.
c
c     define the data grid (assume vectors are defined at the same 
c     points as scalers on the data grid for this simple example) 
c
      do j=1,jmtg
        ytg(j) = (j-1)*del + ytgoff
        yug(j) = ytg(j)
      enddo
      do i=1,imtg
        xtg(i) = (i-1)*del + xtgoff
        xug(i) = xtg(i)
      enddo
c
c     define sample data
c
      do j=1,jmtg
        do i=1,imtg
c     set scaler data to data grid latitude
	  tg(i,j)  = ytg(j)
c     set vector data to point to the geographic north pole
	  ug(i,j,1) = 0.
	  ug(i,j,2) = 1.
        enddo
      enddo
c
c     interpolate scalar data "tg" to "tr" on rotated model grid 
      call rot_intrp_sclr (tg, xtg, ytg, imtg, jmtg, tr, xtr, ytr
     &,                    imtr, jmtr, psir, thetar, phir)
c
c     display scaler data in "tg" on geographic data grid
      scale = 1.
      write (stdout,'(/,a,/)') 
     & ' Geographic latitude on geographic data grid:'
      call matrix (tg, imtg, 1, imtg, -1, -jmtg, scale)
c
c     display scaler data in "tr" on rotated model grid
      scale = 1.
      write (stdout,'(/,a,/)') 
     & ' Interpolated geographic latitude on rotated model grid:'
      call matrix (tr, imtr, 1, imtr, -1, -jmtr, scale)
c
c     interpolate vectors "ug" to "ur" on rotated model grid
      call rot_intrp_vctr (ug, xug, yug, imtg, jmtg, ur, xur, yur
     &,                    imtr, jmtr, psir, thetar, phir)
c
c     display vector data in "ur" on rotated model grid
      scale = 1.
      write (stdout,'(/,a,/,a,/)') 
     &  ' Interpolated zonal component of a geographical northward'
     &, ' pointing unit vector on rotated model grid:'
      call matrix (ur(1,1,1), imtr, 1, imtr, -1, -jmtr, scale)
      write (stdout,'(/,a,/,a,/)') 
     &  ' Interpolated meridional component of a geographical northward'
     &, ' pointing unit vector on rotated model grid:'
      call matrix (ur(1,1,2), imtr, 1, imtr, -1, -jmtr, scale)
c
      return
      end
# endif


      subroutine rot_intrp_vctr (g, xg, yg, ig, jg, r, xr, yr
     &,                          ir, jr, psir, thetar, phir)
c
c=======================================================================
c     interpolate vector data from an geographic data grid to a
c     rotated model grid
c
c     input
c     psir, thetar, phir = Euler angles defining rotation
c     g  = vector on geographic data grid
c     xg = longitude of data points on geographic data grid
c     yg = latitude of data points on geographic data grid
c     ig = number of longitudes in geographic data grid 
c     jg = number of latitudes in geographic data grid 
c     xr = longitude of points on rotated model grid
c     yr = latitude of points on rotated model grid
c     ir = number of longitudes in rotated model grid 
c     jr = number of latitudes in rotated model grid 
c
c     output
c     r  = vector on rotated model grid
c
c     internal
c     (rln,rlt) = (longitude,latitude) in rotated coordinates
c     (gln,glt) = (longitude,latitude) in geographic coordinates
c     xg(iw) = point on the geographic grid to the west of (gln,glt)
c     xg(ie) = point on the geographic grid to the east of (gln,glt)
c     yg(js) = point on the geographic grid to the south of (gln,glt)
c     yg(jn) = point on the geographic grid to the north of (gln,glt)
c
c=======================================================================
c
      dimension g(ig,jg,2), xg(ig), yg(jg)
      dimension r(ir,jr,2), xr(ir), yr(jr)
c
      rad = acos(-1.)/180.
c
c     interpolate vector components as scalers on rotated model grid 
c
      call rot_intrp_sclr (g(1,1,1), xg, yg, ig, jg, r(1,1,1), xr, yr
     &,                          ir, jr, psir, thetar, phir)
      call rot_intrp_sclr (g(1,1,2), xg, yg, ig, jg, r(1,1,2), xr, yr
     &,                          ir, jr, psir, thetar, phir)
c
c     correct vector direction
c
      do j=1,jr
        do i=1,ir
          vmag = sqrt(r(i,j,1)**2 + r(i,j,2)**2)
          if (vmag .gt. 0.) then 
            a = r(i,j,1)/vmag
            a = min(a, 1.)
            a = max(a, -1.)
            a = acos(a)
            if (r(i,j,2) .lt. 0.) a = -a
           call rotate (yr(j), xr(i), -psir, -thetar, -phir, glt, gln)
           call rotvec(glt, gln, phir, thetar, psir, angle)
            a = a + angle*rad
            r(i,j,1) = vmag*cos(a)
            r(i,j,2) = vmag*sin(a)
          else
	    r(i,j,1) = 0.
	    r(i,j,2) = 0.
          endif
        enddo
      enddo
c
      return
      end


      subroutine rot_intrp_sclr (g, xg, yg, ig, jg, r, xr, yr
     &,                          ir, jr, psir, thetar, phir)
c
c=======================================================================
c     interpolate scaler data from an geographic data grid to a 
c     rotated model grid
c
c     input
c     psir, thetar, phir = Euler angles defining rotation
c     g  = scaler on geographic data grid
c     xg = longitude of data points on geographic data grid
c     yg = latitude of data points on geographic data grid
c     ig = number of longitudes in on geographic data grid 
c     jg = number of latitudes in on geographic data grid 
c     xr = longitude of points on rotated model grid
c     yr = latitude of points on rotated model grid
c     ir = number of longitudes in rotated model grid 
c     jr = number of latitudes in rotated model grid 
c
c     output
c     r  = scaler on rotated model grid
c
c     internal
c     (rln,rlt) = (longitude,latitude) in rotated coordinates
c     (gln,glt) = (longitude,latitude) in geographic coordinates
c     xg(iw) = point on the geographic grid to the west of (gln,glt)
c     xg(ie) = point on the geographic grid to the east of (gln,glt)
c     yg(js) = point on the geographic grid to the south of (gln,glt)
c     yg(jn) = point on the geographic grid to the north of (gln,glt)
c
c=======================================================================
c
# include "stdunits.h"
      dimension g(ig,jg), xg(ig), yg(jg), r(ir,jr), xr(ir), yr(jr)
c
      epsln = 1.e-10
      glt_min = 90.
      glt_max = -90.
      gln_min = 360.
      gln_max = -360.
      ln_err = 0
      lt_err = 0
c
c     find longitude points of data within interval [0., 360.]
      istrt = 1
      do i=2,ig
        if (xg(i-1) .lt. 0. .and. xg(i) .ge. 0.) istrt = i
      enddo 
      iend = ig
      do i=2,ig
        if (xg(i-1) .lt. 360. .and. xg(i) .ge. 360.) iend = i
      enddo
c
c     find latitude points of data within interval [-90., 90.]
      jstrt = 1
      do j=2,jg
        if (yg(j-1) .lt. -90. .and. yg(j) .ge. -90.) jstrt = j
      enddo 
      jend = jg
      do j=2,jg
        if (yg(j-1) .lt. 90. .and. yg(j) .ge. 90.) jend = j
      enddo
c
c     interpolate data to model grid 
c
      do j=1,jr
        do i=1,ir
          call rotate (yr(j), xr(i), -psir, -thetar, -phir, glt, gln)
          if (gln .lt. 0.) gln = gln + 360.
          if (gln .ge. 360.) gln = gln - 360.
          glt_min = min(glt,glt_min)
          glt_max = max(glt,glt_max)
          gln_min = min(gln,gln_min)
          gln_max = max(gln,gln_max)
c
          iw = indp (gln, xg(istrt), iend - istrt + 1) + istrt - 1
          if (xg(iw) .gt. gln) iw = iw - 1
	  ie = iw + 1
	  if (iw .ge. istrt .and. ie .le. iend) then
	    del = xg(ie) - xg(iw)
  	    wtw = (xg(ie) - gln)/del
	  else
c     east or west of the last data value. this could be because a
c     cyclic condition is needed or the dataset is too small. in either 
c     case apply a cyclic condition
            ln_err = 1
	    iw = iend
            ie = istrt
	    del = xg(ie) + 360. + epsln - xg(iw) 
            if (xg(ie) .ge. gln) then
  	      wtw = (xg(ie) - gln)/del
            else
	      wtw = (xg(ie) + 360. + epsln - gln)/del
            endif
	  endif
	  wte = 1. - wtw
c
          js = indp (glt, yg(jstrt), jend - jstrt + 1) + jstrt - 1
	  if (yg(js) .gt. glt) js = max(js - 1,jstrt)
	  jn = min(js + 1,jend)
	  if (yg(jn) .ne. yg(js) .and. yg(js) .le. glt) then
            wts = (yg(jn) - glt)/(yg(jn) - yg(js))
          else
c     north or south of the last data value. this could be because a
c     pole is not included in the data set or the dataset is too small.
c     in either case extrapolate north or south
            lt_err = 1
            wts = 1.
          endif
	  wtn = 1. - wts
c
          r(i,j) = g(ie,jn)*wte*wtn + g(ie,js)*wte*wts
     &           + g(iw,jn)*wtw*wtn + g(iw,js)*wtw*wts
c
        enddo
      enddo
c
      if (ln_err .eq. 1) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: the geographic data set does not extend far   '
     &,   '             enough east or west - a cyclic boundary       '
     &,   '             condition was applied. check if appropriate   '
        write (stdout,'(/,(1x,a,2f8.2))')
     &    '    data required between longitudes:', gln_min, gln_max
     &,   '      data set is between longitudes:', xg(istrt), xg(iend)
      endif
c
      if (lt_err .eq. 1) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: the geographic data set does not extend far   '
     &,   '             enough north or south - extrapolation from    '
     &,   '             the nearest data was applied. this may create '
     &,   '             artificial gradients near a geographic pole   ' 
        write (stdout,'(/,(1x,a,2f8.2))')
     &    '    data required between latitudes:', glt_min, glt_max
     &,   '      data set is between latitudes:', yg(jstrt), yg(jend)
      endif
c
      return
      end


      subroutine rotate (glt, gln, phir, thetar, psir, rlt, rln)
c
c=======================================================================
c     subroutine rotate takes a geographic latitude and longitude and 
c     finds the the equivalent latitude and longitude on a rotated grid.
c     when going from a geographic grid to a rotated grid, all of the 
c     defined rotation angles given to rotate by the calling program 
c     are positive, but when going from a rotated grid back to the 
c     geographic grid, the calling program must reverse the angle order 
c     (phir and psir are switched) and all of the angles made negative.
c
c     the first rotation angle phir is defined as a rotation about the
c     original z axis. the second rotation angle thetar is defined as a
c     rotation about the new x axis. the final rotation angle psir is
c     defined as a rotation about the new z axis. these rotation angles
c     are just the Euler angles as defined in "classical mechanics"
c     Goldstein (1951).
c
c     author:   M. Eby            e-mail eby@uvic.ca
c=======================================================================
c
# include "stdunits.h"
c     g...  = geographic value
c     r...  = rotated value
c     ...lt = latitude (or equivalent spherical coordinate)
c     ...ln = longitude (or equivalent spherical coordinate)
c     ...x  = x coordinate
c     ...y  = y coordinate
c     ...z  = z coordinate
c     psir, thetar, phir = Euler angles defining rotation
c
c     define rad for conversion to radians.
      rad = acos(-1.)/180.
c
c     convert latitude and longitude to spherical coordinates
      thetas = gln
      if (thetas .gt. 180.) thetas = thetas - 360.
      if (thetas .lt. -180.) thetas = thetas + 360.
      phis = (90. - glt)*rad
      thetas = thetas*rad
c
c     translate point into Cartesian coordinates for rotation.
      gx = sin(phis)*cos(thetas)
      gy = sin(phis)*sin(thetas)
      gz = cos(phis)
c
c     rotate the point (gx, gy, gz) about the z axis by phir then the x
c     axis by thetar and finally about the z axis by psir.
c 
      rx = gx*(cos(psir)*cos(phir) - cos(thetar)*sin(phir)*sin(psir)) +
     &     gy*(cos(psir)*sin(phir) + cos(thetar)*cos(phir)*sin(psir)) +
     &     gz*sin(psir)*sin(thetar)
c
      ry = gx*(-sin(psir)*cos(phir) - cos(thetar)*sin(phir)*cos(psir)) +
     &     gy*(-sin(psir)*sin(phir) + cos(thetar)*cos(phir)*cos(psir)) +
     &     gz*(cos(psir)*sin(thetar))
c
      rz = gx*(sin(thetar)*sin(phir)) + gy*(-sin(thetar)*cos(phir)) +
     &     gz*(cos(thetar))
c
c     convert rotated point back to spherical coordinates
c
c     check for rounding error (arccos(x): abs(x) must be .le. 1)
      rz = min(rz, 1.)
      rz = max(rz, -1.)
      rlt = acos(rz)
c     if point is at a pole set rotated longitude equal to initial.
      if (rlt .le. 0. .or. rlt .ge. 180.*rad) then
        rln = thetas
      else
c     if rln lies between -135 and -45 or between 45 and 135 degrees
c     it is more accurate to use an arccos calculation.
        if (abs(rx/sin(rlt)) .lt. cos(45.*rad)) then
          rln = rx/sin(rlt)
c     check for rounding error (arccos(x): abs(x) must be .le. 1)
          rln = min(rln, 1.)
          rln = max(rln, -1.)
          rln = acos(rln)
c     arccos will give rln between 0 and 180 degrees.  if the point
c     is negative in y, rln must be equal to negative rln.
          if (ry .lt. 0.) rln = -rln
        else
c     if rln lies between -45 and 45 or between 135 and -135 degrees
c     it is more accurate to use an arcsin calculation.
          rln = ry/sin(rlt)
c     check for rounding error (arcsin(x): abs(x) must be .le. 1)
          rln = min(rln, 1.)
          rln = max(rln, -1.)
          rln = asin(rln)
c     arcsin will give rln between -90 and 90 degrees. if the point
c     is negative in x, rln must be equal to 180 degrees minus rln.
          if (rx .lt. 0.) rln = 180.*rad - rln
        endif
      endif
c
c     convert back to degrees of latitude and longitude.
      rlt = 90. - rlt/rad
      rln = rln/rad
      if (rln .gt. 180.) rln = rln - 360.
      if (rln .le. -180.) rln = rln + 360.
c
      return
      end


      subroutine rotvec (glt, gln, phir, thetar, psir, angle)
c
c=======================================================================
c     subroutine rotvec takes a geographic latitude and longitude and 
c     finds the the vector rotation angle angle (in degrees) for a 
c     vector on the rotated grid. when going from the geographic to a 
c     rotated grid, all of the defined rotation angles given to rotvec 
c     by the calling program are positive, but when going from a 
c     rotated grid back to the geographic grid, the calling program 
c     must reverse the angle order (phir and psir are switched) and all 
c     of the angles made negative. if a pole is detected then an angle 
c     of zero is returned.
c
c     rotvec rotates the point defining the head of a very short north
c     or south pointing geographic vector. the angle between this 
c     vector and a similar direction vector defined in the new grid is 
c     calculated using the law of cosines. the accuracy of this 
c     calculation depends on the size of the direction vector (delta) 
c     and the precision of the computation. the smaller the vector the 
c     more accurate the calculation but then the more precision 
c     required. double precision is strongly recommended.
c
c     author:   M. Eby            e-mail eby@uvic.ca
c=======================================================================
c
# include "stdunits.h"
c     if gx and gy are vector components at glt and gln, the corrected 
c     vector components rx and ry at rlt and rln can be calculated as 
c     follows:
c
c      rad = acos(-1.)/180.
c      call rotate (rlt, rln, -psir, -thetar, -phir, glt, gln)
c      call get_vector (glt, gln, gx, gy)   ! some routine to get vector
c      r = sqrt(gx**2 + gy**2)
c      if (r .gt. 0.) then 
c        a = gx/r
c        a = min(a, 1.)
c        a = max(a, -1.)
c        a = acos(a)
c        if (gy .lt. 0.) a = -a
c        call rotvec(glt, gln, phir, thetar, psir, angle)
c        a = a + angle*rad
c        rx = r*cos(a)
c        ry = r*sin(a)
c      endif
c
c     g...   = geographic value
c     r...   = rotated value
c     ...lt  = latitude
c     ...ln  = longitude
c     ...lth = latitude of head of vector
c     ...lnh = longitude of head of vector
c     dst    = distance between heads of vectors
c     delta  = length of vector
c     angle  = angle to rotate vectors back to original orientation
c     psir, thetar, phir = Euler angles defining rotation
c
c     define multiplier rad for conversion to radians.
      rad = acos(-1.)/180.
c
c     define length of direction vector (single precision may require a 
c     longer and thus less accurate vector length).
      delta = 0.001
c
c     if the base is in the north of the geographic grid use a south 
c     pointing vector to avoid any possibility of going over the pole.
      if (glt .ge. 0.) delta = -delta
c
c     find the base of the direction vectors in the rotated grid.
      call rotate (glt, gln, phir, thetar, psir, rlt, rln)
c
c     if base in the rotated grid is near a pole return an angle of zero
      if (abs(rlt) .ge. 90.-abs(delta)) then
        angle = 0.
        return
      endif
c
c     find the head of the geographic grid direction vector in the 
c     rotated grid.
      call rotate (glt+delta, gln, phir, thetar, psir, glth, glnh)
c
c     if the base is in opposite hemispheres switch the vector
c     direction for better accuracy.
      if (glt*rlt .lt. 0) delta = -delta
c
c     find the head of the rotated grid direction vector.
      rlth = rlt + delta
      rlnh = rln
c
c     find the distance between the heads of the direction vectors.
      call dist (glth, glnh, rlth, rlnh, dst)
c
c     find the angle between direction vectors with the law of cosines.
      delta = abs(delta)
      angle = (cos(dst)-cos(delta)**2)/(sin(delta)**2)
      angle = min(angle, 1.)
      angle = max(angle, -1.)
      angle = acos(angle)/rad
      t = abs(delta)
c
c     adjust the angle if the direction vectors are opposite.
      if (glt*rlt .lt. 0) angle = 180. - angle
c     determine the sign of the angle by checking the offset longitudes.
      if (glnh+360. .gt. rlnh+360.) angle = -angle
c     change sign if the original direction vector was pointing south.
      if (glt .ge. 0.) angle = -angle
c
      return
      end


      subroutine dist (lat1, lng1, lat2, lng2, dst)
c
c=======================================================================
c     subroutine dist calculates the arc distance between two 
c     points given their latitudes and longitudes
c=======================================================================
c
      real lat1, lng1, lat2, lng2, dst, rad
c
c     define multiplier rad for conversion to radians.
      rad = acos(-1.)/180.
c
c     check input.
      lat1 = min(lat1, 90.)
      lat1 = max(lat1, -90.)
      if (lng1 .lt. 0.) lng1 = lng1 + 360.
      if (lng1 .gt. 360.) lng1 = lng1 - 360.
      lat2 = min(lat2, 90.)
      lat2 = max(lat2, -90.)
      if (lng2 .lt. 0.) lng2 = lng2 + 360.
      if (lng2 .gt .360.) lng2 = lng2 - 360.
c
      dst = sin(lat1*rad)*sin(lat2*rad)+cos(lat1*rad)*
     &      cos(lat2*rad)*cos((lng1-lng2)*rad)
      dst = min(dst, 1.)
      dst = max(dst, -1.)
      dst = (acos(dst)/rad)
      return
      end
#else
      subroutine rotation
      return
      end
#endif

