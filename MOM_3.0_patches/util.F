#ifdef test_util
      program driver
c
c=======================================================================
c     To run a test of the utilities in this util.F module, run the
c     "run_util" script
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
# include "stdunits.h"
      parameter (im=50, jm=30, km=25, imf=70, jmf=45, imc=20, jmc=12)
      parameter (lenw=9*imf)
      dimension a(im,jm), ia(im,jm), z(km), x(im), y(jm)
      dimension f(imf,jmf), xf(imf), yf(jmf)
      dimension c(imc,jmc), xc(imc), yc(jmc)
      dimension f1(imf,jmf), c1(imc,jmc)
      dimension mask(im,jm), scr1(im,jm), scr2(im,jm)
      dimension work1(lenw), work2(lenw)
c
      stdout = 6
c
c     define grid coordinates for a data grid
c
      dx = 360.0/im
      do i=1,im
        x(i) = (i-0.5)*dx
      enddo
      dy = 180.0/jm
      do j=1,jm
        y(j) = (j-0.5)*dy - 90.0
      enddo
c
c     define some data on the data grid
c
      pi    = 4.0*atan(1.0)
      twopi = 2.0*pi
      rim   = 1.0/im
      rjm   = 1.0/jm
      do j=1,jm
        do i=1,im
          a(i,j)  = 100.0*sin(twopi*i*rim)*cos(twopi*j*rjm)
          ia(i,j) = a(i,j)
        enddo
      enddo
c
c     define "fine" (high resolution) grid coordinates
c
      dx = 360.0/imf
      do i=1,imf
        xf(i) = (i-0.5)*dx
      enddo
      dy = 180.0/jmf
      do j=1,jmf
        yf(j) = (j-0.5)*dy - 90.0
      enddo
c
c     define "coarse" (low resolution) grid coordinates
c
      dx = 360.0/imc
      do i=1,imc
        xc(i) = (i-0.5)*dx
      enddo
      dy = 180.0/jmc
      do j=1,jmc
        yc(j) = (j-0.5)*dy - 90.0
      enddo
c
      write (stdout,'(/,/,a,/,/)')
     &' ============= Testing the util.F module ============='
c
c-----------------------------------------------------------------------
c     testing "checksum": printing a checksum of a real array
c-----------------------------------------------------------------------
c
      write (stdout,'(/,/,a,/)')
     & ' testing "checksum": printing a checksum of a real array'
      call print_checksum (a, im, jm, 'checksum for "a" = ')
c
c-----------------------------------------------------------------------
c     testing "matrix": printing a real array as a matrix
c-----------------------------------------------------------------------
c
      write (stdout,'(/,/,a,/)')
     & ' testing "matrix": printing a real array as a matrix'
      is = im/5
      ie = (3*im)/4
      js = jm/5
      je = (3*jm)/4
      scale = 0.0
      call matrix (a, im, is, ie, -js, -je, scale)
c
c-----------------------------------------------------------------------
c     testing "imatrx": formatted print of an integer field
c-----------------------------------------------------------------------
c
      write (stdout,'(/,/,a,/)')
     & ' testing "imatrx": formatted print of an integer field'
      is = im/5
      ie = (3*im)/4
      js = jm/5
      je = (3*jm)/4
      iform = 4
      call imatrx (ia, im, is, ie, js, je, iform)
c
c-----------------------------------------------------------------------
c     testing extrap:  extrapolation data into masked out areas.
c-----------------------------------------------------------------------
c
      write (stdout,'(/,/,a,/)') ' testing "extrap": filling land areas'
      is = im/3
      ie = im/2+2
      js = jm/3
      je = jm/2+2
c
c     mask = 0 for land, 1 for ocean
c
      land = 0
      do j=1,jm
        do i=1,im
          if (i .lt. is .or. i .gt. ie .or.
     &      j .lt. js .or. j .gt. je) then
            mask(i,j) = 1
          else
            mask(i,j) = land
            a(i,j) = 0
          endif
        enddo
      enddo
      is = is-2
      ie = ie+2
      js = js-2
      je = je+2
      scale = 1.0
c
c     show original data
c
      call matrix (a, im, is, ie, -js, -je, scale)
c
c     extrapolate original data into land areas
c
      call extrap (a, mask, scr1, scr2, im, jm, 100, 0.001, 'extrap'
     &, 1, land)
c
c     show extrapolated data
c
      call matrix (a, im, is, ie, -js, -je, scale)
c
c-----------------------------------------------------------------------
c     testing "ctf": coarse to fine grid interpolation. do it twice:
c     first to store the weights into the "work" area and secondly
c     to use the "work" area without recalculating the weights. verify
c     that the results are the same.
c-----------------------------------------------------------------------
c
      write (stdout,'(/,/,a,/)')
     & ' testing "ctf": interpolating coarse to fine resolution grids'
      n = 1
      call ctf (a, im, jm, x, y, f, imf, jmf, 5, imf-4, 5, jmf-4
     &,         xf, yf, n, work1, lenw)
      n = 2
      call ctf (a, im, jm, x, y, f1, imf, jmf, 5, imf-4, 5, jmf-4
     &,         xf, yf, n, work1, lenw)
c
      do j=1,jmf
        do i=1,imf
          if (f1(i,j) .ne. f(i,j)) then
            write (stdout,*) 'Error f1 <> f at i=',i,', j=',j
     &,                      ' in "ctf" test'
            call abort()
          endif 
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     testing "ftc": fine to coarse grid interpolation. do it twice:
c     first to store the weights into the "work" area and secondly
c     to use the "work" area without recalculating the weights
c-----------------------------------------------------------------------
c
      write (stdout,'(/,/,a,/)')
     & ' testing "ftc": interpolating fine to coarse resolution grids'
      n = 1
      call ftc (a, im, jm, x, y, c, imc, jmc, 1, imc, 1, jmc
     &,         xc, yc, n, work2, lenw)
      n = 2
      call ftc (a, im, jm, x, y, c1, imc, jmc, 1, imc, 1, jmc
     &,         xc, yc, n, work2, lenw)
c
      do j=1,jmc
        do i=1,imc
          if (c1(i,j) .ne. c(i,j)) then
            write (stdout,*) 'Error c1 <> c at i=',i,', j=',j
     &,                      ' in "ftc" test'
            call abort()
          endif 
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     testing "indp": index to nearest data point
c-----------------------------------------------------------------------
c
      write (stdout,'(/,/,a,/)')
     & ' testing "indp": locating index of nearest data point'
      do k=1,km
        z(k) = 0.1*k**2
      enddo
      value = 5.6
      kk = indp (value, z, km)
      do k=1,km
        write (stdout,*) ' k=',k,' z(k)=',z(k)
      enddo
      write (stdout,*) value,' is closest to z(',kk,') =', z(kk) 
c
c-----------------------------------------------------------------------
c     testing "scope": finding min, max, and avg of an array
c-----------------------------------------------------------------------
c
      write (stdout,'(/,/,a,/)')
     & ' testing "scope": interrogating an array for min, max, and avg'
      call scope (a, im, im, jm, ' testing')
c
c-----------------------------------------------------------------------
c     testing "plot": simple contour of a field using characters
c-----------------------------------------------------------------------
c
      write (stdout,'(/,/,a,/)')
     & ' testing "plot": simple contour of a field using characters'
      is = 10
      ie = 19
      js = 15
      je = 25
      nbins = 50
      zmin  = 0.0
      zmax  = 0.0
      call plot (a, im, 1, im, 1, jm, zmin, zmax, nbins, 'full array')
c
c-----------------------------------------------------------------------
c     testing "iplot": contour of an integer field using characters
c-----------------------------------------------------------------------
c
      write (stdout,'(/,/,a,/)')
     & ' testing "iplot": contour of integer field using characters'
      call iplot (ia(1,1), im, im, jm)
c
      call abort()
      end     
#endif


      function indp (value, array, ia)
c
c=======================================================================
c
c     indp = index of nearest data point within "array" corresponding to
c            "value".
c
c     inputs:
c
c     value  = arbitrary data...same units as elements in "array"
c     array  = array of data points  (must be monotonically increasing)
c     ia     = dimension of "array"
c
c     output:
c
c     indp =  index of nearest data point to "value"
c             if "value" is outside the domain of "array" then indp = 1
c             or "ia" depending on whether array(1) or array(ia) is
c             closest to "value"
c
c             note: if "array" is dimensioned array(0:ia) in the calling
c                   program, then the returned index should be reduced
c                   by one to account for the zero base.
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c     example:
c
c     let model depths be defined by the following:
c     parameter (km=5)
c     dimension z(km)
c     data z /5.0, 10.0, 50.0, 100.0, 250.0/
c
c     k1 = indp (12.5, z, km)
c     k2 = indp (0.0, z, km)
c
c     k1 would be set to 2, & k2 would be set to 1 so that
c     z(k1) would be the nearest data point to 12.5 and z(k2) would
c     be the nearest data point to 0.0
c
c=======================================================================
c
#include "stdunits.h"
      dimension array(ia)
c
      logical keep_going
c
      do i=2,ia
        if (array(i) .lt. array(i-1)) then
          write (stdout,*)
     &   ' => Error: array must be monotonically increasing in "indp"' 
     &,  '           when searching for nearest element to value=',value
          write (stdout,*) '           array(i) < array(i-1) for i=',i 
          write (stdout,*) '           array(i) for i=1..ia follows:'
          do ii=1,ia
            write (stdout,*) 'i=',ii, ' array(i)=',array(ii)
          enddo
          call abort()
        endif
      enddo
      if (value .lt. array(1) .or. value .gt. array(ia)) then
        if (value .lt. array(1))  indp = 1
        if (value .gt. array(ia)) indp = ia
      else
        i=1
        keep_going = .true.
        do while (i .le. ia .and. keep_going)
          i = i+1
          if (value .le. array(i)) then
            indp = i
            if (array(i)-value .gt. value-array(i-1)) indp = i-1
            keep_going = .false.
          endif
        enddo
      endif
      return
      end

      subroutine ftc (f, if1, jf, xf, yf, c, ic, jc, istart, iend
     &,               jstart, jend, xc, yc, init, work, lenw)
c
c=======================================================================
c
c     "ftc" is a acronym for "fine to coarse".
c
c     obtain a coarse grid representation of a fine grid dataset by area
c     averaging grid boxe values on the fine grid which overlay coarse 
c     grid boxes. note: the coarse grid boxes do not have to contain an
c     integral number of fine grid boxes.
c
c
c     inputs:
c
c     f      = data on fine grid
c     if1     = inner dimension of "f"
c     jf     = outer dimension of "f"
c     xf     = coordinates for inner dimension of "f" (eg: longitudes)
c     yf     = coordinates for outer dimension of "f" (eg: latitudes)
c
c     ic     = inner dimension of coarse grid "c"
c     jc     = outer dimension of coarse grid "c"
c     istart = starting index along inner dimension of "c" for which
c              averaged values are desired
c     iend   = ending index along inner dimension of "c" for which
c              averaged values are desired
c     jstart = starting index along outer dimension of "c" for which
c              averaged values are desired
c     jend   = ending index along outer dimension of "c" for which
c              averaged values are desired
c     xc     = coordinates for inner dimension of "c" (eg: longitudes)
c     yc     = coordinates for outer dimension of "c" (eg: latitudes)
c     init   = initialize the averaging factors
c              "init" should be set = 1 on the first call. 
c              "init" <> 1 uses the previously computed factors stored
c              in "work" array.
c     work   = work array of averaging factors when "init" <> 1
c              (previously calculated by "ftc" when "init" = 1)
c     lenw   = size of work array. lenw should be >= 9*max(if1,jf)
c       
c     output:
c
c     c  = coarse grid average of "f" defined over
c          ((c(i,j),i=istart,iend),j=jstart,jend)
c     work   = work array of averaging factors when "init" = 1
c
c     restrictions:
c
c     fine and coarse grids are assumed rectangular with "xf" and "xc"
c     having the same units. "yf" and "yc" must also have the same units
c     the coarse domain xc(istart)...xc(iend) must be within 
c     the fine domain xf(1)...xf(if1). similarly,
c     yc(jstart)...yc(jend) must be within yf(1)...yf(jf). all 
c     coordinates must be strictly monotonically increasing.
c 
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
      logical error, show_coord
      parameter (len=10000, p5=0.5, c0=0.0)
      dimension iso(0:len), ieo(0:len), jso(0:len), jeo(0:len)
     &,         dx(0:len,2), dy(0:len,2), edgecx(0:len), edgecy(0:len)
     &,         edgefx(0:len), edgefy(0:len)
      dimension f(if1,jf), xf(if1), yf(jf)
      dimension c(ic,jc), xc(ic), yc(jc)
      dimension work(lenw)
c
        write (stdout,*) ' '
        write (stdout,*)
     & '            Averaging data from "fine" to "coarse" grid'
c
c-----------------------------------------------------------------------
c     initialize weights or use previously calculated weights
c-----------------------------------------------------------------------
c
      if (init .eq. 1) then
        error = .false.
        write (stdout,*)
     & '              (initializing the averaging weights)'
        write (stdout,*) ' '
c
c       test to verify that array sizes do not exceed limits
c
        if (if1 .gt. len .or. jf .gt. len) then
          i = max(if1,jf)
          write (stdout,*) '=>Error: increase "len" in "ftc" to ',i
          call abort()
        endif
        if (lenw .lt. 9*max(if1,jf)) then
          write (stdout,*) '=>Error: increase size of "work" array',
     &      ' to at least ',9*max(if1,jf),' for calls to "ftc"'
          error = .true.
        endif
c
c       verify that the "coarse" grid lies within the "fine" grid
c
        if (xf(1) .gt. xc(istart) .or. xf(if1) .lt. xc(iend)) then
          write (stdout,*)
     &     '=>Error: Coarse grid "xc" is outside "fine" grid "xf".'
          if (xf(1) .gt. xc(istart)) then
            write (stdout,*) '  xc(',istart,')  .lt.  xf(1)'  
          endif
          if (xf(if1) .lt. xc(iend)) then
            write (stdout,*) '  xc(',iend,')  .gt.  xf(',if1,')'  
          endif
	  stop 'in ftc'
        endif 
        if (yf(1) .gt. yc(jstart) .or. yf(jf) .lt. yc(jend)) then
          write (stdout,*)
     &     '=>Error: Coarse grid "yc" is outside "fine" grid "yf".'
          if (yf(1) .gt. yc(jstart)) then
            write (stdout,*) '  yc(',jstart,') .lt. yf(1)' 
          endif
          if (yf(jf) .lt. yc(jend)) then
            write (stdout,*) '  yc(',jend,')  .gt.  yf(',jf,')'  
          endif
	  stop 'in ftc'
        endif
c
c       construct edges of "coarse" grid boxes
c
        do i=1,ic-1
          edgecx(i) = p5*(xc(i) + xc(i+1))
        enddo
        edgecx(0)  = xc(1) - (edgecx(1) - xc(1))
        edgecx(ic) = xc(ic) + (xc(ic) - edgecx(ic-1))
c
        do j=1,jc-1
          edgecy(j) = p5*(yc(j) + yc(j+1))
        enddo
        edgecy(0)  = yc(1) - (edgecy(1) - yc(1))
        edgecy(jc) = yc(jc) + (yc(jc) - edgecy(jc-1))
c
c       construct edges of "fine" grid boxes
c
        do i=1,if1-1
          edgefx(i) = p5*(xf(i) + xf(i+1))
        enddo
        edgefx(0)  = xf(1) - (edgefx(1) - xf(1))
        edgefx(if1) = xf(if1) + (xf(if1) - edgefx(if1-1))
c
        do j=1,jf-1
          edgefy(j) = p5*(yf(j) + yf(j+1))
        enddo
        edgefy(0)  = yf(1) - (edgefy(1) - yf(1))
        edgefy(jf) = yf(jf) + (yf(jf) - edgefy(jf-1))
c
c       calculate "dx" and "dy" for the "fine" grid boxes
c
        do i=1,if1
          dx(i,1) = edgefx(i) - edgefx(i-1)       
          dx(i,2) = dx(i,1)
        enddo
        dx(0,1) = dx(1,1)
        dx(0,2) = dx(1,2)
c
        do j=1,jf
          dy(j,1) = edgefy(j) - edgefy(j-1)
          dy(j,2) = dy(j,1)
        enddo
        dy(0,1) = dy(1,1)
        dy(0,2) = dy(1,2)
c
c       modify "dx" and "dy" for possibly partial "fine" grid boxes 
c       near the edges of each coarse grid box.
c       "ii" is the index of the fine grid box which contains the
c       eastern edge of coarse grid box with index "i".
c       dx(ii,1) is the portion of the fine grid box to the west of the
c       edge and dx(ii,2) is the portion to the east. similarly,
c       dy(jj,1) is to the south and dy(jj,2) is to the north of the
c       northern edge of coarse box with index "j".
c        
c       note: edgefx and edgefy are zero based and need the -1 when
c       using "indp" 
c
        do i=0,ic
          ii   = indp (edgecx(i), edgefx, if1+1) - 1
          frac = abs(edgecx(i) - edgefx(ii))
          if (edgefx(ii) .lt. edgecx(i)) then
            ii = ii + 1
            dx(ii,2) = (edgefx(min(if1,ii)) - edgefx(ii-1)) - frac
            dx(ii,1) = frac
          else
            dx(ii,2) = frac
            dx(ii,1) = (edgefx(ii) - edgefx(max(ii-1,0))) - frac
          endif
          ieo(i) = min(if1,max(1,ii))
        enddo
        do i=1,ic
          iso(i) = max(1,ieo(i-1))
        enddo
        iso(0) = ieo(0)
c
        do j=0,jc
          jj   = indp (edgecy(j), edgefy, jf+1) - 1
          frac = abs(edgecy(j) - edgefy(jj))
          if (edgefy(jj) .lt. edgecy(j)) then
            jj = jj + 1
            dy(jj,2) = (edgefy(min(jf,jj)) - edgefy(jj-1)) - frac
            dy(jj,1) = frac
          else
            dy(jj,2) = frac
            dy(jj,1) = (edgefy(jj) - edgefy(max(jj-1,0))) - frac
          endif
          jeo(j) = min(jf,max(1,jj))
        enddo
        do j=1,jc
          jso(j) = max(1,jeo(j-1))
        enddo
        jso(0) = jeo(0)
c
c       store the weights into the "work" array
c
        indx = 1
        do j=0,jc
          work(indx)   = jso(j)
          work(indx+1) = jeo(j)
          indx         = indx + 2
        enddo
c
        do i=0,ic
          work(indx)   = iso(i)
          work(indx+1) = ieo(i)
          indx         = indx + 2
        enddo
c
        do j=0,jf
          work(indx)   = dy(j,1)
          work(indx+1) = dy(j,2)
          indx         = indx + 2
        enddo
c
        do i=0,if1
          work(indx)   = dx(i,1)
          work(indx+1) = dx(i,2)
          indx         = indx + 2
        enddo
c
c       verify that coarse grid is coarser than the fine grid
c
        do j=jstart,jend
          if ((jso(j) .eq. jso(j+1)) .and. (yc(j) .ge. yf(1))
     &      .and. (yc(j) .le. yf(jf))) then
            write (stdout,*)
     &         '=>Warning: "Coarse" grid is finer than "fine" grid'
     &,     ' near yf(',jso(j),') =',yf(jso(j))
     &,     ' (average may not be accurate)'
          endif
        enddo
c
        do i=istart,iend
          if ((iso(i) .eq. iso(i+1)) .and. (xc(i) .ge. xf(1))
     &      .and. (xc(i) .le. xf(if1))) then
            write (stdout,*)
     &      '=>Warning: "Coarse" grid is finer than "fine" grid'
     &,     ' near xf(',iso(i),') = ',xf(iso(i))
     &,     ' (average may not be accurate)'
          endif
        enddo
#ifdef skip_interp_details
        show_coord = .false.
        write (stdout,'(/a/)')
     &   'Remove -Dskip_interp_details to show the indices'
#else
        show_coord = .true.
        write (stdout,'(/a/)')
     &   'Use -Dskip_interp_details to not show the following'
#endif
        if (error .or. show_coord) then
          write (stdout,*)
     & ' Indices for averaging fine grid to coarse grid:'
          write (stdout,*)
     & ' (fractional grid boxes are accounted for)'
          write (stdout,8700)
          write (stdout,9000) (m,iso(m),ieo(m),m=istart,iend)
          write (stdout,*) ' '
          write (stdout,*) ' Coordinates for coarse grid points "xc" ='
          write (stdout,8500) xc
          write (stdout,*) ' Coordinates for fine grid points "xf" ='
          write (stdout,8500) xf
c
          write (stdout,8800)
          write (stdout,9000) (m,jso(m),jeo(m),m=jstart,jend)
          write (stdout,*) ' '
          write (stdout,*) ' Coordinates for coarse grid points "yc" ='
          write (stdout,8500) yc
          write (stdout,*) ' Coordinates for fine grid points "yf" ='
          write (stdout,8500) yf
        endif
        if (error) call abort()
      else
        write (stdout,*)
     & '              (using previously initialized averaging weights)'
        write (stdout,*) ' '
c
c       extract the weights from the "work" array
c
        indx = 1
        do j=0,jc
          jso(j) = nint(work(indx))   
          jeo(j) = nint(work(indx+1))
          indx   = indx + 2
        enddo
c
        do i=0,ic
          iso(i) = nint(work(indx))
          ieo(i) = nint(work(indx+1))
          indx   = indx + 2
        enddo
c
        do j=0,jf
          dy(j,1) = work(indx)
          dy(j,2) = work(indx+1)
          indx    = indx + 2
        enddo
c
        do i=0,if1
          dx(i,1) = work(indx)
          dx(i,2) = work(indx+1)
          indx    = indx + 2
        enddo
      endif
c
c-----------------------------------------------------------------------
c     average the "fine" grid to the "coarse" grid
c-----------------------------------------------------------------------
c
      do m=jstart,jend
        do i=istart,iend
          weight = c0
          sum    = c0
          do j=jso(m),jeo(m)
            indy = 2
            if (j .eq. jeo(m)) indy = 1
            wty = dy(j,indy)
            do ii=iso(i),ieo(i)
              indx = 2
              if (ii .eq. ieo(i)) indx = 1
              area   = dx(ii,indx)*wty
              weight = weight + area
              sum    = sum + f(ii,j)*area
            enddo
          enddo
          c(i,m) = sum/weight
        enddo
      enddo
      return
8500  format (1x,10g11.4)
8700  format (/' Along the 1st dimension, the form is (Coarse grid'
     &,' point "xc": range of fine grid points "xf" to average)'/)
8800  format (/' Along the 2nd dimension, the form is (Coarse grid'
     &,' point "yc": range of fine grid points "yf" to average)'/)
9000  format (5(1x,'(',i4,': ',i4,' to ',i4,')'),/)
      end


      subroutine ctf (c, ic, jc, xc, yc, f, if1, jf, istart, iend
     &,               jstart, jend, xf, yf, init, work, lenw)
c
c=======================================================================
c
c     "ctf" is a acronym for "coarse to fine".
c
c     obtain a fine grid representation of a coarse grid dataset by 
c     linear interpolation of grid box values on the coarse grid to grid
c     boxes on the fine grid.
c
c
c     inputs:
c
c     c      = coarse grid data
c     ic     = inner dimension of coarse grid "c"
c     jc     = outer dimension of coarse grid "c"
c     xc     = coordinates for inner dimension of "c" (eg: longitudes)
c     yc     = coordinates for outer dimension of "c" (eg: latitudes)
c
c     if1     = inner dimension of "f"
c     jf     = outer dimension of "f"
c     xf     = coordinates for inner dimension of "f" (eg: longitudes)
c     yf     = coordinates for outer dimension of "f" (eg: latitudes)
c
c     istart = starting index along inner dimension of "f" for which
c              interpolated values are desired
c     iend   = ending index along inner dimension of "f" for which
c              interpolated values are desired
c     jstart = starting index along outer dimension of "f" for which
c              interpolated values are desired
c     jend   = ending index along outer dimension of "f" for which
c              interpolated values are desired
c     init   = initialize the interpolation factors
c              "init" should be set = 1 on the first call. 
c              "init" <> 1 uses the previously computed factors stored
c              in "work" array.
c     work   = work array of interpolation factors when "init" <> 1
c              (previously calculated by "ctf" when "init" = 1)
c     lenw   = size of work array. lenw should be >= 8*max(if1,jf)
c       
c     output:
c
c     f      = interpolated data on fine grid defined over
c              ((f(i,j),i=istart,iend),j=jstart,jend)
c     work   = work array of interpolation factors when "init" = 1
c
c     restrictions:
c
c     fine and coarse grids are assumed rectangular with "xf" and "xc"
c     having the same units. "yf" and "yc" must also have the same units
c     the fine domain xf(istart)...xf(iend) must be within 
c     the coarse domain xc(1)...xc(ic). Similarly,   
c     and yc(js)...yc(je) must be within yf(1)...yf(jf). all coordinates
c     must be strictly monotonically increasing.
c 
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
      logical error, show_coord
      parameter (len=10000, p5=0.5, c0=0.0)
      dimension indxi(len), indxj(len), dnorth(len), dsouth(len)
     &,         deast(len), dwest(len), width(len), height(len)
      dimension f(if1,jf), xf(if1), yf(jf)
      dimension c(ic,jc), xc(ic), yc(jc)
      dimension work(lenw)
c
        write (stdout,*) ' '
        write (stdout,*) 
     & '            Interpolating data from "coarse" to "fine" grid'
c
c-----------------------------------------------------------------------
c     initialize weights or use previously calculated weights
c-----------------------------------------------------------------------
c
      if (init .eq. 1) then
        error = .false.
        write (stdout,*) 
     & '              (initializing interpolation weights)'
        write (stdout,*) ' '
c
c       test to verify that array sizes do not exceed limits
c
        if (if1 .gt. len .or. jf .gt. len) then
          i = max(if1,jf)
          write (stdout,*) '=>Error: increase "len" in "ctf" to ',i
          error = .true.
        endif
        if (lenw .lt. 8*max(if1,jf)) then
          write (stdout,*) '=>Error: increase size of "work" array',
     &      ' to at least ',8*max(if1,jf),' for calls to "ctf"'
          error = .true.
        endif
c
c       verify that the "fine" grid lies within the "coarse" grid
c
        epsilon = 1.e-5
        xcminus = xc(1) - epsilon*(xc(2)-xc(1))
        xcplus = xc(ic) + epsilon*(xc(ic)-xc(ic-1))
        if (xf(istart) .lt. xcminus .or. xf(iend) .gt. xcplus) then
          error = .true.
          write (stdout,*)
     &    '=>Warning: "fine" grid outside "coarse" grid in "ctf".'
          if (xf(istart) .lt. xc(1)) then
            write (stdout,*) '  xf(',istart,')  .lt.  xc(1)'  
          endif
          if (xc(ic) .lt. xf(iend)) then
            write (stdout,*) '  xf(',iend,')  .gt.  xc(',ic,')'  
          endif
        endif 
        ycminus = yc(1) - epsilon*(yc(2)-yc(1))
        ycplus = yc(jc) + epsilon*(yc(jc)-yc(jc-1))
        if (ycminus .gt. yf(jstart) .or. ycplus .lt. yf(jend)) then
          error = .true.
          write (stdout,*)
     &    '=>Warning: "fine" grid outside "coarse" grid in "ctf".'
          if (yc(1) .gt. yf(jstart)) then
            write (stdout,*) '  yf(',jstart,') .lt. yc(1)' 
          endif
          if (yc(jc) .lt. yf(jend)) then
            write (stdout,*) '  yf(',jend,')  .gt.  yc(',jc,')'  
          endif
        endif
c
c       find interpolation factors 
c
        indx = 1
        do j=jstart,jend
          jj = indp (yf(j), yc, jc)
          if (yc(jj) .gt. yf(j) .or. jj .eq. jc) jj = jj - 1
          indxj(j) = jj
          dnorth(j) = yc(jj+1) - yf(j)
          dsouth(j) = yf(j) - yc(jj)
          height(j) = yc(jj+1) - yc(jj)
c
c         store into "work" array for future use (when "init" <> 1)
c
          work(indx)   = indxj(j)
          work(indx+1) = dnorth(j)
          work(indx+2) = dsouth(j)
          work(indx+3) = height(j)
          indx         = indx + 4
        enddo
c
        do i=istart,iend
          ii = indp (xf(i), xc, ic)
          if (xc(ii) .gt. xf(i) .or. ii .eq. ic) ii = ii - 1
          indxi(i) = ii
          deast(i) = xc(ii+1) - xf(i)
          dwest(i) = xf(i) - xc(ii)
          width(i) = xc(ii+1) - xc(ii)
c
c         store into "work" array for future use (when "init" <> 1)
c
          work(indx)   = indxi(i)
          work(indx+1) = deast(i)
          work(indx+2) = dwest(i)
          work(indx+3) = width(i)
          indx         = indx + 4
        enddo
#ifdef skip_interp_details
        show_coord = .false.
        write (stdout,'(/a/)')
     &    'Remove -Dskip_interp_details to show the coordinates'
#else
        show_coord = .true.
        write (stdout,'(/a/)')
     &   'Use -Dskip_interp_details to not show the following'
#endif
        if (error .or. show_coord) then
          write (stdout,*) ' '
          write (stdout,*) ' Coordinates for coarse grid points "xc" ='
          write (stdout,8500) xc
          write (stdout,*) ' Coordinates for fine grid points "xf" ='
          write (stdout,8500) xf

          write (stdout,*) ' Coordinates for coarse grid points "yc" ='
          write (stdout,8500) yc
          write (stdout,*) ' Coordinates for fine grid points "yf" ='
          write (stdout,8500) yf
        endif
        if (error) call abort()
      else
c
c       extract previously calculated interpolation weights from "work"
c
        write (stdout,*) 
     &'      (using previously initialized interpolation weights)'
c
        indx = 1
        do j=jstart,jend
          indxj(j)  = nint(work(indx))
          dnorth(j) = work(indx+1)
          dsouth(j) = work(indx+2)
          height(j) = work(indx+3)
          indx      = indx + 4
        enddo
c
        do i=istart,iend 
          indxi(i) = nint(work(indx))
          deast(i) = work(indx+1)
          dwest(i) = work(indx+2)
          width(i) = work(indx+3)
          indx     = indx + 4
        enddo
      endif
c
c-----------------------------------------------------------------------
c     interpolate data from "coarse" to "fine" grid
c-----------------------------------------------------------------------
c
      do jj=jstart,jend
        j = indxj(jj)
        do ii=istart,iend
          i = indxi(ii)
          f(ii,jj) = (c(i,j)    *deast(ii)*dnorth(jj)
     &              + c(i+1,j)  *dwest(ii)*dnorth(jj)
     &              + c(i,j+1)  *deast(ii)*dsouth(jj)
     &              + c(i+1,j+1)*dwest(ii)*dsouth(jj)) / 
     &               (width(ii)*height(jj))
        enddo
      enddo
c
      return
8500  format (1x,10g11.3)
      end


      subroutine extrap (a, land, sor, res, il, jl, maxscn, crit, text
     &,                  gtype, fill)
c
c=======================================================================
c
c     utility to extrapolate values into specified areas neglecting
c     non-uniformity or asymmetry in the grid by solving a simple
c     heat eqn: del**2(a) = 0 over cells to be "extrapolated" using 
c     values over cells which are not to be extrapolated as boundaries.
c     this alleviates the problem of mismatched land/sea areas due to
c     different geometries or resolutions when interpolating between
c     atmospheric and ocean model grids.
c     the intent is to force reasonable values into land areas
c     near coastlines. far from coasts, the extrapolations may not be
c     reasonable.
c
c     note: incomming values over "extrapolated areas" are used as  
c           an initial guess and need to be specified
c
c
c     inputs:
c
c     a       = array with cells to be extrapolated. cells to
c               be extrapolated need to be set with an initial guess.
c     land    = mask to control which cells are to be extrapolated.
c     fill    = cells with mask=fill will get extrapolated values
c               other cells remain unchanged.
c     il      = number of points along 1st dimension of arrays
c     jl      = number of points along 2nd dimension of arrays
c     maxscn  = maximum number of passes allowed in relaxation
c     crit    = criterion for ending relaxation before "maxscn" limit
c     text    = character string (up to 15 chars) to identify data
c     gtype   = grid type = (1,2) to identify (ocean, atmosphere) grid
c     sor     = scratch area
c     res     = scratch area
c
c
c     outputs:
c
c     a       = array with extrapolated values in land areas.
c               non land areas remain unchanged.
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      logical done
#include "stdunits.h"
      integer gtype, fill
      character*(*) text
      parameter (c0=0.0, p25=0.25)
      dimension a(il,jl), land(il,jl), res(il,jl), sor(il,jl)
c
c-----------------------------------------------------------------------
c
c     solve a simple poisson eqn by relaxation to extrapolate data into
c     land areas using values over non land areas as boundary values.
c
c     note: successive calls to extrap will require fewer scans because
c           the initial guess field over land areas gets better with
c           each call.
c-----------------------------------------------------------------------
c
c
c     check on the grid type: atmosphere or ocean
c
      if (gtype .ne. 1 .and. gtype .ne. 2) then
        write (stdout,98) gtype
        call abort()
      endif     
c
c-----------------------------------------------------------------------
c     set the relaxation coefficient to zero over cells to be filled
c     relc is somewhat arbitrary
c-----------------------------------------------------------------------
c
      relc = 0.6
      do j=1,jl
        do i=1,il
          if (land(i,j) .eq. fill) then
            sor(i,j) = relc
          else
            sor(i,j) = c0
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     iterate until errors are acceptable.
c-----------------------------------------------------------------------
c     
      n = 0
      done = .false.
      do while (.not. done .and. n .le. maxscn)
        resmax = c0
        done   = .true.
        n    = n + 1
        do j=2,jl-1
          do i=2,il-1
            res(i,j) = p25*(a(i-1,j) + a(i+1,j) + a(i,j-1) + a(i,j+1)) 
     &                 - a(i,j)
          enddo
        enddo
        do j=2,jl-1
          do i=2,il-1
            res(i,j) = res(i,j)*sor(i,j)
            a(i,j) = a(i,j) + res(i,j)
            absres = abs(res(i,j))
            if (absres .gt. crit) done = .false.
            resmax = max(absres,resmax)
          enddo
        enddo
c
c-----------------------------------------------------------------------
c       set conditions at edge of grid
c-----------------------------------------------------------------------
c
        if (gtype .eq. 1) then
c
c         use cyclic or no flux conditions on ocean grids
c
          do j=1,jl
#ifdef cyclic
            a(1,j)  = a(il-1,j)
            a(il,j) = a(2,j)
#else
            a(1,j)  = a(2,j)
            a(il,j) = a(il-1,j)
#endif
          enddo
        elseif (gtype .eq. 2) then
c
c         always put cyclic conditions on atmosphere grids
c
          do j=1,jl
            a(1,j)  = a(il-1,j)
            a(il,j) = a(2,j)
          enddo
        endif
c
c       no flux condition at northern and southern boundaries
c
        do i=1,il
          a(i,1)  = a(i,2)
          a(i,jl) = a(i,jl-1)
        enddo
c
      enddo
c
      write (stdout,99) text, n, resmax
99    format (1x,'==> Extrapolated ',a15,' into land using ',i4
     &,       ' scans.  max residual=', g14.7)
98    format (1x,'==> Error:   gtype =',i6,' in extrap')
      return
      end



      subroutine setbcx (a, imt, jmtorkm)
c
c=======================================================================
c     set zonal boundary condition on the first index of array "a" for
c     every second index. the first index corresponds to the "x"
c     or longitude direction. 
c
c     input:
c      a = array in need of setting the zonal b.c.
c     output
c      a = array with zonal b.c. set 
c=======================================================================
c
      dimension a(imt,jmtorkm)
      do k=1,jmtorkm
#ifdef cyclic
        a(1,k)   = a(imt-1,k)
        a(imt,k) = a(2,k)
#else
# ifdef obc_west
        a(1,k)   = a(2,k)
# else
        a(1,k)   = 0.0
# endif
# ifdef obc_east
        a(imt,k) = a(imt-1,k)
# else
        a(imt,k) = 0.0
# endif
#endif
      enddo
      return
      end


      subroutine iplot (iarray, im, il, jl)
c
c=======================================================================
c      map integer array "iarray" into characters for printing with
c      format (a1) to provide a contour map of the integer field.
c      note: max number of unique characters = 80
c
c     inputs:
c
c     iarray = integer array to be plotted
c     im     = inner dimension of "iarray"
c     il     = number of points along inner dimension to plot (along x)
c     jl     = number of points along outer dimension to plot (along y)
c
c     output: prints contour map of "iarray"
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
      dimension iarray(im,jl)
      character*80 levels
      character*80 lev1
      save levels
      write (stdout,*) ' '
c
c     set character markers
c
      lev1(1:51) = '.+*ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuv'
      levels = lev1(1:51)//'wxyz0123456789-=!@#$%<>[]{}()'
c
c     find range of integers
c
      maxint = iarray(1,1)
      minint = iarray(1,1)
      do j=1,jl
        do i=1,il
          maxint = max(maxint,iarray(i,j))
          minint = min(minint,iarray(i,j))
        enddo
      enddo
c
c     show mapping of integers into characters
c
      write (stdout,*) ' '
      write (stdout,*)
     & ' "iplot" mapping of integers to characters is as follows:'
      inc  = 3
      last = min(minint+80-1,maxint) 
      do i=minint,last,inc
        ii = i-minint+1
        if (i+inc .le. last) then
          jinc = inc
        else
          jinc = last-i+1
        endif
        write (stdout,'(6(1x,i6,a,a,3x))') 
     &  (j+minint-1,' is printed as ',levels(j:j),j=ii,ii+jinc-1)
      enddo
      write (stdout,*) ' '
c
      if (maxint - minint + 1 .gt. 80) then
        write (stdout,*)
     & ' => Note: there are ',maxint-minint+1,' integers in the field'
        write (stdout,*) '          "iplot" cannot uniquely assign '
     &,' more than 80 characters for plotting symbols.'
        write (stdout,*) '          therefore integers are'
     &, ' represented by cyclicly reusing the list of plotting symbols'
        write (stdout,*) '   '
      endif
c
c     print character representation of integers
c
      inc=124
      do l=0,il,inc
        incr = min(inc,il-l)
        write (stdout,8800) (l+i,i=1,incr,4)
        do  jj=1,jl
          j  = jl+1-jj
          write (stdout,8900) j, 
     &    (levels(mod(iarray(l+i,j)-minint+1-1,80)+1:
     &            mod(iarray(l+i,j)-minint+1-1,80)+1),i=1,incr) 
        enddo
      enddo
      return
8800  format (/, 2x, 31i4)
8900  format (1x,i3,1x, 124a1)
      end


      subroutine imatrx (iarray, im, istrt, iend, jstrt, jend, iform)
c
c=======================================================================
c      integer matrix print with various formats
c
c     inputs:
c
c     iarray = integer array to be printed
c     im     = the 1st dimension of array
c     istrt  = starting index along 1st dimension for plot
c     iend   = ending index along 1st dimension for plot
c     jstrt  = starting index along 2nd dimension for plot
c     jend   = ending index along 2nd dimension for plot
c              note: if jstrt and jend are negative then the vertical
c                   y-axis is inverted.
c     iform  = format designator: 1 => format (i1) 
c                                 2 => format (i2)
c                                 3 => format (i3)
c                                 4 => format (i4)
c
c
c     output: print integer "iarray" with user format control
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
      dimension iarray(im,1000)
c
c     choose the plotting domain
c
      js = min(abs(jstrt),abs(jend))
      je = max(abs(jstrt),abs(jend))
      is = min(abs(istrt),abs(iend))
      ie = max(abs(istrt),abs(iend))
      il = ie-is+1
c
      write (stdout,*) ' '
      if (iform .eq. 1) then
c
c       use I1 format to print the integer array
c
        inc=120
        do l=0,il,inc
          incr = min(inc,il-l)
          write (stdout,8800) (l+i+is-1,i=1,incr,4)
          do jj=js,je
            if (jstrt .lt. 0 .or. jend .lt. 0) then
              j = je - (jj-js)
            else
              j = jj
            endif
            write (stdout,8900) j, (iarray(l+i+is-1,j),i=1,incr)
          enddo
        enddo
      elseif (iform .eq. 2) then
c
c       use I2 format to print the integer array
c
        inc=60
        do l=0,il,inc
          incr = min(inc,il-l)
          write (stdout,9000) (l+i+is-1,i=1,incr,2)
          do jj=js,je
            if (jstrt .lt. 0 .or. jend .lt. 0) then
              j = je - (jj-js)
            else
              j = jj
            endif
            write (stdout,9200) j, (iarray(l+i+is-1,j),i=1,incr)
          enddo
        enddo
      elseif (iform .eq. 3) then
c
c       use I3 format to print the integer array
c
        inc=40
        do l=0,il,inc
          incr = min(inc,il-l)
          write (stdout,9400) (l+i+is-1,i=1,incr,2)
          do jj=js,je
            if (jstrt .lt. 0 .or. jend .lt. 0) then
              j = je - (jj-js)
            else
              j = jj
            endif
            write (stdout,9500) j, (iarray(l+i+is-1,j),i=1,incr)
          enddo
        enddo
      elseif (iform .eq. 4) then
c
c       use I4 format to print the integer array
c
        inc=30
        do l=0,il,inc
          incr = min(inc,il-l)
          write (stdout,9600) (l+i+is-1,i=1,incr,2)
          do jj=js,je
            if (jstrt .lt. 0 .or. jend .lt. 0) then
              j = je - (jj-js)
            else
              j = jj
            endif
            write (stdout,9700) j, (iarray(l+i+is-1,j),i=1,incr)
          enddo
        enddo
      endif
      return
8800  format (/, 2x, 30i4)
8900  format (1x,i3,1x, 120i1)
9000  format (/, 3x, 30i4)
9200  format (1x,i3,1x, 60i2)
9400  format (/,/,/,2x,20i6/)
9500  format (1x,i3,1x,40i3)
9600  format (/,/,/,1x,20i8/)
9700  format (1x,i3,1x,30i4)
      end


      subroutine matrix (array, irdim, istrt, im, jstrt, jm, scale)
c
c=======================================================================
c
c     matrix is a general two-dimensional array printing routine,
c     input:
c     array = the array to be printed
c     irdim = the 1st dimension of array
c     istrt = the 1st element of the 1st dimension to be printed
c     im    = the last element of the 1st dimension to be printed
c     jstrt = the 1st element of the 2nd dimension to be printed
c     jm    = the last element of the 2nd dimension to be printed
c             the 2nd dimension is printed in reverse order if both
c             jstrt & jm are negative
c     scale = a scaling factor by which array is divided before
c             printing.  (if this is zero, no scaling is done.)
c             if scale=0, 10 columns are printed across in e format
c             if scale>0, 20 columns are printed across in f format
c
c     output: print "array" as a matrix
c=======================================================================
c
#include "stdunits.h"
      parameter (c0=0.0, c1=1.0)
      dimension array(irdim,abs(jstrt):abs(jm))
c
      if (jstrt*jm .lt. 0) then
        write (stdout,999)  jstrt, jm
        call abort()
      endif
c
c     allow for inversion of 2nd dimension
c
      if (jm .lt. 0) then
        js   = -jm
        je   = -jstrt
        jinc = -1
      else
        js   = jstrt
        je   = jm
        jinc = 1
      endif
c
      if (scale .eq. c0) then
c
        do is=istrt,im,10
          ie = min(is + 9,im)
          write (stdout,9001) (i, i=is,ie)
          do l=js,je,jinc
            write (stdout,9002) l, (array(i,l),i=is,ie)
          enddo
          write (stdout,'(/)')
        enddo
      else
        scaler = c1/scale
        do is=istrt,im,20
          ie = min(is + 19,im)
          write (stdout,9003) (i, i=is,ie)
          do l=js,je,jinc
            write (stdout,9004) l, (array(i,l)*scaler,i=is,ie)
          enddo
          write (stdout,'(/)')
        enddo
      endif
      return
999   format (1x,'jstrt=',i5,' jm=',i5,' in matrix')
9001  format(10i13)
9002  format(1x,i2,10(1pe13.5))
9003  format(3x,20i6)
9004  format(1x,i3,1x,20f6.2)
      end


      subroutine scope (array, im, il, jl, text)
c
c=======================================================================
c     scope interrogates "array" for the min, max (with respective
c            locations), and simple unweighted average
c
c     inputs:
c     array = the array to be interrogated
c     im    = the inner dimension of "array"
c     il    = the number of points along the inner dimension to consider
c     jl    = the number of points along the inner dimension to consider
c     text  = descriptive text (up to 15 chars) to be printed
c
c     output: prints min, max (with locations), and average of "array"
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
# include "stdunits.h"
      character*(*) text
      dimension array(im,jl)
      umax  = array(1,1)
      umin  = array(1,1)
      iumax = 1
      jumax = 1
      iumin = 1
      jumin = 1
      sum    = 0.0
      do j=1,jl
        do i=1,il
          sum = sum + array(i,j)
          if (array(i,j) .gt. umax) then
            umax = array(i,j)
            iumax = i
            jumax = j
          endif
          if (array(i,j) .lt. umin) then
            umin = array(i,j)
            iumin = i
            jumin = j
          endif
        enddo
      enddo
      avg = sum/(il*jl)
      write (stdout,9200) text, il, jl, iumin, jumin, umin, iumax, jumax
     &,              umax, avg
      return
9200  format (1x,'Scope: ',a15,': il=',i4,', jl=',i4
     &,', min(',i4,',',i4,')=',g10.3,', max(',i4,',',i4,')=',g10.3
     &,', avg=',g10.3/)
      end


      subroutine sum1st (a, imt, jmt, text)
c
c=======================================================================
c     inputs:
c     a     = the array to be interrogated
c     imt   = the 1st dimension of "a"
c     imt   = the 2nd dimension of "a"
c     text  = descriptive text to be printed
c
c     output: sums the first index of array "a" for each 2nd index
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
# include "stdunits.h"
      character*(*) text
      dimension a(imt,jmt)
      write (stdout,*) ' '
      write (stdout,*) text
      big = abs(a(1,1))
      imax = 0
      jmax = 0
      do j=1,jmt
        sum = 0.0
        do i=1,imt
          sum = sum + a(i,j)
          if (abs(a(i,j)) .gt. big) then
            big = abs(a(i,j))
            imax = i
            jmax = j
          endif
        enddo
        if (sum .ne. 0.0) then
          write (stdout,'(a,i4,a,e14.7)')
     &    '2nd index=',j,'. sum over 1st index =',sum
        endif
      enddo
      if (imax .ne. 0 .and. jmax .ne. 0) then
        write (stdout,*) ' biggest=',a(imax,jmax),' at i=',imax
     &,                  ' j=',jmax
      else
        write (stdout,*) ' ---field is a constant =', a(1,1)
      endif
      return
      end


      subroutine plot (array, im, istrt, iend, jstrt,  jend 
     &,                  zmin, zmax, nbin, title)
c
c=======================================================================
c     "plot" contours "array" by dividing the array values into bins, 
c     assigning a character to each bin, and printing the characters to
c     produce a contour map.
c     The inner (1st) dimension is plotted horizontally (x-axis) and the
c     outer (2nd) dimension is plotted vertically (y-axis) and can be
c     inverted.
c
c     inputs:
c     array = the array to be contoured
c     im    = the inner dimension of "array"
c     istrt = starting index along 1st dimension for plot
c     iend  = ending index along 1st dimension for plot
c     jstrt = starting index along 2nd dimension for plot
c     jend  = ending index along 2nd dimension for plot
c             note: if jstrt and jend are negative then the vertical
c                   y-axis is inverted.
c     zmin  = the minimum value to plot
c     zmax  = the maximum value to plot
c             note: if zmin=zmax then then min and max of "array"
c                   is used.
c     nbin  = the number of bins between zmin and zmax
c     title = descriptive text string
c
c     output: contours "array"
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
# include "stdunits.h"
      parameter (maxbin=52, ncols=124)
      character*(*) title
      dimension array(im,1000), bin(maxbin,3)
      character*(maxbin) c, cc
      character*1 line(ncols)
      data c/'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'/
      save c
c
c     choose the plotting domain
c
      js = min(abs(jstrt),abs(jend))
      je = max(abs(jstrt),abs(jend))
      is = min(abs(istrt),abs(iend))
      ie = max(abs(istrt),abs(iend))
      il = ie-is+1
c
c     set the max and min
c
      if (zmax .eq. zmin) then
        big = array(is,js)
        sml = array(is,js)
        do j=js,je
          do i=is,ie
            if (array(i,j) .gt. big) big = array(i,j)
            if (array(i,j) .lt. sml) sml = array(i,j)
          enddo
        enddo
      else
        big = max(zmin,zmax)
        sml = min(zmin,zmax)
      endif
c
c     set up the discretization into bins
c
      write (stdout,*) ' '
      write (stdout,*) ' => contour bins for: ', title
      delta = (big-sml)/max(1,abs(nbin))
      limit = min(maxbin,max(1,abs(nbin)))
      do n=1,limit
        bin(n,1) = sml + (n-1)*delta
        bin(n,2) = sml + n*delta
        bin(n,3) = 0.5*(bin(n,1)+bin(n,2))
        cc(n:n)  = c(n:n)
        if (n .eq. 1) cc(n:n) = '-'
        if (n .ge. max(1,abs(nbin)) .or. n .eq. maxbin) cc(n:n) = '+'
        if (bin(n,1) .le. 0.0 .and. bin(n,2) .gt. 0.0) cc(n:n) = '.'
        if (n .lt. limit) then
          write (stdout,'(1x,a,a,g10.3,a,g10.3,a,g10.3,a)')  cc(n:n)
     &,   ' = ', bin(n,3) 
     &,   ', (', bin(n,1),' <= x < ', bin(n,2),')'
        else
          write (stdout,'(1x,a,a,g10.3,a,g10.3,a,g10.3,a)') cc(n:n)
     &,   ' = ', bin(n,3) 
     &,   ', (', bin(n,1),' <= x <=', bin(n,2),')'
        endif
      enddo
      if (abs(nbin) .gt. maxbin) then
        write (stdout,*)
     &    ' => Note: numbers larger than this are also plotted as +'
      endif
c
c     plot the character representation of the bins
c
      write (stdout,*) ' '
      do l=0,il,ncols
        incr = min(ncols,il-l)
        write (stdout,8800) (l+i+is-1,i=1,incr,4)
        do jj=js,je
          if (jstrt .lt. 0 .or. jend .lt. 0) then
            j = jj
          else
            j = je - (jj-js)
          endif
          do i=1,incr
            k = min(ifix((array(l+i+is-1,j) - sml)/delta) + 1, limit)
            line(i) = cc(k:k)
          enddo
          write (stdout,8900) j, (line(i),i=1,incr)
        enddo
      enddo
      return
8800  format (/, 2x, 31i4)
8900  format (1x,i3,1x, 124a1)
      end


      subroutine print_checksum (a, im, jm, text)
# include "stdunits.h"
      dimension a(im,jm)
      character*(*) text
      sum = checksum (a, im, jm)
      write (stdout,*) text, sum
      return
      end

      function checksum (a, im, jm)
      dimension a(im,jm)
      sumj = 0.0
      do j=1,jm
        sumi = 0.0
        do i=1,im
          sumi = sumi + abs(a(i,j))
        enddo
	sumj = sumj + sumi
      enddo
      checksum = sumj
      return
      end


      subroutine wrufio (iounit, array, len)
c
c=======================================================================
c     write unformatted fortran i/o
c
c     input:
c
c     iounit = fortran unit number
c     array  = array to be written to "iounit"
c     len    = length of array
c
c     output: writes to unit "iounit"
c=======================================================================
c
      dimension array(len)
      write (iounit) array
      return
      end


      subroutine rdufio (iounit, array, len)
c
c=======================================================================
c     read unformatted fortran i/o
c
c     input:
c
c     iounit = fortran unit number
c     array  = array to be read from "iounit"
c     len    = length of array
c
c     output: none
c=======================================================================
c
      dimension array(len)
      read (iounit) array
      return
      end

      subroutine tranlon (c, ic, il, jl, t, cx, fx, ifl, tx)
c
c-----------------------------------------------------------------------
c     grid interpolators in MOM require the grid coordinates to be
c     monotonic. when the prime meridian lies within the model grid,
c     global datasets (e.g. Scripps topography) must be translated in
c     longitude to remove the jump in
c     longitudes (358.5 359.5, 0.5, 1.5) across the meridian before
c     interpolating to the model grid. This is only of concern in
c     limited domain grids (e.g. atlantic basin) that contain the 
c     prime meridian.
c     
c     translate longitudes  "cx" to "tx" so that tx(i) i=1..ic
c     completely encloses model longitudes  fx(i) i=1..ifl
c     note that "tx" may extend beyond 360 degrees to contain "fx".
c     use same mapping to translate data in "c"
c
c     input:
c     c  = original data array
c     t  = temp array for translating data
c     cx = original data longitudes
c     tx = translated data longitudes
c     fx = model longitudes
c
c     output
c     c  = translated data array
c-----------------------------------------------------------------------
c
# include "stdunits.h"
      dimension c(ic,jl), t(ic), tx(ic), cx(ic), fx(ifl)
c
c-----------------------------------------------------------------------
c     find the index of the 1st model grid point on the data grid
c-----------------------------------------------------------------------
c
      iw = indp (fx(1), cx, ic)
      if (cx(iw) .gt. fx(1)) iw = max(1,iw-1)
c
c-----------------------------------------------------------------------
c     translate data longitudes so that tx(1) = cx(iw), tx(2) = cx(iw+1)
c-----------------------------------------------------------------------
c
      do i=1,ic
        tx(i) = cx(mod(i+iw-2,il) + 1)
        im1   = max(1,i-1)
        if (tx(i) .lt. tx(im1)) tx(i) = tx(i) + 360.0
      enddo
      if (fx(ifl) .gt. tx(ic)) then
        write (stdout,997) iw, ic, ifl
        write (stdout,998) 'tx= ',(tx(i),i=1,ic)
        write (stdout,998) 'fx= ',(fx(i),i=1,ifl)
        call abort()
      endif
c
c-----------------------------------------------------------------------
c     translate data to match translated longitudes
c-----------------------------------------------------------------------
c
      do j=1,jl
        do i=1,ic
          t(i) = c(mod(i+iw-2,il) + 1,j)
        enddo
        do i=1,ic
          c(i,j) = t(i)
        enddo
      enddo
c
      return
997   format (1x, ' ===>  tx(ic) < fx(ifl) in tranlon. iw=',i6,
     1         ' il=',i6,' ifl=',i6)
998   format (1x,a4,(5x,10e11.4))
      end



      subroutine fir2d (A,npass,nx,ny)
      dimension A(nx,ny), B(0:nx+1,0:ny+1)
      dimension fir(-1:1,-1:1)
c     
c     Two dimensional finite impulse response filter.
c     Will be applied npass number of times. 
c     No account taken for land, and does not take
c     account of the spherical lat/lon grid.
c     
      if(npass .gt. 0) then
        fir(-1,-1) = 1.0/16.0
        fir(0,-1)  = 1.0/8.0
        fir(1,-1)  = 1.0/16.0
        fir(-1,0)  = 1.0/8.0
        fir(0,0)   = 1.0/4.0
        fir(1,0)   = 1.0/8.0
        fir(-1,1)  = 1.0/16.0
        fir(0,1)   = 1.0/8.0
        fir(1,1)   = 1.0/16.0         
        do n = 1,npass,1
          do i=1,nx
            do j=1,ny
              B(i,j) = A(i,j)
            enddo
          enddo
          do i=1,nx
            B(i,0)    = B(i,1)
            B(i,ny+1) = B(i,ny)
          enddo
          do j=1,ny
            B(0,j)    = B(1,j)
            B(nx+1,j) = B(nx,j)
          enddo
          do i=1,nx
            do j=1,ny
              A(i,j) = 0.0
              do ip = -1,1
                do jq = -1,1 
                  A(i,j) = A(i,j) + fir(ip,jq)*B(i+ip,j+jq)
                enddo
              enddo           
            enddo
          enddo
        enddo
      endif
c     
      return
      end         

      
      subroutine edge_maker (it, edges, xt, xu, imt)
# include "stdunits.h"
c
c     make edges for grid cells
c
      dimension edges(0:imt), xt(imt), xu(imt)
      if (it .eq. 1) then
c
c       make edges for T cells
c
        edges(0) = xt(1) - (xu(1)-xt(1))
        do i=1,imt
	  edges(i) = xu(i)
        enddo
      elseif (it .eq. 2) then
c
c       make edges for U cells
c
        edges(imt) = xu(imt) + (xu(imt)-xt(imt))
        do i=1,imt
	  edges(i-1) = xt(i)
        enddo
      else
        write (stdout,*) 'Error:  it = ',it, ' in edge_maker'
	call abort()
      endif
      return
      end
   
