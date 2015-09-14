#ifdef drive_grids

# define driver_only
# include "util.F"
# include "iomngr.F"

      program driver
c
c=======================================================================
c
c              G R I D   G E N E R A T O R   M O D U L E
c
c     To experiment with designing grids for MOM: 
c     
c        1) inputs: specify the grid in the USER INPUT section below
c                          
c        2) compile and run this module by choosing the appropriate 
c           options in the "run_grids" script and executing the script.
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
# include "stdunits.h"
      parameter (maxlen=10000)
      dimension dxtdeg(maxlen), dytdeg(maxlen), dzt(maxlen)
      dimension dxudeg(maxlen), dyudeg(maxlen), dzw(0:maxlen)
      dimension xt(maxlen), xu(maxlen), yt(maxlen), yu(maxlen)
      dimension zt(maxlen), zw(maxlen)
c
      stdout = 6
c
# if defined read_my_grid || defined write_my_grid
      call ioinit
# endif
c
c-----------------------------------------------------------------------
c     calculate resolution and coordinates for T, U, W grid cells
c-----------------------------------------------------------------------
c
      write (stdout,'(//,39x,a,/)') 'Driving the grid module:'
c
      call gcoord (maxlen, imt, jmt, km, dxtdeg, dytdeg, dxudeg, dyudeg
     &,            dzt, dzw, xt, xu, yt, yu, zt, zw)
c
      if (imt .gt. maxlen .or. jmt .gt. maxlen .or. km .gt. maxlen) then
        write (stdout,*) ' => increase maxlen in the grids driver'
        call abort()
      endif
c
c-----------------------------------------------------------------------
c     prompt user to change imt, jmt, and km in "size.h"
c-----------------------------------------------------------------------
c
      write (stdout,'(//1x,a/)') '==>GRID INSTALLATION DIRECTIONS'
      write (stdout,'(6x,a,/,6x,a,i4,a,i4,a,i3,a,/)')
     & 'To use this grid in MOM, change parameter in "size.h" to:'
     &,'parameter (imt=',imt,', jmt=',jmt,', km=',km,')'
c
      write (stdout,'(/a/)')
     &' (To change resolution, read the manual)'
c
# ifdef write_my_grid
      write (stdout,'(//a/)')
     & ' =>Writing the grid definition to file "grid.dta.out" '
      call getunit (io, 'grid.dta.out'
     &, 'unformatted sequential rewind ieee')
      write (io) imt, jmt, km
      write (io) (dxtdeg(i),i=1,imt)
     &,         (dytdeg(j),j=1,jmt)
     &,         (dxudeg(i),i=1,imt)
     &,         (dyudeg(j),j=1,jmt)
     &,         (dzt(k),k=1,km)
     &,         (dzw(k),k=0,km)
     &,         (xt(i),i=1,imt)
     &,         (xu(i),i=1,imt)
     &,         (yt(j),j=1,jmt)
     &,         (yu(j),j=1,jmt)
     &,         (zt(k),k=1,km)
     &,         (zw(k),k=1,km)
      call relunit (io)
# endif
      end
#endif



      subroutine gcoord (maxlen, imt, jmt, km, dxtdeg, dytdeg, dxudeg
     &,                  dyudeg, dzt, dzw, xt, xu, yt, yu, zt, zw)
c
c=======================================================================
c
c                      G R I D    C O O R D I N A T E S
c
c     Construct grid point coordinates and resolution
c
c     input:
c
c     maxlen = maximum number of grid cells in latitude, longitude,
c              and depth
c
c              set grid specifications in USER INPUT section.
c
c     output:
c
c     imt    = number of longitudes
c     jmt    = number of latitudes
c     km     = number of depths
c     dxtdeg = width of T grid cells (degrees)
c     dytdeg = height of T grid cells (degrees)
c     dxudeg = width of U grid cells (degrees)
c     dyudeg = height of U grid cells (degrees)
c     dzt    = thickness of T grid cells (cm)
c     dzw    = thickness of W grid cells (cm)
c     xt     = longitude at centers of T grid cells (degrees)
c     xu     = longitude at centers of U grid cells (degrees)
c     yt     = latitude at centers of T grid cells (degrees)
c     yu     = latitude at centers of U grid cells (degrees)
c     zt     = depth at centers of T grid cells (centimeters)
c     zw     = depth at centers of U grid cells (centimeters)
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
      parameter (maxbounds=11)
      dimension xt(maxlen), yt(maxlen), xu(maxlen), yu(maxlen)
      dimension zw(maxlen), zt(maxlen)
c
      dimension dxtdeg(maxlen), dytdeg(maxlen), dzt(maxlen)
      dimension dxudeg(maxlen), dyudeg(maxlen), dzw(0:maxlen)
c
      dimension x_lon(maxbounds),   dx_lon(maxbounds)
      dimension y_lat(maxbounds),   dy_lat(maxbounds)
      dimension z_depth(maxbounds), dz_depth(maxbounds)
c
c-----------------------------------------------------------------------
c     USER INPUT ===>  read this section and specify the grid.
c-----------------------------------------------------------------------
c
c     How to define a grid resolution: Read the manual.
c
c            
c-----------------------------------------------------------------------
c     The following specifies the global test grid resolution: 
c-----------------------------------------------------------------------
c
c     constant 4.0 degree resolution in longitude
c
      parameter (nxlons=2)
      data (x_lon(i), i=1,nxlons) / 0.0, 360.0/
      data (dx_lon(i),i=1,nxlons) / 3.0,   3.0/        
#if !defined isotropic_grid
c
c     constant 3.0 degree resolution in latitude
c
      parameter (nylats=2)
      data (y_lat(j), j=1,nylats) / -90.0,    90.0/
      data (dy_lat(j),j=1,nylats) /   2.8125,    2.8125/        
#endif
#if defined isotropic_grid
c
c     define square grid cells
c
      parameter (nylats=2)
      data (y_lat(j), j=1,nylats) / -65.0,   65.0/
      data (dy_lat(j),j=1,nylats) /   3.0,    3.0/        
#endif
c
c     two vertical regions: constant 25m in the upper 100m, and variable
c     resolution in the lower region. Note that by setting the stretch
c     factor "stretch_z" > 1.0, the last region in the vertical can
c     have the grid cells stretched further. To see the effect of
c     stretching on vertical cell thickness and number of vertical
c     cells, try using "stretch_z"=1.1 and gradually increase it.  
c
      parameter (nzdepths=3)
      data (z_depth(k), k=1,nzdepths) /  0.0e2,   100.0e2,  5600.0e2/
      data (dz_depth(k),k=1,nzdepths) / 25.0e2,    25.0e2,   975.0e2/
      data stretch_z /1.0/
c
c     Some other examples:
c
c      4.0 degrees at the poles down to 2.0 deg at the equator
c
c      parameter (nylats=3)
c      data (y_lat(j), j=1,nylats) / -90.0,  0.0,  90.0/
c      data (dy_lat(j),j=1,nylats) /   4.0,  2.0,   4.0/        
c
c      set 1 deg at -30 deg lat down to 1/3 deg at -10 deg lat
c      constant 1/3 deg resolution from -10 deg to +10 deg lat and
c      1/3 deg at +10 deg lat to 1.0 deg at +30 deg lat 
c
c      parameter (nylats=4)
c      data (y_lat(j), j=1,nylats) / -30.0,  -10.0,    10.0,  30.0/
c      data (dy_lat(j),j=1,nylats) /   1.0,  0.3333,  0.3333,  1.0/
c
# ifdef obctest
c-----------------------------------------------------------------------
c     The following specifies the box for the open boundary test:
c     see Stevens (1990), Chap. 5: A Passive Open Boundary Condition
c-----------------------------------------------------------------------
c
c     constant 1.0 degree resolution
c
      parameter (nxlons=2)
      data (x_lon(i), i=1,nxlons) / -0.5,  18.5/
      data (dx_lon(i),i=1,nxlons) /  1.0,   1.0/
c
c     constant 0.5 degree resolution
c
      parameter (nylats=2)
# if defined obc_south || defined small_basin
      parameter (bound_south=56.75)
# else
      parameter (bound_south=52.75)
# endif
# if defined obc_north || defined small_basin
      parameter (bound_north=67.75)
# else
      parameter (bound_north=72.25)
# endif
      data (y_lat(j), j=1,nylats)  / bound_south,  bound_north/
      data (dy_lat(j),j=1,nylats) /   0.5,    0.5 /
c
      parameter (nzdepths=2)
      data (z_depth(k), k=1,nzdepths) /   0.0e2,  2700.0e2/  
      data (dz_depth(k),k=1,nzdepths) /  50.0e2,   850.0e2/ 
      data stretch_z /1.0/
c
# endif
# ifdef obctest2
c-----------------------------------------------------------------------
c     The following specifies the box for the second open boundary test:
c     see Stevens (1990), 6. An Active Open Boundary Condition
c-----------------------------------------------------------------------
c
c     constant 2.0 degree resolution
c
      parameter (nxlons=2)
      data (x_lon(i), i=1,nxlons) / -1.0,  37.0/
      data (dx_lon(i),i=1,nxlons) /  2.0,   2.0/
c
c     constant 2.0 degree resolution
c
      parameter (nylats=2)
# ifdef obc_south
      data (y_lat(j), j=1,nylats)  / 32.00,  54.00/
# else
      data (y_lat(j), j=1,nylats)  / 16.00,  54.00/
# endif
      data (dy_lat(j),j=1,nylats) /   2.0,    2.0 /
c
      parameter (nzdepths=2)
      data (z_depth(k), k=1,nzdepths) /   0.0e2,  3750.0e2/  
      data (dz_depth(k),k=1,nzdepths) /  50.0e2,  1200.0e2/ 
      data stretch_z /1.0/
c
# endif
c-----------------------------------------------------------------------
c     ==> end of USER INPUT
c-----------------------------------------------------------------------
c
c     set some constants
c
      p5 = 0.5
c
      write (stdout,'(/,a,/(/a))')
     &  'The following grid options have been selected:'
#ifdef generate_a_grid
     &,'-Dgenerate_a_grid'
#endif
#ifdef read_my_grid
     &,'-Dread_my_grid'
#endif
#ifdef write_my_grid
     &,'-Dwrite_my_grid'
#endif
#ifdef centered_t
     &,'-Dcentered_t (T points centered in T cells as in MOM 1)'
#endif
#ifdef isotropic_grid
     &,'-Disotropic_grid (dy = dx * cos(lat) )'
#endif
#ifdef extend_isotropic_grid
     &,'-Dextend_isotropic_grid (to poles)'
#endif
c
      ncase = 0
#ifdef read_my_grid
c
      ncase = ncase + 1
      write (stdout,'(//a/a/)')
     & ' =>Reading the grid definition from file "grid.dta" '
     &,'   (this assumes it was created elsewhere and is available)'
      call getunit (io, 'grid.dta'
     &, 'unformatted sequential rewind ieee')
      read (io) imt, jmt, km
      read (io) (dxtdeg(i),i=1,imt)
     &,         (dytdeg(j),j=1,jmt)
     &,         (dxudeg(i),i=1,imt)
     &,         (dyudeg(j),j=1,jmt)
     &,         (dzt(k),k=1,km)
     &,         (dzw(k),k=0,km)
     &,         (xt(i),i=1,imt)
     &,         (xu(i),i=1,imt)
     &,         (yt(j),j=1,jmt)
     &,         (yu(j),j=1,jmt)
     &,         (zt(k),k=1,km)
     &,         (zw(k),k=1,km)
      call relunit (io)
# ifdef centered_t
      write (stdout,*)
     &'=>Error: option -Dcentered_t not compatible with -Dread_my_grid'
     &,'         ...remove one of them.'
      call abort()
# endif
#endif
#ifdef generate_a_grid
      ncase = ncase + 1
c
      write (stdout,'(//,36x,a,/)') 'G R I D   G E N E R A T I O N'
      write (stdout,'(/,43x,a,/)') 'Grid resolution:'
c
      if (nxlons .gt. maxbounds .or. nylats .gt. maxbounds .or.
     &    nzdepths .gt. maxbounds) then
        write (stdout,'(/a/)')
     &  ' Increase parameter "maxbounds" to contain all regions.'
        call abort()
      endif
c
c-----------------------------------------------------------------------
c     Calculate resolution in longitude. Add one boundary cell at the
c     start and end of the domain so that calculations are meaningful
c     for grid cells i=2,imt-1
c-----------------------------------------------------------------------
c
      nbpts = 2
      write (stdout,'(/a,i1,a)')
     & ' Generating the longitudinal resolution ( ', nbpts
     &, ' extra boundary points will be used).'
c
c     convert negative longitudes (degrees west of Greenwich) to
c     positive longitudes
c
      offset = 0.0
      do n=1,nxlons
        if (x_lon(n) .lt. 0.0) then
          offset = 360.0
        endif
      enddo
      if (offset .gt. 0.0) then
          write (stdout,'(/,a/,a/,a)')
     &    '=>Warning: adding 360.0 degrees to all longitudes to convert'
     &,   '           negative values (degrees west of Greenwich) to'
     &,   '           positive values (degrees east of Greenwich)'
      endif
      do n=1,nxlons
        x_lon(n) = x_lon(n) + offset
      enddo
c
c
c     if a region contains the greenwich meridian, compensate by adding
c     360.0 degrees to all remaining regions
c
      do n=2,nxlons
        if (x_lon(n-1) .gt. x_lon(n)) then
          x_lon(n) = x_lon(n) + 360.0
          write (stdout,'(/,a,i1,a,a/)')
     &    ' Warning: adding 360.0 degrees to x_lon(',n,')'
     &,   ' to insure region boundaries increase monotonically'
        endif
      enddo
c
c     if width of domain exceeds 360.0 deg... limit width to 360.0
c
      if (x_lon(nxlons) - x_lon(1) .gt. 360.0) then
        write (stdout,'(/a,a,g14.7/)')
     &  '=>Warning: Domain width exceeds 360 deg. Restricting last'
     &, ' x_lon to ',x_lon(1)+360.0
          do m=1,nxlons
            write (stdout,'(i3,f10.5)') m, x_lon(m)
          end do
      endif
      if (x_lon(nxlons) - x_lon(1) .eq. 360.0) then
        if (dx_lon(nxlons) .ne. dx_lon(1)) then
          write (stdout,'(/a,a)')
     &    '=>Error: dx_lon(1) must equal dx_lon(last) when domain'
     &,   ' width = 360.0 degrees'
        endif
      endif

      stretch_x = 1.0
      call gcell (maxlen, nxlons, x_lon,  dx_lon, nbpts
     &,           imt, dxtdeg, dxudeg, stretch_x)
c
c-----------------------------------------------------------------------
c     Build the longitudinal grid points for a "B" grid
c     (account for an extra boundary point at the start)
c-----------------------------------------------------------------------
c
      xt(1) = x_lon(1) - p5*dx_lon(1)
      xu(1) = x_lon(1)
c
      do i=2,imt
        xu(i) = xu(i-1) +  dxtdeg(i)
        xt(i) = xt(i-1) +  dxudeg(i-1)
      enddo
c
c-----------------------------------------------------------------------
c     Calculate resolution in latitude. Add one boundary cell at the
c     start and end of the domain so that calculations are meaningful
c     for grid cells jrow=2,jmt-1
c-----------------------------------------------------------------------
c
      nbpts = 2
      write (stdout,'(/a,i1,a)')
     & ' Generating the latitudinal resolution ( ', nbpts
     &, ' extra boundary points will be used).'
c
c     if width of domain exceeds 180.0 deg... limit width to 180.0
c
      if (y_lat(nylats) - y_lat(1) .gt. 180.0) then
        write (stdout,'(/a,a,g14.7/)')
     &  '=>Warning: Latitudinal domain width exceeds 180 deg.'
     &, ' Restricting last y_lat to ',y_lat(1)+180.0
      endif
# if defined isotropic_grid
      call iso_grid (maxlen, nylats, y_lat, dy_lat, jmt, dytdeg, dyudeg
     &,              yt, yu)
# else
      do n=2,nylats
        if (y_lat(n-1) .gt. y_lat(n)) then
          write (stdout,'(/,a,/a/)')
     &    ' =>Error: latitude boundaries must increase monotonically'
     &,   '          check the specifications in the USER INPUT section'
          do m=1,nylats
            write (stdout,'(i3,f10.5)') m, y_lat(m)
          end do
          call abort()
        endif
      enddo
      stretch_y = 1.0
      call gcell (maxlen, nylats, y_lat, dy_lat, nbpts
     &,           jmt, dytdeg, dyudeg, stretch_y)
c
c-----------------------------------------------------------------------
c     Build the latitudinal grid points on a "B" grid
c     (account for an extra boundary point at the start)
c-----------------------------------------------------------------------
c
      yt(1) = y_lat(1) - p5*dy_lat(1)
      yu(1) = y_lat(1)
      do jrow=2,jmt
        yu(jrow) = yu(jrow-1) +  dytdeg(jrow)
        yt(jrow) = yt(jrow-1) +  dyudeg(jrow-1)
      enddo
# endif
c
c-----------------------------------------------------------------------
c     Calculate resolution in depth. No boundary cells added here so
c     calculations are meaningful for k=1,km. Allow the bottom region 
c     to be stretched further if desired.
c-----------------------------------------------------------------------
c
      nbpts = 0
      write (stdout,'(/a,i1,a)')
     & ' Generating the vertical resolution ( ', nbpts
     &, ' extra boundary points will be used).'
      do n=2,nzdepths
        if (z_depth(n-1) .gt. z_depth(n)) then
          write (stdout,'(/,a,/a/)')
     &    ' =>Error: depth boundaries must increase monotonically'
     &,   '          check the specifications in the USER INPUT section'
          call abort()
        endif
      enddo
      call gcell (maxlen, nzdepths, z_depth, dz_depth, nbpts
     &,           km, dzt, dzw(1), stretch_z)
c
c-----------------------------------------------------------------------
c     Build the vertical grid points on a "B" grid. The T and U
c     cells are staggered in the horizontal but at the same level in
c     the vertical. However, the W cells here refer to the vertical
c     advection velocities at the bottoms of the U and T cells.
c     (no extra boundary point at the start)
c-----------------------------------------------------------------------
c
      zt(1) = z_depth(1) + p5*dz_depth(1)
      zw(1) = z_depth(1) + dzt(1)
      do k=2,km
        zw(k) = zw(k-1) + dzt(k)
        zt(k) = zt(k-1) + dzw(k-1)
      enddo
c
c     set W cell thickness at surface and bottom to ocean part of cell
c
      dzw(0)  = zt(1)
      dzw(km) = zw(km) - zt(km)
#endif
c
#ifdef bbl_ag
c
c-----------------------------------------------------------------------
c     add an extra level for the bottom boundary layer
c-----------------------------------------------------------------------
c
      bbl_thickness = 50.0e2
      km            = km + 1
      zt(km)        = zw(km-1) + 0.5*bbl_thickness
      zw(km)        = zw(km-1) + bbl_thickness
      dzt(km)       = bbl_thickness
      dzw(km-1)     = zt(km) - zt(km-1)
      dzw(km)       = 0.5*bbl_thickness
      write (stdout,'(/a,i4,/)') 'Adding one level for BBL: new km=',km
#endif
c
c-----------------------------------------------------------------------
c     Print grid T cell resolution in longitude, latitude and depth
c     (also for W cells in depth)
c-----------------------------------------------------------------------
c
      write (stdout,9101) 
     &'Vertical resolution of cells "dzw(k)" k=0,', km, 'cm'
      write (stdout,9002) (dzw(k),k=0,km)
      write (stdout,9101) 
     &'Vertical resolution of cells "dzt(k)" k=1,', km, 'cm'
      write (stdout,9002) (dzt(k),k=1,km)
      write (stdout,9101) 
     &'Longitudinal resolution of cells "dxtdeg(i)" i=1,', imt, 'deg'
      write (stdout,9001) (dxtdeg(i),i=1,imt)
      write (stdout,9101) 
     &'Latitudinal resolution of cells "dytdeg(jrow)" jrow=1,',jmt,'deg'
      write (stdout,9001) (dytdeg(jrow),jrow=1,jmt)
      if (ncase .eq. 0) then
        write (stdout,'(/a/a/)')
     & '=>Error: One of the following options must be enabled:'
     &,'         generate_a_grid    read_my_grid'
        call abort()
      elseif (ncase .gt. 1) then
        write (stdout,'(/a/a/)')
     & '=>Error: Only one of the following options may be enabled:'
     &,'         generate_a_grid    read_my_grid'
        call abort()
      endif
c
c-----------------------------------------------------------------------
c     Check if the T grid resolution is an average of the
c     U cell resolution. This insures more accurate advection of
c     tracers within a stretched grid.
c-----------------------------------------------------------------------
c
      num  = 0
      tolr = 1.e-5
      write (stdout,'(/)')
      do i=2,imt-1
        dxubar = p5*(dxudeg(i) + dxudeg(i-1))
        if (abs(dxubar-dxtdeg(i)) .gt. tolr) then
          num = num + 1
          write (stdout,'(a,i5,a)')
     &    '=>Warning: T cell delta x at i=',i
     &,   ' is not an average of adjacent U cell delta x`s'     
        endif     
      enddo
c
      do jrow=2,jmt-1
        dyubar = p5*(dyudeg(jrow) + dyudeg(jrow-1))
        if (abs(dyubar-dytdeg(jrow)) .gt. tolr) then
          num = num + 1
          write (stdout,'(a,i5,a)')
     &    '=>Warning: T cell delta y at jrow=',jrow
     &,   ' is not an average of adjacent U cell delta y`s' 
        endif     
      enddo
c
      tolr = 1.e0
      do k=2,km-1
        dzwbar = p5*(dzw(k) + dzw(k-1))
        if (abs(dzwbar-dzt(k)) .gt. tolr) then
          num = num + 1
          write (stdout,'(a,i5,a)')
     &    '=>Warning: T cell delta z at k=',k
     &,   ' is not an average of adjacent W cell delta z`s'     
        endif     
      enddo
c
c-----------------------------------------------------------------------
c     Print all grid coordinates
c-----------------------------------------------------------------------
c
      zz = 0.0
      write (stdout
     &,'(//,40x,a,//,a,g14.7,a,/a/,a,g14.7,a/a,/,a,g14.7,a/)') 
     &  ' Grid Point Coordinate details: '
     &, ' The western edge of the T-cell at i=2 is at:'
     &,  xu(1),' (deg) longitude',' (the 1st T-cell is a boundary cell)'
     &, ' The southern edge of the T-cell at jrow=2 is at:'
     &,  yu(1),' (deg) latitude',' (the 1st T-cell is a boundary cell)'
     &,' The top edge of the T-cell at k=1 is at z =',zz,' (cm)'
      write (stdout,'(a,i4,a,g14.7,a/a,i4,a)') 
     & ' The eastern edge of the T-cell at i=',imt-1,' is at:'
     &,  xu(imt-1),' (deg) longitude'
     &,' (the last T-cell (imt=',imt,') is a boundary cell)'
     &, ' The northern edge of the T-cell at jrow=',jmt-1,' is at:'
     &,  yu(jmt-1),' (deg) latitude'
     &,' (the last T-cell (jmt=',jmt,') is a boundary cell)'
     &,' The bottom edge of T-cell at k=',km,' is at z ='
     &,  zw(km),' (cm)'
      write (stdout,'(/,a,g14.7,a/a/,a,g14.7,a/a/,a,g14.7,a/)') 
     &  ' The western edge of the U-cell at i=2 is at:', xt(1)
     &, ' (deg) longitude',' (the 1st U-cell is a boundary point)'
     &, ' The southern edge of the U cell at jrow=1 is at:', yt(1)
     &, ' (deg) latitude',' (the 1st U cell is a boundary point)'
     &, ' The top edge of the W-cell at k=1 is at z =',zt(1),' (cm)'
      write (stdout,'(a,i4,a,g14.7,a/a,i4,a)') 
     & ' The eastern edge of the U-cell at i=',imt-1,' is at:'
     &,  xt(imt),' (deg) longitude'
     &,' (the last U-cell (imt=',imt,') is a boundary cell)'
     &, ' The northern edge of the U-cell at jrow=',jmt-1,' is at:'
     &,  yt(jmt),' (deg) latitude'
     &,' (the last U-cell (jmt=',jmt,') is a boundary cell)'
     &,' The bottom edge of W-cell at k=',km,' is at z ='
     &,  zw(km),' (cm)'
      write (stdout,9103) km
      write (stdout,9002) (zt(k),k=1,km)
      write (stdout,9104) km
      write (stdout,9002) (zw(k),k=1,km)
      write (stdout,9105) jmt
      write (stdout,9001) (yt(jrow),jrow=1,jmt)
      write (stdout,9106) jmt
      write (stdout,9001) (yu(jrow),jrow=1,jmt)
      write (stdout,9107) imt
      write (stdout,9001) (xt(i),i=1,imt)
      write (stdout,9108) imt
      write (stdout,9001) (xu(i),i=1,imt)
c
#ifdef symmetry
c
c-----------------------------------------------------------------------
c     insure that yu(jmt-1) = 0.0 (equator) when using symmetry     
c-----------------------------------------------------------------------
c
      if (yu(jmt-1) .ne. 0.0) then
        write (stdout,*) '=> Error: yu(jmt-1) must = 0.0 for symmetry'
        call abort()
      endif
#endif
c
c---------------------------------------------------------------------
c     compute a grid checksum
c---------------------------------------------------------------------
c
      cksum = 0.0
      cksum = cksum + checksum (xt, imt, 1)
      cksum = cksum + checksum (yt, jmt, 1)
      cksum = cksum + checksum (zt, km, 1)
      cksum = cksum + checksum (xu, imt, 1)
      cksum = cksum + checksum (yu, jmt, 1)
      cksum = cksum + checksum (zw, km, 1)
      cksum = cksum + checksum (dxtdeg, imt, 1)
      cksum = cksum + checksum (dytdeg, jmt, 1)
      cksum = cksum + checksum (dxudeg, imt, 1)
      cksum = cksum + checksum (dyudeg, jmt, 1)
      cksum = cksum + checksum (dzt, km, 1)
      cksum = cksum + checksum (dzw, km+1, 1)
      write (stdout,'(/)')
      write (stdout,*) 'Grid checksum = ',cksum
      write (stdout,'(/)')
      return
9001  format (1x,10f10.4)
9002  format (1x,10f10.2)
9101  format (/,  a,i4,' in units of ',a,' as follows:')
9103  format (/,' Depth to T & U grid points (cm): zt(k) k=1,',i3)
9104  format (/,' Depth to W grid points (cm): zw(k) k=1,',i3)
9105  format (/,' Latitude of T points (deg): yt(j) j=1,',i4)
9106  format (/,' Latitude of U points (deg): yu(j) j=1,',i4)
9107  format (/,' Longitude of T points (deg): xt(i) i=1,',i4)
9108  format (/,' Longitude of U points (deg): xu(i) i=1,',i4)
      end





      subroutine gcell (maxlen, n_bounds, bounds, d_bounds, nbpts
     &,                 num, deltat, deltau, stretch)
c
c=======================================================================
c
c              G R I D   C E L L   C O N S T R U C T I O N
c
c     A domain is composed of one or more regions:
c     Build "num" T  cells with resolution "deltat(n) n=1,num" 
c     within the domain composed of regions bounded by "bounds".
c     Also construct "num" U  cells of resolution "deltau(n) n=1,num"
c     with the relation between T and U cells given by:
#ifdef centered_t
c     deltau(n) = 0.5*(deltat(n+1) + deltat(n))
#else
c     deltat(n) = 0.5*(deltau(n-1) + deltau(n))
#endif
c     Resolution may be constant or smoothly varying within each
c     region AND there must be an integral number of grid cells within 
c     each region. The domain is the sum of all regions.
c
c     inputs:
c
c     maxlen   = maximum length of "deltat" and "deltau"
c     n_bounds = number of bounds needed to define the regions
c     bounds   = latitude, longitude, or depth at each bound
c     d_bounds = delta (resolution) at each of the "bounds"
c     nbpts    = number of extra boundary cells to add to the domain.
c                (usually one at the beginning and end)
c     stretch  = stretching factor for last region (should only be used
c                in the vertical to provide increased stretching of grid
c                points. "stretch" = 1.0 gives no increased stretching.
c                "stretch" = 1.2 gives increased stretching...etc
c             
c     outputs:
c
c     num    = total number of grid cells within the domain
c     deltau = resolution of U grid cells: n=1,num
c     deltat = resolution of T grid cells: n=1,num
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
      dimension deltat(maxlen), deltau(maxlen)
      dimension d_bounds(n_bounds), bounds(n_bounds)
      logical keep_going
c
c     Set some constants
c
      p5 = 0.5
      pi = 4.0*atan(1.0)
c
c     Do all regions, one at a time, to construct the domain
c
      num  = 1
      do l=1,n_bounds-1
c
        write (stdout,'(2x,a,i2,a,g14.7,a,g14.7,a,g14.7,a,g14.7,a)')
     & ' region # ',l,'  going from ',bounds(l),' (res=',d_bounds(l)
     &,') to ',  bounds(l+1),' (res=',d_bounds(l+1),')' 
c
c       avg_res = average resolution of T cells within region
c       chg_res = change in resolution across the region
c       wid     = width of region
c       tol     = tolerance for fitting T cels within region width
c
c       provide for stretching last region if needed
c
        if (l .eq. n_bounds-1) then
          avg_res = p5*(d_bounds(l) + stretch*d_bounds(l+1))
          chg_res = (stretch*d_bounds(l+1) - d_bounds(l))
        else
          avg_res = p5*(d_bounds(l) + d_bounds(l+1))
          chg_res = (d_bounds(l+1) - d_bounds(l))
        endif
c
        tol = 1.e-5
        wid = abs(bounds(l+1) - bounds(l))
        an  = wid/avg_res
        m   = nint(an)
c
#ifdef centered_t
c
c       Calculate resolution of T cells: "deltat". Note that
c       T grid points are centered in these cells (as in MOM 1)
c       n = number of T cells fitting within the region boundaries
c
        sum = 0.0
        n   = 0
        i = 0
        keep_going = .true.
        do while (i .le. 100000 .and. keep_going)
          i = i + 1
          del = avg_res - p5*chg_res*cos((pi/m)*(i-0.5))
          if (sum + del .le. wid*(1.0 + tol)) then
            sum = sum + del
            if (num+i-1 .gt. maxlen) then
              write (stdout,*) "=>Error: maxlen exceeded in gcell. "
     &,                        " ...increase size of maxlen"
              call abort()
            endif
            deltat(num+i-1) = del
            n = n + 1
          else
            keep_going = .false.
          endif
        enddo 
#else
c
c       Calculate resolution of U cells: "deltau"
c       U grid points will be centered in these cells
c       n = number of T cells fitting within the region boundaries
c       note: "sum" initially discounts half of the U cells widths
c       at the boundaries
c
        sum = 0.5*d_bounds(l) - 0.5*d_bounds(l+1)
        n   = 0
        i = 0
        keep_going = .true.
        do while (i .le. 100000 .and. keep_going)
          i = i + 1
          del = avg_res - p5*chg_res*cos((pi/m)*i)
          if (sum + del .le. wid*(1.0 + tol)) then
            sum = sum + del
            if (num+i-1 .gt. maxlen) then
              write (stdout,*) "=>Error: maxlen exceeded in gcell. "
     &,                        " ...increase size of maxlen"
              call abort()
            endif
            deltau(num+i-1) = del
            n = n + 1
          else
            keep_going = .false.
          endif
        enddo 
#endif
        if (l .eq. n_bounds-1 .and. stretch .ne. 1.0) then
          write (stdout,'(a,i3,a,f5.2)')
     & '    constructed ',n,' cells with a stretch factor of ', stretch
          write (stdout,'(/2(a,g14.7),/2(a,g14.7),/a,a/)')
     &    'Note: you specified the ocean bottom at ',bounds(l+1)
     &,   ' cm with a bottom cell thickness of ',d_bounds(l+1)
     &,   '      The stretch factor puts the bottom at ',bounds(l)+sum
#ifdef centered_t
     &,   ' cm with a bottom cell thickness of ',deltat(num+n-1)
#else
     &,   ' cm with a bottom cell thickness of '
     &,   0.5*(deltau(num+n-1) + deltau(num+n-2))
#endif
     &,   '      Adjust "stretch_z" in subroutine "gcoord" to get'
     &,   ' closer to the desired specifications if needed.'

        else
          write (stdout,'(a,g14.7,a)')
     &   '    constructed ',an,' grid cells for this region'
          if (abs(an-n) .gt. 0.01) then
            write (stdout, '(/,a,i2,/,a,g14.7/,a,g14.7,a//a/a)')
     & '==>Error: non integral number of cells in region #',l
     &,'          average resolution within region =',avg_res
     &,'          this implies ',an,' grid cells'
     &,'          Change grid specifications within USER INPUT section'
     &,'          Here is some help...'
            d_new = (2.0*wid)/(n-1) - d_bounds(l)
            write (stdout,'(/a,i4,a,i2,a,1pe14.7,a,1pe14.7/)')
     &      ' Note: to get ',n-1,' grid cells within region ',l
     &,     ', change resolution from ', d_bounds(l+1), ' to ', d_new
            d_new = (2.0*wid)/n - d_bounds(l)
            write (stdout,'(/a,i4,a,i2,a,1pe14.7,a,1pe14.7/)')
     &      ' Note: to get ',n,' grid cells within region ',l
     &,     ', change resolution from ', d_bounds(l+1), ' to ', d_new
            d_new = (2.0*wid)/(n+1) - d_bounds(l)
            write (stdout,'(/a,i4,a,i2,a,1pe14.7,a,1pe14.7/)')
     &      ' Note: to get ',n+1,' grid cells within region ',l
     &,     ', change resolution from ', d_bounds(l+1), ' to ', d_new
            call abort()
          endif
        endif
        num = num + n
      enddo
c
c     adjust "num" to reflect the total number of cells contained in
c     all regions
c
      num = num - 1
c
      do i=1,num
#ifdef centered_t
c
c       build resolution for U cells: "deltau". Note that
c       variable resolution (stretched grid) implies U points are
c       off center as in MOM 1
c
        if (i .eq. num) then
          deltau(i) = d_bounds(n_bounds)
        else
          deltau(i) = p5*(deltat(i+1) + deltat(i))
        endif
#else
c
c       build resolution for T cells: "deltat". Note that
c       variable resolution (stretched grid) implies T points are
c       off center
c
        if (i .eq. 1) then
          deltat(i) = p5*(d_bounds(1) + deltau(i))
        else
          deltat(i) = p5*(deltau(i) + deltau(i-1))
        endif
#endif
      enddo     
c
c     add boundary points if needed
c
      if (nbpts .ne. 0) then
        do i=num,1,-1
          deltat(i+1) = deltat(i) 
          deltau(i+1) = deltau(i)
        enddo
        deltat(1)     = deltat(2)
        deltau(1)     = d_bounds(1)
        deltat(num+2) = deltat(num+1) 
        deltau(num+2) = deltau(num+1)
        num           = num + 2 
      endif
      return
      end

#if defined isotropic_grid
      subroutine iso_grid (maxlen, nylats, y_lat, dy_lat, jmt
     &,                    dytdeg, dyudeg, yt, yu)
# include "stdunits.h"
c
c=======================================================================
c     compute latitude resolution of grid cells to match the convergence
c     of meridians.
c     
c     inputs:
c
c     maxlen   = maximum length of "dytdeg" and "dyudeg"
c     nylats   = should equal 2 (defines one region)
c     dy_lat   = latitudinal resolution of grid cell on equator
c     y_lat(1) = southern boundary of the domain (it will be adjusted
c                to fit an integral number of cells)
c     y_lat(2) = northern boundary of the domain (it will be adjusted
c                to fit an integral number of cells)
c
c     outputs:
c
c     jmt  = number of grid cells
c     dytdeg = resolution of T cells (degrees)
c     dyudeg = resolution of U cells (degrees)
c     yt     = latitude of point within T cell (degrees)
c     yu     - latitude of point within U cell (degrees)
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      dimension y_lat(nylats), dy_lat(nylats)
      dimension dytdeg(maxlen), dyudeg(maxlen)
      dimension yt(maxlen), yu(maxlen)
      parameter (jmaxlat = 1000, lenjmax = 2*jmaxlat+1)
      dimension dusq(-jmaxlat-1:jmaxlat+1), dtsq(-jmaxlat-1:jmaxlat+1)
      dimension usq(-jmaxlat-1:jmaxlat+1), tsq(-jmaxlat-1:jmaxlat+1)
      dimension y_bound(2), dy_bound(2)
c
      if (nylats .gt. 2) then
        write (stdout,'(/a,a/)')
     &  '=>Error: Must only have nylats=2 when option isotropic_grid'
     &,' is enabled'
        call abort()
      endif
      if (dy_lat(1) .ne. dy_lat(2)) then
        write (stdout,'(/a,a/)')
     &  '=>Error: dy_lat(1) must equal dy_lat(2) when option '
     &,' isotropic_grid is enabled'
        call abort()
      endif
      pi   = 4.0*atan(1.0)
      rrad = pi/180.0
c
c     build an isotropic grid
c
      usq(0)  = 0.0
      dusq(0) = dy_lat(1)
      tsq(1)  = 0.5*dusq(0)
      tsq(0)  = -tsq(1)
      dtsq(1) = dusq(0)*cos(tsq(1)*rrad)
      dtsq(0) = dtsq(1)
      do n=1,jmaxlat
        dusq(n)   = 2.0*dtsq(n) - dusq(n-1)
	usq(n)    = tsq(n) + 0.5*dusq(n)
	if (tsq(n) .lt. 90) then
	  tsq(n+1)  = tsq(n) + dusq(n)
	  dtsq(n+1) = dusq(0)*cos(tsq(n+1)*rrad)
	else
	  tsq(n+1)  = tsq(n)
	  dtsq(n+1) = 0.0
	endif
      enddo
      do n=-1,-jmaxlat,-1
        dusq(n) = dusq(-n)
	usq(n)  = -usq(-n)
	dtsq(n) = dtsq(-(n-1)) 
	tsq(n)  = -tsq(-(n-1))
      enddo
c
c     pick out cells between bounding latitudes
c
      n1 = -jmaxlat
      n2 = jmaxlat
      do n=-jmaxlat,jmaxlat
        if (usq(n) .lt. y_lat(1)) n1 = n
	if (usq(n) .lt. y_lat(2)) n2 = n+1
      enddo
# if defined symmetry
      n2 = n2 + 1
# endif
c
c     re-define bounding latitudes to match isotropic boundaries
c
      y_lat(1) = usq(n1)
      y_lat(2) = usq(n2) 
c
      if (n1 .eq. -jmaxlat .or. n2 .eq. jmaxlat) then
        write (stdout,'(/a,a/)')
     &  '=>Error: Need to increase jmaxlat to reach the max latitude'
     &,' of the domain specified for option isotropic_grid.'
        call abort()
      endif
# if defined extend_isotropic_grid
c
c     extend square grid in southern hemisphere to South pole
c
c      pole1 = -90.0
c     set pole1 to skip over land in antarctica
      pole1 = -80.0
      nps = nint((y_lat(1) - pole1)/dusq(n1)) - 1
      wid = (y_lat(1) - pole1)
      y_bound(1) = pole1
      y_bound(2) = y_lat(1)
      dy_bound(1) = 2.0*wid/nps - dusq(n1)
      dy_bound(2) = dusq(n1)
      j1 = n1-nps+1
      call gcell (maxlen, nylats, y_bound, dy_bound, 0
     &,           jmts, dtsq(j1), dusq(j1), 1.0)
c
c     add boundary point
c
      dtsq(j1-1) = dtsq(j1)
      dusq(j1-1) = dusq(j1)
      jmts = jmts + 1
      n1   = n1 - jmts + 1
c
c     construct latitude of T and U grid points
c
      usq(n1) = pole1
      tsq(n1) = usq(n1) - 0.5*dusq(n1)
      do n=n1+1,n1+jmts-1
	usq(n) = usq(n-1) + dtsq(n)
	tsq(n) = tsq(n-1) + dusq(n-1)
      enddo
c
c     extend square grid to North pole
c
      pole2 = 90.0
      npn = nint((pole2 - y_lat(2))/dusq(n2)) - 1
      wid = (pole2 - y_lat(2))
      y_bound(1) = y_lat(2)
      y_bound(2) = pole2
      dy_bound(1) = dusq(n2)
      dy_bound(2) = 2.0*wid/npn - dusq(n2)
      call gcell (maxlen, nylats, y_bound, dy_bound, 0
     &,         jmtn, dtsq(n2+1), dusq(n2+1), 1.0)
c
c     add boundary point
c
      j2 = (n2+1)+jmtn-1
      dtsq(j2+1) = dtsq(j2)
      dusq(j2+1) = dusq(j2)
      jmtn = jmtn + 1
      n2   = (n2) + jmtn - 1
c
c     construct latitude of T and U grid points
c
      do n=n2-jmtn+1,n2
        usq(n) = usq(n-1) + dtsq(n)
	tsq(n) = tsq(n-1) + dusq(n-1)
      enddo
c
# endif
      jmt = n2 - n1 + 1
      do jrow=1,jmt
        dytdeg(jrow) = dtsq(jrow+n1-1)
	dyudeg(jrow) = dusq(jrow+n1-1)
        yu(jrow) = usq(jrow+n1-1)
        yt(jrow) = tsq(jrow+n1-1)
      enddo
      return
      end
#endif

#if !defined driver_only
      subroutine grids
c
c=======================================================================
c     set up a staggered "B" grid for MOM and compute grid related
c     variables
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#if !defined implicitvmix || defined isoneutralmix
# include "accel.h"
#endif
#include "coord.h"
#include "grdvar.h"
#include "hmixc.h"
#include "scalar.h"
#include "vmixc.h"
c
c-----------------------------------------------------------------------
c     set some constants
c-----------------------------------------------------------------------
c
      pi     = c4*atan(c1)
      radian = c360/(c2*pi)
      degtcm = radius/radian
c
c-----------------------------------------------------------------------
c     calculate coordinates for T and U grid cells.
c-----------------------------------------------------------------------
c
      maxlen = max(imt,jmt,km)
      call gcoord (maxlen, imt2, jmt2, km2, dxtdeg, dytdeg, dxudeg
     &,            dyudeg, dzt, dzw, xt, xu, yt, yu, zt, zw)
c
c-----------------------------------------------------------------------
c     verify that the number of grid points match the number set in
c     the parameter statement in "size.h".
c-----------------------------------------------------------------------
c
      call size_check (imt2, jmt2, km2, 'sub grids', 'stop')
c
c-----------------------------------------------------------------------
c     convert grid resolution to cm
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        dyt(jrow) = dytdeg(jrow)*degtcm
        dyu(jrow) = dyudeg(jrow)*degtcm
      enddo
c
      do i=1,imt
        dxt(i) = dxtdeg(i)*degtcm
        dxu(i) = dxudeg(i)*degtcm
      enddo
#ifdef cyclic
      dxt(1)   = dxt(imt-1)
      dxt(imt) = dxt(2)
      dxu(1)   = dxu(imt-1)
      dxu(imt) = dxu(2)
#endif
c
c-----------------------------------------------------------------------
c     compute all quantities derived from the grid spacings
c-----------------------------------------------------------------------
c
      do k=1,km
        c2dzt(k) = c2*dzt(k)
        dzt2r(k) = c1/c2dzt(k)
      enddo
c
      dzwr(km)  = c1/dzw(km)
      dzw2r(km) = p5/dzw(km)
c
      do k=1,km
        dzwr(k-1)    = c1/dzw(k-1)
        dzw2r(k-1)   = p5/dzw(k-1)
      enddo
c
      do k=1,km
#if defined implicitvmix || defined isoneutralmix
        dztur(k) = c1/(dzw(k-1)*dzt(k))
        dztlr(k) = c1/(dzw(k)*dzt(k))
#endif
#ifdef tcvmix
        dzwur(k) = c1/(dzt(k)*dzw(k))
        if (k .lt. km) dzwlr(k) = c1/(dzt(k+1)*dzw(k))
#endif
        dztr(k)  = c1/dzt(k)
      enddo
c
      tiny = 1.e-20
      do jrow=1,jmt
        dytr(jrow)  = c1/dyt(jrow)
        dyt2r(jrow) = p5/dyt(jrow)
        dyt4r(jrow) = p25/dyt(jrow)
        dyur(jrow)  = c1/dyu(jrow)
        dyu2r(jrow) = p5/dyu(jrow)
        dyu4r(jrow) = p25/dyu(jrow)
#if defined beta_plane || defined f_plane
        phi(jrow)   = yu(1)/radian
        phit(jrow)  = phi(jrow)
#else
        phi(jrow)   = yu(jrow)/radian
        phit(jrow)  = yt(jrow)/radian
#endif
        cst(jrow)   = cos(phit(jrow))
        csu(jrow)   = cos(phi (jrow))
        sine(jrow)  = sin(phi(jrow))
        if (cst(jrow)  .eq. 0.0) then
          write (stdout,'(/a,e14.7,a,i4,/a)')
     &   ' Warning: setting cst(jrow) = ',tiny, ' for jrow =',jrow
     &,  '          to prevent division by zero at the pole'
          cst(jrow)  = tiny
        endif
        if (csu(jrow)  .eq. 0.0) then
          write (stdout,'(/a,e14.7,a,i4,/a)')
     &   ' Warning: setting cst(jrow) = ',tiny, ' for jrow =',jrow
     &,  '          to prevent division by zero at the pole'
          csu(jrow)  = tiny
        endif
        cstr(jrow)     = c1/cst(jrow)
        csur(jrow)     = c1/csu(jrow)
        tng(jrow)      = sine(jrow)/csu(jrow)
        cstdytr(jrow)  = c1/(cst(jrow)*dyt(jrow))
        cstdyt2r(jrow) = cstdytr(jrow)*p5
        csudyur(jrow)  = c1/(csu(jrow)*dyu(jrow))
        csudyu2r(jrow) = p5/(csu(jrow)*dyu(jrow))
        cst_dytr(jrow) = cst(jrow)/dyt(jrow)
        csu_dyur(jrow) = csu(jrow)/dyu(jrow)
c
      enddo
c
      do i=1,imt
        dxtr(i)  = c1/dxt(i)
        dxt2r(i) = p5/dxt(i)
        dxt4r(i) = p25/dxt(i)
        dxur(i)  = c1/dxu(i)
        dxu2r(i) = p5/dxu(i)
        dxu4r(i) = p25/dxu(i)
      enddo
c
      do i=2,imtm1
        dxmetr(i) = c1/(dxt(i) + dxt(i+1))
      enddo
c
      do i=1,imt
        duw(i) = (xu(i) - xt(i))*degtcm
      enddo
      do i=1,imtm1
        due(i) = (xt(i+1) - xu(i))*degtcm
      enddo
# ifdef cyclic
      due(imt) = due(2)
# else
      due(imt) = due(imtm1)
# endif
c
      do jrow=1,jmt
        dus(jrow) = (yu(jrow) - yt(jrow))*degtcm
      enddo
c
      do jrow=1,jmtm1
        dun(jrow) = (yt(jrow+1) - yu(jrow))*degtcm
      enddo
      dun(jmt) = dun(jmtm1)
c
      do k=1,km
        if (k .gt. 1) then
          fracdz(k,0) = (zt(k) - zw(k-1))/dzt(k)
        else
          fracdz(k,0) = (zt(k) - c0)/dzt(k)
        endif
        fracdz(k,1) = (zw(k) - zt(k))/dzt(k)
      enddo
c
# if !defined implicitvmix || defined isoneutralmix
c
c     for convection code, compute values needed to include effects
c     of tracer timestep acceleration on effective layer thicknesses
c
      do k=1,km
        dztxcl(k) = dzt(k)/dtxcel(k)
      enddo
c
      do k=1,kmm1
        dzwxcl(k) = c1/(dztxcl(k)+dztxcl(k+1))
      enddo
      dzwxcl(km) = c0
c
# endif
c
# ifdef quicker
c
c     coefficients for 3rd order advective scheme
c
      do i=2,imt-1
        if (i .lt. imt-1) then
          ip2 = i+2
        else
          ip2 = 3
        endif
        quick_x(i,1) = c2*dxt(i+1)/(dxt(i+1)+dxt(i))
        quick_x(i,2) = c2*dxt(i  )/(dxt(i+1)+dxt(i))
c
        curv_xp(i,1) = c2*dxt(i)*dxt(i+1)/
     &               ((dxt(i-1)+2.0*dxt(i)+dxt(i+1))*(dxt(i)+dxt(i+1)))
        curv_xp(i,2) =-c2*dxt(i)*dxt(i+1)/((dxt(i)+dxt(i+1))
     &                                  *(dxt(i-1)+dxt(i)))
        curv_xp(i,3) = c2*dxt(i)*dxt(i+1)/
     &               ((dxt(i-1)+2.0*dxt(i)+dxt(i+1))*(dxt(i-1)+dxt(i)))
c
        curv_xn(i,1) = c2*dxt(i)*dxt(i+1)/
     &             ((dxt(i)+2.0*dxt(i+1)+dxt(ip2))*(dxt(i+1)+dxt(ip2)))
        curv_xn(i,2) =-c2*dxt(i)*dxt(i+1)/((dxt(i+1)+dxt(ip2))
     &                                  *(dxt(i)+dxt(i+1)))
        curv_xn(i,3) = c2*dxt(i)*dxt(i+1)/
     &               ((dxt(i)+2.0*dxt(i+1)+dxt(ip2))*(dxt(i)+dxt(i+1)))
      enddo
      call setbcx (quick_x, imt, 2)
      call setbcx (curv_xp, imt, 3)
      call setbcx (curv_xn, imt, 3)
c
      do j=1,jmt
        jp2 = min(j+2,jmt)
        jp1 = min(j+1,jmt)
        jm1 = max(j-1,1)
        quick_y(j,1) = c2*dyt(jp1)/(dyt(jp1)+dyt(j))
        quick_y(j,2) = c2*dyt(j  )/(dyt(jp1)+dyt(j))
c
        curv_yp(j,1) = c2*dyt(j)*dyt(jp1)/
     &               ((dyt(jm1)+2.0*dyt(j)+dyt(jp1))*(dyt(j)+dyt(jp1)))
        curv_yp(j,2) =-c2*dyt(j)*dyt(jp1)/((dyt(j)+dyt(jp1))
     &                                  *(dyt(jm1)+dyt(j)))
        curv_yp(j,3) = c2*dyt(j)*dyt(jp1)/
     &               ((dyt(jm1)+2.0*dyt(j)+dyt(jp1))*(dyt(jm1)+dyt(j)))
c
        curv_yn(j,1) = c2*dyt(j)*dyt(jp1)/
     &             ((dyt(j)+2.0*dyt(jp1)+dyt(jp2))*(dyt(jp1)+dyt(jp2)))
        curv_yn(j,2) =-c2*dyt(j)*dyt(jp1)/((dyt(jp1)+dyt(jp2))
     &                                  *(dyt(j)+dyt(jp1)))
        curv_yn(j,3) = c2*dyt(j)*dyt(jp1)/
     &               ((dyt(j)+2.0*dyt(jp1)+dyt(jp2))*(dyt(j)+dyt(jp1)))
      enddo
c
      do k=1,km
        kp2 = min(k+2,km)
        kp1 = min(k+1,km)
        km1 = max(k-1,1)
        quick_z(k,1) = c2*dzt(kp1)/(dzt(kp1)+dzt(k))
        quick_z(k,2) = c2*dzt(k  )/(dzt(kp1)+dzt(k))
c
        curv_zp(k,1) = c2*dzt(k)*dzt(kp1)/
     &               ((dzt(km1)+2.0*dzt(k)+dzt(kp1))*(dzt(k)+dzt(kp1)))
        curv_zp(k,2) =-c2*dzt(k)*dzt(kp1)/((dzt(k)+dzt(kp1))
     &                                  *(dzt(km1)+dzt(k)))
        curv_zp(k,3) = c2*dzt(k)*dzt(kp1)/
     &               ((dzt(km1)+2.0*dzt(k)+dzt(kp1))*(dzt(km1)+dzt(k)))
c
        curv_zn(k,1) = c2*dzt(k)*dzt(kp1)/
     &             ((dzt(k)+2.0*dzt(kp1)+dzt(kp2))*(dzt(kp1)+dzt(kp2)))
        curv_zn(k,2) =-c2*dzt(k)*dzt(kp1)/((dzt(kp1)+dzt(kp2))
     &                                  *(dzt(k)+dzt(kp1)))
        curv_zn(k,3) = c2*dzt(k)*dzt(kp1)/
     &               ((dzt(k)+2.0*dzt(kp1)+dzt(kp2))*(dzt(k)+dzt(kp1)))
      enddo
c
# endif
      return
      end
#endif


