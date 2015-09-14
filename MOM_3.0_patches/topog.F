#ifdef drive_topog
# define driver_only
# include "grids.F"
# include "util.F"
# include "iomngr.F"
# include "size_check.F"

      program driver
c
c=======================================================================
c
c              T O P O G R A P H Y   M O D U L E
c
c     To generate a geometry and topography, use the "run_topog" script
c     Consult the manual for available options.
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use mpp_io_mod
      use mom3_mpp_mod
# include "size.h"
# include "stdunits.h"
# include "coord.h"
      dimension kmt(imt,jmt), map(imt,jmt)
      dimension ht(imt,jmt)
c
      stdout = 6
      call mpp_io_init
      call ioinit
c
      write (stdout,'(//1x,a,//)') 'Construct the MOM grid'
c
c-----------------------------------------------------------------------
c     construct the grid.
c-----------------------------------------------------------------------
c
      maxlen = max(imt, jmt, km)
      call gcoord (maxlen, imt2, jmt2, km2
     &,            dxtdeg, dytdeg, dxudeg, dyudeg
     &,            dzt, dzw, xt, xu, yt, yu, zt, zw)
c
c     check that returned grid sizes match those in file "size.h"
c
      call size_check (imt2, jmt2, km2, 'driver', 'stop')
c
c-----------------------------------------------------------------------
c     generate the kmt field on the grid defined by grids.F
c-----------------------------------------------------------------------
c
      write (stdout,'(/,20x,a/)') ' Driving the topog module...'
      call topog (kmt, map, xt, yt, zt, xu, yu, zw, dzt
     &,           imt2, jmt2, km2, ht)
c
      call mpp_io_exit
      end
#endif



      subroutine topog (kmt, map, xt, yt, zt, xu, yu, zw, dzt
     &,                 imt2, jmt2, km2, ht)
c
c-----------------------------------------------------------------------
c     construct the topographic mask "kmt" which determines the geometry
c     and topography by defining the number of vertical levels on model
c     T points.
c
c     OPTION:
c      -Dpartial_cell: allow bottom cells to be partially filled with
c                      land to give accurate estimate of bottom. 
c                      Otherwise, full-cells are assumed with
c                      inaccurate bottom.
c
c
c     A base kmt field is constructed using one of the following:
c
c      (A) Construct a kmt field using one of the following:
c        -Dscripps_kmt : interpolate Scripps 1 degree topography to
c                        "kmt" on MOM grid
c        -Detopo_kmt : interpolate Etopo 1/12 degree topography to
c                        "kmt" on MOM grid
c        -Dread_my_kmt : read in your own "kmt" from a file "kmt.dta"
c        -Didealized_kmt : generate an "idealized" not very realistic
c                          topography without reading the database
c        -Drectangular_box : sets "kmt" to a flat bottom rectangular box
c        -Dgaussian_kmt : sets "kmt" to a gaussian bump on a sloping
c                          bottom
c
c       Options:
c        -Dcyclic        : sets cyclic conditions in zonal direction
c                          (otherwise solid_walls is assumed)
c        -Dflat_bottom   : sets "kmt" over ocean to maximum levels (km)
c
c     Modifications to the kmt field:
c
c      1) Fill isolated cells (optional)
c
c       Option:
c        -Dfill_isolated_cells : These are places where "T" cells are
c                                 deeper than nearest neighbors. Also
c                                 fills in "U" cells with land at all
c                                 four corners.
c
c      2) Add USER CHANGES to kmt field (if any)
c
c      3) Limit minimum levels based on parameter "kmt_min"
c          "kmt_min" is the minimum number of vertical levels
c           permitted at any ocean "T" cell.  kmt_min must be at >= 2.
c          
c           Options:
c           -Dfill_shallow: makes all T cells for which
c                           kmt(i,jrow) < kmt_min into land.
c           -Ddeepen_shallow: deepens such cell to depth kmt_min.
c           -Dround_shallow: changes depth to the nearer of 0 and
c                            zw(kmt_min)
c
c-----------------------------------------------------------------------
c     input:
c
c     xt  = longitudes of model grid T points (degrees)
c     yt  = latitudes of model grid T points (degrees)
c     zt  = depths of box "u/v" points on model grid (cm)
c     xu  = longitudes of model grid U points (degrees)
c     yu  = latitudes of model grid U points (degrees)
c     zw  = depths of box bottoms on model grid T boxes (cm)
c     imt = number of model longitudes
c     jmt = number of model latitudes
c     km  = maximum number of levels
c
c     output:
c
c     kmt = number of vertical levels on T points
c           0 is a land point.
c           n implies n boxes between the surface and bottom
c             where (n >= kmt_min) & (n <= km)
c     map = array indicating land masses and their perimeters
c           0 is a mid ocean point
c           n > 0 is a land mass
c           n < 0 is the ocean perimeter of land mass |n|
c
c     ht  = depth to ocean bottom at T-cells (cm)
# ifdef bbl_ag
c           not including the thickness of the BBL (the actual ocean
c           depth will be deeper than "ht" by the BBL thickness).
# endif
c-----------------------------------------------------------------------
c
#if defined sponges && !defined driver_only
      use xsponge_mod
#endif
#include "stdunits.h"
#include "size.h"
#if defined partial_cell && ( defined firfil || defined fourfil )
# if !defined driver_only
#  include "index.h"
# endif
#endif
# ifdef bbl_ag
#  include "bbl.h"
# endif
c
c     kmt_min is the minimum permitted number of levels at an ocean
c     T cell.  kmt_min must be at least 2.
c
      parameter (kmt_min = 2)
c
      dimension kmt(imt,jmt), xt(imt), yt(jmt), zt(km)
      dimension xu(imt), yu(jmt), zw(km), dzt(km)
      dimension hu(imt,jmt), slope(7)
      dimension ruff(imt,jmt)
      dimension ht(imt,jmt)
      dimension kmt_old(imt,jmt)
#if defined fill_isolated_dont_change_landmask
      dimension kmt_dont_change_landmask(imt,jmt)
#endif
#if defined sink_isolated_land_cells
      dimension kmt_isolated_land_cells(imt,jmt)
#endif
#if defined smooth_topog_after_user_changes
      dimension depth_after_user_changes(imt,jmt)
#endif
c
      dimension map(imt,jmt)
      dimension iperm(maxipp)
      dimension jperm(maxipp)
      dimension nippts(mnisle)
      dimension iofs(mnisle)
c
      real(kind=4) ruff4(imt,jmt),ht4(imt,jmt),spdpm4,spng4(jmt)
      
      character*80 iotext
      character*60 expnam
      character*32 stamp
c
      kmu(i,j) = min(kmt(i,j), kmt(i+1,j), kmt(i,j+1), kmt(i+1,j+1))
      hum(i,j) = min(ht(i,j), ht(i+1,j), ht(i,j+1), ht(i+1,j+1))
c
#if defined partial_cell && ( defined firfil || defined fourfil )
# if defined driver_only
      write (stdout,*) '=>Error: polar filering (firfil or fourfil)'
     &,                '         is not allowed with partial_cell  '
     &,                '         when executing topog.F as a stand '
     &,                '         alone module.'
      call abort()
# endif
#endif
c
c-----------------------------------------------------------------------
c     check that grid sizes in argument list and file "size.h" agree
c-----------------------------------------------------------------------
c
      call size_check (imt2, jmt2, km2, 'topog', 'stop')
c
      ncase = 0
c
      write (stdout,'(//,t37,a,/)')
     & 'T O P O G R A P H Y    G E N E R A T I O N'
      write (stdout,'(a,i4,/)') 'kmt_min = ', kmt_min
c
c     show options
c
      call topog_options
c
c     initialize bottom depth, roughness, and number of vertical levels
c     to T cell bottoms
c
      do jrow=1,jmt
        do i=1,imt
	  ht(i,jrow)   = 0.0
	  ruff(i,jrow) = 0.0
	  kmt(i,jrow)  = 0
	enddo
      enddo
c
c-----------------------------------------------------------------------
c     discretized ocean depth to number of T cells "kmt" and construct
c     corresponding depth to bottom T cell "ht"
c-----------------------------------------------------------------------
c
#ifdef scripps_kmt
c
c-----------------------------------------------------------------------
c     construct model "kmt" based on Scripps 1 deg topography data.
c-----------------------------------------------------------------------
c
      ncase = ncase + 1
      write (stdout,'(/a/)')
     & ' Interpolating Scripps 1 degree topography to the MOM grid:'
c
      call scripp (kmt, xt, yt, zt, xu, yu, zw, dzt, imt, jmt, km, ht)
#endif
c
#ifdef etopo_kmt
c
c-----------------------------------------------------------------------
c     construct model "kmt" based on NGDC ETOPO5 1/12 deg topography
c-----------------------------------------------------------------------
c
      ncase = ncase + 1
      write (stdout,'(/a/)')
     & ' Interpolating ETOPO 1/12 degree topography to the MOM grid:'
c
      call etopo (kmt, xt, yt, zt, xu, yu, zw, dzt, imt, jmt, km, ht
     &,           ruff)
#endif
#ifdef iowxyz_kmt
c
c-----------------------------------------------------------------------
c     construct model "kmt" based on IOW-topography
c-----------------------------------------------------------------------
c
      ncase = ncase + 1
      write (stdout,'(/a/)')
     & ' Interpolating iow-database topography to the MOM grid:'
c
      call iowxyz (kmt, xt, yt, zt, xu, yu, zw, dzt, imt, jmt, km, ht
     &,           ruff)
#endif
#if !defined etopo_kmt && !defined read_my_roughness
      write (stdout,'(/a,a,/)')
     & '=>Warning: roughness will be zero unless option "etopo_kmt"'
     &,' or option read_my_roughness is used.'
#endif
#if defined read_my_roughness
c
c---------------------------------------------------------------------
c     import the estimated roughness field
c---------------------------------------------------------------------
c
      write (stdout,'(/a/)')
     & ' =>Reading a "roughness" field from file "rough.dta" '
      call getunit (io, 'rough.dta'
     &,             'unformatted sequential rewind ieee')
c
      read (io) stamp, iotext, expnam
      read (io) imax, jmax, kmax
      call size_check (imax, jmax, kmax, 'reading roughness', 'stop')
      read (io) stamp, iotext, expnam
# if defined _CRAY
      read (io) ruff
# else      
      read (io) ruff4
      ruff=ruff4
# endif      
      call relunit (io)
      call print_checksum (ruff, imt, jmt,' imported ruff checksum = ')
      write (stdout,'(/a/)')'==> The roughness field has been read.'
#endif
#ifdef read_my_kmt
c
c-----------------------------------------------------------------------
c     import model "kmt" field
# ifdef partial_cell
c     import "ht" as the total depth on T cells
# endif
c-----------------------------------------------------------------------
c
      ncase = ncase + 1
      write (stdout,'(/a/)') ' =>Reading a "kmt" field from "kmt.dta" '
      call shuttle_kmt (kmt, imt, jmt, km, 'import', 'kmt.dta')
      call print_checksumi (kmt, imt, jmt, 'Imported "kmt" checksum =')
# ifdef partial_cell
      write (stdout,'(/a/)')' =>Reading a "ht" field from "ht.dta"'
      call getunit (io, 'ht.dta'
     &,             'unformatted sequential rewind ieee')
c
      read (io) 
      read (io) imax, jmax, kmax
      call size_check (imax, jmax, kmax, 'reading ht', 'stop')
      read (io)
#if defined _CRAY 
      read (io) ht
#else      
      read (io) ht4
      ht=ht4
#endif
      call relunit (io)
      call print_checksum (ht, imt, jmt,' imported ht checksum = ')
      write (stdout,'(5(/a)/)')  '==> The "ht" field has been read.'
     &,                       '    Make sure it is consitent with kmt.'
     &,                       '    (i.e. Any kmt(i,jrow) that was     '
     &,                       '     changed should have ht(i,jrow) set'
     &,                       '     to the depth of the full cell.)   '
# endif
c
      nerror = 0
      do i=1,imt
        do j=1,jmt
          if (kmt(i,j) .lt. 0 .or. kmt(i,j) .gt. km) then
            nerror = nerror + 1
          end if
        end do
      end do
      if (nerror .gt. 0) then
        write (stdout,'(a,i4,a,a,i4,/)')
     &    '==>ERROR: There are ', nerror, ' depth values in the file '
     &,   '"kmt.dta" that are outside the range 0 <= kmt(i,j) <= km = '
     &,   km
        call abort()
      end if
      write (stdout,'(/a/)') '==> The kmt field has been imported.'
#endif
#ifdef idealized_kmt
c
c-----------------------------------------------------------------------
c     construct a highly "idealized" world ... piece by piece
c
c     note: the purpose of this geometry/topography is to automatically
c           map into arbitrary resolution as grid dimensions "imt" and
c           "jmt" are changed, thereby facilitating the implementation
c           and verification of the model on various computer platforms
c           without referencing the topographic data base.  Although it
c           somewhat resembles the real world, it is NOT realistic.
c-----------------------------------------------------------------------
c
      ncase = ncase + 1
      write (stdout,'(/a/a/)')
     &' Constructing an idealized test case "kmt" which only crudely'
     &,' approximates world geometry. Topography is unrealistic.'
c
      call idealized (kmt, xt, yt, zw, imt, jmt, kbot, ht)
c
# ifdef partial_cell
      write (stdout,'(/a/a/)')
     &' ==> Warning: Partial-cells are being set to full-cell thickness'
     &,'            for the option idealized_kmt'
# endif
#endif
c
#ifdef gaussian_kmt
c
c-----------------------------------------------------------------------
c     construct a topographic bump with a gaussian structure on a
c     gently sloping ocean floor.
c     USER INPUT: set the slope of the floor and scale of the bump
c                 as desired
c-----------------------------------------------------------------------
c
      ncase = ncase + 1
      bump_height = 0.5*zw(km)
      bump_scale = 0.25*min(xt(imt-1)-xt(2), abs(yt(jmt-1)-yt(2)))
      slope_x = 0.e2
      slope_y = 0.e2
      write (stdout,'(/a,f8.2,a/a,f8.2,a/a,f8.2,a/a,f8.2,a/)')
     &' Constructing a gaussian bump of height =', bump_height/100.0
     &,' meters',' and a scale width of ', bump_scale,' degrees.'
     &,' The ocean bottom rises with a slope of ',slope_x/100.0
     &,' meters/deg towards the east',' and ', slope_y/100.0
     &,' meters/deg to the north.'
      xcent = xu(imt/2)
      ycent = yu(jmt/2) 
      do jrow=2,jmt-1
        do i=2,imt-1
	  arg = ((xt(i)-xcent)**2 + (yt(jrow) - ycent)**2)
	  bottom = zw(km) - bump_height*exp(-arg/bump_scale**2)
c
          bottom = bottom - slope_x*(xt(i)-xu(1))
     &                    - slope_y*(yt(jrow)-yu(1))
c
	  bottom = max(bottom,zw(2))	  
# ifdef bbl_ag
          if (bottom .gt. dzt(1)+dzt(km)) bottom = bottom - dzt(km)
# endif
c
	  kmt(i,jrow) = indp (bottom, zw, kbot) 
# ifdef partial_cell
c
c         keep physical "bottom" <= zw(kmt(i,jrow))
c
          ht(i,jrow) = bottom
          if (zw(kmt(i,jrow)) .lt. ht(i,jrow)) then
	    if (kmt(i,jrow) .eq. kbot) then
              ht(i,jrow) = zw(kmt(i,jrow))
#  if defined bbl_ag
     &                     - dzt(km)
#  endif
	    else
	      kmt(i,jrow) = kmt(i,jrow) + 1
	    endif
	  endif
# else
          if (kmt(i,jrow) .eq. 0) then
            ht(i,jrow) = 0.0
	  else
            ht(i,jrow) = zw(kmt(i,jrow))
	  endif
# endif
	enddo
      enddo
#endif
c
#ifdef rectangular_box
c
c-----------------------------------------------------------------------
c     construct a flat bottomed rectangular box
c-----------------------------------------------------------------------
c
      ncase = ncase + 1
      write (stdout,'(/a,i4,a/)')
     & 'Constructing a flat bottomed rectangular box with ',km,' levels'
c
      do jrow=2,jmt-1
        do i=2,imt-1
          kmt(i,jrow) = kbot
# ifdef partial_cell
	  ht(i,jrow) = zw(kmt(i,jrow))
# endif
        enddo
      enddo
#endif
c
c-----------------------------------------------------------------------
c     check that only one method of generating topography was chosen
c-----------------------------------------------------------------------
c
      if (ncase .eq. 0) then
        write (stdout,'(/a/a/a/)')
     &  '=>Error: One of the following options must be enabled:'
     &, '      idealized_kmt, scripps_kmt, read_my_kmt, rectangular_box'
     &, '      etopo_kmt, gaussian_kmt, iowxyz_kmt'
        call abort()
      elseif (ncase .gt. 1) then
        write (stdout,'(/a/a/a/)')
     &  '=>Error: Only one of the following options may be enabled:'
     &, '      idealized_kmt, scripps_kmt, read_my_kmt, rectangular_box'
     &, '      etopo_kmt, gaussian_kmt, iowxyz_kmt'
        call abort()
      endif
c
#ifdef flat_bottom
c
c-----------------------------------------------------------------------
c     set the bottom to the deepest level
c-----------------------------------------------------------------------
c
      do jrow=2,jmt-1
        do i=2,imt-1
          if (kmt(i,jrow) .ne. 0) then
            kmt(i,jrow) = kbot
	    ht(i,jrow) = zw(kbot)
          endif
        enddo
      enddo
#endif
c
#ifdef bbl_ag
      nbbl = 0
      do jrow=2,jmt-1
        do i=2,imt-1
          if (kmt(i,jrow) .gt. km-1) then
            nbbl = nbbl + 1
            write (stdout,'(/a/)')
     &      '=>Error: kmt(i,jrow) cannot equal km with bbl_ag enabled.'
            write(stdout,*) "location coordinates: i,jrow=",i,jrow
          end if
        enddo
      enddo
      if (nbbl .gt. 0) then
        write (stdout,'(/a/)')
     &  '=>Error: kmt(i,jrow) cannot equal km with bbl_ag enabled.'
        call abort()
      endif
#endif
c
c-----------------------------------------------------------------------
c     apply lateral boundary conditions
c-----------------------------------------------------------------------
c
      call kmtbc (kmt, imt, jmt, ht)
c
#ifdef obc
c
c-----------------------------------------------------------------------
c     over-ride boundary conditions with open boundary conditions
c-----------------------------------------------------------------------
c
# ifdef obc_south
      do jrow=1,jmt
        do i=1,imt
          if (jrow .eq. 2) then
	    write(stdout,*)' Changing Topography at southern OB. '
	    kmt(i,jrow)=min(min(kmt(i,jrow),kmt(i,jrow+1)),
     &                  min(kmt(i,jrow),kmt(i,jrow+2)))
            kmt(i,jrow-1)=kmt(i,jrow)
            kmt(i,jrow+1)=kmt(i,jrow)
	    kmt(i,jrow+2)=kmt(i,jrow+1)
	  endif
        enddo
      enddo
# endif
# ifdef obc_north
      do jrow=1,jmt
        do i=1,imt
          if (jrow .eq. jmt-1) then
	    write(stdout,*)' Changing Topography at northern OB. '
	    kmt(i,jrow)=min(min(kmt(i,jrow),kmt(i,jrow-1)),
     &                  min(kmt(i,jrow),kmt(i,jrow-2)))
            kmt(i,jrow+1)=kmt(i,jrow)
            kmt(i,jrow-1)=kmt(i,jrow)
	    kmt(i,jrow-2)=kmt(i,jrow-1)
	  endif
        enddo
      enddo
# endif
# ifdef obc_east
      do jrow=1,jmt
        do i=1,imt
          if (i .eq. 2) then
	    write(stdout,*)' Changing Topography at eastern OB. '
	    kmt(i,jrow)=min(min(kmt(i,jrow),kmt(i+1,jrow)),
     &                  min(kmt(i,jrow),kmt(i+2,jrow)))
            kmt(i-1,jrow)=kmt(i,jrow)
            kmt(i+1,jrow)=kmt(i,jrow)
	    kmt(i+2,jrow)=kmt(i+1,jrow)
	  endif
        enddo
      enddo
# endif
# ifdef obc_west
      do jrow=1,jmt
        do i=1,imt
          if (i .eq. imt-1) then
	    write(stdout,*)' Changing Topography at western OB. '
	    kmt(i,jrow)=min(min(kmt(i,jrow),kmt(i-1,jrow)),
     &                  min(kmt(i,jrow),kmt(i-2,jrow)))
            kmt(i+1,jrow)=kmt(i,jrow)
            kmt(i-1,jrow)=kmt(i,jrow)
	    kmt(i-2,jrow)=kmt(i-1,jrow)
	  endif
        enddo
      enddo
# endif
#endif
#ifdef debug_kmt
c
c-----------------------------------------------------------------------
c     show the original (unmodified) kmt field
c-----------------------------------------------------------------------
c
      write (stdout,'(/a/)')
     & 'Remove -Ddebug_kmt to eliminate the following map'
      write (stdout,'(/,20x,a/)') ' The original "kmt" field follows:'
      call iplot (kmt, imt, imt, jmt)
#else
      write (stdout,'(/a/)') 'Use -Ddebug_kmt to show the first kmt map'
#endif
c
      call print_checksumi (kmt, imt, jmt, 'Original "kmt" checksum =')
      call area_volume (kmt, xu, yu, zw, ht)
c
c     save a copy of kmt before making changes to kmt
c
      do jrow=1,jmt
        do i=1,imt
	  kmt_old(i,jrow) = kmt(i,jrow)
	enddo
      enddo
c
# ifdef bbl_mask
c-----------------------------------------------------------------------
c     Use mask_t_bbl to specify ocean locations where the bbl is not
c     active (i.e. is set to land).  Do this by setting mask_t_bbl to
c     zero at these locations.  Note: bbl-umask is updated based on
c     tmask after the user input section.
c-----------------------------------------------------------------------
c     initialize the bbl t-mask:
      mask_t_bbl(1:imt,1:jmt) = 1
      write(stdout,*)" Note: ifdef bbl_mask is used.  Can specify ocean"
      write(stdout,*)" grid points where the bbl is inactive by setting"
      write(stdout,*)" mask_t_bbl(i,j) to zero (land) at these points."
# endif

c-----------------------------------------------------------------------
c     USER INPUT
c     other changes to "kmt" may be made here. 
c-----------------------------------------------------------------------
c
      write (stdout,'(/a/)') 'USER MODIFICATIONS to "kmt" follow:'
c
#ifdef test_case_A      
# include "topog_test_case_A.inc"
#endif
c-----------------------------------------------------------------------
c     END USER INPUT
c-----------------------------------------------------------------------
c
# ifdef bbl_mask
c     calculate bbl-umask based on tmask after user input:
      mask_u_bbl(1:imt,1:jmt) = 0
      do i=2,imt-1
        do j=2,jmt-1
          mask_u_bbl(i,j)=min(
     &         mask_t_bbl(i,j), mask_t_bbl(i+1,j)
     &         ,mask_t_bbl(i,j+1), mask_t_bbl(i+1,j+1))
        end do
      end do
# endif

#if defined smooth_topog_after_user_changes
c
c     Smooth topog after user changes:
c
c     1) Convert kmt field back to depth field
c
      do jrow=1,jmt
        do i=1,imt
	  if (kmt(i,jrow) .ne. 0) then
# ifdef partial_cell
            if (kmt(i,jrow).ne.kmt_old(i,jrow)) then
              depth_after_user_changes(i,jrow) = zw(kmt(i,jrow))
            else
              depth_after_user_changes(i,jrow) = ht(i,jrow)
            end if
# else
	    depth_after_user_changes(i,jrow) = zw(kmt(i,jrow))
# endif
	  else
	    depth_after_user_changes(i,jrow) = 0.0
	  endif
	enddo
      enddo
      call kmtbc (kmt, imt, jmt, depth_after_user_changes)
c
c     2) Smooth depth field
c
      crit_lat = -85.0
      num_pass = 2
      j_crit = indp (crit_lat, yt, jmt)
c
      write (stdout,'(/a,f5.1,a,i3,a/)')
     & '=>Filtering topography after user changes northward of lat='
     &,yt(j_crit),' with ',num_pass,' passes.'
c
      call filter_topo (depth_after_user_changes(1,1), imt, jmt
     &,                 j_crit, jmt-1, num_pass)
c
c     3) Discretize depth back to kmt
c
# if defined partial_cell && ( defined fourfil || defined firfil )
c
c     make sure there are no partial cells where polar filtering occurs:
c     calculate filtering rows (copied from setocn):
c
#  ifdef driver_only
c     set up model indices for filtering polar latitudes
      call set_polar_filtering_indices
#  endif
      j_ft1   = jft1
      j_ft2   = jft2
# else
      j_ft1   = 0
      j_ft2   = jmt+1
# endif

      do jrow=1,jmt
        do i=2,imt-1
          kmt(i,jrow) = 0
          if (depth_after_user_changes(i,jrow) .ne. 0.0) then
            kmt(i,jrow) = indp(depth_after_user_changes(i,jrow),zw,kbot)
# ifdef partial_cell
c
c           adjust "kmt" so that "a(i,jrow)" <= zw(kmt(i,jrow))
c
            if (jrow .gt. j_ft1 .and. jrow .lt. j_ft2) then
              if (zw(kmt(i,jrow)) .lt. depth_after_user_changes(i,jrow))
     &             then
                if (kmt(i,jrow) .eq. kbot) then
                  depth_after_user_changes(i,jrow) = zw(kmt(i,jrow))
                else
                  kmt(i,jrow) =kmt(i,jrow) + 1
                endif
              endif
            else
              depth_after_user_changes(i,jrow) = zw(kmt(i,jrow))
            endif
# endif
          endif
        enddo
# ifdef cyclic
        kmt(1,jrow)   = kmt(imt-1,jrow)
        kmt(imt,jrow) = kmt(2,jrow)
# else
        kmt(1,jrow)   = 0
        kmt(imt,jrow) = 0
# endif
      enddo
#endif
c
#if defined sponges && !defined driver_only
c
c     retrieve the number of sponge rows and damping timescales
c
      allocate (spng(jmt))
      call getunit (io, 'sponge.mom'
     &,             'unformatted sequential rewind ieee')
      read (io)
#if defined _CRAY      
      read (io) sstamp, spdpm, im, kk, jm, max_rows, mm, spng
#else
      read (io) sstamp, spdpm4, im, kk, jm, max_rows, mm, spng4
      spdpm=spdmp4
      spng=spng4
#endif      
      call relunit (io)
c
c     allocate storage for sponge arrays
c
      allocate (bufmb(imt,km,max_rows,2,2))
      allocate (bufmb1(imt,km,max_rows,2))
      allocate (sponge_row(jmt))
c
# if !defined sponges_old
c
c     reset kmt within sponge zones
c
      write (stdout,*) '=> kmt changes for sponges'
      j_south_spng = 1
      do jrow=jmt/2,2,-1
        if (spng(jrow) .ne. 0.0) then
	  j_south_spng = jrow
	  exit
	endif
      enddo
      do jrow = 2,j_south_spng-1
        write (stdout,*) 
     &'=> Resetting kmt on jrow=',jrow,' to kmt on jrow=',j_south_spng
        do i=1,imt
	  kmt(i,jrow) = kmt(i,j_south_spng)
	enddo
      enddo
c
      j_north_spng = jmt
      do jrow=jmt/2,jmt-1
        if (spng(jrow) .ne. 0.0) then
	  j_north_spng = jrow
	  exit
	endif
      enddo
      do jrow = j_north_spng+1,jmt-1
        write (stdout,*)
     &'=> Resetting kmt on jrow=',jrow,' to kmt on jrow=',j_north_spng
        do i=1,imt
	  kmt(i,jrow) = kmt(i,j_north_spng)
	enddo
      enddo
# endif
#endif
c
c-----------------------------------------------------------------------
c     reset ht to be consistant with kmt
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        do i=1,imt
	  if (kmt(i,jrow) .ne. kmt_old(i,jrow)) then
	    if (kmt(i,jrow) .ne. 0) then
	      ht(i,jrow) = zw(kmt(i,jrow))
	    else
	      ht(i,jrow) = 0.0
	    endif
	  endif
	enddo
      enddo
      call kmtbc (kmt, imt, jmt, ht)
c
c-----------------------------------------------------------------------
c     limit the minimum number of levels. kmt_min should be >= 2
c-----------------------------------------------------------------------
c
      call min_depth (kmt, kmt_min, xu, yu, zw, ht)
c
#ifdef fill_isolated_cells
c
c-----------------------------------------------------------------------
c     fill isolated cells (potholes and trenches) at all levels in kmt
c     filled kmt array is the maximum of the surrounding kmu levels
c     symmetry conditions are automatic
c-----------------------------------------------------------------------
c
# if defined fill_isolated_dont_change_landmask
c
c     Save kmt before filling isolated cells
c
      do jrow=1,jmt
        do i=1,imt
	  kmt_dont_change_landmask(i,jrow) = kmt(i,jrow)
	enddo
      enddo
# endif
      n = 0
      write (stdout,'(/a)') 'Searching for isolated T cells...'
      do j=2,jmt-1
        do i=2,imt-1
          k = max (kmu(i,j), kmu(i-1,j), kmu(i,j-1), kmu(i-1,j-1))
          if (kmt(i,j) .ne. k) then
            n = n + 1
            kmt(i,j) = k
	    ht(i,j) = max(hum(i,j), hum(i-1,j), hum(i,j-1),hum(i-1,j-1))
# ifdef cyclic
            kmt(1,j)   = kmt(imt-1,j)
	    kmt(imt,j) = kmt(2,j)
            ht(1,j)   = ht(imt-1,j)
	    ht(imt,j) = ht(2,j)
# endif
          endif 
        enddo
      enddo
      if (n .gt. 0) then
        write (stdout,*) '-> Found ',n,' and filled them in.'
      else
        write (stdout,*) '-> None Found.'
      endif 
c
# if defined fill_isolated_dont_change_landmask
c
c     restore kmt to previous value where landmask has changed:
c
      do jrow=1,jmt
        do i=1,imt
          if (kmt_dont_change_landmask(i,jrow).ne.0 .and.
     &        kmt(i,jrow).eq.0) then
            iii=kmt_dont_change_landmask(i,jrow)-kmt(i,jrow)
            kmt(i,jrow) = kmt_dont_change_landmask(i,jrow)
            kmt_dont_change_landmask(i,jrow)=iii
            write (stdout,*) "Note: not filling isolated cell "
     &,     "which changed land mask at i,jrow =",i,jrow
          else
            kmt_dont_change_landmask(i,jrow) = 0
          end if
	enddo
      enddo
      write (stdout,'(/a/)')
     & 'locations where fill_isolated_cells was overridden:'
      call iplot (kmt_dont_change_landmask, imt, imt, jmt)
# endif
      call kmtbc (kmt, imt, jmt, ht)
#else
      write (stdout,'(/a)') 'Note=> Not filling isolated T cells...'
#endif
      call kmtbc (kmt, imt, jmt, ht)

#ifdef sink_isolated_land_cells
c     iterate on sinking isolated cells...: 
      num_sunk_isolated_cells=-1
      do while(num_sunk_isolated_cells.ne.0)
        call sink_isolated_cells(kmt,kmt_isolated_land_cells
     &       ,imt,jmt,num_sunk_isolated_cells)
      end do
      call kmtbc (kmt, imt, jmt, ht)
#endif

c
c-----------------------------------------------------------------------
c     construct the map of land masses and perimeters
c     Note: a meta land mass is composed of distinct land masses with
c           a minimum separation of one ocean T cell. No net flow exists
c           between members of the meta land mass
c-----------------------------------------------------------------------
c
      call isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle
     &,                    imt, jmt, km, mnisle, maxipp)
      call kmtbc (kmt, imt, jmt, ht)
#ifdef partial_cell
c
c     don`t allow partial cells thickness less than min(dzt(k))
c     (arbitrarily 50m minimum in coarse models ... otherwise dzt(1))
c
      p_cell_min = min(50.0e2,zw(1))
      do k=2,km
        if (p_cell_min .gt. (zw(k)-zw(k-1))) p_cell_min = zw(k)-zw(k-1) 
      enddo
      do jrow=2,jmt-1
        do i=2,imt-1
	  k = kmt(i,jrow)
	  if (k .lt. 2) then
	    ht(i,jrow) = 0.0
	  else
	    if ((ht(i,jrow)-zw(k-1)) .lt. p_cell_min) then
	      ht(i,jrow) = zw(k-1) + p_cell_min
	    endif
	  endif
        enddo
      enddo 
      call kmtbc (kmt, imt, jmt, ht)
#endif
      call area_volume (kmt, xu, yu, zw, ht)
c
#  ifndef skip_island_map
      write (stdout,'(/a/)')
     &   'Use -Dskip_island_map to eliminate the following map'
      call showmap (map, imt, jmt, linewidth)
#  else
      write (stdout,'(/a/)')
     &   'Remove -Dskip_island_map to show perimeter map'
#  endif
#  ifdef show_perimeter_details
      write (stdout,'(/a/)')
     & 'Remove -Dshow_perimeter_details to hide the following indices'
      do isle=1,nisle
         write (stdout,'(/a,i3,/)')
     &   'Perimeter indices for land mass #',isle
         write (stdout,'(10(a,i4,a,i4,a,2x))')
     &   ('(',iperm(iofs(isle)+n),',',jperm(iofs(isle)+n),')'
     &,  n=1,nippts(isle))
      enddo
#  else
      write (stdout,'(/a/)')
     &   'Use -Dshow_perimeter_details to show perimeter indices'
#  endif
c
c-----------------------------------------------------------------------
c     show the final kmt field
c-----------------------------------------------------------------------
c
#  ifndef skip_kmt_map
      write (stdout,'(/a/)')
     & 'Use -Dskip_kmt_map to eliminate the following map'
      write (stdout,'(/,20x,a/)') ' The final "kmt" field follows:'
      call iplot (kmt, imt, imt, jmt)
#  else
      write (stdout,'(/a/)') 'Remove -Dskip_kmt_map to show the kmt map'
#  endif
      call area_volume (kmt, xu, yu, zw, ht)
!
!-----------------------------------------------------------------------
!     check the position of rivers
!-----------------------------------------------------------------------
!
#  ifdef river_inflow
      call check_river(kmt, xu, yu)
#  endif
c
c-----------------------------------------------------------------------
c     construct hu as the depth to bottom on U cells
c-----------------------------------------------------------------------
c
      do i=1,imt
        hu(i,jmt) = 0.0
      enddo
      do j=1,jmt-1
        do i=2,imt-1
# ifdef partial_cell
          hu(i,j) = min( ht(i,j), ht(i+1,j), ht(i,j+1), ht(i+1,j+1) )
# else
          k = kmu(i,j)
	  if (k .gt. 0) then
            hu(i,j) = zw(k)
	  else
            hu(i,j) = 0.0
	  endif
# endif
	enddo
# ifdef cyclic
        hu(1,j)   = hu(imt-1,j)
	hu(imt,j) = hu(2,j)
# else
        hu(1,j)   = 0
	hu(imt,j) = 0
# endif
      enddo
c
c-----------------------------------------------------------------------
c      analysis of topographic slopes
c-----------------------------------------------------------------------
c
      do n=1,7
        slope(n) = 0
      enddo
      c2rad = 180.0/(4.0*atan(1.0))
      c2cm  = 1.1132387e7
      do j=2,jmt-1
        cs = cos(yu(j)*c2rad)
	dy = (yu(j+1)-yu(j))*c2cm
        do i=2,imt-1
	  dx = cs*(xu(i+1)-xu(i))*c2cm
	  do n=1,2
	    if (n .eq. 1) then
              if (hu(i+1,j) .ge. zw(1) .and. hu(i,j) .ge. zw(1)) then
	        grad = abs((hu(i+1,j)-hu(i,j))/dx)
              else
	        grad = -1.0
              endif
	    else
              if (hu(i,j+1) .ge. zw(1) .and. hu(i,j) .ge. zw(1)) then
	        grad = abs((hu(i,j+1)-hu(i,j))/dy)
              else
	        grad = -1.0
              endif
	    endif
	    nbin = 0
	    if (grad .le. 0.00001) then
	      nbin = 1
	    elseif (grad .le. 0.0001) then
	      nbin = 2
	    elseif (grad .le. 0.001) then
	      nbin = 3
	    elseif (grad .le. 0.01) then
	      nbin = 4
	    elseif (grad .le. 0.1) then
	      nbin = 5
	    elseif (grad .le. 1.0) then
	      nbin = 6
	    elseif (grad .gt. 1.0) then
	      nbin = 7
	    endif
	    if (nbin .ne. 0) slope(nbin) = slope(nbin) + 1
	  enddo
	enddo
      enddo
      write (stdout,'(//,10x,a/)') "Topographic slope analysis"
      sum = 0
      do n=1,7
        sum = sum + slope(n)
      enddo
      write (stdout,'(1x,f7.3,a)') 100.0*slope(1)/sum
     &,"% of slopes are GE 0 and LE 0.00001"
      write (stdout,'(1x,f7.3,a)') 100.0*slope(2)/sum
     &,"% of slopes are GT 0.00001 and LE 0.0001"
      write (stdout,'(1x,f7.3,a)') 100.0*slope(3)/sum
     &,"% of slopes are GT 0.0001 and LE 0.001"
      write (stdout,'(1x,f7.3,a)') 100.0*slope(4)/sum
     &,"% of slopes are GT 0.001 and LE 0.01"
      write (stdout,'(1x,f7.3,a)') 100.0*slope(5)/sum
     &,"% of slopes are GT 0.01 and LE 0.1"
      write (stdout,'(1x,f7.3,a)') 100.0*slope(6)/sum
     &,"% of slopes are GT 0.1 and LE 1.0"
      write (stdout,'(1x,f7.3,a)') 100.0*slope(7)/sum
     &,"% of slopes are GT 1.0 "
c
#if defined topog_diagnostic || defined hl_diffusivity 
c
c-----------------------------------------------------------------------
c     save the depth at T cells, f/H, and kmt field as output File
c-----------------------------------------------------------------------
c
      call topog_data (kmt, ht, hu, ruff)
#endif
c
c---------------------------------------------------------------------
c     compute a topography checksum
c---------------------------------------------------------------------
c
      call print_checksumi (kmt, imt, jmt, 'Final "kmt" checksum =')
c
#ifdef write_my_kmt
c
c-----------------------------------------------------------------------
c     export the final "kmt" field 
# ifdef partial_cell
c     export "ht" as the total depth on T cells
# endif
c-----------------------------------------------------------------------
c
      write (stdout,'(/a/)')' =>Writing a "kmt" field to file "kmt.dta"'
      call shuttle_kmt (kmt, imt, jmt, km, 'export', 'kmt.dta')
      write (stdout,'(/a/)')  '==> The "kmt" field has been written.'
# ifdef partial_cell
      write (stdout,'(/a/)')' =>Writing a "ht" field to file "ht.dta" '
      call getunit (io, 'ht.dta'
     &,             'unformatted sequential rewind ieee')
c
      stamp = ' no stamp'
      iotext = 'read (io) imt, jmt, km'
      expnam = '  '
      write (io) stamp, iotext, expnam
      write (io) imt, jmt, km
      iotext = 'read (io) ((ht(i,j),i=1,imt),j=1,jmt)'
      write (io) stamp, iotext, expnam
      write (io) ht
      call relunit (io)
      call print_checksum (ht, imt, jmt,' exported ht checksum = ')
      write (stdout,'(/a/)')  '==> The "ht" field has been written.'
# endif
#endif
#if defined write_my_roughness 
c
c---------------------------------------------------------------------
c     export the estimated roughness field
c---------------------------------------------------------------------
c
      write (stdout,'(/a/)')
     & ' =>Writing a "roughness" field to file "rough.dta" '
      call getunit (io, 'rough.dta'
     &,             'unformatted sequential rewind ieee')
c
      stamp = ' no stamp'
      iotext = 'read (io) imax, jmax, kmax'
      expnam = '  '
      write (io) stamp, iotext, expnam
      write (io) imt, jmt, km
      iotext = 'read (io) ((ruff(i,j),i=1,imax),j=1,jmax)'
      write (io) stamp, iotext, expnam
      write (io) ruff
      call relunit (io)
      call print_checksum (ruff, imt, jmt,' exported ruff checksum = ')
      write (stdout,'(/a/)')'==> The roughness field has been written.'
#endif
      return
      end

      subroutine print_checksumi (kmt, im, jm, text)
# include "stdunits.h"
      dimension kmt(im,jm)
      character*(*) text
      cksum = 0.0
      do jrow=1,jm
        do i=1,im
          cksum = cksum + i*jrow*kmt(i,jrow)
        enddo
      enddo
      write (stdout,*) text, cksum
      return
      end


      subroutine topog_options
      parameter (maxifdefs=40)
      character*(30) ifdefs(maxifdefs)
#include "stdunits.h"
        write (stdout,'(/,10x,a,/,10x,a,/)')
     &  ' The following "ifdef" options have been enabled:'
     &, '  (Consult the manual for their meaning)'
        n = 0
#ifdef drive_topog
        n = n + 1
        ifdefs(n) = 'drive_topog'
#endif
#ifdef smooth_topo
        n = n + 1
        ifdefs(n) = 'smooth_topo'
#endif
#ifdef fill_isolated_dont_change_landmask
        n = n + 1
        ifdefs(n) = 'fill_isolated_dont_change_landmask'
#endif
#ifdef smooth_topog_after_user_changes
        n = n + 1
        ifdefs(n) = 'smooth_topog_after_user_changes'
#endif
#ifdef smooth_topo_allow_deepening
        n = n + 1
        ifdefs(n) = 'smooth_topo_allow_deepening'
#endif
#ifdef scripps_kmt
        n = n + 1
        ifdefs(n) = 'scripps_kmt'
#endif
#ifdef etopo_kmt
        n = n + 1
        ifdefs(n) = 'etopo_kmt'
#endif
#ifdef iowxyz_kmt
        n = n + 1
        ifdefs(n) = 'iowxyz_kmt'
#endif
#ifdef rough_mixing
        n = n + 1
        ifdefs(n) = 'rough_mixing'
#endif
#ifdef write_my_roughness
        n = n + 1
        ifdefs(n) = 'write_my_roughness'
#endif
#ifdef read_my_roughness
        n = n + 1
        ifdefs(n) = 'read_my_roughness'
#endif
#ifdef read_my_kmt
        n = n + 1
        ifdefs(n) = 'read_my_kmt'
#endif
#ifdef write_my_kmt
        n = n + 1
        ifdefs(n) = 'write_my_kmt'
#endif
#ifdef read_unformatted_kmt
        n = n + 1
        ifdefs(n) = 'read_unformatted_kmt'
#endif
#ifdef idealized_kmt
        n = n + 1
        ifdefs(n) = 'idealized_kmt'
#endif
#ifdef gaussian_kmt
        n = n + 1
        ifdefs(n) = 'gaussian_kmt'
#endif
#ifdef rectangular_box
        n = n + 1
        ifdefs(n) = 'rectangular_box'
#endif
#ifdef flat_bottom
        n = n + 1
        ifdefs(n) = 'flat_bottom'
#endif
#ifdef solid_walls
        n = n + 1
        ifdefs(n) = 'solid_walls'
#endif
#ifdef cyclic
        n = n + 1
        ifdefs(n) = 'cyclic'
#endif
#ifdef fill_shallow
        n = n + 1
        ifdefs(n) = 'fill_shallow'
#endif
#ifdef deepen_shallow
        n = n + 1
        ifdefs(n) = 'deepen_shallow'
#endif
#ifdef round_shallow
        n = n + 1
        ifdefs(n) = 'round_shallow'
#endif
#ifdef fill_isolated_cells
        n = n + 1
        ifdefs(n) = 'fill_isolated_cells'
#endif
#ifdef partial_cell
        n = n + 1
        ifdefs(n) = 'partial_cell'
#endif
#ifdef bbl_ag
        n = n + 1
        ifdefs(n) = 'bbl_ag'
#endif
      if (n .gt. maxifdefs) then
        write (stdout,*) n
     &, '=>Error: increase "maxifdefs" to ',n,' in "topog_options"'
        call abort()
      endif
      nacros = 4
      nn = n/float(nacros) + 1
      do i=1,nn
        ns = 1+nacros*(i-1)
        ne = ns + nacros-1
        if (ne .gt. n) ne = n
        if (ns .le. n) write (stdout,9101) (ifdefs(n1),n1=ns,ne)
      enddo
      write (stdout,'(/a/)')  ' END OF ENABLED OPTIONS LIST'
c
#if !defined fill_isolated_cells && !defined partial_cell
      write (stdout,*)
     & '=>Warning: option "fill_isolated_cells" is recommended'
      write (stdout,*)'           but is not enabled'
#endif
      n = 0
#if !(defined deepen_shallow || defined fill_shallow)
# define no_deepen_or_fill
#endif
#if defined round_shallow || defined no_deepen_or_fill
      n = n + 1
#endif
#if defined fill_shallow
      n = n + 1
#endif
#if defined deepen_shallow
      n = n + 1
#endif
      if (n .gt. 1) then
        write (stdout,*)
     &    '=>Error: only one of the following is permitted'
	write (stdout,*)
     &    '         round_shallow, fill_shallow, deepen_shallow'
        call abort()
      endif
      return
9101  format(4(1x,a30))
      end





      subroutine area_volume (kmt, xu, yu, zw, ht)
c
c-----------------------------------------------------------------------
c     compute surface area and volume of ocean (T cells and U cells)
c     (note that areas are defined at each level)
c-----------------------------------------------------------------------
c
#include "size.h"
#include "stdunits.h"
      dimension kmt(imt, jmt), xu(imt), yu(jmt), zw(km)
#ifdef partial_cell
      dimension ht(imt,jmt)
#endif
      dimension tcella(km)
      dimension cst(jmt)
      dimension dxt(imt), dyt(jmt)
c
      do k=1,km
        tcella(k) = 0
      enddo
      ocnp   = 0
      tcellv = 0
      pi = 4.0 * atan(1.0)
      degrad = pi / 180.0
      radius = 0.6371e9
      do i=2,imt-1
        dxt(i) = radius * degrad * (xu(i) - xu(i-1))
      end do
      do jrow = 2,jmt-1
        cst(jrow) = cos(0.5*(yu(jrow)+yu(jrow-1))*degrad)
        dyt(jrow) = radius * degrad * (yu(jrow) - yu(jrow-1))
      end do
c
      do jrow=2,jmt-1
        do i=2,imt-1
          if (kmt(i,jrow) .gt. 0) then
            do k=1,kmt(i,jrow)
              tcella(k) = tcella(k) +
     &           cst(jrow)*dxt(i)*dyt(jrow)
            enddo
            tcellv = tcellv + cst(jrow)*dxt(i)*dyt(jrow)
#ifdef partial_cell
     &                        *ht(i,jrow)
#else
     &                        *zw(kmt(i,jrow))
#endif
            ocnp   = ocnp + float(kmt(i,jrow))
          endif
        enddo
      enddo
c
      write (stdout,'(/)')
      write (stdout,'(a,f8.0)')    'number of ocean T cells =', ocnp
      write (stdout,'(a,e15.6,a)') 'surface area (T cells)  =',tcella(1)
     &,                            ' cm**2'
      write (stdout,'(a,e15.6,a)') 'ocean volume (T cells)  =', tcellv
     &,                            ' cm**3'
      write (stdout,'(/)')
c
      return
      end



      subroutine kmtbc (kmt, imt, jmt, ht)
c
c-----------------------------------------------------------------------
c     set lateral boundary conditions on kmt and ht
c     set cyclic, solid wall and symmetry conditions on T grid
c-----------------------------------------------------------------------
c
      dimension kmt(imt,jmt), ht(imt,jmt)
c
#ifdef cyclic
c
c     set cyclic conditions on eastern and western boundary
c
      do jrow=1,jmt
        kmt(1,jrow)   = kmt(imt-1,jrow)
        kmt(imt,jrow) = kmt(2,jrow)
        ht(1,jrow)    = ht(imt-1,jrow)
	ht(imt,jrow)  = ht(2,jrow)
      enddo
#else
c
c     set solid wall conditions on eastern and western boundary
c     (if "cyclic" is not enabled then option "solid_walls" is assumed)
c
      do jrow=1,jmt
# ifdef obc_west
        kmt(1,jrow) = kmt(2,jrow)
        ht(1,jrow)  = ht(2,jrow)
# else	
        kmt(1,jrow) = 0
        ht(1,jrow)  = 0
# endif
# ifdef obc_east
        kmt(imt,jrow) = kmt(imt-1,jrow)
        ht(imt,jrow)  = ht(imt-1,jrow)
# else	
        kmt(imt,jrow) = 0
        ht(imt,jrow)  = 0
# endif
      enddo
#endif
#ifdef symmetry
c
c     set symmetry conditions at U row jmt-1
c
      do i=1,imt
        kmt(i,jmt) = kmt(i,jmt-1)
        ht(i,jmt)  = ht(i,jmt-1)
      enddo
#else
# ifdef obc_north
c
c     extrapolate the domain to the northernmost row
c
      do i=1,imt
        kmt(i,jmt) = kmt(i,jmt-1)
        ht(i,jmt)  = ht(i,jmt-1)
      enddo
# else
c
c     close the domain at the northernmost row
c
      do i=1,imt
        kmt(i,jmt) = 0
        ht(i,jmt)  = 0
      enddo
# endif
#endif
# ifdef obc_south
c
c     extrapolate the domain to the southernmost row
c
      do i=1,imt
        kmt(i,1) = kmt(i,2)
        ht(i,1)  = ht(i,2)
      enddo
# else
c
c     close the domain at the southernmost row
c
      do i=1,imt
        kmt(i,1) = 0
        ht(i,1)  = 0
      enddo
# endif
      return
      end



#ifdef scripps_kmt
      subroutine scripp (kmt, xt, yt, zt, xu, yu, zw, dzt
     &,                  imt2, jmt2, km2, a)
c
c-----------------------------------------------------------------------
c     construct "kmt" from scripps 1 deg topography.
c-----------------------------------------------------------------------
c
# include "stdunits.h"
      parameter (iw=362, jw=182)
      parameter (lenw=20*iw)
# include "size.h"
# if defined partial_cell && ( defined fourfil || defined firfil )
#  if !defined driver_only
#   include "index.h"
#  endif
# endif
c
      dimension kmt(imt2,jmt2), xt(imt2), yt(jmt2), zt(km2)
      dimension xu(imt2), yu(jmt2), zw(km2), dzt(km2)
      dimension a(imt,jmt)
c
c     file 'scripps.top' has real *4 ieee format data
c
      real*4 lons4, lats4, cx4, cy4, d4

      dimension d(iw,jw), cx(iw), cy(jw), dt(iw), cxt(iw), work(lenw)
      dimension d4(iw,jw), cx4(iw), cy4(jw)
c
c     check that grid sizes in argument list and file "size.h" agree
c
      call size_check (imt2, jmt2, km2, 'scripps', 'stop')
c
c-----------------------------------------------------------------------
c     read  Scripps 1 deg topography data
c-----------------------------------------------------------------------
c
      call getunit (io_scripps, 'scripps.top'
     &,             'unformatted sequential rewind ieee')
      read (io_scripps)
#if defined _CRAY
      read (io_scripps) lons, lats, cx, cy, d      
#else
      read (io_scripps) lons4, lats4, cx4, cy4, d4
      do i=1,iw
        cx(i) = cx4(i)
      end do
      do j=1,jw
        cy(j) = cy4(j)
      end do
      do i=1,iw
        do j=1,jw
          d(i,j) = d4(i,j)
        end do
      end do
#endif
      call relunit (io_scripps)
      
      write (stdout,'(/,a,/)')
     & ' => Preparing model "kmt" from Scripps 1 deg Topography'
c
c-----------------------------------------------------------------------
c     translate Scripps grid longitudes to eliminate non-monotonic
c     coordinates (359.5, 0.5, 1.5) across prime meridian from being in
c     the interior of the model grid. (only for limited domain grids
c     that contain the prime meridian)
c-----------------------------------------------------------------------
c
      call tranlon (d, iw, iw-2, jw, dt, cx, xt(2), imt-2, cxt)
c
#ifndef skip_translation_details
      write (stdout,'(/a/)')
     & 'Use -Dskip_translation_details to not show the following'
      write(stdout,'(///,50x, a30,/)') 'Scripps 1 deg topography grid'
      write (stdout,'(//,1x,a30,/)') '            latitudes:'
      write (stdout,'(1x,10f10.5)') cy
      write (stdout,'(//,1x,a30,/)') '  original longitudes:'
      write (stdout,'(1x,10f10.5)') cx
      write (stdout,'(//,1x,a30,/)') 'translated longitudes:'
      write (stdout,'(1x,10f10.5)') cxt
#else
      write (stdout,'(/a/)')
     & 'Remove -Dskip_translation_details to show translation details'
#endif
c
c-----------------------------------------------------------------------
c     interpolate Scripps topography to MOM T grid points
c     intrp = 1 => model grid is coarser than 1 deg Scripps data
c     intrp = 2 => model grid is finer than 1 deg Scripps data
c-----------------------------------------------------------------------
c
c     USER INPUT: set "intrp" appropriately if this test fails
c
      avgdx = (xu(imt-1) - xu(1))/(imt-2)
      avgdy = (yu(jmt-1) - yu(1))/(jmt-2)
      if (avgdx .gt. 1.0 .and. avgdy .gt. 1.0) then
        intrp = 1
      else
        intrp = 2
      end if
c     END USER INPUT
      if (intrp .eq. 1) then
        write (stdout,'(/a/a/)')
     & ' "intrp=1" assumes MOM resolution is coarser than Scripps 1 deg'
     &,' if not ... then reset it in the USER INPUT section of topog.F'
        if (cy(1) .gt. yt(1)) cy(1) = yt(1)
        if (cy(jw) .lt. yt(jmt)) cy(jw) = yt(jmt)
        call ftc (d, iw, jw, cxt, cy, a, imt, jmt, 2, imt-1, 1, jmt
     &,         xt, yt, 1, work, lenw)
      else if (intrp .eq. 2) then
        write (stdout,'(/a/a/)')
     & ' "intrp=2" assumes MOM resolution is finer than Scripps 1 deg'
     &,' if not ... then reset it in the USER INPUT section of topog.F'
        call ctf (d, iw, jw, cxt, cy, a, imt, jmt, 2, imt-1, 1, jmt
     &,         xt, yt, 1, work, lenw)
      endif
c
      do jrow=1,jmt
# ifdef cyclic
        a(1,jrow)   = a(imt-1,jrow)
        a(imt,jrow) = a(2,jrow)
# else
        a(1,jrow)   = a(2,jrow)
        a(imt,jrow) = a(imt-1,jrow)
# endif
      enddo
c
c-----------------------------------------------------------------------
c     convert depths from meters to cm and change sign.
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        do i=1,imt
          if (a(i,jrow) .lt. 0.0) then
            a(i,jrow) = -a(i,jrow)*100.0
# if defined bbl_ag
            if (a(i,jrow) .gt. dzt(1)+dzt(km)) then
	      a(i,jrow) = a(i,jrow) - dzt(km)
	    endif
# endif
          else
            a(i,jrow) = 0.0
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     set topographic modifications on model grid here
c-----------------------------------------------------------------------
c
c     set artificial walls at northern & southern boundaries
c
      do i=1,imt
        a(i,1 )  = 0.0
        a(i,jmt) = 0.0
      enddo
c
# if defined smooth_topo
c
c-----------------------------------------------------------------------
c     filter topography northward of latitude "crit_lat" by repeated
c     applications of a 2D finite impulse filter.
c     USER INPUT: set crit_lat and num_pass as needed
c-----------------------------------------------------------------------
c
#ifdef test_case_A
      crit_lat = -85.0
      num_pass = 2
#else
      crit_lat = 85.0
      num_pass = 10
#endif
      j_crit = indp (crit_lat, yt, jmt)
c
      write (stdout,'(/a,f5.1,a,i3,a/)')
     & '=>Filtering topography poleward of lat='
     &,yt(j_crit),' with ',num_pass,' passes.'
c
      call filter_topo (a(1,1), imt, jmt, j_crit, jmt-1, num_pass)
# endif
c
c-----------------------------------------------------------------------
c     discretize topography "kmt" to nearest "zw" levels
c-----------------------------------------------------------------------
c
# if defined partial_cell && ( defined fourfil || defined firfil )
c
c     make sure there are no partial cells where polar filtering occurs:
c     calculate filtering rows (copied from setocn):
c
#  ifdef driver_only
c     set up model indices for filtering polar latitudes
      call set_polar_filtering_indices
#  endif
      j_ft1   = jft1
      j_ft2   = jft2
# else
      j_ft1   = 0
      j_ft2   = jmt+1
# endif

      do jrow=1,jmt
        do i=2,imt-1
          kmt(i,jrow) = 0
          if (a(i,jrow) .ne. 0.0) then
            kmt(i,jrow) = indp (a(i,jrow), zw, kbot)
# ifdef partial_cell
c
c           adjust "kmt" so that "a(i,jrow)" <= zw(kmt(i,jrow))
c
            if (jrow .gt. j_ft1 .and. jrow .lt. j_ft2) then
              if (zw(kmt(i,jrow)) .lt. a(i,jrow)) then
                if (kmt(i,jrow) .eq. kbot) then
                  a(i,jrow) = zw(kmt(i,jrow))
                else
                  kmt(i,jrow) =kmt(i,jrow) + 1
                endif
              endif
            else
              a(i,jrow) = zw(kmt(i,jrow))
            endif
# endif
          endif
        enddo
# ifdef cyclic
        kmt(1,jrow)   = kmt(imt-1,jrow)
        kmt(imt,jrow) = kmt(2,jrow)
# else
        kmt(1,jrow)   = 0
        kmt(imt,jrow) = 0
# endif
      enddo
c
      return
      end
#endif


#ifdef etopo_kmt
      subroutine etopo (kmt, xt, yt, zt, xu, yu, zw, dzt, imt2, jmt2
     &,                 km2, a, ruff)
c
c-----------------------------------------------------------------------
c     construct "kmt" from NGDC ETOPO5 1/12 deg topography.
c     This dataset is not available from GFDL.
c     The ETOPO5 dataset can be purchased from the
c     Marine Geology and Geophysics Division of the
c     National Geophysical Data Center.
c
c     author:      Bill Hurlin      e-mail=> wh@gfdl.gov
c-----------------------------------------------------------------------
c
#include "stdunits.h"
      parameter (iw=4322, jw=2160, iwm2=iw-2)
      parameter (lenw=20*iw)
#include "size.h"
# if defined partial_cell && ( defined fourfil || defined firfil )
#  if !defined driver_only
#   include "index.h"
#  endif
# endif
c
      dimension kmt(imt2,jmt2), xt(imt2), yt(jmt2), zt(km2)
      dimension xu(imt2), yu(jmt2), zw(km2), dzt(km2)
      dimension a(imt,jmt), ruff(imt,jmt)
c
c     file 'etopo' has real*4 ieee format data
c
      dimension lon(iwm2)
      dimension d(iw,jw), cx(iw), cy(jw), dt(iw), cxt(iw), work(lenw)
      dimension d_buf(iw) 
      
      real(kind=4) lon4,d_buf4(iw)
c
c     check that grid sizes in argument list and file "size.h" agree
c
      call size_check (imt2, jmt2, km2, 'etopo', 'stop')
c
c-----------------------------------------------------------------------
c     read  etopo 1/12 deg topography
c-----------------------------------------------------------------------
c
      call getunit (io_etopo, '/net/rcp/ETOPO5.NGDCunformat_ieee'
     &,             'unformatted sequential rewind ieee')
c
      twelvdg = 1.0/12.0
      cx(1)=-twelvdg
      do i = 2,iw         
        cx(i) = cx(i-1) + twelvdg         
      enddo
c
      write (stdout,*) ' Reading Etopo 1/12 deg data...'
      rlat = 90.0
      do j=1,jw
        jj = jw - j + 1
#if defined _CRAY        
        read (io_etopo) lon
#else        
        read (io_etopo) lon4
        lon=lon4
#endif        
        cy(jj) = rlat
	rlat = rlat - twelvdg
        do i=1,iwm2
          d(i+1,jj) = lon(i)
        enddo
        d(1,jj)  = d(iw-1,jj)            
        d(iw,jj) = d(2,jj)
# ifdef debug_etopo
        write(stdout,105) jj, cy(jj), (d(ii,jj),ii=1,2)
     &,              (d(ii,jj),ii=iw-1,iw)
 105    format(1x,'jj=',i5,' cy(jj)=',f8.2
     &,         ' dpt(1..2)=',2f10.4,' dpt(iw-1..iw)=',2f10.4)
# endif
      enddo
      call relunit (io_etopo)
      write (stdout,'(/,a,/)')
     & ' => Preparing model "kmt" from etopo 1/12 deg Topography'
c
c-----------------------------------------------------------------------
c     translate Etopo grid longitudes to eliminate non-monotonic
c     coordinates (360-1/12, 0, 1/12) across prime meridian from being
c     in the interior of the model grid. (only for limited domain grids
c     that contain the prime meridian)
c-----------------------------------------------------------------------
c
      call tranlon (d, iw, iw-2, jw, dt, cx, xt(2), imt-2, cxt)
c
      write(stdout,'(///,50x, a30,/)') 'Etopo 1/12 deg topography grid'
      write (stdout,'(//,1x,a30,/)') '            latitudes:'
      write (stdout,'(1x,10f10.5)') cy
      write (stdout,'(//,1x,a30,/)') '  original longitudes:'
      write (stdout,'(1x,10f10.5)') cx
      write (stdout,'(//,1x,a30,/)') 'translated longitudes:'
      write (stdout,'(1x,10f10.5)') cxt
c
c-----------------------------------------------------------------------
c     interpolate Etopo topography to MOM T grid points
c     intrp = 1 => model grid is coarser than 1 deg Etopo data
c     intrp = 2 => model grid is finer than 1 deg Etopo data
c-----------------------------------------------------------------------
c
c     USER INPUT: set "intrp" appropriately if this test fails
c
      avgdx = (xu(imt-1) - xu(1))/(imt-2)
      avgdy = (yu(jmt-1) - yu(1))/(jmt-2)
      if (avgdx .ge. twelvdg .and. avgdy .ge. twelvdg) then
        intrp = 1
      else
        intrp = 2
      end if
      intrp=1
c     END USER INPUT
      if (intrp .eq. 1) then
        write (stdout,'(/a/a/)')
     & '"intrp=1" assumes MOM resolution is coarser than Etopo 1/12 deg'
     &,'if not ... then reset it in the USER INPUT section of topog.F'
        if (cy(1) .gt. yt(1)) cy(1) = yt(1)
        if (cy(jw) .lt. yt(jmt)) cy(jw) = yt(jmt)
        call ftc (d, iw, jw, cxt, cy, a, imt, jmt, 2, imt-1, 1, jmt
     &,         xt, yt, 1, work, lenw)
      else if (intrp .eq. 2) then
        write (stdout,'(/a/a/)')
     & '"intrp=2" assumes MOM resolution is finer than Etopo 1/12 deg'
     &,' if not ... then reset it in the USER INPUT section of topog.F'
        call ctf (d, iw, jw, cxt, cy, a, imt, jmt, 2, imt-1, 1, jmt
     &,         xt, yt, 1, work, lenw)
      endif
#if defined rough_mixing
c
c-----------------------------------------------------------------------
c     calculate roughness (Questions => Ron Pacanowski rcp@gfdl.gov)
c     note: Array "d" is now roughness
c-----------------------------------------------------------------------
c
      call getunit (io_ruff, 'ruff.dta'
     &,             'unformatted sequential rewind ieee')
      ddxr= 1.0e2/(twelvdg*1.1132387e7)**2
      ddyr = ddxr
      do j=1,jw
        jm1 = max(j-1,1)
	jp1 = min(j+1,jw)
        do i=2,iwm2
	  ddx = (d(i+1,j)-2*d(i,j)+d(i-1,j))*ddxr
	  ddy = (d(i,jp1)-2*d(i,j)+d(i,jm1))*ddyr
	  d_buf(i) = sqrt(ddx**2 + ddy**2)
	enddo
	d_buf(1) = d_buf(iw-1)
	d_buf(iw) = d_buf(2)
	write (io_ruff) d_buf
      enddo
      rewind io_ruff
      do j=1,jw
#if defined _CRAY         
        read(io_ruff) d_buf
#else        
        read(io_ruff) d_buf4
	do i=1,iw
	  d(i,j) = d_buf4(i)
       enddo
#endif       
      enddo
      call relunit (io_ruff)
#else
      do j=1,jw
	do i=1,iw
	  d(i,j) = 0.0
	enddo
      enddo
#endif
c
      if (intrp .eq. 1) then
        write (stdout,'(/a/a/)')
     & '"intrp=1" assumes MOM resolution is coarser than Etopo 1/12 deg'
     &,'if not ... then reset it in the USER INPUT section of topog.F'
        call ftc (d, iw, jw, cxt, cy, ruff, imt, jmt, 2, imt-1, 1, jmt
     &,         xt, yt, 1, work, lenw)
      else if (intrp .eq. 2) then
        write (stdout,'(/a/a/)')
     & '"intrp=2" assumes MOM resolution is finer than Etopo 1/12 deg'
     &,' if not ... then reset it in the USER INPUT section of topog.F'
        call ctf (d, iw, jw, cxt, cy, ruff, imt, jmt, 2, imt-1, 1, jmt
     &,         xt, yt, 1, work, lenw)
      endif
c
c-----------------------------------------------------------------------
c     Apply boundary conditions
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
#ifdef cyclic
        a(1,jrow)   = a(imt-1,jrow)
        a(imt,jrow) = a(2,jrow)
#else
        a(1,jrow)   = a(2,jrow)
        a(imt,jrow) = a(imt-1,jrow)
#endif
      enddo
c
c-----------------------------------------------------------------------
c     convert depths from meters to cm and change sign.
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        do i=1,imt
          if (a(i,jrow) .lt. 0.0) then
            a(i,jrow) = -a(i,jrow)*100.0
# if defined bbl_ag
            if (a(i,jrow) .gt. dzt(1)+dzt(km)) then
	      a(i,jrow) = a(i,jrow) - dzt(km)
	    endif
# endif
          else
            a(i,jrow) = 0.0
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     set topographic modifications on model grid here
c-----------------------------------------------------------------------
c
c     set artificial walls at northern & southern boundaries
c
      do i=1,imt
        a(i,1 )  = 0.0
        a(i,jmt) = 0.0
      enddo
c
# if defined smooth_topo
c
c-----------------------------------------------------------------------
c     filter topography northward of latitude "crit_lat" by repeated
c     applications of a 2D finite impulse filter.
c     USER INPUT: set crit_lat and num_pass as needed
c-----------------------------------------------------------------------
c
      crit_lat = 85.0
      num_pass = 10
      j_crit = indp (crit_lat, yt, jmt)
c
      write (stdout,'(/a,f5.1,a,i3,a/)')
     & '=>Filtering topography poleward of lat='
     &,yt(j_crit),' with ',num_pass,' passes.'
c
      call filter_topo (a(1,1), imt, jmt, j_crit, jmt-1, num_pass)
# endif
c
c-----------------------------------------------------------------------
c     discretize topography "kmt" to nearest "zw" levels
c-----------------------------------------------------------------------
c
# if defined partial_cell && ( defined fourfil || defined firfil )
c
c     make sure there are no partial cells where polar filtering occurs:
c     calculate filtering rows (copied from setocn):
c
#  ifdef driver_only
c     set up model indices for filtering polar latitudes
      call set_polar_filtering_indices
#  endif
      j_ft1   = jft1
      j_ft2   = jft2
# else
      j_ft1   = 0
      j_ft2   = jmt+1

# endif

      do jrow=1,jmt
        do i=2,imt-1
          kmt(i,jrow) = 0
          if (a(i,jrow) .ne. 0.0) then
            kmt(i,jrow) = indp (a(i,jrow), zw, kbot)
# ifdef partial_cell
c
c           keep physical bottom "a(i,jrow)" <= zw(kmt(i,jrow))
c
            if (jrow .gt. j_ft1 .and. jrow .lt. j_ft2) then
              if (zw(kmt(i,jrow)) .lt. a(i,jrow)) then
                if (kmt(i,jrow) .eq. kbot) then
                  a(i,jrow) = zw(kmt(i,jrow))
                else
                  kmt(i,jrow) =kmt(i,jrow) + 1
                endif
              endif
            else
              a(i,jrow) = zw(kmt(i,jrow))
            endif
# endif
          endif
        enddo
#ifdef cyclic
        kmt(1,jrow)   = kmt(imt-1,jrow)
        kmt(imt,jrow) = kmt(2,jrow)
#else
        kmt(1,jrow)   = 0
        kmt(imt,jrow) = 0
#endif
      enddo
c
      return
      end

#endif


#ifdef idealized_kmt
      subroutine idealized (kmt, xt, yt, zw, imt, jmt, km, ht)
c
c-----------------------------------------------------------------------
c     construct a highly "idealized" world ... piece by piece
c
c     note: the purpose of this geometry/topography is to automatically
c           map into arbitrary resolution as grid dimensions "imt" and
c           "jmt" are changed, thereby facilitating the implementation
c           and verification of the model on various computer platforms
c           without referencing the topographic data base.  Although it
c           somewhat resembles the real world, it is NOT realistic.
c
c     author: r. c. pacanowski            e-mail==>rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      dimension kmt(imt,jmt), xt(imt), yt(jmt), zw(km)
# ifdef partial_cell
      dimension ht(imt,jmt)
# endif
c
      do jrow=2,jmt-1
        do i=2,imt-1
          kmt(i,jrow) = km
        enddo
      enddo
c
      do jrow=1,jmt
        kmt(imt,jrow) = kmt(2,jrow)
        kmt(1,jrow)   = kmt(imt-1,jrow)
      enddo
      do i=1,imt
        kmt(i,jmt) = 0
        kmt(i,1)   = 0
      enddo
c
c     antarctica
c
      call setkmt (kmt, xt, yt, imt, jmt, -90.0, 0.0, 360.0, -80.0, 0.0
     &,            360.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -80.0, 360.0-25.0, 360.0
     &,            -70.0, 360.0, 360.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -80.0, 0.0, 360.0, -70.0, 0.0
     &,            170.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -80.0, 360.0-135.0, 360.0-60.0
     &,            -68.0, 360.0-75.0, 360.0-60.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -70.0, 0.0, 155.0, -67.0, 50.0
     &,            145.0, 0)
c
c     australia
c
      call setkmt (kmt, xt, yt, imt, jmt, -35.0, 116.0, 120.0, -31.0
     &,            114.0, 130.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -38.0, 140.0, 151.0, -31.0
     &,            130.0, 151.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -31.0, 115.0, 153.0, -20.0
     &,            113.0, 149.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -20.0, 113.0, 149.0, -11.0
     &,            131.0, 143.0, 0)
c
c     south america
c
      call setkmt (kmt, xt, yt, imt, jmt, -50.0, 360.0-74.0, 360.0-68.0
     &,            -40.0, 360.0-73.0, 360.0-62.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -40.0, 360.0-73.0, 360.0-62.0
     &,            -20.0, 360.0-70.0, 360.0-40.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -20.0, 360.0-70.0, 360.0-40.0
     &,            -16.0, 360.0-81.0, 360.0-35.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, -16.0, 360.0-81.0, 360.0-35.0
     &,            0.0, 360.0-80.0, 360.0-50.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 0.0, 360.0-80.0, 360.0-50.0
     &,            11.0, 360.0-75.0, 360.0-60.0, 0)
c
c     central america
c
      call setkmt (kmt, xt, yt, imt, jmt, 6.0, 360.0-78.0, 360.0-75.0
     &,            20.0, 360.0-105.0, 360.0-97.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 20.0, 360.0-105.0, 360.0-97.0
     &,            30.0, 360.0-115.0, 360.0-94.0, 0)
c
c     north america
c
      call setkmt (kmt, xt, yt, imt, jmt, 25.0, 360.0-82.0, 360.0-80.0
     &,            30.0, 360.0-85.0, 360.0-81.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 30.0, 360.0-115.0, 360.0-80.0
     &,            40.0, 360.0-124.0, 360.0-74.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 40.0, 360.0-124.0, 360.0-74.0
     &,            50.0, 360.0-124.0, 360.0-57.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 50.0, 360.0-124.0, 360.0-57.0
     &,            60.0, 360.0-140.0, 360.0-64.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 60.0, 360.0-165.0, 360.0-64.0
     &,            65.0, 360.0-140.0, 360.0-64.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 65.0, 360.0-140.0, 360.0-64.0
     &,            70.0, 360.0-162.0, 360.0-72.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 70.0, 360.0-162.0, 360.0-140.0
     &,            72.0, 360.0-157.0, 360.0-157.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 70.0, 360.0-130.0, 360.0-70.0
     &,            75.0, 360.0-120.0, 360.0-80.0, 0)
c
c     greenland
c
      call setkmt (kmt, xt, yt, imt, jmt, 60.0, 360.0-45.0, 360.0-45.0
     &,            75.0, 360.0-58.0, 360.0-19.0, 0)
c
c     africa
c
      call setkmt (kmt, xt, yt, imt, jmt, -35.0, 19.0, 28.0, 6.0, 8.0
     &,            50.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 6.0, 0.0, 50.0, 18.0, 0.0
     &,            56.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 18.0, 0.0, 56.0, 26.0, 0.0
     &,            59.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 6.0, 360.0-10.0, 360.0, 18.0
     &,            360.0-18.0, 360.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 18.0, 360.0-18.0, 360.0, 26.0
     &,            360.0-15.0, 360.0, 0)
c
c     northern africa &  europe & asia
c
      call setkmt (kmt, xt, yt, imt, jmt, 26.0, 360.0-15.0, 360.0, 40.0
     &,            360.0-7.0, 360.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 40.0, 360.0-7.0, 360.0, 50.0
     &,            360.0, 360.0, 0)
c
      call setkmt (kmt, xt, yt, imt, jmt, 8.0, 77.0, 78.0, 26.0, 65.0
     &,            90.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 4.0, 99.0, 100.0, 26.0, 90.0
     &,            115.0, 0)
c
      call setkmt (kmt, xt, yt, imt, jmt, 26.0, 0.0, 126.0, 40.0, 0.0
     &,            122.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 40.0, 0.0, 130.0, 50.0, 0.0
     &,            140.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 50.0, 0.0, 140.0, 60.0, 8.0
     &,            140.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 60.0, 8.0, 163.0, 65.0, 13.0
     &,            180.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 65.0, 13.0, 188.0, 70.0, 20.0
     &,            180.0, 0)
      call setkmt (kmt, xt, yt, imt, jmt, 70.0, 70.0, 180.0, 75.0, 90.0
     &,            100.0, 0)
c
c     add an "idealized" undulating topography
c
      bot = zw(km)
      pi  = 4.0*atan(1.0)
      do jrow=2,jmt
        do i=2,imt-1
          if (kmt(i,jrow) .ne. 0) then
            arg = bot*(1-0.4*abs(cos((jrow*pi)/jmt)*sin((i*2*pi)/imt)))
            kmt(i,jrow) = indp (arg, zw, km)
          endif
        enddo
      enddo
c
c     add "idealized" ridges
c
      level = indp (0.666*zw(km), zw, km)
c
      call setkmt (kmt, xt, yt, imt, jmt, -20.0, 360.0-20.0, 360.0-10.0
     &,            30.0, 360.0-45.0, 360.0-35.0, level)
      call setkmt (kmt, xt, yt, imt, jmt, 30.0, 360.0-45.0, 360.0-35.0
     &,            60.0, 360.0-20.0,  360.0-30.0, level)
      call setkmt (kmt, xt, yt, imt, jmt, -60.0,360.0-100.0, 360.0-130.0
     &,            40.0, 360.0-160.0, 180.0, level)
c
      level = indp (0.5*zw(km), zw, km)
c
      call setkmt (kmt, xt, yt, imt, jmt, -50.0, 360.0-120.0
     &,            360.0-120.0, 30.0, 190.0, 190.0, level)
c
# ifdef partial_cell
c
c     arbitrarily set partial bottom cells to full depth.
c
      do jrow=1,jmt
        do i=1,imt
	  if (kmt(i,jrow) .ne. 0) then
	    ht(i,jrow) = zw(kmt(i,jrow))
	  else
	    ht(i,jrow) = 0.0
	  endif
	enddo
      enddo
# endif
      return
      end
#endif



      subroutine setkmt (kmt, xt, yt, imt, jmt, alat1, slon1, elon1
     &,                  alat2, slon2, elon2, num)
c
c-----------------------------------------------------------------------
c     set the topography mask "kmt(i,j)" = "num" within the area of
c     the trapezoid bounded by vertices:
c     (alat1,slon1), (alat1,elon1), (alat2,slon2), & (alat2,elon2)
c
c     inputs:
c
c     xt = longitudes of T points in degrees
c     yt = latitudes of T points in degrees
c     imt = number of model longitudes
c     jmt = number of model latitudes
c     alat1 = southern latitude of trapezoid (degrees)
c     slon1 = starting longitude of southern edge of trapezoid (deg)
c     elon1 = ending longitude of southern edge of trapezoid (deg)
c     alat2 = northern latitude of trapezoid (degrees)
c     slon2 = starting longitude of northern edge of trapezoid (deg)
c     elon2 = ending longitude of northern edge of trapezoid (deg)
c     num   = number of vertical levels
c
c     outputs:
c
c     kmt   = mask of vertical levels on model T points
c
c     author: r. c. pacanowski            e-mail==>rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      dimension kmt(imt,jmt), xt(imt), yt(jmt)
c
c     convert the four vertices into model indices
c     (js,is1), (js,ie1), (je,is2), (je,ie2)
c
      j1 = indp (alat1, yt, jmt)
      j2 = indp (alat2, yt, jmt)
      js = min (j1,j2)
      je = max (j1,j2)
c
      i1  = indp (slon1, xt, imt)
      i2  = indp (elon1, xt, imt)
      is1 = min (i1,i2)
      ie1 = max (i1,i2)
c
      i1  = indp (slon2, xt, imt)
      i2  = indp (elon2, xt, imt)
      is2 = min (i1,i2)
      ie2 = max (i1,i2)
c
      is = is1
      ie = ie1
c
c     fill in the area bounded by (js,is1), (js,ie1), (je,is2), (je,ie2)
c     the nudging of 1.e-5 is to insure that the test case resolution
c     generates the same topography and geometry on various computers.
c
      c1 = 1.0
      if (js .eq. je) then
        rdj = c1
      else
        rdj = c1/(je-js)
      endif
      do jrow=js,je
        do i=is,ie
          kmt(i,jrow) = num
        enddo
        is = nint(rdj*((jrow-js)*is2 + (je-jrow)*is1) + 1.0e-5)
        ie = nint(rdj*((jrow-js)*ie2 + (je-jrow)*ie1) + 1.0e-5)
      enddo
      return
      end


      subroutine min_depth (kmt, kmt_min, xu, yu, zw, ht)
c
c     limit the minimum number of levels. kmt_min should be >= 2
c
#include "size.h"
#include "stdunits.h"
c
      dimension kmt(imt,jmt), ht(imt,jmt)
      dimension xu(imt), yu(jmt), zw(km)
      write (stdout,'(/a/)')
     & 'Searching for and correcting minimum level violations'
      n = 0
      do i=2,imt-1
        do jrow=jmt-1,2,-1
          if (kmt(i,jrow) .ne. 0 .and. kmt(i,jrow) .lt. kmt_min) then
            n = n + 1
	    kmt_shallow = kmt(i,jrow)
#if !(defined deepen_shallow || defined fill_shallow)
# define no_deepen_or_fill
#endif
#if defined round_shallow || defined no_deepen_or_fill
            if (zw(kmt(i,jrow)) .lt. 0.5*zw(kmt_min)) then
              kmt(i,jrow) = 0
              ht(i,jrow)  = 0.0
            else
              kmt(i,jrow) = kmt_min
              ht(i,jrow)  = zw(kmt_min)
            end if
#endif
#if defined fill_shallow
            kmt(i,jrow) = 0
            ht(i,jrow)  = 0.0
#endif
#if defined deepen_shallow
            kmt(i,jrow) = kmt_min
            ht(i,jrow)  = zw(kmt_min)
#endif
#ifdef debug_topog
            write (stdout,'(a,i4,a,i4,a,i4,a,i4)')
     &      'a) changed shallow kmt(',i,',',jrow,') from ', kmt_shallow
     &,     ' to ', kmt(i,jrow)
#endif
          endif
        enddo
      enddo
      if (n .gt. 0) then
        write (stdout,'(a,i5,a/)')  '->Modified', n,' shallow cells'
      else
        write (stdout,'(a/)') '->No modifications needed'
      endif
      return
      end

#if defined smooth_topo
      subroutine filter_topo (d, im, jm, js, je, num_pass)
c
c=======================================================================
c     smooth topographic depth "d" with "num_pass" applications of a 2D
c     version of a (1/4, 1/2, 1/4) filter. 
c     allow filtering to decrease the bottom depth but not increase it.
c     do not allow original geometry to change.
c     note: depth "d" should be on a grid of uniformly constant spacing
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      dimension d(im,jm)
      dimension rmask(im,js-1:je+1), f(-1:1,-1:1), s(im,js:je)
c
c     2D symmetric filter weights
c
      f(-1,-1) = 1.0/16.0
      f( 0,-1) = 1.0/8.0
      f( 1,-1) = 1.0/16.0
      f(-1, 0) = 1.0/8.0
      f( 0, 0) = 1.0/4.0
      f( 1, 0) = 1.0/8.0
      f(-1, 1) = 1.0/16.0
      f( 0, 1) = 1.0/8.0
      f( 1, 1) = 1.0/16.0
c
c     1D symmetric filter weights
c
c      f( 0,-1) = 1.0/4.0
c      f( 0, 0) = 1.0/2.0
c      f( 0, 1) = 1.0/4.0
c
c     geometry mask
c
      do j=js-1,je+1
        do i=1,im
	  if (d(i,j) .eq. 0.0) then
	    rmask(i,j) = 0.0
	  else
	    rmask(i,j) = 1.0
	  endif
	enddo
      enddo
c
      do n=1,num_pass
        do j=js,je
	  do i=2,im-1
	    s(i,j) = 0.0
	    d_old  = d(i,j)
	    do ip=-1,1
	      do jp=-1,1
	        if (rmask(i+ip,j+jp) .eq. 0.0) then
	          s(i,j) = s(i,j) + d(i,j)*f(ip,jp)
		else
	          s(i,j) = s(i,j) + d(i+ip,j+jp)*f(ip,jp)
		endif
	      enddo
	    enddo
#if !defined smooth_topo_allow_deepening
	    if (s(i,j) .gt. d_old) then
	      s(i,j) = d_old
	    endif
#endif
	  enddo
c
          do i=2,im-1
	    s(i,j) = s(i,j)*rmask(i,j)
	  enddo
# if defined cyclic
	  s(1,j)  = s(im-1,j)
	  s(im,j) = s(2,j)
# else
	  s(1,j)  = d(1,j)
	  s(im,j) = d(im,j)
# endif
	enddo
c
        do j=js,je
	  do i=1,im
	    d(i,j) = s(i,j)
	  enddo
	enddo
      enddo
c
      return
      end
#endif

#if defined topog_diagnostic || defined hl_diffusivity 
      subroutine topog_data (kmt, htp, hu, ruff)

c-----------------------------------------------------------------------
c     Saves depth over T points, f/H, and real(kmt)
c
c     author: r. c. pacanowski            e-mail==>rcp@gfdl.gov
c
c     updates to mpp_io made August 1999 by 
c                  s. m. griffies        e-mail=> smg@gfdl.gov
c     Note: no domains used here; only written for single processor.
c-----------------------------------------------------------------------
      use mom3_mpp_mod
      use mpp_io_mod
# include "param.h"
# include "coord.h"
# include "diag.h"      
# include "vers.h"
# if defined hl_diffusivity
#  include "isopyc.h"
# else 
      dimension betaeff(imt,jmt)
# endif
      character*120 file_name, gvtitle
      dimension bufzt(km)

      dimension buf_kmt(imt,jmt),buf_kmu(imt,jmt)
      dimension buf_ruff(imt,jmt)  
      dimension buf_betaeff(imt,jmt),buf_ht(imt,jmt)
# ifdef partial_cell
      dimension buf_hup(imt,jmt),buf_htp(imt,jmt)      
      dimension buf_dht(imt,jmt,km),buf_dhu(imt,jmt,km)      
# endif
c-----------------------------------------------------------------------
c     field types
c-----------------------------------------------------------------------
     
      type(fieldtype), save :: field_ruff
      type(fieldtype), save :: field_ht
      type(fieldtype), save :: field_foverht
      type(fieldtype), save :: field_kmt
      type(fieldtype), save :: field_kmu
      type(fieldtype), save :: field_betaeff

      type(fieldtype), save :: field_dxt
      type(fieldtype), save :: field_dyt
      type(fieldtype), save :: field_dzt
      type(fieldtype), save :: field_dxu
      type(fieldtype), save :: field_dyu
      type(fieldtype), save :: field_dzw

# ifdef partial_cell
      type(fieldtype), save :: field_htp
      type(fieldtype), save :: field_hup
      type(fieldtype), save :: field_foverhtp 
      type(fieldtype), save :: field_dht
      type(fieldtype), save :: field_dhu

      dimension htp(imt,jmt), foverhtp(imt,jmt), hu(imt,jmt)
      dimension dht(imt,km), dhu(imt,km)
# endif
      dimension ruff(imt,jmt)
      dimension ht(imt,jmt), foverht(imt,jmt)
      dimension kmt(imt,jmt), kkmu(imt,jmt)
      dimension dxt(imt), dyt(jmt), dxu(imt), dyu(jmt)
      dimension zw0(0:km), dzw0(0:km)
      dimension htbeta(imt,jmt)

      type (axistype), save :: axis_xt
      type (axistype), save :: axis_yt
      type (axistype), save :: axis_zt
      type (axistype), save :: axis_xu
      type (axistype), save :: axis_yu
      type (axistype), save :: axis_zw

      file_name = 'topog.dta'
      period = 0.0
      gvtitle = 'Topography and geometry' // momver
      bufzt(:)  = zt(:)*0.01
      zw0(1:km) = zw(1:km)*0.01      
      zw0(0)    = 0.0

      call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')

c-----------------------------------------------------------------------
c     axes 
c-----------------------------------------------------------------------

      ibeg = 1
      iend = imt
      jbeg = 1
      jend = jmt

      call mpp_write_meta(
     &       io, axis_xt, 'xt_i', 'degrees_E', 'Longitude of T points',
     &       cartesian='X', data=xt(ibeg:iend))

      call mpp_write_meta(
     &         io, axis_yt, 'yt_j', 'degrees_N', 'Latitude of T points',
     &         cartesian='Y', data=yt(jbeg:jend))

      call mpp_write_meta(
     &         io, axis_zt, 'zt_k', 'm', 'Depth of T grid point',
     &         cartesian='Z', sense=-1, data=bufzt)
 
      call mpp_write_meta(
     &       io, axis_xu, 'xu_i', 'degrees_E', 'Longitude of U points',
     &       cartesian='X',  data=xu(ibeg:iend))

      call mpp_write_meta(
     &         io, axis_yu, 'yu_j', 'degrees_N', 'Latitude of U points',
     &         cartesian='Y', data=yu(jbeg:jend))

      call mpp_write_meta(
     &         io, axis_zw, 'zw_k', 'm', 'Depth of T cell bottom',
     &         cartesian='Z', sense=-1, data=zw0)

c-----------------------------------------------------------------------
c       data attributes 
c-----------------------------------------------------------------------

      do i=1,imt
        dxt(i) = dxtdeg(i)*1.1132387e7
        dxu(i) = dxudeg(i)*1.1132387e7
      enddo
      do j=1,jmt
        dyt(j) = dytdeg(j)*1.1132387e7
        dyu(j) = dyudeg(j)*1.1132387e7
      enddo
      dzw0(0) = 0.5*dzt(1)
      do k=1,km
        dzw0(k) = dzw(k)
      enddo

      call mpp_write_meta(
     &      io, field_dxt, (/axis_xt/), 
     &      'dxt', 'cm', 'T-Cell delta x',
     &      0.0, 1.e12, -1.0E+34, +1.0E+34)

      call mpp_write_meta(
     &      io, field_dyt, (/axis_yt/), 
     &      'dyt', 'cm', 'T-Cell delta y',
     &      0.0, 1.e12, -1.0E+34, +1.0E+34)

      call mpp_write_meta(
     &      io, field_dzt, (/axis_zt/), 
     &      'dzt', 'cm', 'Full T-Cell thickness',
     &      0.0, 1.e12, -1.0E+34, +1.0E+34)

      call mpp_write_meta(
     &      io, field_dxu, (/axis_xu/), 
     &      'dxu', 'cm', 'U-Cell delta x',
     &      0.0, 1.e12, -1.0E+34, +1.0E+34)

      call mpp_write_meta(
     &      io, field_dyu, (/axis_yu/), 
     &      'dyu', 'cm', 'U-Cell delta y',
     &      0.0, 1.e12, -1.0E+34, +1.0E+34)

      call mpp_write_meta(
     &      io, field_dzw, (/axis_zw/), 
     &      'dzw', 'cm', 'Full T-Cell thickness',
     &      0.0, 1.e6, -1.0E+34, +1.0E+34)

      call mpp_write_meta(
     &      io, field_betaeff, (/axis_xt,axis_yt/), 
     &      'betaeff', '1/cm/sec', 'H|grad(f/H)|',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

      call mpp_write_meta(
     &      io, field_ht, (/axis_xt,axis_yt/), 
     &      'ht', 'cm', 'H(kmt(i,jrow))',
     &      0.0, 1.e6, -1.0E+34, +1.0E+34)

      call mpp_write_meta(
     &      io, field_ruff, (/axis_xt,axis_yt/), 
     &      'ruff', '1/cm', 'ruffness',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

      call mpp_write_meta(
     &      io, field_foverht, (/axis_xt,axis_yt/), 
     &      'foverht', '1/cm/sec', 'f/H(kmt(i,jrow))',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

      call mpp_write_meta(
     &      io, field_kmt, (/axis_xt,axis_yt/), 
     &      'kmt', 'cells', 'kmt(i,jrow)',
     &      0.0, 1.e6, -1.0E+34, +1.0E+34)

      call mpp_write_meta(
     &      io, field_kmu, (/axis_xu,axis_yu/), 
     &      'kmu', 'cells', 'kmu(i,jrow)',
     &      0.0, 1.e6, -1.0E+34, +1.0E+34)

# ifdef partial_cell
      call mpp_write_meta(
     &      io, field_htp, (/axis_xt,axis_yt/), 
     &      'htp', 'cm', 'Depth of T-cell Partial bottom',
     &      0.0, 1.e6, -1.0E+34, +1.0E+34)

      call mpp_write_meta(
     &      io, field_hup, (/axis_xu,axis_yu/), 
     &      'hup', 'cm', 'Depth of U-cell Partial bottom',
     &      0.0, 1.e6, -1.0E+34, +1.0E+34)

      call mpp_write_meta(
     &      io, field_foverhtp, (/axis_xu,axis_yu/), 
     &      'foverhtp', '1/cm/sec', 'f/htp(kmt(i,jrow))',
     &      -1.e6, 1.e6, -1.0E+34, +1.0E+34)

      call mpp_write_meta(
     &      io, field_dht, (/axis_xt,axis_yt,axis_zt/), 
     &      'dht', 'cm', 'T Cell thickness',
     &      0.0, 1.e6, -1.0E+34, +1.0E+34)

      call mpp_write_meta(
     &      io, field_dhu, (/axis_xu,axis_yu,axis_zt/), 
     &      'dhu', 'cm', 'U Cell thickness',
     &      0.0, 1.e6, -1.0E+34, +1.0E+34)
# endif
c
c     figure title 
c
      call mpp_write_meta(io, 'title', cval=gvtitle )

c-----------------------------------------------------------------------
c       write axes now that all meta data has been written
c-----------------------------------------------------------------------
      
      call mpp_write (io,axis_xt)
      call mpp_write (io,axis_yt)
      call mpp_write (io,axis_zt)
      call mpp_write (io,axis_xu)
      call mpp_write (io,axis_yu)
      call mpp_write (io,axis_zw)

c-----------------------------------------------------------------------
c     compute the kmu field
c-----------------------------------------------------------------------
      
      do j=1,jmt
        jp1 = min(j+1,jmt)
        do i=2,imt-1
	  kkmu(i,j) = min( kmt(i,j),kmt(i+1,j),kmt(i,jp1),kmt(i+1,jp1))
	enddo
#  ifdef cyclic
        kkmu(1,j)   = kkmu(imt-1,j)
	kkmu(imt,j) = kkmu(2,j)
#  else
        kkmu(1,j)   = 0
	kkmu(imt,j) = 0
#  endif
      enddo
   
c-----------------------------------------------------------------------
c     Compute H and f/H
c-----------------------------------------------------------------------
  
      pi      = 4.0*atan(1.0)
      omega   = pi/43082.0
      crad    = pi/180.0
      do j=1,jmt
        f = 2.0*omega*sin(yt(j)*crad)
        do i=1,imt
	  if (kmt(i,j) .eq. 0) then
	    ht(i,j) =  0.0
	    foverht(i,j) = 0.0
# ifdef partial_cell
	    foverhtp(i,j) = 0.0
# endif
	  else
	    ht(i,j) = zw(kmt(i,j))
	    foverht(i,j) = f/ht(i,j)
# ifdef partial_cell
	    foverhtp(i,j) = f/htp(i,j)
# endif
	  endif
        enddo
      enddo

c-----------------------------------------------------------------------
c     Compute betaeff
c     Use a filtered topography if employing hl_diffusivity
c-----------------------------------------------------------------------

      do i = 1,imt
        do j = 1,jmt
# if defined partial_cell
            htbeta(i,j) = htp(i,j)
# else
            htbeta(i,j) = ht(i,j)
# endif
        enddo
      enddo
# if defined hl_diffusivity
      numfltrtopog = 1
      call fir2d(htbeta,numfltrtopog,imt,jmt)
# endif
c
      radius  = 6371.0e5
      do j=1,jmt
        f    = 2.0*omega*sin(yt(j)*crad)
        beta = 2.0*omega*cos(yt(j)*crad)/radius
        jm1  = max(1,j-1)
        jp1  = min(jmt,j+1)
        do i=1,imt
          im1  = max(1,i-1)
          ip1  = min(imt,i+1)
	  if (kmt(i,j) .eq. 0) then
            betaeff(i,j) = 0.0
	  else
            hx=0.5*((htbeta(ip1,j)-htbeta(i,j))/dxu(ip1)
     &             +(htbeta(i,j)-htbeta(im1,j))/dxu(im1))
     &             /cos(yt(j)*crad)
            hy=0.5*((htbeta(i,jp1)-htbeta(i,j))/dyu(jp1)
     &             +(htbeta(i,j)-htbeta(im1,j))/dyu(jm1))
            betaeff(i,j) = sqrt((beta - f*hy/htbeta(i,j))**2
     &             + (f*hx/htbeta(i,j))**2)
	  endif
        enddo
        call setbcx(betaeff(1,j),imt,1)
      enddo

c-----------------------------------------------------------------------
c     add missing point values 
c-----------------------------------------------------------------------

      buf_kmt(:,:)     = kmt(:,:)
      buf_kmu(:,:)     = kkmu(:,:)      
      buf_ruff(:,:)    = ruff(:,:)
      buf_betaeff(:,:) = betaeff(:,:)
      buf_ht(:,:)      = ht(:,:)
# ifdef partial_cell
      buf_hup(:,:)     = hu(:,:)
      buf_htp(:,:)     = htp(:,:)
# endif
      do j=1,jmt
        do i=1,imt
          if (kmt(i,j) .eq. 0) then
             buf_kmt(i,j)     = field_kmt%missing
             buf_ruff(i,j)    = field_ruff%missing
             buf_betaeff(i,j) = field_betaeff%missing
             buf_ht(i,j)      = field_ht%missing
             foverht(i,j)     = field_foverht%missing 
# ifdef partial_cell
             buf_htp(i,j)     = field_htp%missing
             foverhtp(i,j)    = field_foverhtp%missing
# endif         
           endif
           if (kkmu(i,j) .eq. 0) then
            buf_kmu(i,j)      = field_kmu%missing
# ifdef partial_cell
            buf_hup(i,j)      = field_hup%missing
# endif
          endif 
        enddo
      enddo   

c-----------------------------------------------------------------------
c     write out the diagnostics 
c-----------------------------------------------------------------------

      call mpp_write(io, field_dxt, dxt)
      call mpp_write(io, field_dyt, dyt)
      call mpp_write(io, field_dzt, dzt)
      call mpp_write(io, field_dxu, dxu)
      call mpp_write(io, field_dyu, dyu)
      call mpp_write(io, field_dzw, dzw)

      call mpp_write(io, field_betaeff, buf_betaeff)
      call mpp_write(io, field_ht,  buf_ht)
      call mpp_write(io, field_ruff, buf_ruff)
      call mpp_write(io, field_foverht, foverht)
      call mpp_write(io, field_kmt, buf_kmt)
      call mpp_write(io, field_kmu, buf_kmu)

# ifdef partial_cell
      call mpp_write(io, field_htp, buf_htp)
      call mpp_write(io, field_hup, buf_hup)
      call mpp_write(io, field_foverhtp, foverhtp)
      do jrow=1,jmt
	do i=2,imt-1
	  do k=1,km
	    dht(i,k) = dzt(k)
	    dhu(i,k) = dzt(k)
          enddo
	  kz = kmt(i,jrow)
          if (kz .gt. 1) then
            dht(i,kz) = htp(i,jrow) - zw(kz-1)
          endif
          if (kz .ne. km) then
            do k=kz+1,km
              dht(i,k) = 0.0
            enddo
	  endif
c
          jp1 = min(jrow+1,jmt)          
	  kz = min(kmt(i,jrow), kmt(i+1,jrow), kmt(i,jp1)
     &,            kmt(i+1,jp1))
          if (kz .gt. 1) then
            dhu(i,kz) = min (htp(i,jrow), htp(i+1,jrow)
     &,                        htp(i,jp1), htp(i+1,jp1)) - zw(kz-1)
          endif
          if (kz .ne. km) then
            do k=kz+1,km
              dhu(i,k) = 0.0
            enddo
	  endif
        enddo
	do k=1,km
#  ifdef cyclic
	  dht(1,k)   = dht(imt-1,k)
	  dht(imt,k) = dht(2,k)
	  dhu(1,k)   = dhu(imt-1,k)
	  dhu(imt,k) = dhu(2,k)
#  else
	  dht(1,k)   = dzt(k)
	  dht(imt,k) = dzt(k)
	  dhu(1,k)   = dzt(k)
	  dhu(imt,k) = dzt(k)
#  endif
        enddo   
        do i=1,imt
          do k=1,km 
            if(dht(i,k) .eq. 0.0 ) then
              buf_dht(i,jrow,k) = field_dht%missing
            else
              buf_dht(i,jrow,k) = dht(i,k)
            endif
            if(dhu(i,k) .eq. 0.0 ) then
              buf_dhu(i,jrow,k) = field_dhu%missing
            else
              buf_dhu(i,jrow,k) = dhu(i,k)
            endif
          enddo
        enddo
      enddo  
      call mpp_write(io, field_dht, buf_dht)
      call mpp_write(io, field_dhu, buf_dhu)
# endif

      write (stdout,'(a,a)')
     &   ' => Topography fields written to file ', file_name
      return
      end
#endif





      subroutine shuttle_kmt (kmt, im, jm, km, action, filename)
c
c=======================================================================
c
c     export a "kmt" field under formatted control to file "filename" 
c     to allow for editing "filename" with a text based editor. Also,
c     import file "filename" under formatted control to "kmt" for use
c     in MOM
c      
c     inputs:
c
c     kmt      = integer array
c     im       = the 1st dimension of array
c     jm       = the 2nd dimension of array
c     action   = "import" to read the kmt field
c                "export" to write the kmt field
c     filename = "kmt.dta" (if changed, a stop is executed. If the stop
c                is removed, proceed with caution.
c
c     output: "kmt" written/read as formatted array to/from "filename"
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
      dimension kmt(im,jm)
      character*20 cell_type
      character*20 form1, form2
      character*(*) action, filename
      
      real(kind=4) kmt4
c
      if (filename .ne. "kmt.dta") then
        write (stdout,*) '=>Error: filename=',filename,' in shuttle_kmt'
	write (stdout,*) '         remove the "stop" to continue'
        call abort()
      endif
c
      if (action .ne. 'import' .and. action .ne. 'export') then
        write (stdout,*) 'action =',action
     &,  '=>Error: Only "import" or "export" are allowed in shuttle_kmt'
        call abort()
      endif
c
c     read or write "kmt"
c
#ifdef read_unformatted_kmt
c
c     read older kmt.dta files written unformatted
c
      if (action .eq. 'import') then
        call getunit (io, 'kmt.dta'
     &,             'unformatted sequential rewind ieee')
        read (io) 
        read (io) imax, jmax, kmax
        call size_check (imax, jmax, kmax, 'shuttle_kmt unform', 'stop')
        read (io)
#if defined _CRAY     
        read (io) kmt
#else        
        read (io) kmt4
        kmt=kmt4
#endif
        call relunit (io)
	if (imax .eq. imax) return
      endif
#endif
c
c     determine format to use
c
      if (km .le. 9) then
        inc=120
c       bug: istep should equal 4 but original code used istep=2
        write(stdout,*) '=>Warning: shuttle_kmt. istep=4 may not be '
     &,'compatible with older kmt.dta files. If so, change to istep=2'
	istep=4
        write (form1,'(a)') '(/, 2x, 30i4)'
        write (form2,'(a)') '(1x,i3,1x, 120i1)'
      elseif (km .le. 99) then
        inc=40
	istep=2
        write (form1,'(a)') '(/,/,/,2x,20i6/)'
        write (form2,'(a)') '(1x,i3,1x,40i3)'
      elseif (km .le. 999) then
        inc=30
	istep=2
        write (form1,'(a)') '(/,/,/,1x,20i8/)'
        write (form2,'(a)') '(1x,i3,1x,30i4)'
      else
        write (stdout,*)'=>Error: Can`t handle km > 999 in shuttle_kmt'
        call abort()
      endif
c
c     read/write kmt.dta files which are formatted for editing purposes
c
      call getunit (io, filename, 'formatted sequential rewind')
      if (action .eq. 'import') then
        read (io,'(i6,i6,i6,a20)') imax, jmax, kmax, cell_type
        call size_check (imax, jmax, kmax, 'shuttle_kmt', 'stop')
#ifdef partial_cell
        if (cell_type .ne. 'partial cells') then
          write (stdout,*)
     &    '=>Error: the kmt.dta file is not for partial cells'
          write (stdout,*) 'cell_type=',cell_type
          call abort()
        endif
#else
        if (cell_type .ne. 'full cells') then
          write (stdout,*)
     &    '=>Error: the kmt.dta file is not for full cells'
          write (stdout,*) 'cell_type=',cell_type
          call abort()
        endif
#endif
      else
#ifdef partial_cell
        cell_type = 'partial cells'
#else
        cell_type = 'full cells'
#endif
        write (io,'(i6,i6,i6,a20)') im, jm, km, cell_type
      endif
      do l=0,im,inc
        incr = min(inc,im-l)
        if (action .eq. 'import') then
          read (io,form1) (ii,i=1,incr,istep)
        else
          write (io,form1) (l+i,i=1,incr,istep)
	endif
        do jrow=jm,1,-1
          if (action .eq. 'import') then
            read (io,form2) jrowin, (kmt(l+i,jrow),i=1,incr)
          else
            write (io,form2) jrow, (kmt(l+i,jrow),i=1,incr)
	  endif
        enddo
      enddo
      call relunit (io)
      return
      end



      subroutine isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle
     &,                    imt, jmt, km, mnisle, maxipp)
c
c=======================================================================
c          Computes map of land masses and island perimeters
c
c          The main computational subroutine, expand, uses a "floodfill"
c          algorithm to expand one previously unmarked land
c          point to its entire connected land mass and its perimeter
c          ocean points. Disconnected land masses separated by a minimum
c          of one ocean T-cell are considered a meta land mass. There
c          is zero net flow between the members of a meta land mass.
c
c          The subroutine expand uses a queue of size maxq of
c          coordinate pairs of candidate points.  Suggested
c          size for maxq is 4*(imt+jmt).  Queue overflow stops
c          execution with a message to increase the size of maxq.
c          Similarly a map with more that maxipp island perimeter
c          points or more than mnisle land masses stops execution
c          with an appropriate error message.
c
c          Input:
c                  kmt = array of depths.  (0 for land) (>0 for ocean)
c          Outputs:
c                  map = map of land masses and their ocean perimeters
c                           mid-ocean cells are labeled 0
c                           land masses are labeled 1, 2, 3, ...,
c                           their perimeter ocean cells -1, -2, -3, ...,
c                  iperm = i coordinates of perimeter points
c                  jperm = j coordinates of perimeter points
c                  iofs = offset of each land mass in iperm, jperm
c                  nippts = number of island perimeter points by isle
c                  nisle = number of land masses
c          Array size inputs:
c                  imt = east/west array extents
c                  jmt = north/south array extents
c                  mnisle = maximum number of land masses
c                  maxipp = maximum number of island perimeter points
c
c
c         author: Charles Goldberg 
c         December 1993/revised February 1995
c
c         modifications:  R.C.Pacanowski   rcp@gfdl.gov
c                         eliminate perimeter violations by allowing
c                         disconnected land masses separated by one
c                         ocean T-cell to be thought of as one
c                         "meta" land mass. The implication is zero
c                         flow between masses within the meta mass.
c                         Also, remove misc useless code.
c=======================================================================
c
#include "stdunits.h"
c                                                                       
c     a "floodfill" algorithm is used to compute land mass ocean 
c     perimeter points. The maximum suggested size for the queue
c     "maxq" is 4*(imt+jmt). 
c
      parameter (maxq=10000)
      common /qsize/ maxqsize
c
      dimension kmt(imt,jmt)
      dimension map(imt,jmt)
      dimension iperm(maxipp)
      dimension jperm(maxipp)
      dimension nippts(mnisle)
      dimension iofs(mnisle)
c
      dimension iq(maxq), jq(maxq)
      integer qfront, qback, ocean
      parameter (land=1, ocean=0)
c
      write (stdout,'(/,a,/)') 'Finding perimeters of all land masses'
c
c-----------------------------------------------------------------------
c     copy kmt to map changing notation
c     initially, 0 means ocean and 1 means unassigned land in map
c     as land masses are found, they are labeled 2, 3, 4, ...,
c     and their perimeter ocean cells -2, -3, -4, ...,
c     when no land points remain unassigned, land mass numbers are
c     reduced by 1 and their perimeter ocean cells relabelled
c-----------------------------------------------------------------------
c
      do i=1,imt
        do j=1,jmt
          if (kmt(i,j) .gt. 0) then
            map(i,j) = ocean
          else
            map(i,j) = land
          end if
        end do
#ifdef obc_south
c
c       first row is land for calculating poisson equation
c
        map(i,1)     = land
#endif
#ifdef obc_north
c
c       last row is land for calculating poisson equation
c
        map(i,jmt)   = land
#endif
      end do
#if defined obc_west || defined obc_east      
      do j=1,jmt
# ifdef obc_west
c
c       first column is land for calculating poisson equation
c
        map(1,j)     = land
# endif
# ifdef obc_east
c
c       last column is land for calculating poisson equation
c
        map(imt,j)   = land
# endif
      enddo
#endif
c
c-----------------------------------------------------------------------
c     find unassigned land points and expand them to continents
c-----------------------------------------------------------------------
c
      maxqsize = 0
      call qinit (iq, jq, qfront, qback)
      label = 2
      iofs(label) = 0
      nippts(label) = 0
      nerror = 0
#ifdef symmetry
      jnorth = jmt-1
#else
      jnorth = jmt
#endif
#ifdef cyclic
      iwest = 2
      ieast = imt-1
#else
      iwest = 1
      ieast = imt
#endif
      do j=jnorth,1,-1
        do i=iwest,ieast
          if (map(i,j) .eq. land) then
            call qpush (i, j, iq, jq, qfront, qback)
            call expand (map, label, iq, jq, qfront, qback, nerror
     &,                  iperm, jperm, iofs, nippts
     &,                  imt, jmt, km, mnisle, maxipp, kmt)
            write (stdout,'(a,i2,a,i4)')
     &        'number of island perimeter points: nippts(',label-1,')=',
     &         nippts(label)
            label = label + 1
            if (label .gt. mnisle) then
              write (stdout,'(a,i3,a)') 'ERROR==> mnisle=',mnisle
     &,                                 ' is too small'
              call abort()
            end if
            iofs(label) = iofs(label-1) + nippts(label-1)
            nippts(label) = 0
          end if
        end do
      end do
      nisle = label - 1
c
c-----------------------------------------------------------------------
c     relabel land masses and their ocean perimeters
c-----------------------------------------------------------------------
c
      do i=iwest,ieast
        do j=1,jnorth
          if (map(i,j) .ne. 0) then
            map(i,j) = map(i,j) - sign(1, map(i,j))
          end if
        end do
      end do
      do isle=2,nisle
        iofs(isle-1) = iofs(isle)
        nippts(isle-1) = nippts(isle)
      end do
      nisle = nisle - 1

#ifdef symmetry
      do i=iwest,ieast
        map(i,jmt) = map(i,jmt-1)
      end do
#endif
#ifdef cyclic
      do j=1,jmt
        map(1,j) = map(imt-1,j)
        map(imt,j) = map(2,j)
      end do
#endif
c
#ifdef debug_island_perimeters
      call showmap (map, imt, jmt, linewidth)
#endif
c
      write (stdout,*) ' Island perimeter statistics:'
      write (stdout,*) 'maximum queue size was ',maxqsize
      write (stdout,*) 'number of land masses is ', nisle
      write (stdout,*) 'number of island perimeter points is ',
     &        nippts(nisle) + iofs(nisle)
c
      if (nerror .gt. 0) then
        write (stdout,*)  ' '
        write (stdout,*) 
     &  '==>Yipes. nerror =',nerror,' in isleperim'
        call abort()
      end if
c
      return
      end



      subroutine showmap (map, imt, jmt, linewidth)
#include "stdunits.h"
      dimension map(imt,jmt)
#if !defined narrow_map
      linewidth = 125
#else
      linewidth = 70
#endif
      write (stdout,'(/,38x,132a)') 'Land Masses and Perimeters'
      istart = 0
      iremain = imt
      do isweep=1,imt/linewidth + 1
        iline = min(iremain, linewidth)
        iremain = iremain - iline
        if (iline .gt. 0) then
          write (stdout,*)  ' '
          write (stdout,'(t6,32i5)') (istart+i+4,i=1,iline,5)
          do j=jmt,1,-1
            write (stdout,'(i4,t6,160i1)') j
     &           ,(mmm(map(istart+i,j)),i=1,iline)
          end do
          write (stdout,'(t6,32i5)') (istart+i+4,i=1,iline,5)
          istart = istart + iline
        end if
      end do
      write (stdout,*)  ' '
      return
      end


      subroutine expand (map, label, iq, jq, qfront, qback, nerror
     &,                  iperm, jperm, iofs, nippts
     &,                  imt, jmt, km, mnisle, maxipp, kmt)
c
c-----------------------------------------------------------------------
c          The subroutine expand uses a "flood fill" algorithm
c          to expand one previously unmarked land
c          point to its entire connected land mass and its perimeter
c          ocean points.   Diagonally adjacent land points are
c          considered connected. 
c
c          The subroutine expand uses a queue of size maxq of
c          coordinate pairs of candidate points.  Suggested
c          size for maxq is 4*(imt+jmt).  Queue overflow stops
c          execution with a message to increase the size of maxq.
c          Similarly a map with more that maxipp island perimeter
c          points or more than mnisle land masses stops execution
c          with an appropriate error message.
c
c          modification:  eliminate perimeter violations by allowing
c                         disconnected land masses separated by one
c                         ocean T-cell to be thought of as one
c                         "meta" land mass. The implication is zero 
c                         flow between masses within the meta mass.
c                         Also, remove other useless code.
c                         R.C.Pacanowski   rcp@gfdl.gov
c                         
c-----------------------------------------------------------------------
c
      dimension map(imt,jmt), kmt(imt,jmt)

      dimension iperm(maxipp)
      dimension jperm(maxipp)
      dimension nippts(mnisle)
      dimension iofs(mnisle)
#include "stdunits.h"
      character * 32 problem
c                                                                       
c     a "floodfill" algorithm is used to compute land mass ocean 
c     perimeter points. The maximum suggested size for the queue
c     "maxq" is 4*(imt+jmt). 
c
      parameter (maxq=10000)
      dimension iq(maxq)
      dimension jq(maxq)
      integer qfront, qback
      logical qempty
c
      integer offmap, ocean
      parameter (offmap = -1)
      parameter (land = 1, ocean = 0)
      logical keep_going
      keep_going=.true.
c
      write (stdout,'(a,i3)') 'Exploring land mass ',label-1
c
c-----------------------------------------------------------------------
c     main loop:
c        Pop a T cell off the queue and process it.
c-----------------------------------------------------------------------
c
 1000 continue

      do while (keep_going)
      keep_going=.false.
c
      if (qempty (iq, jq, qfront, qback)) then
        call qinit (iq, jq, qfront, qback)
      keep_going=.false.
      else
        call qpop (i, j, iq, jq, qfront, qback)
        if (i .eq. offmap .or. j .eq. offmap) then
c
c         case: (i,j) is off the map
c
           keep_going=.true.
        else if (map(i,j) .eq. label) then
c
c         case: map(i,j) is already labeled for this land mass
c
           keep_going=.true.
        else if (map(i,j) .eq. -label) then
c
c         case: map(i,j) is an ocean perimeter point of this land mass
c
           keep_going=.true.
        else if (map(i,j) .eq. land) then
c
c         case: map(i,j) is an unassigned land point
c
          map(i,j) = label
c         write (stdout,*)  'labeling T point ',i,j,' as ',label
c
c         force nearest neighbors onto queue
c
          jp1 = jpn(j)
	  jm1 = jmn(j)
	  ip1 = ipn(i)
	  im1 = imn(i)
c
          call qpush (i,   jp1, iq, jq, qfront, qback)
          call qpush (ip1, jp1, iq, jq, qfront, qback)
          call qpush (ip1, j,   iq, jq, qfront, qback)
          call qpush (ip1, jm1, iq, jq, qfront, qback)
          call qpush (i,   jm1, iq, jq, qfront, qback)
          call qpush (im1, jm1, iq, jq, qfront, qback)
          call qpush (im1, j,   iq, jq, qfront, qback)
          call qpush (im1, jp1, iq, jq, qfront, qback)
c
c         force land masses disconnected by one ocean T cell onto queue
c
          keep_going=.true.
        else if (map(i,j) .eq. ocean .or. map(i,j) .lt. 0) then
c
c         case: map(i,j) is an ocean point adjacent to this land mass
c
          if (map(i,j) .lt. 0) then
c
c           map(i,j) is a perimeter ocean point of another mass
c
            nerror = nerror + 1
            write (stdout,'(a,a,i3,a,i3,a,a,i3,a,i3)')
     &            'PERIMETER VIOLATION==> ',
     &            'map(',i,',',j,') is in the perimeter of both ',
     &            'land masses ', -map(i,j)-1, ' and ', label-1
            keep_going=.true.
          end if
c
c         case: map(i,j) is an ocean point--label it for current mass
c
          map(i,j) = -label
          nippts(label) = nippts(label) + 1
          if (iofs(label) + nippts(label) .gt. maxipp) then
            write (stdout,*)  'ERROR==>  maxipp=',maxipp
     &,                       ' is not large enough'
            call abort()
          end if
          iperm(iofs(label) + nippts(label)) = i
          jperm(iofs(label) + nippts(label)) = j
c
c         force nearest ummarked land neighbors onto queue
c
          jp1 = jpn(j)
	  jm1 = jmn(j)
	  ip1 = ipn(i)
	  im1 = imn(i)
c
          if (map(i,jp1) .eq. land) then
	    call qpush (i,jp1, iq, jq, qfront, qback)
	  endif
          if (map(ip1,jp1) .eq. land) then
            call qpush (ip1, jp1, iq, jq, qfront, qback)
	  endif
          if (map(ip1,j) .eq. land) then
            call qpush (ip1, j,   iq, jq, qfront, qback)
	  endif
          if (map(ip1,jm1) .eq. land) then
            call qpush (ip1, jm1, iq, jq, qfront, qback)
	  endif
          if (map(i,jm1) .eq. land) then
            call qpush (i,   jm1, iq, jq, qfront, qback)
	  endif
          if (map(im1,jm1) .eq. land) then
            call qpush (im1, jm1, iq, jq, qfront, qback)
	  endif
          if (map(im1,j) .eq. land) then
            call qpush (im1, j,   iq, jq, qfront, qback)
	  endif
          if (map(im1,jp1) .eq. land) then
            call qpush (im1, jp1, iq, jq, qfront, qback)
	  endif
          keep_going=.true.
        else
c
c         case: problems
c
          nerror = nerror + 1
          write (stdout,'(a,a,i3,a,i3,a,a,i3,a,i3)')
     &          'ERROR ==>  ',
     &          'map(',i,',',j,') is labeled for both ',
     &          'land masses ', map(i,j)-1,' and ',label-1
        end if
        keep_going=.true.
c
      end if
      end do
      return
      end


      subroutine qinit (iq, jq, qfront, qback)
c                                                                       
c     a "floodfill" algorithm is used to compute land mass ocean 
c     perimeter points. The maximum suggested size for the queue
c     "maxq" is 4*(imt+jmt). 
c
      parameter (maxq=10000)
      dimension iq(maxq)
      dimension jq(maxq)
      integer qfront, qback
c
      qfront = 1
      qback = 0
c
c     fake assignments to iq and jq to avoid "flint" warning
c
      iq(qfront) = 0
      jq(qfront) = 0
c
      return
      end


      subroutine qpush (i, j, iq, jq, qfront, qback)
c
#include "stdunits.h"
c                                                                       
c     a "floodfill" algorithm is used to compute land mass ocean 
c     perimeter points. The maximum suggested size for the queue
c     "maxq" is 4*(imt+jmt). 
c
      parameter (maxq=10000)
      common /qsize/ maxqsize
c
      dimension iq(maxq)
      dimension jq(maxq)
      integer qfront, qback
c
      qback = qback + 1
      if (qback .gt. maxq) then
        if (qfront .ne. 1) then
c         shift queue left to make room
          ishift = qfront - 1
          do ip=qfront,qback-1
            iq(ip-ishift) = iq(ip)
            jq(ip-ishift) = jq(ip)
          end do
          qfront = 1
          qback = qback - ishift
        else
          write (stdout,*)  'queue fault in qpush'
          call abort()
        end if
      end if
      iq(qback) = i
      jq(qback) = j
c
      maxqsize = max(maxqsize, (qback-qfront))
c
      return
      end



      subroutine qpop (i, j, iq, jq, qfront, qback)
c
c                                                                       
c     a "floodfill" algorithm is used to compute land mass ocean 
c     perimeter points. The maximum suggested size for the queue
c     "maxq" is 4*(imt+jmt). 
c
      parameter (maxq=10000)
      dimension iq(maxq)
      dimension jq(maxq)
      integer qfront, qback
c
      i = iq(qfront)
      j = jq(qfront)
      qfront = qfront + 1
      return
      end


      function qempty (iq, jq, qfront, qback)
c
c                                                                       
c     a "floodfill" algorithm is used to compute land mass ocean 
c     perimeter points. The maximum suggested size for the queue
c     "maxq" is 4*(imt+jmt). 
c
      parameter (maxq=10000)
      dimension iq(maxq)
      dimension jq(maxq)
      integer qfront, qback
      logical qempty
c
      qempty = (qfront .gt. qback)
c      
      return
      end

      function mmm(m)
      if (m .eq. 0) then
        mmm = 0
      else if (m .gt. 0) then
        mmm = mod(m,10)
      else
        mmm = m
      end if
      return
      end


      function jpn(j)
#include "size.h"
      integer offmap
      parameter (offmap = -1)
c
      jpn = j + 1
c
#ifdef symmetry
      if (jpn .gt. jmt) then
        jpn = offmap
      elseif (jpn .eq. jmt) then
        jpn = jmt-2
      endif
#else
      if (jpn .gt. jmt) jpn = offmap
#endif
      return
      end


      function jmn(j)
      integer offmap
      parameter (offmap = -1)
c
      jmn = j - 1
c
      if (jmn .lt. 1) jmn = offmap
      return
      end



      function ipn(i)
#include "size.h"
      integer offmap
      parameter (offmap = -1)
c
      ipn = i + 1
c
#ifdef cyclic
      if (ipn .eq. imt) ipn = 2
      if (ipn .eq. imt+1) ipn = 3
      if (ipn .eq. imt+2) ipn = 4
#else
      if (ipn .gt. imt) ipn = offmap
#endif
      return
      end



      function imn(i)
#include "size.h"
      integer offmap
      parameter (offmap = -1)
c
      imn = i - 1
c
#ifdef cyclic
      if (imn .eq. 1) imn = imt-1
      if (imn .eq. 0) imn = imt-2
      if (imn .eq. -1) imn = imt-3
#else
      if (imn .lt. 1) imn = offmap
#endif
      return
      end

#if defined sink_isolated_land_cells
      subroutine sink_isolated_cells(kmt,kmt_isolated_land_cells
     &     ,imt,jmt,num_sunk_isolated_cells)
# include "stdunits.h"
      dimension kmt(imt,jmt),kmt_isolated_land_cells(imt,jmt)

c     Find one-cell islands and land points that stick out of a land
c     mass, and sink them to a depth that is the average of their
c     neighbors.  An example of a sticking land point is:

c This mask:    00000000    is changed into this one:    00000000
c               00010000			         00000000
c               11111111			         11111111
c               11111111			         11111111

      do j=1,jmt
        do i=1,imt
          kmt_isolated_land_cells(i,j)=kmt(i,j)
        end do
      end do

      do j=2,jmt-1
        do i=2,imt-1
          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          if (kmt(i,j).eq.0) then
c           Is this a one-cell islands?: i_one_cell will be zero if
c           there is at least one land neighbor cell:
            i_one_cell =
     &           kmt(ip,j)*kmt(im,j)*kmt(i,jp)*kmt(i,jm)
     &           *kmt(ip,jp)*kmt(im,jp)*kmt(ip,jm)*kmt(im,jm)
c           Is this an isolated land point sticking out of a larger land mass?:
            isticking_land_point = 0
            if (
c           A land point is a "sticking land point" if there is ocean in the points
c           marked by 1.  points marked by * could be either land or ocean:
c           111 
c           101 <= 1: ocean, 0: land, *: ocn/land, middle point is (i,j):
c           000  
     &(kmt(ip,jp)*kmt(i,jp)*kmt(im,jp)*kmt(im,j)*kmt(ip,j).ne.0
     &           .and.(kmt(im,jm)+kmt(i,jm)+kmt(ip,jm).eq.0)) .or.
c           000  
c           101  <= 1: ocean, 0: land, *: ocn/land, middle point is (i,j):
c           111  
     &(kmt(ip,jm)*kmt(i,jm)*kmt(im,jm)*kmt(im,j)*kmt(ip,j).ne.0
     &           .and.(kmt(im,jp)+kmt(i,jp)+kmt(ip,jp).eq.0)) .or.
c           011  
c           001  <= 1: ocean, 0: land, *: ocn/land, middle point is (i,j):
c           011  
     &(kmt(ip,jm)*kmt(ip,j)*kmt(ip,jp)*kmt(i,jm)*kmt(i,jp).ne.0
     &           .and.(kmt(im,jm)+kmt(im,j)+kmt(im,jp).eq.0)) .or.
c           110  
c           100  <= 1: ocean, 0: land, *: ocn/land, middle point is (i,j):
c           110  
     &(kmt(im,jm)*kmt(im,j)*kmt(im,jp)*kmt(i,jm)*kmt(i,jp).ne.0
     &           .and.(kmt(ip,jm)+kmt(ip,j)+kmt(ip,jp).eq.0))
     & ) isticking_land_point = 1
            if (i_one_cell.gt.0  .or. isticking_land_point.eq.1) then
c             Sink this grid point:
              ikmt_avg=0
              kmt_avg=0
              do i_p=-1,1
                do j_p=-1,1
                  if (kmt(i+i_p,j+j_p).ne.0) then
                    if (.not.(i_p.eq.0.and.j_p.eq.0)) then
                      kmt_avg=kmt_avg+kmt(i+i_p,j+j_p)
                      ikmt_avg=ikmt_avg+1
                    end if
                  end if
                end do
              end do
              kmt_isolated_land_cells(i,j)=
     &             nint((1.0*kmt_avg)/(1.0*ikmt_avg))
              write(stdout,*) " sink_isolated_land_cells: at point i,j="
     &             ,i,j,"; kmt_prev,kmt(i,j)="
     &             ,kmt(i,j),kmt_isolated_land_cells(i,j)
     &             ,"; i_one_cell, isticking_land_point="
     &             ,i_one_cell, isticking_land_point
            end if
          end if
        end do
      end do

      num_sunk_isolated_cells=0
      do j=2,jmt-1
        do i=2,imt-1
          if (kmt(i,j).ne.kmt_isolated_land_cells(i,j)) then
            num_sunk_isolated_cells=num_sunk_isolated_cells+1
            kmt(i,j)=kmt_isolated_land_cells(i,j)
          end if
        end do
      end do
      write(stdout,*) " sink_isolated_cells: sank "
     &     ,num_sunk_isolated_cells," isolated land grid points."

      return
      end
#endif
