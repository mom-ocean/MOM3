!==================== include file "readic.h" ==========================
# include "netcdf.inc"
# include "grdvar.h"
# ifdef free_surf_on
#  include "scalar.h"
# endif  
# include "derived_options.h"
!     ioic - io-unit number
!     start, count
!          - fields to specify the four dimensional vectors of start
!            of field indicees and increments of field indicees. For 
!            details consult the netcdf-manual
!     filemap - stores information on netcdf file structure
!               (see netcdfmod.F90)
!     ic_file - Filename with initial data
!
!               offsets and file-grid
!               snapshots are written from i=2..imt-1, j=2..jmt-1.
!               Hence, the is an offset of snapshot-data to the model grid.
!               The sign of the offset is defined as
!               offset = index_model - index_file
!     ioff_ic - zonal offset of the data in the file to the model grid
!     joff_ic - meridional offset of the data in the file to the model grid
!               The default is ioff_ic=1, joff_ic=1 for snapshots. 
!               A negative offset can be used, to extract 
!               subsets from larger files.
!     idim_ic   zonal grid dimension in the file-grid.
!     jdim_ic   meridional grid dimension in the file-grid.
!     kdim_ic   vertical grid dimension in the file-grid.
!               The default is idim_ic=imt-2, jdim_ic=jmt-2 kdim_ic=kbot.
!               Smaller file grids than the default are not allowed.  
!
!     nic2d   - number of 2-dimensional variables to be read
!     nic3d   - number of 3-dimensional variables to be read
!     icmap.d(1,:) - names of variables to be read from the netcdf-file
!     icmap.d(2,:) - Corresponding variable name for MOM. This is usually
!               the tracer name as defined in trname or for baroclinic
!               velocity components 'u' and 'v'. For the surface elevation
!               use 'eta'. The prescribed setup shows the defaults.
!               Go through subroutine rowi in setocn.F to see more
!               details. 
!               Since trname may not be well defined, the alternative
!               syntax "tracer_n" can be used, to put data to tracer
!               number n. 
!    trmap    - A field of maximum dimension nt+2. It is defined for 
!               m=1, nic3d and contains the model tracer index, where
!               the data read from file correspond to. If the file variables
!               are 'u' and 'v' the index is set to nt+1 and nt+2 respectively. 
!    idmap    - The netcdf-id of the variables m=1, nic3d which are read
!               from file.               
!    max2d    - The maximum number of 2-D variables which can be read from file.
!    max3d    - The maximum number of 3-D variables which can be read from file.
!               Increase these parameters if more variables are required.
!                  
      character*12 trname
      common /docnam/ trname(nt)
      integer                         :: ioic, ictimestep
      integer                         :: ioff_ic, joff_ic
      integer                         :: idim_ic, jdim_ic, kdim_ic 
      data ioff_ic, joff_ic              /1, 1/
      data idim_ic, jdim_ic, kdim_ic     /imtm2, jmtm2, kbot/
      integer, dimension(4)           :: start, count
      type(filemap_type)              :: filemap      
      character(len=40)               :: ic_file, filename
      logical                         :: error
      integer                         :: nic2d, nic3d
      integer, dimension(nt+2)        :: trmap, idmap
      data nic2d, nic3d                  /1, 2/
      parameter (max2d = 6, max3d = 15)
      character(len=24), dimension(2,max2d) :: icmap2d
      character(len=24), dimension(2,max3d) :: icmap3d
      data ((icmap2d(n,m),n=1,2), m=1, max2d) 
     &         /
#if defined stream_function
     &           'psi', 'psi'
     &,          'none','none'
     &,          'none','none'
#else
# if defined explicit_free_surface
     &           'eta', 'etat'
# else
     &,          'eta' ,'ps' 
# endif
     &,          'ubar','ubar'
     &,          'vbar','vbar'
#endif
     &,          'none','none'
     &,          'none','none'
     &,          'none','none'
     &          /
!     &,          'taux' ,'taux' 
!     &,          'tauy' ,'tauy' 
!     &,          'sedi' ,'sedi' 
      data ((icmap3d(n,m),n=1,2), m=1, max3d) 
     &         /
     &           'temp','temp' 
     &,          'salinity','salinity'
     &,          'u'   ,'u'
     &,          'v'   ,'v'
     &,          'none','none'
     &,          'none','none'
     &,          'none','none'
     &,          'none','none'
     &,          'none','none'
     &,          'none','none'
     &,          'none','none'
     &,          'none','none'
     &,          'none','none'
     &,          'none','none'
     &,          'none','none'
     &          /
!     &,         'nitr','tracer_3',
!     &,         'phyt','tracer_4',
      data ic_file /'ic.dta.nc'/
      data ictimestep /1/	    
      namelist /icdata/ 
     &          ic_file, ictimestep, nic2d, nic3d, icmap2d, icmap3d ! specify user input there
     &         ,idim_ic, jdim_ic, kdim_ic, ioff_ic, joff_ic         ! specify user input there
! data buffers
      real(kind=8), allocatable :: buff8 (:,:,:,:)
!
#if defined use_netCDF
      integer, parameter :: INPUT_FORM = MPP_NETCDF
#else
      integer, parameter :: INPUT_FORM = MPP_IEEE32
#endif
      real               :: missing
      real, allocatable  :: csum(:)
      data missing /-1.e33/
