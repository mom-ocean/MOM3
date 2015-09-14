  module river_mod
! call this module after setocn!
  use arrays_2d 
  use mpp_mod
  use time_manager_mod
  private
!
!-----------------------------------------------------------------------
! For parallel treatment the rivers are split in river boxes.
! One dataset belong to one box. One river has normally several
! boxes. However, the dataset of one box contains some information of the
! corresponding river, i.e., index, name and area.
!
! Variables:
!   rivbox_in_pe  - number of grid cells with river runoff in the domain
!                   of the PE
!   river_data    - field of data structures, will be allocated with 
!                   dimension rivbox_in_pe. Each element stores all 
!                   information on the river runoff of one grid cell.
! 
!   river_rec     - field of data structures, will be allocated with 
!                   dimension rivbox_in_pe. Each element contains all 
!                   information on the river runoff of one grid cell 
!                   interpolated to model time.
! 
!  iprevmriv, inextmriv
!                 - index pointers in memory for timeinterp, not used
!  iprevdriv, inextdriv
!                 - index pointers determined with timeinterp. Point to
!                   the elements in river_data used for time interpolation.
!  wprevriv       - weight of the data record iprevdriv in the time
!                   interpolation
!  methodriv      - interpolation method - see timeinterp
!
!     author:      M. Schmidt   e-mail=> martin.schmidt@io-warnemuende.de
!
!-----------------------------------------------------------------------
  type, public :: river_data_type     
     integer               :: index    ! identifies river
     character*16          :: name     ! name of the river 
     integer               :: ir,jr    ! zonal and meridional grid index
     real                  :: area     ! area of all boxes of the river
     integer               :: mrecriv  ! number of data records
     logical               :: perriv   ! true for periodic data
     integer               :: ninterp  ! data index for timeinterp
     real, pointer         :: aprec(:) ! period lengths
     real, pointer         :: tdrec(:) ! times of data records
     type(time_type)       :: start_time, end_time
     character*32, pointer :: dstamp(:)! time stamps marking the end of each record
     real, pointer         :: rff(:)   ! river fresh water flux(time)
     real, pointer         :: rtf(:,:) ! river tracer flux(tracer index, time)
     character*20, pointer :: trname(:) ! Name of the tracers ( tracer index)
  end type river_data_type

  type, public  :: river_type     
     integer        :: index   ! identifies river
     character*16   :: name    ! name of the river 
     integer        :: ir,jr   ! zonal and meridional grid index
     real           :: area    ! area of all boxes of the river
     character*32   :: dstamp  ! time stamps marking the end of each record
     real, pointer  :: rff     ! river fresh water flux
     real, pointer  :: rtf(:)  ! river tracer flux(tracer index)
     character*20, pointer :: trname(:) ! Name of the tracers ( tracer index)
  end type river_type
      
  integer, public, allocatable:: iprevdriv(:), inextdriv(:)
  integer, public, allocatable:: iprevmriv(:), inextmriv(:)
  integer, public             :: methodriv
  real, public, allocatable   :: wprevriv(:) 
  
  integer, public:: rivbox_in_pe
  type(time_type), public :: river_time_offset
  
  type(river_data_type), public, allocatable::  river_data(:)
  type(river_type)     , public, allocatable::  river_rec(:)
    
  end module river_mod
