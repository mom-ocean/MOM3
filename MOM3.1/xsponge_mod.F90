      module xsponge_mod
      use time_manager_mod
!
      type(time_type) :: sp_start_time, sp_end_time
      type(time_type) :: sp_time
!
!
!     newtonian damping variables for sponge regions adjacent
!     to artificial southern and northern boundaries for use in limited 
!     domain basins. data must be prepared using run_sponge 
!     included in PREP_DATA.
!
!     spng(jrow)  = coefficient for damping T & S back to sponge row
!     sponge_row(jrow)  = maps each jrow=1..jmt to a sponge row (1 or 2)
!     j_south_spng = northern most jrow of southern sponge
!     j_north_spng = southern most jrow of northern sponge
!
!     iprevd = pointer to previous month on disk
!              (the  month whose mid point has just been passed)
!     inextd = pointer to next month on disk
!              (the next month whose mid point hasn`t been reached yet)
!     iprev  = pointer to memory buffer for the previous month data
!              (the  month whose mid point has just been passed)
!     inext  = pointer to memory buffer for the next month data
!              (the next month whose mid point hasn`t been reached yet)
!     bufmb  = buffer for holding previous and next month disk data
!     bufmb  = sponge buffer data (imt,km,max_rows,2,2) as follows:
!              max_rows = 2 (1 for southern sponge and 1 in the north)
!              (,,,1,) for temperature data
!              (,,,2,) for salinity data
!              (,,,,1..2) space to hold previous and next months data
!              
!     annlev = (t,f) = (replace seasonal data by annual means, use 
!                       seasonal data)
!     spgdpm = period in days for each monthly record
!     tspng  = time at midpoints of monthly records (days)
!     indxsp = index of dataset (needed for the interpolator)
!     readsp = (true,false) = (read, do not read) sponge data
!     wprev  = interpolation weight for previous month`s data
!
!     symbolically:
!     spdata(at time step) = (1-wprev)*bufmb(inext) + wprev*bufmb(iprev)
!
# if defined sponges
      integer :: inext, iprev, inextd, iprevd, indxsp
      integer :: max_rows, j_south_spng, j_north_spng
      
      real ::  wprev
      real, dimension(12) :: spgdpm, tspng
      real, dimension(:,:,:,:,:), allocatable :: bufmb
!                                          bufmb(imt,km,max_rows,2,2)
      integer, dimension(:), allocatable :: sponge_row
      real, dimension(:), allocatable :: spng
!                                        spng(jmt), sponge_row(jmt)
!
!     array bufmb1 is only used locally to increase I/O efficiency
!     when reading and writing sponge data from disk.
!
      real, dimension(:,:,:,:), allocatable :: bufmb1
!                             bufmb1(imt,km,max_rows,2)
!
      character(len=32) :: sstamp, stprev, stnext
      character(len=128) :: opt_sponge
      logical :: annlev, readsp
# endif
      end module
