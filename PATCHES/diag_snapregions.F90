  module diag_mod
#ifdef snapshots_regions
  use mpp_mod
  use datamod
  use mpp_io_mod
  use mpp_domains_mod
  use time_manager_mod
#include "derived_options.h"
  private
  
  type, public :: region2D
     sequence
     character(len=32)       :: name
     type(domain2D), pointer :: domains(:)
     type(domain2D), pointer :: domain
  end type region2D
 
!-----------------------------------------------------------------------
!     field types
!-----------------------------------------------------------------------

   type, public :: region_field
     sequence
!-----------------------------------------------------------------------
!	axis types
!-----------------------------------------------------------------------

     type (axistype) :: axis_xt
     type (axistype) :: axis_xu
     type (axistype) :: axis_yt
     type (axistype) :: axis_yu
     type (axistype) :: axis_zt
     type (axistype) :: axis_zw
     type (axistype) :: axis_surf
     type (axistype) :: axis_t

     type(fieldtype) :: u
     type(fieldtype) :: v
     type(fieldtype) :: w
     type(fieldtype) :: convU

     type(fieldtype) :: u_gm
     type(fieldtype) :: v_gm
     type(fieldtype) :: w_gm

     type(fieldtype) :: sfft
     type(fieldtype) :: taux
     type(fieldtype) :: tauy
#  ifdef air_pressure
     type(fieldtype) :: pair
#  endif

     type(fieldtype) :: eta
     type(fieldtype) :: psiu
     type(fieldtype) :: psiv
     type(fieldtype) :: psi
 
     type(fieldtype) :: hblt
     type(fieldtype) :: hmxl
     type(fieldtype) :: oppsignx
     type(fieldtype) :: oppsigny
 
     type(fieldtype) :: period
     type(fieldtype), pointer :: trcr(:)
     type(fieldtype), pointer :: stf(:)
 
   end type region_field
      
!-----------------------------------------------------------------------
!     axis edge arrays
!-----------------------------------------------------------------------

   type (axistype), save, public :: axis_xt_edges
   type (axistype), save, public :: axis_xu_edges
   type (axistype), save, public :: axis_yt_edges
   type (axistype), save, public :: axis_yu_edges
   type (axistype), save, public :: axis_zt_edges
   type (axistype), save, public :: axis_zw_edges

# ifdef diagnostics_himem

   type, public :: region_buffer 
     sequence
     real, pointer :: trcr(:,:,:,:)

     real, pointer :: u(:,:,:)
     real, pointer :: v(:,:,:)
     real, pointer :: w(:,:,:)
     real, pointer :: u_gm(:,:,:)
     real, pointer :: v_gm(:,:,:)
     real, pointer :: w_gm(:,:,:)

     real, pointer :: stf(:,:,:)
     real, pointer :: smf(:,:,:)

     real, pointer :: eta(:,:)
     real, pointer :: oppsignx(:,:)
     real, pointer :: oppsigny(:,:)
     real, pointer :: convU(:,:)
     real, pointer :: psiu(:,:)
     real, pointer :: psiv(:,:)
     real, pointer :: psi(:,:)
     real, pointer :: sfft(:,:)
     real, pointer :: hmxl(:,:)
     real, pointer :: hblt(:,:)
   end type region_buffer
   
   type, public :: region_output
     type(region_field)  :: fields
     type(region_buffer) :: buffs
   end type region_output
# endif

#endif    



end module diag_mod
 
