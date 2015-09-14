  module posdef_mod
#ifdef posdef
! call this module after setocn!
  private
! This module specifies regions where overshootings may crash the model.
! This applies especially for river mouths where salinity is small.
! In the grid cells specified here, vertical and horizontal diffusion is
! set to the minimum value A = 0.5 u ds to prevent overshooting.
!     author:      M. Schmidt   e-mail=> martin.schmidt@io-warnemuende.de
!
!-----------------------------------------------------------------------
!
  type, public :: posdef_row_type     
     integer          :: jrow           ! This should be jrow
     logical          :: posdef_north   ! If true the row is considered
     logical          :: posdef_east    ! If true the row is considered
     integer          :: iblock_east         ! The number of blocks in the row
     integer          :: iblock_north         ! The number of blocks in the row
     integer, pointer :: is_north(:), ie_north(:)   ! Start-stop indicees of the blocks
     integer, pointer :: is_east(:) , ie_east (:)   ! Start-stop indicees of the blocks
  end type posdef_row_type  
  
  type(posdef_row_type) , public, allocatable :: prow(:)
  logical, public                             :: posdefcell_in_pe
  data posdefcell_in_pe /.false./      
  integer, parameter   :: bufflen=1000
  integer, allocatable :: buff(:)
  logical              :: inblock
  public posdefi, posdef_tracer
  contains 
  
  subroutine posdefi
# include "stdunits.h"
# include "taskrows.h"
# include "posdef.h"
! Check, if there are cells in the domain of the PE which must be treated 
! with a positive definite scheme
!-----------------------------------------------------------------------
!  
  call checks_posdef
!
!-----------------------------------------------------------------------
  do jr=1,jpos
    if (jline(jr).ge.jscomp.and.jline(jr).le.jecomp) then
      posdefcell_in_pe = .true.
    endif
  enddo
  if (posdefcell_in_pe) then
    allocate (prow(jstask:jetask))
    prow(:)%posdef_north = .false.
    prow(:)%posdef_east  = .false.
    do jr=1, jpos
      jrow = jline(jr)
      if (jstask.le.jrow.and.jrow.le.jetask) then
        prow(jrow)%posdef_north   = .true.
        prow(jrow)%posdef_east    = .true.
        iblock = ibpos(jr)
        prow(jrow)%iblock_east  = iblock
        allocate (prow(jrow)%is_east(iblock))
        allocate (prow(jrow)%ie_east(iblock))
!----------------------------------------------------------------------
!     Find the start and stop index for modified diffusion on eastern
!     cell boundaries. is_east = is - 1, since diffusion on the western 
!     boundary of cell i is calculated with cell index i-1.
!----------------------------------------------------------------------
        do ib=1,iblock
          prow(jrow)%is_east(ib) = is(jr,ib) - 1
          prow(jrow)%ie_east(ib) = ie(jr,ib)
          if(is(jr,ib).gt.ie(jr,ib)) then
            write(stdout,'(/a,i2/a)') &
            '=>Error: is > ie for jrow = ',jrow,'Check posdef.h!'
            call abort()
      	  endif
        enddo
!----------------------------------------------------------------------
!     Find the start and stop index for modified diffusion on northern
!     cell boundaries. 
!----------------------------------------------------------------------
        prow(jrow)%iblock_north  = iblock
        allocate (prow(jrow)%is_north(iblock))
        allocate (prow(jrow)%ie_north(iblock))
        do ib=1,iblock
          prow(jrow)%is_north(ib) = is(jr,ib)
          prow(jrow)%ie_north(ib) = ie(jr,ib)
        enddo
      endif
    enddo
!----------------------------------------------------------------------
!   Find the start and stop index for modified diffusion on southern
!   cell boundaries which are northern boundaries at row jrow -1. 
!----------------------------------------------------------------------
!   
    jst = max(jstask+1,3)
    do jrow = jst, jetask
!     The case that nothing is defined at jrow-1. 
!     Define the cells at row jrow-1 like in jrow
      if(prow(jrow)%posdef_north.and.(.not.prow(jrow-1)%posdef_north)) then
        prow(jrow-1)%posdef_north = .true.
        iblock = prow(jrow)%iblock_north 
        prow(jrow-1)%iblock_north    = iblock
        allocate (prow(jrow-1)%is_north(iblock))
        allocate (prow(jrow-1)%ie_north(iblock))
        prow(jrow-1)%is_north(:) = prow(jrow)%is_north(:)
        prow(jrow-1)%ie_north(:) = prow(jrow)%ie_north(:)
      endif
!     The case that something is defined at jrow-1. 
      if(prow(jrow)%posdef_north.and.prow(jrow-1)%posdef_north) then
!       The block in jrow must be added eventually to jrow-1
!       Allocate a buffer row. Set buff=0 for normal cells and buff=1
!       for cells which must be changed. 
        allocate(buff(bufflen))
	buff(:)=0
!       Set the buffer for cells, which are still defined in jrow - 1
 	do ib=1,prow(jrow-1)%iblock_north
	 
	  do i=prow(jrow-1)%is_north(ib), prow(jrow-1)%ie_north(ib)
	    buff(i) = 1
	  enddo
	enddo
!       Set the buffer for cells, which provide the southern flux for 
!       cells in jrow 
	do ib=1,prow(jrow)%iblock_north 
	  do i=prow(jrow)%is_north(ib), prow(jrow)%ie_north(ib)
	    buff(i) = 1
	  enddo
	enddo
!       Find the new start stop blocks
        iblock  = 0
	inblock = .false.
	do i=1, bufflen
!         A new block starts
	  if(buff(i).eq.1.and.(.not.inblock)) then
	    iblock = iblock + 1
	    inblock = .true.
	  endif	
!         A block ends
	  if(buff(i).eq.0.and.inblock) inblock = .false.
	enddo
!       Bail out, if the buffer ends within a cell block.	
	if(inblock) then
	  write(stdout,*)'=>Error! Cannot generate new start-stop index'&
	  ,'for jrow = ',jrow
	  write(stdout,*)'  Is imt > ',bufflen,'?'
	  call abort()
	endif
!       Reallocate the data structure for the reconfigured jrow-1
        deallocate (prow(jrow-1)%is_north)
        deallocate (prow(jrow-1)%ie_north)
        allocate (prow(jrow-1)%is_north(iblock))
        allocate (prow(jrow-1)%ie_north(iblock))
	prow(jrow-1)%iblock_north = iblock
        iblock  = 0
	inblock = .false.
!       Put the information from the buffer to prow(jrow-1)
	do i=1, bufflen
!         A new block starts
	  if(buff(i).eq.1.and.(.not.inblock)) then
	    iblock = iblock + 1
	    inblock = .true.
	    prow(jrow-1)%is_north(iblock) = i
	  endif	
!         A block ends
	  if(buff(i).eq.0.and.inblock) then
	    inblock = .false.
	    prow(jrow-1)%ie_north(iblock) = i - 1
	  endif
	enddo
        deallocate(buff)
      endif
    enddo
!   Write some diagnostics
    write(stdout,'(/a,i2)')&
     'The following tracer cells are specified for posdef at PE ', pn-1
    write(stdout,'(a)') 'North faces: '
    write(stdout,'(a)') '   jrow   is    ie'
    do jrow = jstask, jetask
      if (prow(jrow)%posdef_north) &
      write(stdout,'(2x,100(i4,2x))') jrow,&
       (prow(jrow)%is_north(ib),prow(jrow)%ie_north(ib)&
       ,ib=1,prow(jrow)%iblock_north)
    enddo
    write(stdout,'(/a)') 'East faces: '
    write(stdout,'(a)') '   jrow   is    ie'
    do jrow = jscomp, jecomp
      if (prow(jrow)%posdef_east) &
      write(stdout,'(2x,100(i4,2x))') jrow,&
       (prow(jrow)%is_east(ib),prow(jrow)%ie_east(ib)&
       ,ib=1,prow(jrow)%iblock_east)
    enddo
  else  
    write(stdout,'(/a,i2/)')&
     'No tracer cells are specified for posdef at PE ', pn-1
  endif  
  
  end subroutine posdefi
  
  subroutine checks_posdef
# include "derived_options.h"
# include "stdunits.h"
# include "posdef.h"
  logical :: error
  error = .false.
  if(jpos.le.0) then
    write(stdout,'(/a,a,i2/a)') &
        '=>Error: Option posdef is enabled, but the number of rows jpos ',&
        '  is ',jpos,'Check posdef.h!'
    error = .true.
  endif
# if !defined tracer_horz_mix_var
    write(stdout,'(/a/a)') &
        '=>Error: Option posdef does not work with constant horizontal '//&
        'mixing coefficients.', ' Enable tracer_horz_mix_var!'
    error = .true.
# endif  
  if (error) call abort()
  end subroutine checks_posdef
    
  subroutine posdef_tracer (joff, js, je, is, ie, n)
# if defined tracer_horz_mix_var
!=======================================================================
!     compute tracers at "tau+1" for rows js through je in the MW.
!
!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!
!     author:   m. schmidt
!=======================================================================
!
      use memory_window
#  include "param.h"
#  include "coord.h"
#  include "grdvar.h"
#  include "hmixc.h"
#  include "mw.h"
#  include "vmixc.h"
!
!----------------------------------------------------------------------
!     make the scheme positive definity to avoid negative salinities in
!     rivers
!----------------------------------------------------------------------
!
     if (posdefcell_in_pe) then
       do j=js,je
     	 jrow = j + joff
     	 if(prow(jrow)%posdef_east) then
     	   do k=1,kbot
     	     do ib=1,prow(jrow)%iblock_east
     	       do i=prow(jrow)%is_east(ib),prow(jrow)%ie_east(ib)
#  ifdef partial_cell
		 diff_pos = abs(0.5*adv_vet(i,k,j)*dxu(i)*csu(jrow)) &
                  	     /dht(i,k,j)
#  else
		 diff_pos = abs(0.5*adv_vet(i,k,j)*dxu(i)*csu(jrow))
#  endif
        	 diff_cet(i,k,j) = max(diff_cet(i,k,j),diff_pos)
	       enddo
	     enddo
	   enddo
         endif
       enddo
       do j=js-1,je
         jrow = j + joff
	 if(prow(jrow)%posdef_north) then
           do k=1,kbot
             do ib=1,prow(jrow)%iblock_north
               do i=prow(jrow)%is_north(ib),prow(jrow)%ie_north(ib)
#  ifdef partial_cell
		 diff_pos = abs(0.5*adv_vnt(i,k,j)*dyu(jrow)&
     	            *csur(jrow))/dht(i,k,j)
#  else
	         diff_pos = abs(0.5*adv_vnt(i,k,j)*dyu(jrow)*csur(jrow))
#  endif
             	 diff_cnt(i,k,j) = max(diff_cnt(i,k,j),diff_pos)
	       enddo
	     enddo
	   enddo
         endif
       enddo
       nt2 = min(n,nsmix)
       do j=js,je
         jrow = j + joff
	 if(prow(jrow)%posdef_east) then
           do k=1,kbot-1
             do ib=1,prow(jrow)%iblock_east
               do i=prow(jrow)%is_east(ib),prow(jrow)%ie_east(ib)
		 diff_pos = abs(0.5*adv_vbt(i,k,j)*dzt(k))
        	 diff_cbt(i,k,j,nt2) = max(diff_cbt(i,k,j,nt2),diff_pos)
  	       enddo
  	     enddo
  	   enddo
  	 endif
       enddo
     endif
# endif  

   return
   end subroutine posdef_tracer



#endif  
  end module posdef_mod
