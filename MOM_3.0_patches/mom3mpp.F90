module mom3_mpp_mod

!-----------------------------------------------------------------------
!     MOM3 interface to mpp_mod
!
!     author: V. Balaji      e-mail => vb@gfdl.gov
!-----------------------------------------------------------------------
  use mpp_mod
  use datamod
  use mpp_domains_mod
  public
!
! for each latitude row, list PE on which it is computed
!     (1 and jmt are added for convenience)
!
  integer, dimension(:), allocatable ::  pe_of_lat
  
  private :: ramdisk3d, ramdisk2d
! currently one domain per PE
! defines execution domains (jscomp:jecomp) and domains with buffer regions (jstask:jetask)
  type(domain2D), public, allocatable, dimension(:), target :: mom_domains
  type(domain2D), public, pointer                           :: mom_domain
! defines execution domains (jscomp:jecomp) and domains with buffer regions (jstask-1:jetask+1)
  type(domain2D), public, allocatable, dimension(:), target :: mom_etat_domains
  type(domain2D), public, pointer                           :: mom_etat_domain
! execution domains with global data domains (2:jmt-1)
  type(domain2D), public, allocatable, dimension(:), target :: mom_global_domains
  type(domain2D), public, pointer                           :: mom_global_domain
#ifdef ramdrive
! domains for ramdisk3d
  type(domain2D), public, allocatable, dimension(:), target :: mom_ramdisk_domains
  type(domain2D), public, pointer                           :: mom_ramdisk_domain
#endif ramdrive
#ifdef max_window
! domains for T & U
  type(domain2D), public, allocatable, dimension(:), target :: mom_prognostic_domains
  type(domain2D), public, pointer                           :: mom_prognostic_domain
#endif


  contains

      subroutine define_tasks
!
!=======================================================================
!     output:
!
!     1) define starting and ending "jrows" for each processor
!        for n=1 to num_processors
!           
!        jstask3(n)  = southernmost buffer jrow for processor "n"
!        jetask3(n)  = northernmost buffer jrow for processor "n"
!        jscomp3(n)  = southernmost computed jrow for processor "n"
!        jecomp3(n)  = northernmost computed jrow for processor "n"
!        num_loads3(n)= number of times processor "n`s" MW is loaded
!                       
!     2) set starting and ending "jrows" for this processor (pn). each
!        processor has exactly one MW (memory window) and it is 
!        initially filled with data from "jrows" starting with
!        jrow = jstask. Subsequent northward movement of the MW and
!        loading of new "jrows" occurs until jrow = jetask. Equations
!        are solved from jrow = "jscomp" to jrow = "jecomp".
!     
!        jstask  = southernmost buffer jrow for this processor (pn)
!        jetask  = northernmost buffer jrow for this processor (pn)
!        jscomp  = southernmost computed jrow for this processor (pn)
!        jecomp  = northernmost computed jrow for this processor (pn)
!        num_loads= number of times processor "pn`s" MW is loaded
!
!     author: R. C. Pacanowski   e-mail -> rcp@gfdl.gov
!=======================================================================
!
#include "param.h"
#include "taskrows.h"
!
!-----------------------------------------------------------------------
!     1-D decomposition in latitude:
!
!     divide the physical domain into bands of latitude and assign each
!     band to a processor. The task for each processor
!     is to solve prognostic equations on interior rows within
!     the band (computed rows) but not on edge rows (buffer rows)
!
!     define latitude row limits for all tasks on each processor.
!-----------------------------------------------------------------------
!
!     these are the domains for model fields with halo regions jbuf wide

      allocate( mom_domains(0:npes-1) )
      mom_domain => mom_domains(pe)
#ifdef cyclic 
      call mpp_define_domains( (/2,imt-1,2,jmt-1/), mom_domains, xflags=GLOBAL_COMPUTE_DOMAIN + CYCLIC_GLOBAL_DOMAIN, &
			      xhalo=jbuf, yhalo=jbuf )
#else
      call mpp_define_domains( (/2,imt-1,2,jmt-1/), mom_domains, xflags=GLOBAL_COMPUTE_DOMAIN, &
                              xhalo=jbuf, yhalo=jbuf )
#endif
      allocate( mom_etat_domains(0:npes-1) )
      mom_etat_domain => mom_etat_domains(pe)
#ifdef cyclic 
      call mpp_define_domains( (/2,imt-1,2,jmt-1/), mom_etat_domains, xflags=GLOBAL_COMPUTE_DOMAIN + CYCLIC_GLOBAL_DOMAIN, &
                                    xhalo=jbuf, yhalo=jbuf+1 )
#else
      call mpp_define_domains( (/2,imt-1,2,jmt-1/), mom_etat_domains, xflags=GLOBAL_COMPUTE_DOMAIN , &
                              xhalo=jbuf, yhalo=jbuf+1 )
#endif

      allocate( mom_global_domains(0:npes-1) )
      mom_global_domain => mom_global_domains(pe)
#ifdef cyclic 
      call mpp_define_domains( (/2,imt-1,2,jmt-1/), mom_global_domains, xflags=GLOBAL_COMPUTE_DOMAIN + CYCLIC_GLOBAL_DOMAIN, &
                              yflags=GLOBAL_DATA_DOMAIN, xhalo=jbuf, yhalo=jbuf )
#else
      call mpp_define_domains( (/2,imt-1,2,jmt-1/), mom_global_domains, xflags=GLOBAL_COMPUTE_DOMAIN , &
                              yflags=GLOBAL_DATA_DOMAIN, xhalo=jbuf, yhalo=jbuf )
#endif


# ifdef ramdrive
!     these are domains for updating ramdisk3d
      allocate( mom_ramdisk_domains(0:npes-1) )
      mom_ramdisk_domain => mom_ramdisk_domains(pe)
      call mpp_define_domains( (/1,nslab,2,jmt-1/), mom_ramdisk_domains, xflags=GLOBAL_COMPUTE_DOMAIN, xhalo=0, yhalo=jbuf )
# endif
!     limit all domains to lie within 1:jmt
      if (jbuf .gt. 1)then
        mom_domains(:)%y%data%start_index = max(mom_domains(:)%y%data%start_index,1)
        mom_domains(:)%y%data%end_index = min(mom_domains(:)%y%data%end_index,jmt)
        mom_domains(:)%x%data%start_index = max(mom_domains(:)%x%data%start_index,1)
        mom_domains(:)%x%data%end_index = min(mom_domains(:)%x%data%end_index,imt)

        mom_global_domains(:)%y%data%start_index = max(mom_global_domains(:)%y%data%start_index,1)
        mom_global_domains(:)%y%data%end_index = min(mom_global_domains(:)%y%data%end_index,jmt)
        mom_global_domains(:)%x%data%start_index = max(mom_global_domains(:)%x%data%start_index,1)
        mom_global_domains(:)%x%data%end_index = min(mom_global_domains(:)%x%data%end_index,imt)

        mom_etat_domains(:)%y%data%start_index = max(mom_etat_domains(:)%y%data%start_index,1)
        mom_etat_domains(:)%y%data%end_index = min(mom_etat_domains(:)%y%data%end_index,jmt)
        mom_etat_domains(:)%x%data%start_index = max(mom_etat_domains(:)%x%data%start_index,1)
        mom_etat_domains(:)%x%data%end_index = min(mom_etat_domains(:)%x%data%end_index,imt)

# ifdef ramdrive
        mom_ramdisk_domains(:)%y%data%start_index = max(mom_ramdisk_domains(:)%y%data%start_index,1)
        mom_ramdisk_domains(:)%y%data%end_index = min(mom_ramdisk_domains(:)%y%data%end_index,jmt)
# endif
      endif

# ifdef max_window
!     these are domains for updating prognostic variables
      allocate( mom_prognostic_domains(0:npes-1) )
      mom_prognostic_domain => mom_prognostic_domains(pe)
      call mpp_define_domains( (/1,imt*km,2,jmt-1/), mom_prognostic_domains, xflags=GLOBAL_COMPUTE_DOMAIN, xhalo=0, yhalo=jbuf )
# endif
!
      calculated_rows = float(jmt-2)/num_processors
!
      do n = 1,num_processors
!
!       define task for processor "n"
!
        jstask3(n) = mom_domains(n-1)%y%data%start_index
        jetask3(n) = mom_domains(n-1)%y%data%end_index
        jscomp3(n) = mom_domains(n-1)%y%compute%start_index
        jecomp3(n) = mom_domains(n-1)%y%compute%end_index
!
# if defined max_window
        window_size3(n) = jecomp3(n)-jscomp3(n)+1 + 2*jbuf
# else
        window_size3(n) = jmw
# endif
        nr            = jecomp3(n) - jscomp3(n) + 1
	ncrows        = window_size3(n) - 2*jbuf
	num_loads3(n) = nr/ncrows
	if (mod(nr,ncrows) .ne. 0) num_loads3(n) = num_loads3(n)+1
!
!       clip the tasks at the domain boundaries
!
        jstask3(n) = max(jstask3(n),1)
	jetask3(n) = min(jetask3(n),jmt)
      enddo
!
      allocate (pe_of_lat(jmt))
!
      do j = 2,jmt-1
        do n = 1,num_processors
          if (j .ge. jscomp3(n) .and. j .le. jecomp3(n)) then
	    pe_of_lat(j) = n-1
	  endif
        enddo
      enddo
      pe_of_lat(1)   = pe_of_lat(2)
      pe_of_lat(jmt) = pe_of_lat(jmt-1)
!
      write (stdout,'(/,10x,a,i3,a/8(a13,2x))') &
      ' Assigning tasks for ',num_processors &
      ,' processor(s) as follows:'&
      , ' Processor # ', 'jstask   ', 'jetask   ',' jscomp   '&
      ,'jecomp   ','computed rows','  jmw    ','MW loads'

      do n = 1,num_processors
        write (stdout,'(8(4x,i4,7x))') &
        n, jstask3(n), jetask3(n), jscomp3(n), jecomp3(n)&
        ,jecomp3(n)-jscomp3(n)+1, window_size3(n), num_loads3(n)

        if (window_size3(n) .ne. window_size3(1)) then
	  write (stdout,*) '=>Error: the window size (jmw) must be the'&
          ,' same on all processors.'
          call abort()
	endif
      enddo
!
!-----------------------------------------------------------------------
!     set latitude row limits for this processor (pn)
!-----------------------------------------------------------------------
!
      jstask      = jstask3(pn)
      jetask      = jetask3(pn)
      num_loads   = num_loads3(pn)
      jscomp      = jscomp3(pn)
      jecomp      = jecomp3(pn)
      window_size = window_size3(pn)
!
      if (jstask .lt. 1 .or. jstask .gt. jmt) then
        write (stderr,*) '=>Error in define_tasks: jstask=',jstask
	call abort()
      endif
      if (jetask .lt. 1 .or. jetask .gt. jmt) then
        write (stderr,*) '=>Error in define_tasks: jetask=',jetask
	call abort()
      endif
      if (jscomp .lt. 2 .or. jscomp .gt. jmt-1) then
        write (stderr,*) '=>Error in define_tasks: jscomp=',jscomp
	call abort()
      endif
      if (jecomp .lt. 2 .or. jecomp .gt. jmt-1) then
        write (stderr,*) '=>Error in define_tasks: jecomp=',jecomp
	call abort()
      endif
!
#if defined ramdrive
!
!-----------------------------------------------------------------------
!     allocate space for 3-D and 2-D arrays on a ramdrive. Only the
!     latitude band required by each processor is allocated. 
!-----------------------------------------------------------------------
!
      allocate (ramdisk3d(nslab,jstask:jetask,2))
      allocate (ramdisk2d(imt*(jetask-jstask+1),nkflds))
#endif
      return
      end subroutine
    
end module
