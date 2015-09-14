#if defined test_task_simulator
      program task_simulator
!
!=======================================================================
!     using 1 processor, simulate how tasks are divided among one
!     or more processors using 2nd order, 4th order, and 6th order 
!     memory windows.
!     
!     author: R. C. Pacanowski   e-mail -> rcp@gfdl.gov
!=======================================================================
!
      integer, parameter :: stdout=6
# include "taskrows.h"
      character*100 options

      jmt = 66

      options = "options: "
# if defined second_order_window
      jmw = 3; jbuf=1; options = trim(options) // " second_order_window"
# endif
# if defined fourth_order_window
      jmw = 5; jbuf=2; options = trim(options) // " fourth_order_window"
# endif
# if defined sixth_order_window
      jmw = 7; jbuf=3;  options = trim(options) // " sixth_order_window"
# endif
# if defined max_window
!     jmw is set appropriately for 1 processor here.
      jmw = (jmt-2) + 2*jbuf
      options = trim(options) // " with fully_opened_window"
# endif

# if defined pressure_gradient_average && !defined max_window
      write (stdout,*) '=>Error: option pressure_gradient_average " 
     &,' requires option max_window'
      stop
# endif
      do n=1,4
        num_processors = 2**(n-1) 
	pn = 1

        write (stdout,'(///5x,a,i3,a,a,i3,a/)')
     &  'Testing ',num_processors, ' processor(s). '
     &, ' jmt=',jmt, ', '//trim(options)

        call simulate_define_tasks (jmw, jbuf, jmt, options)
        do pn=1,num_processors
	  jstask      = jstask3(pn)
	  jetask      = jetask3(pn)
	  jscomp      = jscomp3(pn)
	  jecomp      = jecomp3(pn)
	  num_loads   = num_loads3(pn)
	  window_size = window_size3(pn)
          call simulate_memory_window (jmw, jbuf, jmt)
        enddo
      enddo
      stop
      end


      subroutine simulate_memory_window (jmw, jbuf, jmt)
      integer, parameter :: stdout=6
# include "taskrows.h"

      write (stdout,'(/6(a,i3)/)') ' --> Details for processor #',pn

      num_rows_to_copy = 2*jbuf
      do mw=1,num_loads

        ncrows = window_size-2*jbuf
	joff = jscomp - jbuf - 1 + (mw-1)*ncrows
        js = (jstask - joff)*(1/mw) + (1 + num_rows_to_copy)*(1-1/mw)
	je = min(window_size,jetask-joff)

        js_calc = 1+jbuf
        je_calc = ncrows+jbuf

        call show_memory_window_move (mw, jmw, jbuf, jmt)

        write (stdout,'(15(a,i3))')
     & 'MW load #',mw,': js=',js,' to je=',je
     &,' (jrows=',js+joff,' to',je+joff,'), joff=',joff
     &,', solve on MW rows ',js_calc,' to'
     &, je_calc,' (jrows=',js_calc+joff,' to',je_calc+joff,')'
      enddo
      return
      end

      subroutine show_memory_window_move (mw, jmw, jbuf, jmt)
      integer, parameter :: stdout=6
# include "taskrows.h"
      if (mw .eq. 2) then
        num_rows_to_copy = 2*jbuf
        do j=1,num_rows_to_copy
          jfrom = jmw - num_rows_to_copy + j
	  jto   = j
	  write (stdout,'(5x,a,i3,a,i3)')
     &   'Copying data from MW row ',jfrom,' ->',jto
	enddo
      elseif (mw .eq. 3) then
	write (stdout,'(5x,a,a)')
     & 'Copying data within the MW is the same as done previously from '
     &,'here on (not showing it...)'
      endif
      return
      end

      subroutine simulate_define_tasks (jmw, jbuf, jmt)
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
      integer, parameter :: stdout=6
# include "taskrows.h"
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
      calculated_rows = float(jmt-2)/num_processors
!
      do n = 1,num_processors
!
!       define task for processor "n"
!
        jstask3(n)      = nint((n-1)*calculated_rows) + (2-jbuf)
        jetask3(n)      = nint(n*calculated_rows) + (1+jbuf)
	jscomp3(n)      = jstask3(n) + jbuf
	jecomp3(n)      = jetask3(n) - jbuf
# if defined max_window
        window_size3(n) = jecomp3(n)-jscomp3(n)+1 + 2*jbuf
# else
        window_size3(n) = jmw
# endif
        if (window_size3(n) .ge. jecomp3(n)-jscomp3(n)+1) then
          num_loads3(n) = 1
	else
          nr            = jecomp3(n) - jscomp3(n) + 1
	  ncrows        = window_size3(n) - 2*jbuf
	  num_loads3(n) = nr/ncrows
	  if (mod(nr,ncrows) .ne. 0) num_loads3(n) = num_loads3(n)+1
	endif
!
!       limit the tasks at the domain boundaries
!
        jstask3(n) = max(jstask3(n),1)
	jetask3(n) = min(jetask3(n),jmt)
      enddo

      write (stdout,'(/,10x,a,i3,a/8(a13,2x))')
     & ' Assigning tasks for ',num_processors 
     &,' processor(s) as follows:'
     &, ' Processor # ', 'jstask   ', 'jetask   ',' jscomp   '
     &,'jecomp   ','computed rows','  jmw    ','MW loads'

      do n = 1,num_processors
        write (stdout,'(8(4x,i4,7x))') 
     &  n, jstask3(n), jetask3(n), jscomp3(n), jecomp3(n)
     &, jecomp3(n)-jscomp3(n)+1, window_size3(n), num_loads3(n)

        if (window_size3(n) .ne. window_size3(1)) then
	  write (stdout,*) '=>Error: the window size (jmw) must be the'
     &,   ' same on all processors.'
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

      return
      end

#else
      subroutine dummy_task_simulator
      return
      end
#endif
