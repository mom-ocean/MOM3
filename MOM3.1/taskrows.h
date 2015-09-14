!====================== include file "taskrows.h" ======================
!
!     max_tasks      = maximum number of processors
!     num_processors = requested number of processors
!     pn             = this processor number (from 1 to num_processors)
!
!     jstask3(n)  = southernmost buffer jrow for processor "n"
!     jetask3(n)  = northernmost buffer jrow for processor "n"
!     jscomp3(n)  = southernmost computed jrow for processor "n"
!     jecomp3(n)  = northernmost computed jrow for processor "n"
!     num_loads3(n)= number of times processor "n`s" MW is loaded
!     window_size3(n)= the size of the memory window for processor "n"
!
!     jstask  = southernmost buffer jrow for this processor (pn)
!     jetask  = northernmost buffer jrow for this processor (pn)
!     jscomp  = southernmost computed jrow for this processor (pn)
!     jecomp  = northernmost computed jrow for this processor (pn)
!     num_loads= number of times processor "pn`s" MW is loaded
!     window_size= the size of the memory window for this processor (pn)
!                  It is set as follows:
!                  if defined max_window then
!                    jecomp3-jscomp3+1 + 2*jbuf
!                    for second, fourth, and sixth order windows
!                    on processor "pn". The window_size must be the
!                    same on all processors. 
!                  else
!                    3 (for second order windows)
!                    5 (for fourth order windows)
!                    7 (for sixth order windows)
!                    regardless of the number of processors
!                  endif
!
      parameter (max_tasks = 2048)
      integer pn, window_size3, window_size
      common /taskrows/ pn
      common /taskrows/ num_processors
      common /taskrows/ jstask3(max_tasks), jetask3(max_tasks)
      common /taskrows/ jscomp3(max_tasks), jecomp3(max_tasks)
      common /taskrows/ num_loads3(max_tasks)
      common /taskrows/ jstask, jetask, num_loads, jscomp, jecomp
      common /taskrows/ window_size3(max_tasks), window_size
