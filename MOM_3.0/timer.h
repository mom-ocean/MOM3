c====================== include file "timer.h" =========================
c
c     common block for timing routines
c
      parameter (ntimers=64, nsubtimers=45)
      parameter (maxstack=6)
      parameter (namelength=40)
      character*40 timer(0:ntimers)
      character*40 subtimer(0:ntimers, 1:nsubtimers)
      dimension starttime(0:ntimers, 1:nsubtimers)
      dimension elapsedtime(0:ntimers, 1:nsubtimers)
      dimension ncalls(0:ntimers, 1:nsubtimers)
      dimension nresumes(0:ntimers, 1:nsubtimers)
      dimension nested_calls(0:ntimers, 1:nsubtimers)
      integer weighted_calls(0:ntimers, 1:nsubtimers)
      dimension npops(0:ntimers, 1:nsubtimers)
      dimension npushes(0:ntimers, 1:nsubtimers)
      integer active(0:ntimers), subtimer_active(0:ntimers,1:nsubtimers)
      logical group_overflow(0:ntimers), ngroup_overflow
      logical print_timer_warnings, print_tictocs
      dimension jstack(0:ntimers,1:maxstack)
      dimension istackptr(0:ntimers)

      common /timers1/ timer, subtimer
      common /timers2/ starttime, elapsedtime
      common /timers3/ ncalls, active, subtimer_active, ntictoc_errors
      common /timers3/ nested_calls, weighted_calls, npops, nresumes
      common /timers3/ npushes
      common /timers3/ jstack, istackptr
      common /timers3/ max_timer_warnings
      common /timers4/ group_overflow, ngroup_overflow
      common /timers4/ print_timer_warnings, print_tictocs
      common /timers5/ prev_second

c-----------------------------------------------------------------------

