!====================== include file "stdunits.h" ======================
!
!     stdin  = unit number for standard input.
!     stdout = unit number for standard output.
!     stderr = unit number for standard error.
!
      integer stdin, stdout, stderr
      parameter (stdin = 5, stderr = 0)
      common /stdunit/ stdout
