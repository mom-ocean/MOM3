!===================== include file "posdef.h" =========================
#ifdef posdef
!
!     Definition of the grid cells for a positive definite scheme.
!
!     To avoid overshootings in the tracer scheme the Reynolds number
!     must be smaller than 2. This may introduce a large numerical 
!     diffusion in the whole model. Select tracer cells for a separate 
!     treatment below. 
!
!     This file is an example and must be modified.
!     Copy this file to your working directory and edit the data segment.
!     Note: specify the grid cell index, not the index pointing on the 
!     cell boundary!
!
!
!     coded by:     m. schmidt  martin.schmidt@io-warnemuende.de
!
!
!     jpos      -  the number of rows which require posdef
!     ibmax     -  the maximum number of blocks of grid cells in one row
!     ibpos     -  the actual number of blocks of  grid cells in the row
!     jline     -  the meridional row index
!     ie, is    -  the start and stop index for the positive definite scheme.
!                  For ibpos < ibmax specify is, ie = 0.
!
!     Example data statement  /4,  2, 15,25, 35,40, 0,0, 0,0/
!                              | row number 4
!                                  | two grid cell blocks
!                                     | the first one from i=15 to i=25
!                                             | the second from i=35 to i=40
!                                                    | two pairs of dummies
!                                                      (ibmax is 4) 
!
      parameter (jpos  = 5)
      parameter (ibmax = 2)
      integer   jline(jpos), ibpos(jpos), is(jpos,ibmax), ie(jpos,ibmax)
      
      data (jline(jr),ibpos(jr),&
           (is(jr,ib),ie(jr,ib),ib=1,ibmax),jr=1,jpos)&
     / 2, 1,  2,8,  0,0&   
     , 3, 1,  2,8,  0,0&   
     , 4, 1,  2,8,  0,0&   
     , 5, 1,  2,8,  0,0&   
     , 6, 1,  2,8,  0,0&   
     /
#endif
