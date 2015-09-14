!====================== include file "aparam.h" ========================
!
!     parameter file which sets the atmospheric characteristics needed
!     for coupling to MOM
!
!     ima   = number of grid points in the longitudinal direction
!              for the atmospheric S.B.C. grid
!     jma   = number of grid points in the latitudinal direction
!              for the atmospheric S.B.C. grid
!     imap2 = ima + two extra points for the cyclic condition needed
!             for interpolating S.B.C. to MOM
!     
!
      integer, parameter :: ima=36, jma=18, imap2=ima+2
