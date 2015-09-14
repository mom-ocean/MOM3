!====================== include file "pconst.h" ========================
!
!
!     rules for parameter constants
!
!     use prefix of "c" for whole real numbers (eg: c57 for 57.0)
!     use "m" after prefix to designate negative values (minus sign)
!       (eg: cm7 for -7.0)
!     use prefix of "p" for non repeating fractions (eg: p5 for 0.5)
!     use prefix of "r" for reciprocals (eg: r3 for 1/3.0)
!     combine use of prefix above and "e" for scientific notation, with
!       (eg: c5e4 for 5.0e4, c1em10 for 1.0e-10)
!
      parameter (c0=0.0, c1=1.0, c2=2.0, c3=3.0, c4=4.0, c5=5.0, c7=7.0)
      parameter (c8=8.0)
      parameter (c14=14.0, c16=16.0, c360=360.0)
      parameter (p125=0.125, p25=0.25, p5=0.5, p75=0.75)
      parameter (epsln=1.0e-20)
!
      parameter (c24=24.0, c60=60.0, c1440=1440.0)
      parameter (r24=c1/c24, r60=c1/c60, r1440=c1/c1440)
      parameter (secday=c1/(c60*c1440))

!     minimum allowable thickness (cm) of the top cell when using the
!     nonlinear free surface.        
      parameter (thickmin=10.0)      
