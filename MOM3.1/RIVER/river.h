!===================== include file "river.h" ========================
!
!     Definition of the river inflow geometry
!
!     The fresh water flux of the rivers is simulated by an equivalent 
!     surface fresh water flux and is added to the usual surface fresh 
!     water flux. The tracer flux which comes with the fresh water is
!     added to the surface tracer flux. The momentum impact of the 
!     river may be added but is not included yet. 
!      
!     A river cannot be added in an isolated coastal box which is not
!     associated with active velocity boxes. So each river should be
!     distributed at least over two tracer boxes. A single tracer box 
!     river is not forbidden but yields a warning.
!
!
!     coded by:     m. schmidt
!
!     nriv        = total number of rivers to be used
!     ntriv       = number of tracers specified for the rivers
!     nboxmax     = maximum number of t-boxes allocated to a river
!                   set active boxes by indices > 0
!     nboxriv     = number of active river t-boxes (<= nboxmax)
!     iriv        = zonal indices of river t-boxes
!     jriv        = meridional indices of river t-boxes
!     riverindex  = index to identify the river 
!     rivername   = string to identify the river            
!                   river data files are assumed to be "rivername.dat"
!     rivertrname = string to identify the ntriv tracers             

      parameter (nriv = 15)
# ifndef simple_sbc
      parameter (ntriv = 2)
# endif
      parameter (nboxmax = 4)
!
      integer   nboxriv(nriv)
      integer   iriv(nriv,nboxmax), jriv(nriv,nboxmax), riverindex(nriv)
!
      character rivername(nriv) *16
      data (riverindex(nr),nboxriv(nr),
     &     (iriv(nr,nb),jriv(nr,nb),nb=1,nboxmax),
     &	   rivername(nr), nr=1,nriv)		       
     & / 1, 2,  97,98,  98,98, 0,0, 0,0, 'Kemijoki'    
     &,  2, 2,  89,99,  90,99, 0,0, 0,0, 'Lulealv'
     &,  3, 2,  74,78,  74,79, 0,0, 0,0, 'Angermansalv'
     &,  4, 2,  82,86,  83,86, 0,0, 0,0, 'Umealv'
     &,  5, 2,  87,71,  87,72, 0,0, 0,0, 'Kokemaenjoki'
     &,  6, 2, 109,56, 109,57, 0,0, 0,0, 'Narva'
     &,  7, 2, 115,60, 115,61, 0,0, 0,0, 'Neva'
     &,  8, 2,  62,19,  63,19, 0,0, 0,0, 'Oder'
     &,  9, 2,  78,22,  79,22, 0,0, 0,0, 'Weichsel'
     &, 10, 2,  86,31,  86,32, 0,0, 0,0, 'Njemen'
     &, 11, 2,  65,35,  66,35, 0,0, 0,0, 'Helgean'
     &, 12, 2,  69,40,  69,41, 0,0, 0,0, 'Eman-Motal'
     &, 13, 2,  75,55,  75,56, 0,0, 0,0, 'Maelaren'
     &, 14, 2,  93,39,  94,39, 0,0, 0,0, 'Daugawa'
     &, 15, 2,  55,44,  55,45, 0,0, 0,0, 'Goetaaelv'
     & /
# ifndef simple_sbc
      character rivertrname(ntriv) *12
      data (rivertrname(n),n=1,ntriv) /'temperature','salinity'/
# endif
      logical monthly_data
      data monthly_data /.true./
