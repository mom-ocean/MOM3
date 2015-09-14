c===================== include file "river.h" ========================
c
c     Definition of the river inflow geometry
c
c     The fresh water flux of the rivers is simulated by an equivalent 
c     surface fresh water flux and is added to the usual surface fresh 
c     water flux. The tracer flux which comes with the fresh water is
c     added to the surface tracer flux. The momentum impact of the 
c     river may be added but is not included yet. 
c      
c     A river cannot be added in an isolated coastal box which is not
c     associated with active velocity boxes. So each river should be
c     distributed at least over two tracer boxes. A single tracer box 
c     river is not forbidden but yields a warning.
c
c
c     coded by:     m. schmidt
c
c     nriv      = number of rivers
c
c     nboxriv   = number of t-boxes for a river
c     iriv      = zonal indicees of river boxes t-boxes
c     jriv      = meridional indicees of river boxes t-boxes
c
c     areariv   = the horizontal surface of the boxes which simulate
c                 river flux
c     rff       = the equivalent fresh water flux of the river (cm/s)
c     rtf       = the tracer flux of the river
c 
      parameter (nriv = 1)
c
      parameter (nboxmax = 4)
      dimension nboxriv(nriv)
      common /criver/ areariv(nriv),rff(nriv),rtf(nriv,nt)
      data  nboxriv /4/
      dimension iriv(nboxmax, nriv), jriv(nboxmax,nriv)
      data  iriv /58, 58, 59, 59/
      data  jriv /28, 29, 28, 29/

