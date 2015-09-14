c===================== include file "cnep.h" ===========================
c
c     the option "neptune" provides a subgridscale parameterization 
c     for the interaction of eddies and topography
c
c     reference:
c       Holloway, G., 1992: representing topographic stress for large
c     scale ocean models, J. Phys. Oceanogr., 22, 1033-1046
c
c     neptune is calculated as an equilibrium streamfunction given by 
c     pnep=-f*snep*snep*hnep and is applied through eddy viscosity
c
c     hnep = model streamfunction depth
c     snep = spnep + (senep-spnep)*(0.5 + 0.5*cos(2.0*latitude))
c       the neptune length scale snep has a value of senep at the
c       equator and smoothly changes to spnep at the poles
c
c     variables used in applying neptune
c
c     spnep = neptune length scale at the pole
c     senep = neptune length scale at the equator
c
      common /cnep/ spnep, senep
c
