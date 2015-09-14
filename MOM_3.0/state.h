c====================== include file "state.h" =========================
c
c     variables for equation of state
c
c     to = reference temperature for level
c     so = reference salinity for level
c     ro0= reference density for level
c     c  = polynomial coefficients for equation of state
c     tmink = min temperature at level k used for polynomial coeffs
c     tmaxk = max temperature at level k used for polynomial coeffs
c     smink = min salinity at level k used for polynomial coeffs
c     smaxk = max salinity at level k used for polynomial coeffs
c
      common /stater/ ro0(km), to(km), so(km), c(km,9)
      common /stater/ tmink(km), tmaxk(km), smink(km), smaxk(km)
