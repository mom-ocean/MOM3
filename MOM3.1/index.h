c====================== include file "index.h" =========================
c
c    starting & ending indices for controlling where quantities
c    are filtered. this removes the time step restriction due to
c    convergence of meridians.
c
#if defined fourfil || defined firfil
c
c     define latitudinal domain for filtering
c
c     filter_reflat_n = filtering reference latitude in the northern
c                       hemisphere. all latitudes to the north of this
c                       latitude will be filtered to remove length 
c                       scales smaller than dx*cos(filter_reflat_n).
c
c     filter_reflat_s = filtering reference latitude in the southern
c                       hemisphere. all latitudes to the south of this
c                       latitude will be filtered to remove length 
c                       scales smaller than dx*cos(filter_reflat_s).
c     rjfrst          = southern latitude below which no filtering
c                       is to be done. This is used to eliminate
c                       wasteful filtering within antarctica.
c
c    istf  = starting index for filtering T grid points
c    ietf  = ending   index for filtering T grid points
c    isuf  = starting index for filtering U grid points
c    ieuf  = ending   index for filtering U grid points
c    iszf  = starting index for filtering "vorticity" grid points
c                     (T cells at least two cells away from land)
c    iezf  = ending   index for filtering "vorticity" grid points
c                     (T cells at least two cells away from land)
c
c     lsegf  = max number of longitudinal strips per filtering latitude 
c     jmtfil = max number of latitudes to be filtered
c     numflt = specifies the number of filter applications to use
c              for the tracers, vorticity, or divergence when using the
c               finite impulse response filter
c     numflu = specifies the number of filter applications to use
c              for velocities when using the finite impulse response
c               filter
c
c
c
      parameter (lsegf=7, jmtfil=23)
c
      common /indexr/ refcosn, refcoss
      common /indexr/ rjfrst
      common /indexi/  jfrst,  jft1,  jft2,  jfu1,  jfu2
      common /indexi/  jskpt, jskpu, njtbft, njtbfu
# ifdef firfil
      common /indexi/   numflt(jmt), numflu(jmt)
# endif
      common /indexi/ istf(jmtfil,lsegf,km), ietf(jmtfil,lsegf,km)
      common /indexi/ isuf(jmtfil,lsegf,km), ieuf(jmtfil,lsegf,km)
      common /indexi/ iszf(jmtfil,lsegf),    iezf(jmtfil,lsegf)
#endif
      common /indexr/ filter_reflat_n, filter_reflat_s
c
