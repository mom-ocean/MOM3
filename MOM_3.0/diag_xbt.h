c====================== include file "diag_xbt.h" ======================
c
c     Each XBT station is located at latitude "xbtlat" and longitude 
c     "xbtlon". Data is collected at each grid point from the first
c     level down through the nearest model level corresponding to a
c     depth of "xbtdpt" cm. Actually, all coordinates are converted to
c     the nearest model temperature grid point.
c
c     All basic quantities as well as all terms in the momentum,
c     temperature, and salinity equations are averaged over the time
c     period specified by "xbtint".
c     
c     The maximum number of XBTs may be increased by changing parameter
c     "maxxbt" below.
c
c     inputs:
c
c     maxxbt  = maximum number of XBTs allowed.
c     kmxbt   = maximum number of levels from surface downward (<=km)
c               set kmxbt  < km to save space
c     xbtlat  = real latitude of XBTs in degrees
c     xbtlon  = real longitude of XBTs in degrees
c     xbtdpt  = real depths of XBTs in cm
c     items   = number of items in the XBT
c     xname   = character*12 names of XBT quantities
c
c     outputs:
c
c     numxbt  = actual number of XBTs used
c     nxbtts  = current number of time steps in accumulated XTB data 
c
c     ixbt    = longitude index of nearest model temperature grid point
c               corresponding to "xbtlon"
c     jxbt    = latitude index of nearest model temperature grid point
c               corresponding to "xbtlat"
c     kxbt    = depth index of nearest model temperature grid point
c               corresponding to "xbtdpt"
c     nsxbt   = starting number for the XBTs on each latitude
c     nexbt   = ending number for the XBTs on each latitude
c
c     txbt    = accumulator array for time rate of change of
c                tracers. the total time rate of change
c                is broken down into components as follows:
c                the form is d( )/dt = terms (2) ... (10) where each
c                term has the units of "tracer units/sec" using
c                schematic terms for illustration.
c                (1)  = total time rate of change for the tracer
c                (2)  = change due to zonal nonlinear term (UT)x
c                (3)  = change due to meridional nonlinear term (VT)y
c                (4)  = change due to vertical nonlinear term (WT)z
c                (5)  = change due to zonal diffusion: Ah*Txx
c                (6)  = change due to meridional diffusion: Ah*Tyy
c                (7)  = change due to vertical diffusion:  kappa_h*Tzz
c                (8)  = change due to source term
c                (9)  = change due to explicit convection
c                (10) = change due to filtering
c     the nonlinear terms can be broken into two parts: advection and a
c     continuity part: The physically meaningful part is advection. 
c     eg: Zonal advection of tracer "A" is -U(A)x = A(Ux) - (UA)x
c                (11) = zonal advection U(Ax)
c                (12) = meridional advection V(Ay)
c                (13) = vertical advection W(Az)
c                (14) = change of tracer variance
c                (15) = average tracer within volume (tracer units)
c     uxbt    = accumulator array for time rate of change of
c                momentum. the total time rate of change
c                is broken down into components as follows:
c                the form is d( )/dt = terms (2) ... (13) where each
c                term has the units of "cm/sec**2" and "Q" is the
c                momentum component {zonal or meridional} using
c                schematic terms for illustration.
c                (1)  = total time rate of change for the momentum
c                (2)  = change due to the pressure gradient: grad_p
c                       without the surface pressure gradients
c                       (i.e., for computing the internal modes)
c                (3)  = change due to zonal nonlinear term: (UQ)x
c                (4)  = change due to meridional nonlinear term: (VQ)y
c                (5)  = change due to vertical nonlinear term: (wQ)z
c                (6)  = change due to zonal viscosity: Am*Qxx
c                (7)  = change due to meridional viscosity: Am*Qyy
c                (8)  = change due to vertical viscosity: kappa_m*Qzz
c                (9)  = change due to metric terms
c                (10) = change due to coriolis terms: fQ
c                (11) = change due to source terms
c                (12) = change due to surface pressure gradient
c                       this is obtained after solving the external mode
c                       in the stream function technique. It is solved
c                       directly from the elliptic equation for the 
c                       prognostic surface pressure technique
c                (13) = change due to metric advection
c     the nonlinear terms can be broken into two parts: advection and a
c     continuity part: The physically meaningful part is advection. 
c     eg: Zonal advection of vel component "Q" is -U(Q)x = Q(U)x - (UQ)x
c                (14) = zonal advection U(Qx)
c                (15) = meridional advection V(Qy)
c                (16) = vertical advection W(Qz)
c                (17) = average velocity component
c     xbtw    = accumulator array for vertical velocity. (cm/sec)
c               this is the average of adv_vbu at top and bottom of cell
c     txbtsf  = accumulator array for tracer surface flux terms.
c                tracer (#1,#2) units = (cal/cm**2/sec, gm/cm**2/sec)
c     uxbtsf  = accumulator array for wind stress terms. (dynes/cm**2)
c
c     ntxbt   = number of terms for tracers
c     nuxbt   = number of terms for velocity
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
      parameter (maxxbt=3, kmxbt=4)
      parameter (ntxbt=15, nuxbt=17)
c
      common /cxbt0/ nxbtts, numxbt 
c
      character*12 xnamet, xnameu, xnamex
      common /cxbtc/ xnamet(ntxbt), xnameu(nuxbt,2), xnamex(4)
      common /cxbtr/ xbtlat(maxxbt), xbtlon(maxxbt), xbtdpt(maxxbt)
      common /cxbtr/ txbt(kmxbt,ntxbt,nt,maxxbt), txbtsf(nt,maxxbt)
      common /cxbtr/ uxbt(kmxbt,nuxbt,2,maxxbt), uxbtsf(2,maxxbt)
      common /cxbtr/ xbtw(kmxbt,maxxbt)
      common /cxbti/ ixbt(maxxbt), jxbt(maxxbt), kxbt(maxxbt)
      common /cxbti/ nsxbt(jmt), nexbt(jmt)
