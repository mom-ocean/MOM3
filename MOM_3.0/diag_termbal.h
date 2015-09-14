c====================== include file "diag_termbal.h" =================
c
c     term balances are instantaneous breakdowns of all terms in the
c     momentum & tracer equations. They are averaged over ocean volumes
c     defined by horizontal and vertical regional masks:
c
c
c     termbt   = term balance components for time rate of change of
c                tracers within a  volume. the total time rate of change
c                is broken down into components as follows:
c                the form is d( )/dt = terms (2) ... (10) where each
c                term has the units of "tracer units/sec" using
c                schematic terms for illustration.
c                (1)  = total time rate of change for the tracer
c                (2)  = change due to zonal nonlinear term: (UT)x
c                (3)  = change due to meridional nonlinear term: (VT)y
c                (4)  = change due to vertical nonlinear term: (WT)z
c                (5)  = change due to zonal diffusion: Ah*Txx
c                (6)  = change due to meridional diffusion: Ah*Tyy
c                (7)  = change due to vertical diffusion:  kappa_h*Tzz
c                (8)  = change due to source term
c                (9)  = change due to explicit convection
c                (10) = change due to filtering
c     the nonlinear terms can be broken into two parts: advection and a
c     continuity part: The physically meaningful part is advection. 
c     eg: Zonal advection of tracer "A" is -U(A)x = A(U)x - (UA)x
c                (11) = zonal advection U(Ax) 
c                (12) = meridional advection V(Ay)
c                (13) = vertical advection W(Az) 
c                (14) = change of tracer variance (tracer**2 units)
c                (15) = average tracer within volume (tracer units)
c     terr     = error term = (1) - sum (2) ... (10)
c     asst     = average sea surface tracer for regional surface areas
c     stflx    = average surface tracer flux for regional surface areas
c                tracer (#1,#2) units = (cal/cm**2/sec, gm/cm**2/sec)
c
c     termbm   = term balance components for time rate of change of
c                momentum within a volume. the total time rate of change
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
c                (9)  = change due to metric diffusion terms
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
c     smflx    = average surface momentum flux for regional surf areas
c                in dynes/cm**2
c     avgw     = average vertical velocity (cm/sec)
c
c     ustf     = names & units for surface tracer fluxes
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
      parameter (ntterms=15, nuterms=17)
c
      character*15 ustf(nt,2)
      common /termb1/ ustf
      common /termb2/asst(nt,0:nhreg), avgw(numreg)
      common /termb2/termbt(0:km,ntterms,nt,0:numreg)
      common /termb2/termbm(0:km,nuterms,2,numreg)
      common /termb2/smflx(2,0:nhreg)
      common /termb2/stflx(nt,0:nhreg), terr(nt)
