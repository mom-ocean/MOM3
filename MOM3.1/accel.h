c====================== include file "accel.h" =========================
c
c     depth dependent tracer timestep acceleration multipliers used to
c     hasten the convergence to equilibrium of the deeper portions of 
c     ocean-climate models.  
c
c     accelerate abyssal processes by varying the length of the tracer 
c     timestep with depth.  by using longer timesteps at depth, one can
c     in effect reduce the heat capacity of the deeper levels and speed
c     convergence to equilibrium.
c     note:
c     by applying this method, one is assuming that there is a single
c     steady-state solution to the model being considered.
c     also, since the diagnostic timestep calculations of "termbt" do 
c     not attempt to account for depth variant timestep lengths, the 
c     truncation error reported will increase, because it will include 
c     the tracer changes due to variations in "dtxcel".
c
c     reference:
c       Bryan, K., 1984: accelerating the convergence to equilibrium 
c     of ocean climate models, J. Phys. Oceanogr., 14, 666-673. 
c     ("dtxcel" here is the same as 1/gamma in the above reference)
c     set "dtxcel" to 1.0 at the surface and for upper levels not 
c     to be accelerated
c     set "dtxcel" to values greater than 1.0 at deeper levels to
c     accelerate convergence if above requirements are met
c
c     dtxcel   = model level dependent tracer timestep multipliers
#ifdef isoneutralmix
c     dtxsqr   = square root of "dtxcel" (used in computation of
c                maximum slope constraint for isopycnal mixing)
#endif
c#if !defined implicitvmix || defined isoneutralmix
c     dztxcl   = layer thickness divided by the timestep multiplier
c                (needed for convection code)
c     dzwxcl   = multiplication factor relating to the vertical
c                distance between ts points, scaled according
c                to timestep multipliers for use in convection code
c#endif
c     
c     author:      k. dixon      e-mail=> kd@gfdl.gov
c
c
      common /accel/ dtxcel(km)
#ifdef isoneutralmix
      common /accel/ dtxsqr(km)
#endif
c#if !defined implicitvmix || defined isoneutralmix
      common /accel/ dztxcl(km), dzwxcl(km)
c#endif
c
