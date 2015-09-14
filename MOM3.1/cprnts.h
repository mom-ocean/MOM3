c====================== include file "cprnts.h" =======================
c
c     variables used for controlling matrix printouts during diagnostic
c     timesteps:
c
c     nlatpr = maximum number of latitudes for matrix printouts
c              on diagnostic time steps
c     prlat  = latitudes (deg) at which (x,z) printouts are desired
c
c     start & end coordinates for matrix printouts of (x,z) sections
c
c     prslon = starting longitudes (deg)
c     prelon = ending longitudes (deg)
c     prsdpt = starting depths  (cm)
c     predpt = ending depths  (cm)
c
c     start & end coordinates for matrix printouts of (x,y) sections
c
c     slonxy = starting longitude (deg)
c     elonxy = ending longitude (deg)
c     slatxy = starting latitude (deg)
c     elatxy = ending latitude (deg)
c
c     matrix printouts of (y,z) sections will use above coordinates
c
      common /cprnts/ prlat(nlatpr), prslon(nlatpr), prelon(nlatpr)
      common /cprnts/ prsdpt(nlatpr), predpt(nlatpr)
      common /cprnts/ slatxy, elatxy, slonxy, elonxy
c
