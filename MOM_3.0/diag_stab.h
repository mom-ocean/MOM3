c====================== include file "stab.h" ==========================
c
c     CFL and other stability criteria information
c
c     cflons  = starting longitude (degrees) for stability tests
c     cflone  = ending longitude (degrees) for stability tests
c     cflats  = starting latitude (degrees) for stability tests
c     cflate  = ending latitude (degrees) for stability tests
c     cfldps  = starting depth (cm) for stability tests
c     cfldpe  = ending depth (cm) for stability tests
c
c     iscfl   = index corresponding to "cflons"
c     iecfl   = index corresponding to "cflone"
c     jscfl   = index corresponding to "cflats"
c     jecfl   = index corresponding to "cflate"
c     kscfl   = index corresponding to "cfldps"
c     kecfl   = index corresponding to "cfldpe"
c
c     cflcrt  = factor by which the cfl criteria must be exceeded in
c               order to show local values  (see blkdta.F)
c     maxcfl  = maximum number of times the "cflcrt" factor can be
c               exceeded before stopping.
c     numcfl  = counter for number of times the "cflcrt" factor was
c               exceeded.
c     cflum   = zonal velocity which comes closest to its cfl criteria
c     cflup   = percent of cfl criteria reached by "cflum"
c     icflu   = "i" coordinate of "cflum"
c     jcflu   = "j" coordinate of "cflum"
c     kcflu   = "k" coordinate of "cflum"
c
c     cflvm   = meridional velocity which comes closest to its cfl
c               criteria
c     cflvp   = percent of cfl criteria reached by "cflvm"
c     icflv   = "i" coordinate of "cflvm"
c     jcflv   = "j" coordinate of "cflvm"
c     kcflv   = "k" coordinate of "cflvm"
c
c     cflwtm  = vertical velocity at "t" box bottom
c                which comes closest to its cfl criteria
c     cflwtp  = percent of cfl criteria reached by "cflwtm"
c     icflwt  = "i" coordinate of "cflwtm"
c     jcflwt  = "j" coordinate of "cflwtm"
c     kcflwt  = "k" coordinate of "cflwtm"
c
c     cflwum  = vertical velocity at "u,v" box bottom
c                which comes closest to its cfl criteria
c     cflwup  = percent of cfl criteria reached by "cflwum"
c     icflwu  = "i" coordinate of "cflwum"
c     jcflwu  = "j" coordinate of "cflwum"
c     kcflwu  = "k" coordinate of "cflwum"
c
c     reynx   = maximum reynolds number in the zonal direction
c     ireynx  = "i" coordinate of "reynx"
c     jreynx  = "j" coordinate of "reynx"
c     kreynx  = "k" coordinate of "reynx"
c     reynu   = "u" for computing "reynx"
c     reynmu  = zonal mixing of momentum for computing "reynx"
c
c     reyny   = maximum reynolds number in the meridional direction
c     ireyny  = "i" coordinate of "reyny"
c     jreyny  = "j" coordinate of "reyny"
c     kreyny  = "k" coordinate of "reyny"
c     reynv   = "v" for computing "reyny"
c     reynmv  = meridional mixing of momentum for computing "reyny"
c
c     reynz   = maximum reynolds number in the vertical direction
c     ireynz  = "i" coordinate of "reynz"
c     jreynz  = "j" coordinate of "reynz"
c     kreynz  = "k" coordinate of "reynz"
c     reynw   = "w" for computing "reynz"
c     reynmw  = vertical mixing of momentum for computing "reynz"
c
c     peclx   = maximum peclet number in the zonal direction
c     ipeclx  = "i" coordinate of "peclx"
c     jpeclx  = "j" coordinate of "peclx"
c     kpeclx  = "k" coordinate of "peclx"
c     peclu   = "u" for computing "peclx"
c     peclmu  = zonal mixing of tracer for computing "peclx"
c
c     pecly   = maximum peclet number in the meridional direction
c     ipecly  = "i" coordinate of "pecly"
c     jpecly  = "j" coordinate of "pecly"
c     kpecly  = "k" coordinate of "pecly"
c     peclv   = "v" for computing "pecly"
c     peclmv  = meridional mixing of tracer for computing "pecly"
c
c     peclz   = maximum peclet number in the vertical direction
c     ipeclz  = "i" coordinate of "peclz"
c     jpeclz  = "j" coordinate of "peclz"
c     kpeclz  = "k" coordinate of "peclz"
c     peclw   = "w" for computing "peclz"
c     peclmw  = vertical mixing of tracer for computing "peclz"
c
c     tdig    = factor by which local tracer extremum must be exceeded
c               before showing ficticious creation of tracer
c
      common /stabr/ cflons, cflone, cflats, cflate, cfldps, cfldpe
      common /stabii/ iscfl,  iecfl,  jscfl,  jecfl,  kscfl,  kecfl
      common /stabr/ cflup, cflum, cflvp, cflvm
      common /stabii/ icflu, jcflu, kcflu, icflv, jcflv, kcflv
      common /stabr/ cflwtp, cflwtm, cflwup, cflwum 
      common /stabii/ icflwt, jcflwt, kcflwt, icflwu, jcflwu, kcflwu 
      common /stabr/ cflcrt, tdig 
      common /stabii/ numcfl, maxcfl 
      common /stabii/ ireynx, jreynx, kreynx, ireyny, jreyny, kreyny 
      common /stabr/ reynx, reynu, reynmu, reyny, reynv, reynmv
      common /stabii/ ireynz, jreynz, kreynz, ipeclx, jpeclx, kpeclx 
      common /stabr/ reynz, reynw, reynmw, peclx, peclu, peclmu
      common /stabii/ ipecly, jpecly, kpecly, ipeclz, jpeclz, kpeclz 
      common /stabr/ pecly, peclv, peclmv, peclz, peclw, peclmw
c
