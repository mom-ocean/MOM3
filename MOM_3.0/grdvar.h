c====================== include file "grdvar.h" ========================
c
c     variables which are functions of the grid defined by "coord.h"
c
c     dxt     = longitudinal width of "t" grid box at the 
c               equator (in cm)
c     dxtr    = reciprocal of "dxt"
c     dxt2r   = reciprocal of "2*dxt"
c     dxt4r   = reciprocal of "4*dxt"
c     dxu     = longitudinal width of "u,v" grid box at the 
c               equator (in cm)
c     dxur    = reciprocal of "dxu"
c     dxu2r   = reciprocal of "2*dxu"
c     dxu4r   = reciprocal of "4*dxu"
c     dxmetr  = reciprocal of "(dxt(i)+dxt(i+1))"
c     duw     = xu(i) - xt(i)
c     due     = xt(i+1) - xu(i)
c     dus     = yu(jrow) - yt(jrow)
c     dun     = yt(jrow+1) - yu(jrow)
c
c     dyt     = latitudinal height of "t" grid box (in cm)
c     dytr    = reciprocal of "dyt"
c     dyt2r   = reciprocal of "2*dyt"
c     dyt4r   = reciprocal of "4*dyt"
c     dyu     = latitudinal height of "u,v" grid box (in cm)
c     dyur    = reciprocal of "dyu"
c     dyu2r   = reciprocal of "2*dyu"
c     dyu4r   = reciprocal of "4*dyu"
c
c     csu     = cosine of "u,v" grid point latitude
c     csur    = reciprocal of "csu"
c     cst     = cosine of "t" grid point latitude
c     cstr    = reciprocal of "cst"
c     phi     = latitude of "u,v" grid point in radians
c     phit    = latitude of "t" grid point in radians
c     sine    = sine of "u,v" grid point latitude
c     tng     = tan of "u,v" grid point latitude
c     fcor    = 2*omega*sine(j)
c
c     c2dzt(k)= "2*dzt"
c     dztr(k) = reciprocal of dzt ("t" cell vertical resolution)
c     dzt2r(k)= reciprocal of "2*dzt"
c     dzwr(k) = reciprocal of dzw ("w" cell vertical resolution)
c     dzw2r(k)= reciprocal of "2*dzw"
c     dztur(k)= upper diffusion grid factor = 1.0/(dzw(k-1)*dzt(k))
c     dztlr(k)= lower diffusion grid factor = 1.0/(dzw(k)*dzt(k))
c     dzwur(k)= upper diffusion grid factor = 1.0/(dzt(k)*dzw(k))
c     dzwlr(k)= lower diffusion grid factor = 1.0/(dzt(k+1)*dzw(k))
c     fracdz(k,0) = fractional distance between grid point and cell top
c     fracdz(k,1) = fractional distance between grid point and cell bot
c     
      common /grdvar/ dxt(imt),  dxtr(imt),  dxt2r(imt), dxu(imt)
      common /grdvar/ dxur(imt), dxu2r(imt), dxu4r(imt), dxt4r(imt)
      common /grdvar/ dyt(jmt),  dytr(jmt),  dyt2r(jmt), dyu(jmt)
      common /grdvar/ dyur(jmt), dyu2r(jmt), dyu4r(jmt), dyt4r(jmt)
      common /grdvar/ csu(jmt),  csur(jmt),  cst(jmt),   cstr(jmt)
      common /grdvar/ cstdytr(jmt), cstdyt2r(jmt)
      common /grdvar/ csudyur(jmt), csudyu2r(jmt)
      common /grdvar/ cst_dytr(jmt), csu_dyur(jmt)
      common /grdvar/ phi(jmt),  phit(jmt),  sine(jmt),  tng(jmt)
      common /grdvar/ fcor(jmt), c2dzt(km),  dztr(km),   dzt2r(km) 
      common /grdvar/ dzwr(0:km), dzw2r(0:km)
      common /grdvar/ dxmetr(imt), duw(imt), due(imt)
      common /grdvar/ dun(jmt), dus(jmt)
      common /grdvar/ fracdz(km,0:1)
#ifdef tcvmix
      common /grdvar/ dzwur(km), dzwlr(km)
#endif
#if defined implicitvmix || defined isoneutralmix
      common /grdvar/ dztur(km), dztlr(km)
#endif
c
#ifdef quicker
      common /grdvar/ quick_x(imt,2), curv_xp(imt,3), curv_xn(imt,3)
      common /grdvar/ quick_y(jmt,2), curv_yp(jmt,3), curv_yn(jmt,3)
      common /grdvar/ quick_z(km,2),  curv_zp(km,3),  curv_zn(km,3)
#endif
