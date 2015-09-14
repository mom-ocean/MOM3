!====================== include file "grdvar.h" ========================
!
!     variables which are functions of the grid defined by "coord.h"
!
!     dxt     = longitudinal width of "t" grid box at the 
!               equator (in cm)
!     dxtr    = reciprocal of "dxt"
!     dxt2r   = reciprocal of "2*dxt"
!     dxt4r   = reciprocal of "4*dxt"
!     dxu     = longitudinal width of "u,v" grid box at the 
!               equator (in cm)
!     dxur    = reciprocal of "dxu"
!     dxu2r   = reciprocal of "2*dxu"
!     dxu4r   = reciprocal of "4*dxu"
!     dxmetr  = reciprocal of "(dxt(i)+dxt(i+1))"
!     duw     = xu(i) - xt(i)
!     due     = xt(i+1) - xu(i)
!     dus     = yu(jrow) - yt(jrow)
!     dun     = yt(jrow+1) - yu(jrow)
!
!     dyt     = latitudinal height of "t" grid box (in cm)
!     dytr    = reciprocal of "dyt"
!     dyt2r   = reciprocal of "2*dyt"
!     dyt4r   = reciprocal of "4*dyt"
!     dyu     = latitudinal height of "u,v" grid box (in cm)
!     dyur    = reciprocal of "dyu"
!     dyu2r   = reciprocal of "2*dyu"
!     dyu4r   = reciprocal of "4*dyu"
!
!     csu     = cosine of "u,v" grid point latitude
!     csur    = reciprocal of "csu"
!     cst     = cosine of "t" grid point latitude
!     cstr    = reciprocal of "cst"
!     phi     = latitude of "u,v" grid point in radians
!     phit    = latitude of "t" grid point in radians
!     sine    = sine of "u,v" grid point latitude
!     tng     = tan of "u,v" grid point latitude
!     fcor    = 2*omega*sine(j)
!
!     c2dzt(k)= "2*dzt"
!     dztr(k) = reciprocal of dzt ("t" cell vertical resolution)
!     dzt2r(k)= reciprocal of "2*dzt"
!     dzwr(k) = reciprocal of dzw ("w" cell vertical resolution)
!     dzw2r(k)= reciprocal of "2*dzw"
!     dztur(k)= upper diffusion grid factor = 1.0/(dzw(k-1)*dzt(k))
!     dztlr(k)= lower diffusion grid factor = 1.0/(dzw(k)*dzt(k))
!     dzwur(k)= upper diffusion grid factor = 1.0/(dzt(k)*dzw(k))
!     dzwlr(k)= lower diffusion grid factor = 1.0/(dzt(k+1)*dzw(k))
!     fracdz(k,0) = fractional distance between grid point and cell top
!     fracdz(k,1) = fractional distance between grid point and cell bot
!     
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
!
#ifdef quicker
      common /grdvar/ quick_x(imt,2), curv_xp(imt,3), curv_xn(imt,3)
      common /grdvar/ quick_y(jmt,2), curv_yp(jmt,3), curv_yn(jmt,3)
      common /grdvar/ quick_z(km,2),  curv_zp(km,3),  curv_zn(km,3)
#endif
