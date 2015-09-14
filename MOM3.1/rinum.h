c====================== include file "rinum.h" =========================
c
c     variables for computation of richardson number
c
c     rhom1z = rho(k)-rho(k+1)
c     riu    = Richardson number at base of U-cells
c     rit    = Richardson number at base of T-cells
c
      common /rinum0/ rhom1z(imt,km,jmw)
      common /rinum0/ riu(imt,km,1:jmw-1), rit(imt,km,2:jmw-1)
c
