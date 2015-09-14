c====================== include file "cshort.h" ========================
c
c     pen    = double exponential penetration function defined at
c              bottoms of "T" boxes
c
c              note: pen(0) is set 0.0 instead of 1.0 to compensate for
c                    the shortwave part of the total surface flux "stf"
c
c     divpen = divergence of penetration defined at the center of
c              "T" boxes
c
      common /cshort/ pen(0:km), divpen(km)
