c====================== include file "ctavg.h" =========================
c
c     common for tracer averages within horizontal regions
c
c
c     sumbt = volume weighted tracer sum for a given region
c     avgbt = volume weighted tracer average for a given region
c     sumbk = volume weighted tracer sum for a region at a k level
c     avgbk = volume weighted tracer average for a region at a k level
c     sumgt = total global volume weighted tracer sum
c     avggt = total global volume weighted tracer average
c     sumgk = total global volume weighted tracer sum at a k level
c     avggk = total global volume weighted tracer average at a k level
c     sumbf = area weighted tracer flux sum for a given region
c     avgbf = average area weighted tracer flux for a region
c     sumgf = global area weighted tracer flux sum
c     avggf = global average area weighted tracer flux
c
c
      common /tavgr/ sumbk(nhreg,km,nt), sumbt(nhreg,nt), sumgk(km,nt)
      common /tavgr/ sumgt(nt), sumbf(nhreg,nt), sumgf(nt)
      common /tavgr/ avgbk(nhreg,km,nt), avgbt(nhreg,nt), avggk(km,nt)
      common /tavgr/ avggt(nt), avgbf(nhreg,nt), avggf(nt)

