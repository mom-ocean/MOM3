c====================== include file "diag_region.h" ==================
c
c     variables used for computing regional tracer averages (see
c     "diag_region.F" & sub "diagt1" in "tracer.F") and for computing 
c     term  balances for tracer and momentum equations (see 
c     "baroclinic.F",  "tracer.F" and "diag2.F")
c
c
c     mskvr  =  mask field defining regions in the vertical
c              (eg: mskvr(k) = m indicates all points under a horizontal
c               mask at level "k" are in the "mth" vertical region
c               where m=1..nvreg)
c               The "mskvr" masks are used in "diag.F", but not
c               in "diagt1", where tracer averages are calculated
c               for each k-level.
c
c     hregnm = horizontal region name
c     vregnm = vertical region name
c     volbt  = total volume under a given horizontal region
c     volbk  = volume contained in a horizontal region at level "k"
c     volgt  = total ocean volume
c     volgk  = total ocean volume at level "k"
c     areab  = total ocean surface area for a given horizontal region
c     areag  = total ocean surface area
c
c     volt   = ocean volume within a particular horizontal & vertical
c              region (on the "t" grid) for tracer term balances
c              volt(0) represents the sum of all regions
c     rvolt  = 1/volt ( 0.0 if volt = 0.0)
c     areat  = horizontal ocean surface area corresponding to "volt"
c              areat(0) represents the sum of all regions
c     rareat = 1/areat ( 0.0 if areat = 0.0)
c     volu   = ocean volume within a particular horizontal & vertical
c              region (on the "u" grid) for momentum term balances
c              volu(0) represents the sum of all regions
c     rvolu  = 1/volu ( 0.0 if volu = 0.0)
c     areau  = horizontal ocean area corresponding to "volu"
c              areau(0) represents the sum of all regions
c     rareau = 1/areau ( 0.0 if areau = 0.0)
c     llvreg = level limits for defining vertical regions in term
c              balance calculations (not used in computing volume
c              weighted tracer averages)
c              (eg: llvreg(3,1) = 4... means that starting level for
c                  the third region in the vertical is 4. similarly,
c                  llvreg(3,2) = 6 means the ending level is 6 for that
c                  region. note regions should not overlap.)
c
      character*40 hregnm
      character*20 vregnm
c
      common /cregni/ mskvr(km), llvreg(numreg,2)
      common /cregnr/ volbk(nhreg,km), volbt(nhreg), volgk(km)
      common /cregnr/ volgt, areab(nhreg), areag, volt(0:numreg)
      common /cregnr/ volu(0:numreg), areat(0:numreg), areau(0:numreg)
      common /cregnr/ rvolt(0:numreg), rvolu(0:numreg), rareat(0:numreg)
      common /cregnr/ rareau(0:numreg)
      common /cregnc/ hregnm(nhreg), vregnm(nvreg)
c
