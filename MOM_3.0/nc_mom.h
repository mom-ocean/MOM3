c================= include file "nc_mom.h"==============================
c     
c     Interface variables needed to connect MOM to Sheldon`s NETCDF
c     
      parameter(maxncdim = 100)
      parameter(maxncnam = 128)
      parameter (maxlatt=256)
      integer     lendims(ndimsout)
      dimension     dimvals(mxdimout,ndimsout)
      logical       lclobber, lgspval
      character*(maxncnam) cdimnam(ndimsout),   cvarnam(nvarsout)
      character*(maxlatt)  cdimlnam(ndimsout),  cvarlnam(nvarsout)
      character*(maxlatt)  cdimunits(ndimsout), cvarunits(nvarsout)
      character*1          cart_axis(ndimsout)
      integer		   ipositive(ndimsout), nbits(nvarsout)
      character*(maxlatt)  gtitle
      real       gspval
      dimension  validr(2,nvarsout)
      integer  nvdims(nvarsout)
      integer  idimindx(ndimsout,nvarsout)
      integer  maxnitems(nvarsout)
      integer  istart(4),icount(4)
      real*8 time(2000)
c
      character*(maxlatt) time_since
      real*8 timrec
