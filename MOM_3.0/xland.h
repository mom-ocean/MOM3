c====================== include file "xland.h" =========================
#ifdef xlandmix
c
c     this file holds information used to establish communication
c     between bodies of water separated by land.  communication consists
c     of mixing tracers in water columns found on each side of the land
c     barrier using a method that some liken to a kind of partial
c     horizontal convection.
c
c     the vertical extent of the mixing, the volume rate of mixing, and
c     the i,j location of the pairs of points to be mixed, are specified
c     in data statements in subroutine xland1st.
c
c     Note: many of the user specified values are very model dependent.
c
c     user specified variables:
c     ixland = the specified 2 i column locations of the water columns
c              to be mixed for each of the nxland pairs of points
c     jxland = the specified 2 j-row locations of the nxland pairs of
c              points
c     kxland = the specified vertical extent of the mixing by setting
c              the top & bottom model k-levels to be mixed
c     vxland = the specified rate at which of the water columns are to
c              be mixed [cubic cm per second] equivalent to
c              the flow to the east = the flow to the west
c
c     variables calculated from user input:
c     fxland = this logical flag is true for j-rows containing a
c              crossland mixing point, and false for those which do not.
c     bxland = this calculated variable is the fraction of the model
c              grid boxes specified that must be mixed per second
c              in order to meet the requested vxland mixing rate
c
c        Note: for ixland, jxland, and kxland pairs the values of the
c              (n,1) element should be < the corresponding (n,2) value.
c              also, if nxland > 1 performance is enhanced if the
c              jxland(n,1) values are in ascending order.
c
c     author:   k.w.dixon          e-mail  kd@gfdl.gov
c=======================================================================
c
c
c
# ifdef r15_res
      parameter (nxland = 3)
# else
      parameter (nxland = ?)
# endif
c
      logical fxland
      common /xlandl/ fxland(jmt)
c
      common /xlandi/ ixland(nxland,2), jxland(nxland,2)
      common /xlandi/ kxland(nxland,2)
c
      common /xlandr/ vxland(nxland), bxland(nxland,2)

# ifdef xlandmix_eta
      common /xlandr/ source_eta(imt,jmt),bxland_eta(nxland,2)
# endif
c
#endif
