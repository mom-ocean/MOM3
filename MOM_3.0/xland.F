      subroutine xland1st
#ifdef xlandmix
c
c=======================================================================
c
c               Initialization for crossland mixing scheme
c
c     This subroutine handles the initial set up for crossland mixing
c     of tracers.  It neeeds to be called just once per model run from
c     subroutine setocn.
c
c     I,J,K locations of crossland mixing points are set in data
c     statements.
c     Time invariant crossland mixing rates (in units of cm**3/sec) can
c     be set in data statements.
c
c     Checks are performed to ensure that the selected crossland mixing
c     grid locations are valid as the model is currently configured.
c     Logical flags are set indicating which model jrows contain points
c     affected by crossland mixing.
c     A short summary of the locations of crossland mixing points is
c     written out.
c
c     inputs:
c      ixland and jxland = in data statments the user needs to specify
c             nxland pairs of i,j grid locations for the crossland
c             mixing scheme.  (nxland is a parameter set in xland.h)
c      kxland = in data statements the user needs to specify the
c             uppermost (kxland(n,1)) and deepest (kxland(n,2)) model
c             levels over which crossland mixing will be done for each
c             pair of crossland mixing points.
c      vxland = one can set time invariant rates of crossland mixing in
c             data statements. (if vxland changes during the model run,
c             it will need to be reset elsewhere in a manner left to
c             the user)
c
c     output:
c      fxland = logical flags indicating which model j-rows contain
c             points affected by crossland mixing of tracers.
c
c     author:   K.W. Dixon    kd@gfdl.gov
c=======================================================================
c
      use arrays_2d
# include "param.h"
# include "coord.h"
# include "xland.h"
c
c     ===========
c     USER INPUT:
c     ===========
c     set locations for crossland mixing scheme in data statements.
c     values are extremely model dependent (see "xland.h")
c
# ifdef r15_res
c     set specifications for crossland mixing for ocean model matched
c     with the r15 spectral model
c
c     ------------------------------------------------------------------
c
c     set specs for mixing across the unresolved Straits of Gibraltar
c
      data (ixland(1,ix),ix=1,2) /  2,  95 /
c
      data (jxland(1,jx),jx=1,2) / 30,  30 /
c
#  if defined km12_res
c     set Gibraltar crossland mixing k-range for 12-level model
      data (kxland(1,kx),kx=1,2) / 1, 5/
#  endif
#  if defined km18_res
c     set Gibraltar crossland mixing k-range for 18-level model
      data (kxland(1,kx),kx=1,2) / 1, 8/
#  endif
c
c     set Gibraltar inflow and outflow to be 1.75e12 cm**3 per sec
c     (approx. value taken from Knauss, Physical Oceanography, pg. 69)
      data vxland(1) / 1.75e12 /
c
c     ------------------------------------------------------------------
c
c     set specs for crossland mixing across the Bering Strait
c     here it accomplished using two pairs of points
c
      data (ixland(2,ix),ix=1,2) / 52,  52 /
      data (ixland(3,ix),ix=1,2) / 53,  53 /
c
      data (jxland(2,jx),jx=1,2) / 35,  37 /
      data (jxland(3,jx),jx=1,2) / 35,  37 /
c
#  if defined km12_res
c     set Bering Strait crossland mixing k-range for 12-level model
      data (kxland(2,kx),kx=1,2) / 1, 1/
      data (kxland(3,kx),kx=1,2) / 1, 1/
#  endif
#  if defined km18_res
c     set Bering strait crossland mixing k-range for 18-level model
      data (kxland(2,kx),kx=1,2) / 1, 2/
      data (kxland(3,kx),kx=1,2) / 1, 2/
#  endif
c
c     set bering strait inflow and outflow totaling 1 Sverdrup (0.5 for
c     each of the two pairs of points)
      data vxland(2) / 0.400e12 /
      data vxland(3) / 0.400e12 /
c
c     ------------------------------------------------------------------
# else
c
c     specifications not yet set for resolutions other than r15
c
      data vxland(1) = ???
c
# endif
c
c
c=======================================================================
c
c     for each pair of i,j points to be mixed in this manner
c     the following is done:
c
c=======================================================================
c
c     initialize logical flags used to indicate if any crossland mixing
c     points are located in a given j-row
c
      do j=1,jmt
        fxland(j) = .FALSE.
      end do
c
      do nxl=1,nxland
c
c       check for invalid crossland mixing grid locations
c
        do lx = 1,2
          if (ixland(nxl,lx) .lt. 2 .or. ixland(nxl,lx) .gt. imtm1) then
            write(stdout,991) nxl, lx, ixland(nxl,lx)
            call abort()
          endif
          if (jxland(nxl,lx) .lt. 1 .or. jxland(nxl,lx) .gt. jmt ) then
            write(stdout,992) nxl, lx, jxland(nxl,lx)
            call abort()
          else
c           set flag to show that this j-row has a crossland mixing pt
            fxland(jxland(nxl,lx)) = .TRUE.
          endif
          if (kxland(nxl,lx) .lt. 1 .or. kxland(nxl,lx) .gt. km ) then
            write(stdout,993) nxl, lx, kxland(nxl,lx)
            call abort()
          endif
        end do
c
        ktop = kxland(nxl,1)
        kbom = kxland(nxl,2)
c
c       make sure that k-levels are ordered properly and that j-rows
c       of a given crossland mixing pair fit in the memory slab window
c
        if ( ktop .gt. kbom ) then
          write (stdout,994) nxl, ktop, nxl, kbom
          call abort()
        endif
c
        if (abs((jxland(nxl,2) - jxland(nxl,1))) .gt. ((jmw-1)/2) ) then
          write (stdout,995) nxl,nxl,jxland(nxl,2),jxland(nxl,1),jmw
          call abort()
        endif
c
c       check for attempts to mix land rather than sea
c
        do lx = 1,2
          if ( kbom .gt. kmt(ixland(nxl,lx),jxland(nxl,lx))
     &         ) then
            write (stdout,
     &      '(/" Error => problem with crossland tracer mixing:",/,
     &      " improper k-levels requested"/" kxland(",i2,",2) =",i4,
     &      ", and improperly exceeds kmt(",i4,",",i4,")"
     &           ," = ",i4,/)') nxl, kbom 
     &           , ixland(nxl,lx),jxland(nxl,lx)
     &           , kmt(ixland(nxl,lx),jxland(nxl,lx))
          call abort()
          endif
        end do
c
c       write out summary info for this pair of crossland mixing points
c
        if ( ktop .eq. 1 ) then
         ztop = 0.0
        else
         ztop = zw(ktop-1)
        endif
        write(stdout,191) nxl, ixland(nxl,1),     jxland(nxl,1),
     &                      xt(ixland(nxl,1)), yt(jxland(nxl,1)),
     &                         ixland(nxl,2),     jxland(nxl,2),
     &                      xt(ixland(nxl,2)), yt(jxland(nxl,2)),
     &                      ktop, kbom, ztop, zw(kbom)
c
      end do
c
      return
191   format(/' ===== from xland1st ====='/
     &' for crossland sea communication pair number',i2,/
     &' mix  I,J gridpt (',i4,',',i4,') [long=',f8.3,' lat=',f8.3,']',/
     &' with I,J gridpt (',i4,',',i4,') [long=',f8.3,' lat=',f8.3,']',/
     &' from level',i3,' to',i3,' [depths of ',f10.3,' to ',f10.3,'cm]')
991   format(/' Error => problem with crossland tracer mixing:',/,
     &' out of bounds i grid location requested'/
     &'      ixland(',i4,',',i1,') was set incorrectly set to = ',i8,/)
992   format(/' Error => problem with crossland tracer mixing:',/,
     &' out of bounds j-row grid location requested'/
     &'      jxland(',i4,',',i1,') was set incorrectly set to = ',i8,/)
993   format(/' Error => problem with crossland tracer mixing:',/,
     &' out of bounds k-level grid location requested'/
     &'      kxland(',i4,',',i1,') was set incorrectly set to = ',i8,/)
994   format(/' Error => problem with crossland tracer mixing:',/,
     &' improper k-values requested'/' kxland(',i3,',1)=',i5,
     &' and is not less than kxland(',i3,',2)=',i5,/)
995   format(/' Error => problem with crossland tracer mixing:',/,
     &'  memory window not wide enough for j-row pair '/
     &'  abs(jxland(',i3,',2) -  jxland(',i3,',1)) = (',i4,' - ',i4,
     &') and should not be greater than (jmw-1)/2 = (',i4,'-1)/2')
      end



      subroutine xlandvchk (nxl)
c
c=======================================================================
c
c     Mixing rate conversion and check for crossland mixing scheme
c
c     In this routine the crossland mixing rate vxland(nxl) for
c     the pair of points associated with index number nxl is
c     converted into the fraction of the model grid boxes to be mixed
c     per second, and checked.  These checks ensure that the rate of
c     crossland mixing requested is valid in that it can be realized
c     given the timestep length and column volumes involved.
c
c     If a time invariant crossland mixing rate is being used, this
c     routine can be called only once each model run for each pair of
c     crossland mixing points (it can be called from setocn).
c     But, if the crossland mixing rates vary in time, then this routine
c     needs to be called each time the vxland(nxl) mixing rate value
c     changes.
c
c     input:
c      nxl  = crossland mixing pair index number (values from 1 to
c             nxland are valid)
c     output:
c      bxland = one over mixing time scale for cross land mixing (1/sec).
c               calculated as: 
c               rate of mixing (cm^3/sec)/ total volume of mixed boxes (cm^3).
c
c     author:  K.W. Dixon    kd@gfdl.gov
c=======================================================================
c
      use arrays_2d
# include "param.h"
# include "accel.h"
# include "coord.h"
# include "xland.h"
# include "grdvar.h"
# include "scalar.h"
c
c
c     local variables:
c
c     depth  = depth thickness over which crossland mixing is done [cm]
c     colvol = column volume over which crossland mixing is done [cm**3]
c
      dimension colvol(2)
c
c=======================================================================
c
      ktop = kxland(nxl,1)
      kbom = kxland(nxl,2)
c
c     Calculate depth range over which crossland mixing will be done.
c     (This depth range is only approximate if bottom partial cells or
c     nonlinear free surface are used.)
c
      depth = c0
      do k=ktop,kbom
        depth = depth + dzt(k)
      end do
c
c     for each of the two points to be mixed, convert from volume to
c     be mixed per second vxland(nxl) to model box fraction to be
c     mixed per second bxland(nxl,l)      (see xland.h)
c
      do lx=1,2
        colvol(lx) = cst(jxland(nxl,lx))*dxt(ixland(nxl,lx))*
     &               dyt(jxland(nxl,lx))*depth
        bxland(nxl,lx) = vxland(nxl) / colvol(lx)
#ifdef xlandmix_eta
c       For eta-mixing, the volume transfer due to the eta mixing is set
c       to (1/eta_mixing_time_scale)*(eta_1-eta_2)*(area_1+area_2)/2, where:
c       eta_1, eta_2  = eta at the two grid points being mixed
c       area_1,area_2 = area element at the two grid points being mixed
c       The volume flux needs to be normalized by the area element
c       at each point, in order to be used for the source term.

c       First, set the time scale for eta mixing as that required for 
c       the external waves to travel the distance between the two points
c       being mixed:
        eta_mixing_time_scale= 10.0*86400.0
        if (lx.eq.1) then
          write(*,'(" using eta_mixing_time_scale=",g14.6,"(hours)")')
     &         eta_mixing_time_scale/3600.
        end if
c       area of grid box at two mixing points:
        area_1=cst(jxland(nxl,1))*dxt(ixland(nxl,1))*dyt(jxland(nxl,1))
        area_2=cst(jxland(nxl,2))*dxt(ixland(nxl,2))*dyt(jxland(nxl,2))
c       Mixing factor used later to calculate source_eta:
        bxland_eta(nxl,lx)=  
     &       (1.0/eta_mixing_time_scale)
     &       *0.5*(area_1+area_2)
     &       /(cst(jxland(nxl,lx))*dxt(ixland(nxl,lx))*
     &               dyt(jxland(nxl,lx)))
#endif
c
c       check for attempts to mix at a rate so fast that it can not
c       be achieved in a single, leapfrog timestep
c       (for the two boxes to be thoroughly mixed, no more than
c       one-half of the volume of a given grid box can be
c       transported in effect into the other gridbox)
c
        do k=ktop,kbom
          factor = bxland(nxl,lx)*c2*dtts*dtxcel(k)
          if (factor .gt. p5) then
            tslim = 0.25/bxland(nxl,lx)
            write (stdout,997) nxl, factor, ixland(nxl,lx),
     &            jxland(nxl,lx), k ,bxland(nxl,lx), tslim,
     &            vxland(nxl), ktop, kbom, colvol(lx)
            call abort()
          endif
        end do
      end do
c
      write (stdout,197) nxl, ixland(nxl,1), jxland(nxl,1),
     &      ixland(nxl,2), jxland(nxl,2), ktop, kbom, depth,
     &      (colvol(ll),ll=1,2), vxland(nxl),
     &      (bxland(nxl,ll),ll=1,2)
c
      return
c
197   format(/' ===== from xlandvchk ====='/
     &' for crossland sea communication pair number',i2,/
     &' mix I,J gridpoints (',i4,',',i4,') and (',i4,',',i4,')'/
     &' from level',i3,' to',i3,' (a depth range of ',e12.6,' cm)'/
     &' column volumes =',e12.6,' and ',e12.6,' cm**3'/
     &' simulated flow in = flow out = ',e12.6,' cm**3/sec,',/,
     &' so mix ',e12.6,' fraction of 1st column with ',e12.6,
     &' of 2nd column per sec'/)
997   format(/' Error => problem with crossland tracer mixing',
     &' pair #',i3,/,' Error => attempting to mix at too fast a rate:'/
     &' asking for more than half (',g12.6,') the volume of a grid box',
     &' at i, j, k: ',3i5,/,' to be mixed into its neighbor on',
     &' leapfrog timesteps'/
     &' As specified now, fraction of column mixed per sec = ',e12.6,/,
     &' Potential Fixes: shorten dtts*dtxcel(k) to be less than',e13.6,/
     &t18,' reduce mixing rate vxland (now ',e12.6,' cm**3/sec)'/
     &t18,' or mix over a larger depth range (now from k=',i3,' to ',i3
     &,' => water column volume = ',e12.6,' cm**3)'/)
c
      end



      subroutine xland (joff, js, je, is, ie, n)
c
c=======================================================================
c
c
c                       crossland mixing scheme
c
c     allows the user to establish communication between bodies of water
c     separated by land.  communication consists of mixing tracers in
c     water columns found on each side of a land barrier using a method
c     that some users liken to a kind of partial horizontal convection.
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       n    = tracer number
c
c     output:
c       source = source term including effects of crossland mixing
c
c     author:   k.w.dixon            e-mail  kd@gfdl.gov
c=======================================================================
c
      use memory_window
      use arrays_2d
# include "param.h"
# include "mw.h"
# include "xland.h"
# include "scalar.h"
# include "coord.h"
# include "grdvar.h"
c
# ifdef timing
      call tic ('tracer', 'crossland mixing')
# endif
c
      do j=js,je
        jrow = j + joff
        if ( fxland(jrow) ) then
          do nxl = 1,nxland
            do lhere = 1,2
              lthere = mod(lhere,2) + 1
              if (jxland(nxl,lhere) .eq. jrow) then
                ihere  = ixland(nxl,lhere)
                jhere  = jxland(nxl,lhere) - joff
                ithere = ixland(nxl,lthere)
                jthere = jxland(nxl,lthere) - joff
c               Setup a few variables that are needed for inland sea
c               mixing with nonlinear explicit free surface and/or with
c               partial_cell:
                jhere1  = jxland(nxl,lhere)
                jthere1 = jxland(nxl,lthere)
# if defined explicit_free_surface  && defined xlandmix_eta
c               Calculate source for inland sea mixing of free surface height:

                eta_here  = etat(ihere,jhere1,1  )
                eta_there = etat(ithere,jthere1,1)

c               Next, calculate the source to be added to the eta
c               equation: The eta-mixing via
                k=1
c               update 
                source_eta(ihere,jhere1)=source_eta(ihere,jhere1)+
     &               bxland_eta(nxl,lhere)
     &               *(eta_there-eta_here)
     &               *tmask(ihere,k,jhere)*tmask(ithere,k,jthere)

#  if defined explicit_free_surface_linear
c               Note that no source term is added to the tracer
c               equations due to eta mixing, in the case of a **linear**
c               free surface.  That's because we do not add source_eta
c               to the explicit fresh water flux appearing in the tracer
c               equations, no such cancellation is needed.
#  else
c               mixing free surface height should affect the tracer
c               balance, and a corresponding source term should be added
c               to the tracer equations to express the fact that actual
c               sea water is being mixed when we mix eta, and not just
c               fresh water:
#  endif
                d_eta=(eta_there-eta_here)
                dvt_here=cst(jhere1)*dxt(ihere)*dyt(jhere1)
     &               *dht(ihere,k,jhere)
                source(ihere,k,jhere) = source(ihere,k,jhere) + (
     &               0.5*(d_eta+abs(d_eta))*t(ithere,k,jthere,n,taum1) +
     &               0.5*(d_eta-abs(d_eta))*t(ihere ,k,jhere ,n,taum1) )
     &               *tmask(ihere,k,jhere)*tmask(ithere,k,jthere)
     &               /dvt_here
# endif
                do k = kxland(nxl,1),kxland(nxl,2)
# if defined partial_cell || ( defined explicit_free_surface && ! defined explicit_free_surface_linear )
                  dvt_here  = cst(jhere1) *dxt(ihere) *dyt(jhere1)
     &                 *dht(ihere,k,jhere)
                  dvt_there = cst(jthere1)*dxt(ithere)*dyt(jthere1)
     &                 *dht(ithere,k,jthere)
# endif
                  source(ihere,k,jhere) = source(ihere,k,jhere) +
     &             bxland(nxl,lhere) *
     &             (t(ithere,k,jthere,n,taum1)-t(ihere,k,jhere,n,taum1))
     &             *tmask(ihere,k,jhere)*tmask(ithere,k,jthere)
# if defined partial_cell || ( defined explicit_free_surface && ! defined explicit_free_surface_linear )
     &               *(dvt_here+dvt_there)/(2.0*dvt_here)
# endif
                enddo
              endif
            enddo
          enddo
        endif
      enddo
c
# ifdef timing
      call toc ('tracer', 'crossland mixing')
# endif
c
#endif
      return
      end
