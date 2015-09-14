#ifndef test_poisson
      subroutine barotropic (f, itt, js, je)
c
c=======================================================================
c
c        S O L V E   T H E   B A R O T R O P I C   E Q U A T I O N  
c                  (refer to the manual for details)
c=======================================================================
c
      use arrays_2d
      use memory_window
      use mom3_mpp_mod
# include "param.h"
# include "emode.h"
# include "grdvar.h"
# include "iounit.h"
# ifdef obc
#  include "cobc.h"
# endif
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# include "taskrows.h"
      character*8 bc_symm
# ifdef rot_grid
      dimension f(imt,js:je)
# else
      dimension f(js:je)
# endif
      save nminscan
c
# ifdef timing
      call tic ('ocean', 'barotropic')
# endif
c
      call border (zu(1,js,1), 'u even', js, je)
      call border (zu(1,js,2), 'u odd',  js, je)
c
#if defined parallel_1d
      call mpp_update_domains( zu, mom_domain )
# endif
c
# if defined debug_barotropic
      write (stdout,*) 'in barotropic ck2d: js=',js,' je=',je,' pn=',pn
      write (stdout,*) ' ck2d zu1=',ck2d(zu(1,js,1), js, je)
      write (stdout,*) ' ck2d zu2=',ck2d(zu(1,js,2), js, je)
# endif
c
# if defined rigid_lid_surface_pressure || defined free_surf_on
c
c-----------------------------------------------------------------------
c     solve for the "tau+1" barotropic velocities
c-----------------------------------------------------------------------
c
#  if defined rigid_lid_surface_pressure || defined implicit_free_surface
      call impl_freesurf (js, je)
c
c     test accuracy of solving for change in surface pressure
c
      if (.not.converged) then
        write (stdout,'(a,i5,3(a,1pe10.3))')
     &  ' WARNING: SOLVER DID NOT CONVERGE in ',mscan
     &, ' scans. max(ps)='
     &, absmax(ps(1,js,2), js, je), ' max(ptd)=',absmax(ptd, js, je)
     &, ' estimated max(err)=', esterr
      endif
#  endif
c
#  ifdef explicit_free_surface
      if (.not.euler2) then
	call expl_freesurf (js, je)
      endif
#  endif
# endif
c
# ifdef stream_function
c
c-----------------------------------------------------------------------
c     construct the forcing for the stream function equation
c-----------------------------------------------------------------------
c
      call sfforc (zu, dxu, dyu, csu, ztd) 

c
#  if defined fourfil || defined firfil
c
c     filter forcing at high latitudes
c
      call filz (ztd, js, je)
#  endif
c
c-----------------------------------------------------------------------
c     read in solution from 2 previous timesteps for the
c     purpose of computing an initial guess for the present solution.
c-----------------------------------------------------------------------
c
      luptdb = nkflds - mod(itt,2)
      luptd  = nkflds - 1 + mod(itt,2)
      call oget (kflds, nwds, luptdb, res)
      call oget (kflds, nwds, luptd, ptd)
c
      fxa=c1
      if (.not. leapfrog) fxa=p5
      do jrow=1,jmt
        do i=1,imt
          ptd(i,jrow) = fxa*(c2*ptd(i,jrow)-res(i,jrow))
        enddo
      enddo
#  ifdef cyclic
      do jrow=2,jmtm1
        ptd(1,jrow)   = ptd(imtm1,jrow)
        ptd(imt,jrow) = ptd(2,jrow)
      enddo
#  endif
c
c-----------------------------------------------------------------------
c     choose 5 or 9 point numerics
c-----------------------------------------------------------------------
c
#  ifdef sf_5_point
c
c     initialize coefficients using 5 point numerics
c
      call sfc5pt (acor, f, c2dtsf, dxu, dyu, csu, hr, cf)
      npt = 5
#  endif
#  ifdef sf_9_point
c
c     initialize coefficients using 9 point numerics
c
      call sfc9pt (acor, f, c2dtsf, dxu, dyu, csu, hr, cf)
      npt = 9
#  endif
c
c-----------------------------------------------------------------------
c     solve for the "tau+1" stream function change
c-----------------------------------------------------------------------
c
      variable   = 'dpsi'
      bc_symm    = 't odd'
c
#  ifdef timing
      call tic ('congrad', 'stream function (barotropic)')
#  endif
#  ifdef orlanski
#   include "barotropic1_obc.inc"
#  endif
      call congr (npt, variable, bc_symm, ptd, ptd, ztd, res 
     &,           cf, js, je
     &,           mxscan, mscan, tolrsf
     &,           imask, iperm, jperm, iofs, nisle, nippts
     &,           converged, esterr)
#  ifdef timing
      call toc ('congrad', 'stream function (barotropic)')
#  endif
c
c     correct for drifting dpsi on land mass "imain"
c     replace "0" by "epsln" on coast of "imain" so ferret will not
c     treat coastal points as land points.
c
      if (imain .gt. 0) then
        i1 = iperm(iofs(imain)+1)
	j1 = jperm(iofs(imain)+1)
        dpsi1 = ptd(i1,j1)
        call con_adjust (ptd, dpsi1, map)
        do n=1,nippts(imain)
          i = iperm(iofs(imain)+n)
          j = jperm(iofs(imain)+n)
          ptd(i,j) = epsln
        enddo
        do jrow=1,jmt
#  ifdef cyclic
          ptd(1,jrow)   = ptd(imt-1,jrow)
	  ptd(imt,jrow) = ptd(2,jrow)
#  else
#   ifdef obc_west
          ptd(1,jrow)   = ptd(2,jrow)
#   else	
          ptd(1,jrow)   = 0.0
#   endif
#   ifdef obc_east
	  ptd(imt,jrow) = ptd(imt-1,jrow)
#   else	
	  ptd(imt,jrow) = 0.0
#   endif
#  endif
        enddo
      endif
c
c-----------------------------------------------------------------------
c     test accuracy of solving for change in stream function
c-----------------------------------------------------------------------
c
      if (.not.converged) then
        write (stdout,'(a,i5,a,i10,3(a,1pe10.3))')
     &  ' WARNING: SOLVER DID NOT CONVERGE in ',mscan
     &, ' scans on itt=',itt,'. max(psi)='
     &, absmax(psi(1,js,2), js, je), ' max(dpsi)=',absmax(ptd, js, je)
     &, ' estimated max(err)=', esterr
      endif
c
      if (first) then
        nminscan = 0
      endif
      if (mscan .lt. 0.9*max(imt,jmt)) then
        nminscan = nminscan + 1
        if (nminscan .le. 50) then
          write (stdout,'(1x,a,i5,a,i10,a)') 
     &    '=>WARNING: Number of scans (',mscan
     &,   ') may be too low on itt=',itt
     &,   '. Consult manual (under Elliptic Equation Solver)'
	endif
      endif
      if (eorun .and. nminscan .ne. 0) then
          write (stdout,'(/1x,a,i5,a/)') 
     &    '=>WARNING: Small number of scans detected ',nminscan
     &,   ' times. Consult manual (under Elliptic Equation Solver)'
      endif
c
c-----------------------------------------------------------------------
c     update the stream function based upon the solution
c-----------------------------------------------------------------------
c
      if (euler2) then
        do jrow=1,jmt
          do i=1,imt
            psi(i,jrow,1) = psi(i,jrow,2) + ptd(i,jrow)
          enddo
	enddo
      else
        do jrow=1,jmt
          do i=1,imt
            res(i,jrow)    = psi(i,jrow,2) + ptd(i,jrow)
#  if defined robert_time_filter
c
c           time filter streamfunction
c
	    psi(i,jrow,1) = psi(i,jrow,1)
     &      + smooth*(0.5*(res(i,jrow) + psi(i,jrow,2))
     &      - psi(i,jrow,1))
#  endif
            psi(i,jrow,2)  = psi(i,jrow,1)
            psi(i,jrow,1)  = res(i,jrow)
          enddo
	enddo
      endif
c
#  ifdef obc
#   include "barotropic2_obc.inc"
#  endif
c-----------------------------------------------------------------------
c     save ptd to compute 1st guess for next timestep
c     (..note.. on 1st pass of euler backward timestep, bypass this
c            save, since it will be done on the 2nd pass)
c     (..note.. on a mixing timestep, alter ptd to be consistent with
c            normal, leap-frog stepping)
c-----------------------------------------------------------------------
c
      if (.not. euler1) then
c
        if (.not. leapfrog) then
          do jrow=1,jmt
            do i=1,imt
              ptd(i,jrow)=c2*ptd(i,jrow)
            enddo
	  enddo
        endif
c
        call oput (kflds, nwds, luptdb, ptd)
c
      endif
# endif
c
# ifdef timing
      call toc ('ocean', 'barotropic')
# endif
      return
      end
#endif



#if !defined explicit_free_surface
      subroutine sfforc (zu, dxu, dyu, csu, forc)
c
c=======================================================================
c
c           S T R E A M   F U N C T I O N   F O R C I N G
c
c  
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
c
# include "size.h"
      dimension zu(imt,jmt,2), dxu(imt), dyu(jmt), csu(jmt) 
      dimension forc(imt,jmt)
      dimension ustuff(imt,jmt), vstuff(imt,jmt)
      dimension cddxu(0:1,0:1), cddyu(0:1,0:1)
      dimension cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)
      parameter (p5=0.5, c0=0.0)
c
# ifdef timing
      call tic ('sfforc', 'sfforc')
# endif
c
c-----------------------------------------------------------------------
c     initialize the forcing
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        do i=1,imt
          forc(i,jrow) = c0
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     construct coefficients for partial differences. a partial
c     difference in "x" is defined as an "x" difference of a quantity 
c     which is averaged in "y". (and symmetrically for "y" differences).
c     Note that this is an x difference and NOT an x derivative.
c     partial differences of quantities on the "t" grid are defined on
c     the "u" grid and visa versa.
c     therefore partial differences at: 
c     u/v points (i,j), involve nearby t/s points with subscripts: 
c        (i  ,j+1)    (i+1,j+1)
c        (i  ,j  )    (i+1,j  )
c     t/s points (i,j), involve nearby u/v points with subscripts:
c        (i-1,j  )    (i  ,j  )
c        (i-1,j-1)    (i  ,j-1)
c     thus if qu(i,j) is defined on u/v points, its partial
c     difference ddxqt = ddxt(qu) is defined on t/s points and has the
c     value 
c     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
c                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
c-----------------------------------------------------------------------
c
      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5
c
      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5
c
      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5
c
      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5
c
c-----------------------------------------------------------------------
c     multiply the u eqn by dx*cos, the v eqn by dy, then subtract their
c     partial differences to eliminate the unknown surface pressure from
c     the resulting equation
c-----------------------------------------------------------------------
c
      do jrow=1,jmt-1
        do i=1,imt-1
          ustuff(i,jrow) = zu(i,jrow,1)*dxu(i)*csu(jrow)
          vstuff(i,jrow) = zu(i,jrow,2)*dyu(jrow)
        enddo
      enddo
c
      do jrow=2,jmt-1
        do i1=-1,0
          do j1=-1,0
            do i=2,imt-1
              forc(i,jrow) = forc(i,jrow) 
     &                     - cddyt(i1,j1)*ustuff(i+i1,jrow+j1)
     &                     + cddxt(i1,j1)*vstuff(i+i1,jrow+j1)
            enddo
          enddo
        enddo
      enddo
c
# ifdef timing
      call toc ('sfforc', 'sfforc')
# endif
c
      return
      end




      subroutine sfc5pt (acor, f, c2dtsf, dxu, dyu, csu, hr, coef)
c
c=======================================================================
c
c     5  P T    C O E F F I C I E N T   I N I T I A L I A Z A T I O N
c
c     coefficient initialization for 5 point elliptic solvers
c
c     inputs:
c
c     acor   = implicit coriolis factor (0.0 => 1.0)
c     f      = 2*omega*sin(phi(j))
c     c2dtsf = twice the time step (seconds)
c     dxu    = width of "u" grid cell (cm)
c     dyu    = height of "u" grid cell (cm)
c     csu    = cosine of "u" grid cell
c     hr     = 1/depth at "u" cells (cm)
c
c     outputs:
c
c     coeff   = 3 x 3 array of coefficients at each (i,j) point
c
c     author: c.h. goldberg         e-mail => chg@gfdl.gov
c=======================================================================
c
# include "size.h"
      parameter (p5=0.5, c0=0.0)
      dimension csu(jmt), dxu(imt), dyu(jmt), hr(imt,jmt)
# ifdef rot_grid
      dimension f(imt,jmt)
# else
      dimension f(jmt)
# endif
      dimension coef(imt,jmt,-1:1,-1:1)
      dimension ustuff(imt,jmt), vstuff(imt,jmt)
      dimension cddxu(0:1,0:1), cddyu(0:1,0:1)
      dimension cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)
c
# ifdef timing
      call tic ('sfc5pt', 'sfc5pt')
# endif
c
c-----------------------------------------------------------------------
c     initialize the coefficients
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do jj=-1,1
          do ii=-1,1
            do i=1,imt
              coef(i,j,ii,jj) = c0
            enddo
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     construct coefficients for partial differences. a partial
c     difference in "x" is defined as an "x" difference of a quantity 
c     which is averaged in "y". (and symmetrically for "y" differences).
c     Note that this is an x difference and NOT an x derivative.
c     partial differences of quantities on the "t" grid are defined on
c     the "u" grid and visa versa.
c     therefore partial differences at: 
c     u/v points (i,j), involve nearby t/s points with subscripts: 
c        (i  ,j+1)    (i+1,j+1)
c        (i  ,j  )    (i+1,j  )
c     t/s points (i,j), involve nearby u/v points with subscripts:
c        (i-1,j  )    (i  ,j  )
c        (i-1,j-1)    (i  ,j-1)
c     thus if qu(i,j) is defined on u/v points, its partial
c     difference ddxqt = ddxt(qu) is defined on t/s points and has the
c     value 
c     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
c                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
c-----------------------------------------------------------------------
c
      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5
c
      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5
c
      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5
c
      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5
c
c-----------------------------------------------------------------------
c     compute coefficients for all points
c-----------------------------------------------------------------------
c
      do j=1,jmt-1
        do i=1,imt-1
          ustuff(i,j) = dxu(i)*csu(j)*hr(i,j) / (c2dtsf*dyu(j))
          vstuff(i,j) = dyu(j)*hr(i,j) / (c2dtsf*dxu(i)*csu(j))
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     calculate 5 point coefficients
c
c     note that ne [and nw] coefficient adds to n coefficient in 
c     ustuff term, but ne [and se] coefficient adds to e coefficient in
c     vstuff term for the 5 point operator.
c-----------------------------------------------------------------------
c
      do j=2,jmt-1
        do i1=0,1
          do j1=0,1
            do i2=-1,0
              do j2=-1,0
                do i=2,imt-1
                  coef(i,j,0,j1+j2) = coef(i,j,0,j1+j2) +
     &               cddyu(i1,j1)*cddyt(i2,j2)*ustuff(i+i2,j+j2)
                  coef(i,j,i1+i2,0) = coef(i,j,i1+i2,0) +
     &               cddxu(i1,j1)*cddxt(i2,j2)*vstuff(i+i2,j+j2)
                enddo
	      enddo
            enddo
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     augment coefficients for implicit treatment of coriolis term
c     all coefficients are calculated, but corner ones are zero.
c-----------------------------------------------------------------------
c
      if (acor .ne. 0.0) then
        do j=1,jmt-1
          do i=1,imt-1
# ifdef rot_grid
            ustuff(i,j) = acor*hr(i,j)*(-f(i,j))
            vstuff(i,j) = acor*hr(i,j)*( f(i,j))
# else
            ustuff(i,j) = acor*hr(i,j)*(-f(j))
            vstuff(i,j) = acor*hr(i,j)*( f(j))
# endif
          enddo
        enddo

        do j=2,jmt-1
          do i1=0,1
            do j1=0,1
              do i2=-1,0
                do j2=-1,0
                  do i=2,imt-1
                    coef(i,j,i1+i2,j1+j2) = coef(i,j,i1+i2,j1+j2)
     &                 - cddxu(i1,j1)*cddyt(i2,j2)*ustuff(i+i2,j+j2)
                    coef(i,j,i1+i2,j1+j2) = coef(i,j,i1+i2,j1+j2)
     &                 - cddyu(i1,j1)*cddxt(i2,j2)*vstuff(i+i2,j+j2)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
c
# ifdef timing
      call toc ('sfc5pt', 'sfc5pt')
# endif
c
      return
      end



      subroutine sfc9pt  (acor, f, c2dtsf, dxu, dyu, csu, hr, coef)
c
c=======================================================================
c
c     9  P T    C O E F F I C I E N T   I N I T I A L I A Z A T I O N
c
c     coefficient initialization for 9 point elliptic solvers
c
c     inputs:
c
c     acor   = implicit coriolis factor (0.0 => 1.0)
c     f      = 2*omega*sin(phi(j))
c     c2dtsf = twice the time step (seconds)
c     dxu    = width of "u" grid cell (cm)
c     dyu    = height of "u" grid cell (cm)
c     csu    = cosine of "u" grid cell
c     hr     = 1/depth at "u" cells (cm)
c
c     outputs:
c
c     coeff   = 3 x 3 array of coefficients at each (i,j) point
c
c     author: c.h. goldberg         e-mail => chg@gfdl.gov
c=======================================================================
c
# include "size.h"
      parameter (c0=0.0, p5=0.5)
      dimension csu(jmt), dxu(imt), dyu(jmt), hr(imt,jmt)
# ifdef rot_grid
      dimension f(imt,jmt)
# else
      dimension f(jmt)
# endif
      dimension coef(imt,jmt,-1:1,-1:1)
      dimension ustuff(imt,jmt), vstuff(imt,jmt)
      dimension cddxu(0:1,0:1), cddyu(0:1,0:1)
      dimension cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)
c
# ifdef timing
      call tic ('sfc9pt', 'sfc9pt')
# endif
c
c-----------------------------------------------------------------------
c     initialize the work area
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do jj=-1,1
          do ii=-1,1
            do i=1,imt
              coef(i,j,ii,jj) = c0
            enddo
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     generate arrays of coefficients
c
c     construct coefficients for partial differences. a partial
c     difference in "x" is defined as an "x" difference of a quantity 
c     which is averaged in "y". (and symmetrically for "y" differences).
c     Note that this is an x difference and NOT an x derivative.
c     partial differences of quantities on the "t" grid are defined on
c     the "u" grid and visa versa.
c     therefore partial differences at: 
c     u/v points (i,j), involve nearby t/s points with subscripts: 
c        (i  ,j+1)    (i+1,j+1)
c        (i  ,j  )    (i+1,j  )
c     t/s points (i,j), involve nearby u/v points with subscripts:
c        (i-1,j  )    (i  ,j  )
c        (i-1,j-1)    (i  ,j-1)
c     thus if qu(i,j) is defined on u/v points, its partial
c     difference ddxqt = ddxt(qu) is defined on t/s points and has the
c     value 
c     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
c                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
c-----------------------------------------------------------------------
c
      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5
c
      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5
c
      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5
c
      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5
c
c-----------------------------------------------------------------------
c     compute coefficients for all points
c-----------------------------------------------------------------------
c
      do j=1,jmt-1
        do i=1,imt-1
          ustuff(i,j) = dxu(i)*csu(j)*hr(i,j) / (c2dtsf*dyu(j))
          vstuff(i,j) = dyu(j)*hr(i,j) / (c2dtsf*dxu(i)*csu(j))
        enddo
      enddo
c
c---------------------------------------------------------------------
c     calculate 9 point coefficients
c---------------------------------------------------------------------
c
      do j=2,jmt-1
        do i1=0,1
          do j1=0,1
            do i2=-1,0
              do j2=-1,0
                do i=2,imt-1
                  coef(i,j,i1+i2,j1+j2) = coef(i,j,i1+i2,j1+j2) +
     &               cddyu(i1,j1)*cddyt(i2,j2)*ustuff(i+i2,j+j2)  +
     &               cddxu(i1,j1)*cddxt(i2,j2)*vstuff(i+i2,j+j2)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     augment coefficients for implicit treatment of coriolis term
c     all coefficients are calculated, but corner ones are zero.
c-----------------------------------------------------------------------
c
      if (acor .ne. 0.0) then
        do j=1,jmt-1
          do i=1,imt-1
# ifdef rot_grid
            ustuff(i,j) = acor*hr(i,j)*(-f(i,j))
            vstuff(i,j) = acor*hr(i,j)*( f(i,j))
# else
            ustuff(i,j) = acor*hr(i,j)*(-f(j))
            vstuff(i,j) = acor*hr(i,j)*( f(j))
# endif
          enddo
        enddo
        do j=2,jmt-1
          do i1=0,1
            do j1=0,1
              do i2=-1,0
                do j2=-1,0
                  do i=2,imt-1
                    coef(i,j,i1+i2,j1+j2) = coef(i,j,i1+i2,j1+j2)
     &               - cddxu(i1,j1)*cddyt(i2,j2)*ustuff(i+i2,j+j2)
     &               - cddyu(i1,j1)*cddxt(i2,j2)*vstuff(i+i2,j+j2)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
c
# ifdef timing
      call toc ('sfc9pt', 'sfc9pt')
# endif
c
      return
      end


      subroutine spforc (zu, dxu, dyu, csu, h, forc, js, je)
c
c=======================================================================
c
c           S U R F A C E   P R E S S U R E   F O R C I N G
c
c  
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
c
# include "size.h"
# include "taskrows.h"
      dimension zu(imt,js:je,2)
      dimension dxu(imt), dyu(jmt), csu(jmt)
      dimension h(imt,js:je)
      dimension forc(imt,js:je)
      dimension ustuff(imt,js:je), vstuff(imt,js:je)
      dimension cddxu(0:1,0:1), cddyu(0:1,0:1)
      dimension cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)
      parameter (p5=0.5)
c
# ifdef timing
      call tic ('spforc', 'spforc')
# endif
c
c-----------------------------------------------------------------------
c     generate arrays of coefficients
c
c     construct coefficients for partial differences. a partial
c     difference in "x" is defined as an "x" difference of a quantity 
c     which is averaged in "y". (and symmetrically for "y" differences).
c     Note that this is an x difference and NOT an x derivative.
c     partial differences of quantities on the "t" grid are defined on
c     the "u" grid and visa versa.
c     therefore partial differences at: 
c     u/v points (i,j), involve nearby t/s points with subscripts: 
c        (i  ,j+1)    (i+1,j+1)
c        (i  ,j  )    (i+1,j  )
c     t/s points (i,j), involve nearby u/v points with subscripts:
c        (i-1,j  )    (i  ,j  )
c        (i-1,j-1)    (i  ,j-1)
c     thus if qu(i,j) is defined on u/v points, its partial
c     difference ddxqt = ddxt(qu) is defined on t/s points and has the
c     value 
c     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
c                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
c-----------------------------------------------------------------------
c
      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5
c
      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5
c
      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5
c
      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5
c
c     weight "zu" and "zv" by the cell area and take the divergence
c
      do j=js,jecomp
        do i=1,imt-1
	  ustuff(i,j) = h(i,j)*zu(i,j,1)*dyu(j)
	  vstuff(i,j) = h(i,j)*zu(i,j,2)*dxu(i)*csu(j)
	enddo
      enddo
c
      do j=js,je
        do i=1,imt
	  forc(i,j) = 0.0
	enddo
      enddo
c
      do j=jscomp,jecomp
        do i1=-1,0
	  do j1=-1,0
	    do i=2,imt-1
	      forc(i,j) = forc(i,j) + cddxt(i1,j1)*ustuff(i+i1,j+j1)
     &                              + cddyt(i1,j1)*vstuff(i+i1,j+j1)
	    enddo
	  enddo
	enddo
      enddo
c
# ifdef timing
      call toc ('spforc', 'spforc')
# endif
c
      return
      end



      subroutine spc9pt (dxu, dyu, csu, h, coef, js, je)
c
c=======================================================================
c
c     S U R F A C E   P R E S S U R E    C O E F F I C I E N T
c
c                I N I T I A L I A Z A T I O N
c
c     inputs:
c
c     dxu    = width of "u" grid cell (cm)
c     dyu    = height of "u" grid cell (cm)
c     csu    = cosine of "u" grid cell
c     h      = depth at "u,v" cells (cm)
c
c     outputs:
c
c
c     coeff   = 3 x 3 array of coefficients at each (i,j) point
c
c     author: c.h. goldberg         e-mail => chg@gfdl.gov
c=======================================================================
c
# include "size.h"
# include "taskrows.h"
      parameter (c0=0.0, p5=0.5)
      dimension csu(jmt), dxu(imt), dyu(jmt), h(imt,js:je)
      dimension coef(imt,js:je,-1:1,-1:1)
      dimension ustuff(imt,js:je), vstuff(imt,js:je)
      dimension cddxu(0:1,0:1), cddyu(0:1,0:1)
      dimension cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)
c
# ifdef timing
      call tic ('spc9pt', 'spc9pt')
# endif
c
c-----------------------------------------------------------------------
c     generate arrays of coefficients
c
c     construct coefficients for partial differences. a partial
c     difference in "x" is defined as an "x" difference of a quantity 
c     which is averaged in "y". (and symmetrically for "y" differences).
c     Note that this is an x difference and NOT an x derivative.
c     partial differences of quantities on the "t" grid are defined on
c     the "u" grid and visa versa.
c     therefore partial differences at: 
c     u/v points (i,j), involve nearby t/s points with subscripts: 
c        (i  ,j+1)    (i+1,j+1)
c        (i  ,j  )    (i+1,j  )
c     t/s points (i,j), involve nearby u/v points with subscripts:
c        (i-1,j  )    (i  ,j  )
c        (i-1,j-1)    (i  ,j-1)
c     thus if qu(i,j) is defined on u/v points, its partial
c     difference ddxqt = ddxt(qu) is defined on t/s points and has the
c     value 
c     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
c                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
c-----------------------------------------------------------------------
c
      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5
c
      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5
c
      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5
c
      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5
c
c-----------------------------------------------------------------------
c     compute coefficients for all points
c-----------------------------------------------------------------------
c
c     initialize all 9 coefficients to zero
c
      do j=js,je
        do i1=-1,1
          do j1=-1,1
            do i=1,imt
              coef(i,j,i1,j1) = c0
            enddo
          enddo
	enddo
      enddo
c
      do j=js,je
        do i=1,imt
	  ustuff(i,j) = 0.0
	  vstuff(i,j) = 0.0
	enddo
      enddo
      do j=js,jecomp
        do i=1,imt-1
          ustuff(i,j) = h(i,j)*dyu(j)/(dxu(i)*csu(j))
          vstuff(i,j) = h(i,j)*dxu(i)*csu(j)/dyu(j)
        enddo
      enddo
c
c     calculate divergence = ddx (ddx (ustuff)) + ddy( ddy (vstuff))
c
      do j=jscomp,jecomp
        do i1=0,1
          do j1=0,1
            do i2=-1,0
              do j2=-1,0
                do i=2,imt-1
                  coef(i,j,i1+i2,j1+j2) = coef(i,j,i1+i2,j1+j2)
     &                 + cddxu(i1,j1) * cddxt(i2,j2) * ustuff(i+i2,j+j2)
     &                 + cddyu(i1,j1) * cddyt(i2,j2) * vstuff(i+i2,j+j2)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c
# ifdef timing
      call toc ('spc9pt', 'spc9pt')
# endif
      return
      end



      subroutine filz (fext, js, je)
# if defined fourfil || defined firfil
c
c=======================================================================
c     subroutine filz sets up input needed for fourier filtering
c     (when the "fourfil" ifdef is defined) -or- symmetric finite
c     impulse response filtering (when the "firfil" ifdef is defined)
c     of "fext" at the specified high latitudes. "fext" is forcing for
c     the external mode.
c=======================================================================
c
      use arrays_2d
#  include "param.h"
#  include "emode.h"
#  include "grdvar.h"
#  include "index.h"
#  include "taskrows.h"
      dimension fext(imt,js:je)
      dimension temp(imt)
c
      do jrow=max(jfrst,js),jecomp
        if (jrow.le.jft1 .or. jrow.ge.jft2) then
c
c---------------------------------------------------------------------
c          filter fext at high latitudes
c---------------------------------------------------------------------
c
            jj = jrow - jfrst + 1
c
            if (phi(jrow) .gt. c0) then
              jj = jj - jskpt + 1
              refcos = refcosn
            else
              refcos = refcoss
            endif
c
            do l=1,lsegf
              is = iszf(jj,l)
              if (is .ne. 0) then
                ie = iezf(jj,l)
                do ii=is,ie
                  i = mod(ii-2,imtm2) + 2
                  temp(ii+1-is) = fext(i,jrow)
                enddo
                im = ie-is+1
#  ifdef cyclic
                if (im .ne. imtm2) then
                   m = 1
                   n = nint(im*cst(jrow)/refcos)
                else
                   m = 3
                   n = nint(im*cst(jrow)*p5/refcos)
                endif
#  else
                m = 1
                n = nint(im*cst(jrow)/refcos)
#  endif
c
#  ifdef fourfil
                call filtr (temp(1), im, m ,n, 0)
#  endif 
#  ifdef firfil
                call filfir (temp(1), im, m, numflt(jrow))
#  endif
c
                do ii=is,ie
                  i = mod(ii-2,imtm2)+2
                  fext(i,jrow) = temp(ii+1-is)
                enddo
              endif
	    enddo
        endif
      enddo
# endif
      return
      end
#else
      subroutine dumexpl
      return
      end
#endif




