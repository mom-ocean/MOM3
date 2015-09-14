      subroutine impl_freesurf (js, je)
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
      use arrays_2d
      use memory_window
      use mom3_mpp_mod
      use mom3_io_mod
c
c=======================================================================
c
c     compute uncorrected barotropic velocities and their divergence
c     for r.h.s. of surface pressure equation (free surface or rigidlid)
c
c     Parallel Ocean General Circulation Modeling
c     R.D. Smith, J.K. Dukowicz and R.C. Malone
c     Physica D 60 (1992) 38-61
c     North-Holland
c     Elsevier Science Publishers B.V.
c
c     Implicit Free-Surface Method for the Bryan-Cox-Semtner Ocean
c     Model.
c     J.K. Dukowicz, R.D. Smith
c     Submitted to J. Geophysical Research   June 1993
c
c     js = global starting latitude index for the task
c     je = global ending latitude index for the task
c
c     author:       R. Smith             e-mail   rds@acl.lanl.gov
c                   A. Rosati            e-mail   ar@gfdl.gov
c     (mpp items)   v. balaji             e-mail=> vb@gfdl.gov
c=======================================================================
c
# include "param.h"
# include "coord.h"
# include "cprnts.h"
# include "diag.h"
# include "emode.h"
# include "grdvar.h"
# include "index.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# include "taskrows.h"
# include "tmngr.h"
      character *8 bc_symm
c
# ifdef timing
      call tic ('barotropic', 'impl_freesurf')
# endif
c
c-----------------------------------------------------------------------
c      construct uncorrected barotropic velocities
c      based on gradients of surface pressure
c-----------------------------------------------------------------------
c
      do jrow=js,je
        if (jrow .eq. 1 .or. jrow .eq. jmt) then
          do i=1,imt
            ubarm1(i,jrow,1)   = c0
            ubarm1(i,jrow,2)   = c0
          enddo
        endif
      enddo
c
      do jrow=jscomp,jecomp
# ifndef rot_grid
	fx = acor*c2dtsf*cori(jrow,1)
	fy = c1/(c1 + fx**2)
# endif
        factu = p5*c2dtsf*csur(jrow)
        factv = p5*c2dtsf*dyur(jrow)
c
c       leapfrog
c
        if (leapfrog) then 
	  do i=2,imtm1
# ifdef rot_grid
  	    fx = acor*c2dtsf*cori(i,jrow,1)
	    fy = c1/(c1 + fx**2)
# endif
	    d1 = gam*(ps(i+1,jrow+1,1) - ps(i,jrow,1))
     &       + (c1-gam)*(ps(i+1,jrow+1,2) - ps(i,jrow,2))
	    d2 = gam*(ps(i+1,jrow,1) - ps(i,jrow+1,1))
     &       + (c1-gam)*(ps(i+1,jrow,2) - ps(i,jrow+1,2))
            utwid = zu(i,jrow,1)*c2dtsf - factu*(d1 + d2)*dxur(i)
            vtwid = zu(i,jrow,2)*c2dtsf - factv*(d1 - d2)
c
c           calculate uncorrected velocities at time tau+1
c
            uhat(i,jrow,1) = fy*(utwid + fx*vtwid) + ubarm1(i,jrow,1)
# ifdef implicit_free_surface
     &                    + ubar(i,jrow,1)
# endif
            uhat(i,jrow,2) = fy*(vtwid - fx*utwid) + ubarm1(i,jrow,2)  
# ifdef implicit_free_surface
     &                   + ubar(i,jrow,2)
# endif
          enddo
        elseif (euler1 .or. forward) then
c
c         forward or 1st pass of euler
c
	  do i=2,imtm1
# ifdef rot_grid
  	    fx = acor*c2dtsf*cori(i,jrow,1)
	    fy = c1/(c1 + fx**2)
# endif
	    d1 = ps(i+1,jrow+1,1) - ps(i,jrow,1)
	    d2 = ps(i+1,jrow,1) - ps(i,jrow+1,1)
            utwid = zu(i,jrow,1)*c2dtsf - factu*(d1 + d2)*dxur(i)
            vtwid = zu(i,jrow,2)*c2dtsf - factv*(d1 - d2)
c
c           calculate uncorrected velocities at time tau+1
c
            uhat(i,jrow,1) = fy*(utwid + fx*vtwid)
     &	     + ubarm1(i,jrow,1)
            uhat(i,jrow,2) = fy*(vtwid - fx*utwid)
     &       + ubarm1(i,jrow,2)
          enddo
        elseif (euler2) then
c
c         euler 2nd pass
c
	  do i=2,imtm1
# ifdef rot_grid
  	    fx = acor*c2dtsf*cori(i,jrow,1)
	    fy = c1/(c1 + fx**2)
# endif
	    d1 = theta*(pguess(i+1,jrow+1) - pguess(i,jrow))
     &         + (c1-theta)*(ps(i+1,jrow+1,1) - ps(i,jrow,1))
	    d2 = theta*(pguess(i+1,jrow) - pguess(i,jrow+1))
     &       + (c1-theta)*(ps(i+1,jrow,1) - ps(i,jrow+1,1))
            utwid = zu(i,jrow,1)*c2dtsf - factu*(d1 + d2)*dxur(i)
            vtwid = zu(i,jrow,2)*c2dtsf - factv*(d1 - d2)
c
c           calculate uncorrected velocities at time tau+1
c
            uhat(i,jrow,1) = fy*(utwid + fx*vtwid)
     &                    + ubarm1(i,jrow,1)
            uhat(i,jrow,2) = fy*(vtwid - fx*utwid)
     &                    + ubarm1(i,jrow,2)
          enddo
        else
	  write (stdout,*) '=>Error: leapfrog, euler1, forward, euler2='
     &,   leapfrog, euler1, forward, euler2
          call abort()
        endif
c
        do i=2,imtm1
	  if (kmu(i,jrow) .eq. 0 )then
	    uhat(i,jrow,1) = c0
	    uhat(i,jrow,2) = c0
          endif
        enddo
      enddo
c
      call border (uhat(1,js,1), 'u even', js, je)
      call border (uhat(1,js,2), 'u odd', js, je)
c
# if defined fourfil || defined firfil
c
c-----------------------------------------------------------------------
c       filter uhat and vhat at high latitudes
c-----------------------------------------------------------------------
c
      do jrow=jscomp,jecomp
        if ((jrow.le.jfu1 .or. jrow.ge.jfu2) .and. jrow.ge.jfrst) then
          call filhat (uhat(1,jrow,1), uhat(1,jrow,2), jrow)
        endif
      enddo
# endif
# ifdef parallel_1d
      call mpp_update_domains (uhat, mom_domain)
# endif
c
c-----------------------------------------------------------------------
c     form divergence of uncorrected barotropic velocities for
c     r.h.s. of surface pressure eqn
c-----------------------------------------------------------------------
c
c     compute the divergence of the sea surface pressure gradients
c
      call spforc (uhat, dxu, dyu, csu, h, divf, js, je)
c
      fx = c1/(apgr*c2dtsf)
      do jrow=jscomp,jecomp
        do i=2,imtm1
	  divf(i,jrow) = fx*divf(i,jrow)
	enddo
      enddo
c
# ifdef implicit_free_surface
      if (euler2) then
c
c       correct r.h.s on 2nd pass euler backward
c
        fx = c1/(grav*dtsf*c2dtsf*apgr)
        do jrow=jscomp,jecomp
	  fy = dyt(jrow)*cst(jrow)*fx
	  do i=1,imtm1
	    divf(i,jrow) = divf(i,jrow) 
     &                   + fy*dxt(i)*(pguess(i,jrow) - ps(i,jrow,1))
          enddo
        enddo
      endif
# endif
c
      call border (divf, 't even', js, je)
# ifdef parallel_1d
      call mpp_update_domains (divf, mom_domain)
# endif
c     
c-----------------------------------------------------------------------
c     calculate guess for change in pressure from guess for pressure 
c-----------------------------------------------------------------------
c
      if (leapfrog) then
        do jrow=js,je
          do i=1,imt
            ptd(i,jrow) = pguess(i,jrow)-ps(i,jrow,2)
          enddo
        enddo
      elseif (euler1 .or. forward) then
        do jrow=js,je
          do i=1,imt
            ptd(i,jrow) = pguess(i,jrow)-ps(i,jrow,1)
          enddo
        enddo
      else
        do jrow=js,je
          do i=1,imt
            ptd(i,jrow) = c0 
          enddo
        enddo
      endif
c
c-----------------------------------------------------------------------
c     use conjugate gradient 9 point laplacian to solve directly for
c     the change in surface pressure.
c-----------------------------------------------------------------------
c
c     initialize coefficients for the conjugate gradient solver
c
      call spc9pt (dxu, dyu, csu, h, cf, js, je)
c
c-----------------------------------------------------------------------
c     add diag term to central coeff in the implicit free surface
c-----------------------------------------------------------------------
c
# ifdef implicit_free_surface
      do jrow=jscomp,jecomp
	fx = cst(jrow)*dyt(jrow)/(apgr*c2dtsf*dtsf*grav)
	do i=2,imtm1
          if (map(i,jrow) .le. 0) then
	    cf(i,jrow,0,0) = cf(i,jrow,0,0) - fx*dxt(i)
	  endif
        enddo
      enddo
# endif
      npt    = 9
      variable   = 'd(surf press)'
      nislsp = 0
      bc_symm = 't even'
# ifdef implicit_free_surface
      tolr   = tolrfs
# else
      tolr   = tolrsp
# endif
c
# ifdef timing
      call tic ('congrad', 'surface pressure (impl_freesurf)')
# endif
      call congr (npt, variable, bc_symm, ptd, ptd, divf, res 
     &,           cf, js, je
     &,           mxscan, mscan, tolr
     &,           imask, iperm, jperm, iofs, nislsp, nippts
     &,           converged, esterr)
# ifdef timing
      call toc ('congrad', 'surface pressure (impl_freesurf)')
# endif
c
# ifndef implicit_free_surface
c
c     remove null space from ptd (rigid lid only)
c
      call checkerboard (ptd, map, js, je)
      call border (ptd, bc_symm, js, je)
c
c     remove mean
c
      call zero_level (ptd, 'surf press', map, dxt, dyt, cst)
      call border (ptd, bc_symm, js, je)
# endif
c
c-----------------------------------------------------------------------
c     correct barotropic velocities with change in surface pressure
c     gradient. this completes the calculation of barotropic
c     velocity at tau+1
c     also update barotropic velocities for next time step
c-----------------------------------------------------------------------
c
      do jrow=jscomp,jecomp
        factu = p5*apgr*c2dtsf*csur(jrow)
        factv = p5*apgr*c2dtsf*dyur(jrow)
        do i=2,imtm1
          d1 = ptd(i+1,jrow+1) - ptd(i,jrow)
          d2 = ptd(i+1,jrow) - ptd(i,jrow+1) 
c
c         tau + 1
c
	  tempu = uhat(i,jrow,1) - factu*(d1 + d2)*dxur(i)
	  tempv = uhat(i,jrow,2) - factv*(d1 - d2)
	  if (leapfrog) then
# ifdef implicit_free_surface
	    tempu = tempu - ubar(i,jrow,1)
	    tempv = tempv - ubar(i,jrow,2)
# endif
# ifdef robert_time_filter
c
c           time filter ubar
c
	    ubar(i,jrow,1) = ubar(i,jrow,1)
     &      + smooth*(0.5*(tempu + ubarm1(i,jrow,1))
     &      - ubar(i,jrow,1))
            ubar(i,jrow,2) = ubar(i,jrow,2)
     &      + smooth*(0.5*(tempv + ubarm1(i,jrow,2))
     &      - ubar(i,jrow,2))
# endif
c
c           tau - 1  <= tau
c
	    ubarm1(i,jrow,1) = ubar(i,jrow,1)
	    ubarm1(i,jrow,2) = ubar(i,jrow,2)
	  endif
c
c         tau  <= tau + 1
c
          ubar(i,jrow,1) = tempu
          ubar(i,jrow,2) = tempv
        enddo
        do i=2,imtm1
          if (kmu(i,jrow) .eq. 0 ) then
            ubar(i,jrow,1)   = c0
            ubar(i,jrow,2)   = c0
            ubarm1(i,jrow,1) = c0
            ubarm1(i,jrow,2) = c0
          endif
        enddo
      enddo
c
      call border (ubar(1,js,1), 'u even', js, je)
      call border (ubar(1,js,2), 'u odd', js, je)
      call border (ubarm1(1,js,1), 'u even', js, je)
      call border (ubarm1(1,js,2), 'u odd', js, je)
c
# ifdef parallel_1d
      call mpp_update_domains (ubar,   mom_domain)
      call mpp_update_domains (ubarm1, mom_domain)
# endif
c
c-----------------------------------------------------------------------
c     update the surface pressure based upon the relaxation solution
c-----------------------------------------------------------------------
c
      if (leapfrog) then
c
c       leapfrog
c
        do jrow=js,je
	  do i=1,imt
            pnew           = ptd(i,jrow) + ps(i,jrow,2)
            pguess(i,jrow) = c3*(pnew - ps(i,jrow,1)) + ps(i,jrow,2)
            ps(i,jrow,2)   = ps(i,jrow,1)
            ps(i,jrow,1)   = pnew
          enddo
	enddo
      elseif (euler1) then
c
c       eb 1st pass
c
        do jrow=js,je
          do i=1,imt
# ifdef implicit_free_surface
	    pnew = ptd(i,jrow) + ps(i,jrow,1)
	    pguess(i,jrow) = pnew
            ps(i,jrow,1)   = pnew
# else
	    pnew = ptd(i,jrow) + ps(i,jrow,1)
	    pguess(i,jrow) = pnew
# endif
          enddo
	enddo
      elseif (forward) then
c
c       forward
c
        do jrow=js,je
          do i=1,imt
	    pnew           = ptd(i,jrow) + ps(i,jrow,1)
	    pguess(i,jrow) = c3*(pnew - ps(i,jrow,1)) + ps(i,jrow,2)
            ps(i,jrow,2)   = ps(i,jrow,1)
            ps(i,jrow,1)   = pnew
          enddo
	enddo
      elseif (euler2) then
c
c       eb 2nd pass
c
        do jrow=js,je
          do i=1,imt
	    pnew           = ptd(i,jrow) + pguess(i,jrow)
	    pguess(i,jrow) = c3*(pnew - ps(i,jrow,1)) + ps(i,jrow,2)
            ps(i,jrow,2)   = ps(i,jrow,1)
            ps(i,jrow,1)   = pnew
          enddo
	enddo
      endif
c
# ifdef remove_ps_checkerboard
#  if !defined implicit_free_surface
c
c     test accumulation of residual checkerboard 
c
      call checkerboard(ps(1,js,1), map, js, je)
      call border (ps(1,js,1), bc_symm, js, je)
      call zero_level (ps(1,js,1), 'surf press', map, dxt, dyt, cst)
      call border (ps(1,js,1), bc_symm, js, je)
#  endif
# endif
# ifdef timing
      call toc ('barotropic', 'impl_freesurf')
# endif
c
# ifdef parallel_1d
      call mpp_update_domains (ps, mom_domain)
# endif
#endif
      return
      end





#if !defined stream_function
# if defined fourfil || defined firfil
      subroutine filhat (uhati, vhati, jrow)
      use arrays_2d
      use memory_window
c
c=====================================================================
c     filhat sets up input needed for fourier filtering
c     (when the "fourfil" option is defined) -or- symmetric finite
c     impulse response filtering (when ifdef "firfil" is defined) of
c     barotropic velocities at the specifiied high latitude row "jrow".
c
c     author: r.c.pacanowski     e-mail => rcp@gfdl.gov
c=======================================================================
c
#  include "param.h"
#  include "coord.h"
#  include "cpolar.h"
#  include "emode.h"
#  include "grdvar.h"
#  include "index.h"
#  include "mw.h"
#  include "scalar.h"
#  include "switch.h"
      dimension tempik(imt,2), uhati(imt), vhati(imt)
c
#  ifdef timing
      call tic ('free_surf', 'filtering')
#  endif
c
c---------------------------------------------------------------------
c     filter uhati and vhati in polar latitudes
c---------------------------------------------------------------------
c
      jj = jrow - jfrst + 1
c
      if (phi(jrow) .gt. c0) then
        jj = jj - jskpu + 1
        refcos = refcosn
        fx = c1
      else
        refcos = refcoss
        fx = -c1
      endif
c
      isave = 0
      ieave = 0
c
      do l=1,lsegf
        k=1
        if (isuf(jj,l,k) .ne. 0) then
          is = isuf(jj,l,k)
          ie = ieuf(jj,l,k)
          iredo = 1
          if (is.ne.isave .or. ie.ne.ieave) then
            iredo = 0
            im = ie - is + 1
            isave = is
            ieave = ie
#  ifdef cyclic
            if (im .ne. imtm2) then
              m = 2
              n = nint(im*csu(jrow)/refcos)
            else
              m = 3
              n = nint(im*csu(jrow)*p5/refcos)
            endif
#  else
            m = 2
            n = nint(im*csu(jrow)/refcos)
#  endif
          endif
c
c         rotate vector to compensate for polar latitude
c
          ism1 = is - 1
          iea = ie
          if (ie .ge. imt) iea = imtm1
          do i=is,iea
            tempik(i-ism1,1) = -fx*uhati(i)*spsin(i)
     &                             - vhati(i)*spcos(i)
            tempik(i-ism1,2) =  fx*uhati(i)*spcos(i)
     &                             - vhati(i)*spsin(i)
          enddo
          if (ie .ge. imt) then
            ieb = ie - imtm2
            ii  = imtm1 - is
            do i=2,ieb
              tempik(i+ii,1) = -fx*uhati(i)*spsin(i)
     &                             - vhati(i)*spcos(i)
              tempik(i+ii,2) =  fx*uhati(i)*spcos(i)
     &                            -  vhati(i)*spsin(i)
            enddo
          endif
#  ifdef fourfil
          call filtr (tempik(1,1), im, m, n, iredo)
          call filtr (tempik(1,2), im, m, n, 1)
#  endif
#  ifdef firfil
          call filfir (tempik(1,1), im, m, numflu(jrow))
          call filfir (tempik(1,2), im, m, numflu(jrow))
#  endif
c
c         rotate vector back to compensate for polar latitude
c
          do i=is,iea
            uhati(i) = fx*(-tempik(i-ism1,1)*spsin(i)
     &                   + tempik(i-ism1,2)*spcos(i))
            vhati(i) = -tempik(i-ism1,1)*spcos(i)
     &                   - tempik(i-ism1,2)*spsin(i)
          enddo
          if (ie .ge. imt) then
            do i=2,ieb
              uhati(i) = fx*(-tempik(i+ii,1)*spsin(i)
     &                     + tempik(i+ii,2)*spcos(i))
              vhati(i) = -tempik(i+ii,1)*spcos(i)
     &                     - tempik(i+ii,2)*spsin(i)
            enddo
          endif
        endif
      enddo
#  ifdef cyclic
      uhati(1)   = uhati(imt-1)
      uhati(imt) = uhati(2)
      vhati(1)   = vhati(imt-1)
      vhati(imt) = vhati(2)
#  else
      uhati(1)   = c0
      uhati(imt) = c0
      vhati(1)   = c0
      vhati(imt) = c0
#  endif
c
#  ifdef timing
      call toc ('free_surf', 'filtering')
#  endif
      return
      end
# endif
#endif
