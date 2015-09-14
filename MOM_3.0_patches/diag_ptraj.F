#ifdef trajectories
      subroutine ptraji
c
c=======================================================================
c     initialize particle positions
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      character*128 file_name, file_stamp
#include "param.h"
#include "coord.h"
#include "diag_ptraj.h"
#include "iounit.h"
#include "switch.h"
#include "tmngr.h"
c      
c-----------------------------------------------------------------------
c     distribute particles within volume defined using starting and
c     ending longitudes, latitudes, and depths given by 
c     ptslon, ptelon, ptslat, ptelat. ptsdpt, and ptedpt.
c     for example: if nptraj = 500 particles then    
c
c     pacific distribution: drop 100 particles between 180E and 190E
c                           5S to 5N, and down to 200 meters. 
c
c      ptslon = 180.0
c      ptelon = 190.0
c      ptslat = -5.0
c      ptelat = 5.0
c      ptsdpt = 5.1e2
c      ptedpt = 200.e2
c      call pdist (1, 100)    
c
c     deep central atlantic distribution: drop 400 particles between
c     330E and 340E, 5S to 5N and between 200 and 3000 meters depth.
c
c      ptslon = 330.0
c      ptelon = 340.0
c      ptslat = -5.0
c      ptelat =  5.0
c      ptsdpt = 200.1e2
c      ptedpt = 3000.e2
c      call pdist (101, nptraj)    
c-----------------------------------------------------------------------
c
# ifdef timing
      call tic ('diagnostic', 'particle trajectories')
# endif
c
c     distribute all particles within volume defined by the following: 
c
      ptslon = 180.0
      ptelon = 190.0
      ptslat = -5.0
      ptelat = 5.0
      ptsdpt = 5.1e2
      ptedpt = 200.e2
      call pdist (1, nptraj)    
# ifdef lyapunov
c
c-----------------------------------------------------------------------
c     initialize deformation rate matrix
c-----------------------------------------------------------------------
c
      do n=1,nptraj
        em(1,1,n) = c1
	em(1,2,n) = c0
	em(2,1,n) = c0
	em(2,2,n) = c1
      enddo
# endif
c
c-----------------------------------------------------------------------
c     write initial particle positions to file 
c-----------------------------------------------------------------------
c
      reltim = relyr
c
      file_name = file_stamp ('particles_initial', stop_stamp, '.dta')
      call getunit (io, trim(file_name)
     &,             'unformatted sequential append ieee')
c
      iotext = 'read (iotraj) reltim'
      write (io) stamp, iotext, expnam
      write (io) reltim
c
# ifdef lyapunov
      iotext ='read (iotraj) nptraj, pxyz, pijk, em'
      write (io) stamp, iotext, expnam
      write (io) nptraj, pxyz, pijk, em
      write (stdout,*) 
     &' => Initial particle positions (lyapunov) written'
     &,    ' unformatted to file particles.dta on ts = ',itt,' ', stamp
# else
      iotext ='read (iotraj) nptraj, pxyz, pijk'
      write (io) stamp, iotext, expnam
      write (io) nptraj, pxyz, pijk
      write (stdout,*) ' => Initial particle positions written'
     &,    ' unformatted to file particles.dta on ts = ',itt,' ', stamp
# endif
      call relunit (io)
c
# ifdef timing
      call toc ('diagnostic', 'particle trajectories')
# endif
c
      return
      end


      subroutine pdist (ns, ne)
c
c=======================================================================
c     distribute n particles (ne-ns+1) within volume given by
c     ptslon, ptelon, ptslat, ptelat, ptsdpt, ptedpt by uniformly
c     placing approximately n**(1/3) particles along each dimension
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
# include "param.h"
# include "coord.h"
# include "diag_ptraj.h"
c
      if (ns .gt. nptraj .or. ne .gt. nptraj) then
        write (stdout,*) ' => Error: ns=',ns,', ne=',ne
	call abort()
      endif
c
c-----------------------------------------------------------------------
c     constrain the volume (containing starting positions of
c     particles) to lie within the model domain.      
c-----------------------------------------------------------------------
c
      if (ptslon .lt. xu(2))     ptslon = xu(2)
      if (ptslon .gt. xu(imtm1)) ptslon = xu(imtm1)
      if (ptslat .lt. yu(1))     ptslat = yu(1)
      if (ptslat .gt. yu(jmtm1)) ptslat = yu(jmtm1)
      if (ptsdpt .lt. zt(1))     ptsdpt = zt(1)
      if (ptsdpt .gt. zt(km))    ptsdpt = zt(km)
c
      if (ptelon .lt. xu(2))     ptelon = xu(2)
      if (ptelon .gt. xu(imtm1)) ptelon = xu(imtm1)
      if (ptelat .lt. yu(1))     ptelat = yu(1)
      if (ptelat .gt. yu(jmtm1)) ptelat = yu(jmtm1)
      if (ptedpt .lt. zt(1))     ptedpt = zt(1)
      if (ptedpt .gt. zt(km))    ptedpt = zt(km)
c
      if (ptslon .gt. ptelon) then
        t      = ptslon
        ptslon = ptelon
        ptelon = t
      endif
      if (ptslat .gt. ptelat) then
        t      = ptslat
        ptslat = ptelat
        ptelat = t
      endif
      if (ptsdpt .gt. ptedpt) then
        t      = ptsdpt
        ptsdpt = ptedpt
        ptedpt = t
      endif
c
c-----------------------------------------------------------------------
c     distribute the particles throughout the volume
c-----------------------------------------------------------------------
c
      cubr       = (float(ne-ns+1))**0.333333
      distx      = (ptelon - ptslon)/cubr
      disty      = (ptelat - ptslat)/cubr
      distz      = (ptedpt - ptsdpt)/cubr
      pxyz(1,ns) = p5*distx + ptslon
      pxyz(2,ns) = p5*disty + ptslat
      pxyz(3,ns) = p5*distz + ptsdpt
      do n=ns,ne
        if (n .gt. ns) then
	  pxyz(1,n) = pxyz(1,n-1)
	  pxyz(2,n) = pxyz(2,n-1)
	  pxyz(3,n) = pxyz(3,n-1)
	endif
	pxyz(1,n)  = pxyz(1,n) + distx
	if (pxyz(1,n) .gt. ptelon) then
	  pxyz(1,n) = ptslon + (pxyz(1,n)-ptelon)
	  pxyz(2,n) = pxyz(2,n) + disty
	  if (pxyz(2,n) .gt. ptelat) then
	    pxyz(2,n) = ptslat + (pxyz(2,n)-ptelat)
	    pxyz(3,n) = pxyz(3,n) + distz
	  endif
	endif
c
        i = indp (pxyz(1,n), xu, imt)
        if (xu(i) .gt. pxyz(1,n)) then
          pijk(1,n) = i
        else
          pijk(1,n) = i+1
        endif
c
        jrow = indp (pxyz(2,n), yu, jmt)
        if (yu(jrow) .gt. pxyz(2,n)) then
          pijk(2,n) = jrow
        else
          pijk(2,n) = jrow+1
        endif
c
        k = indp (pxyz(3,n), zt, km)
        if (zt(k) .gt. pxyz(3,n)) then
          pijk(3,n) = k
        else
          pijk(3,n) = k+1
        endif
      enddo
      write (stdout,9000)
      num_part = ne-ns+1
      write (stdout,*) num_part,' particles were initialized'
     &,   ' to lie within the volume described by:'
      write (stdout,'(1x,f8.2," < lon <",f8.2)') ptslon, ptelon
      write (stdout,'(1x,f8.2," < lat <",f8.2)') ptslat, ptelat
      write (stdout,'(1x,e10.3," < dpt < ",e10.3)') ptsdpt, ptedpt
      write (stdout,*) ' '
      return
9000  format (/20x,'P A R T I C L E    I N I T I A L I Z A T I O N'/)
      end


      subroutine ptraj (j, jrow)
c
c=======================================================================
c     integrate particle trajectories
c
c     note:
c
c     all indices refer to the "xu","yv" and "zt" grids.
c
c     it may be useful to interpolate other quantities (eg: potential
c     temperature ...) to the particle positions and save them as well.
c     this can help in understanding where and when such quantities
c     are conserved.
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use memory_window
      use mom3_mpp_mod
      use mpp_io_mod
      character*128 file_name, file_stamp, gvtitle
# include "param.h"
# include "coord.h"
# include "diag.h"      
# include "diag_ptraj.h"
# include "grdvar.h"
# include "iounit.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# include "tmngr.h"
# include "vers.h"      
      save file_name
      data num_traj /0/
      save num_traj
      save io

      type (axistype), save :: axis_n
      type (axistype), save :: axis_t
      
      type(fieldtype), save :: field_pos(3)
      type(fieldtype), save :: field_period
# ifdef lyapunov
      type(fieldtype), save :: field_em(2,2)
# endif
      dimension bufpart(nptraj)
      
# ifdef timing
      call tic ('diagnostic', 'particle trajectories')
# endif
c
c-----------------------------------------------------------------------
c     initialize so that every particle needs to be considered
c-----------------------------------------------------------------------
c
      if (jrow .eq. 2) then
        do n=1,nptraj
          ptdone(n) = .false.
        enddo
      endif
c
c-----------------------------------------------------------------------
c    calculate trajectory for all particles between jrow and jrow-1
c-----------------------------------------------------------------------
c
      rrad   = c1/radian
      cmdeg  = 8.982799e-8
c
      do n=1,nptraj
        if (.not. ptdone(n) .and. pijk(2,n) .eq. jrow) then
          ptdone(n) = .true.
c
c-----------------------------------------------------------------------
c         the particle is bounded by the volume with vertices given by
c         the eight nearest surrounding model grid points on the "xu",
c         "yu", and "zt" grids. (i,j,k) is the index of the deepest 
c         northeast corner of this bounding volume.
c-----------------------------------------------------------------------
c         
          i  = pijk(1,n)
c         j in the MW corresponds to  jrow = pijk(2,n)
          k  = pijk(3,n)
c
c-----------------------------------------------------------------------
c         compute volume weights for linear interpolation of velocity
c         at vertices of bounding volume to the particle position.
c
c         distances between particle and bounding volume faces
c
c         xe = distance to the east face
c         xw = distance to the west face
c         yn = distance to the north face
c         ys = distance to the south face
c         za = distance above to the top face
c         zb = distance below to the bottom face
c-----------------------------------------------------------------------
c
          xe = (xu(i) - pxyz(1,n))
          xw = (pxyz(1,n) - xu(i-1))
          yn = (yu(jrow) - pxyz(2,n))
          ys = (pxyz(2,n) - yu(jrow-1))
          za = (pxyz(3,n) - zt(k-1))
          zb = (zt(k) - pxyz(3,n))
	  dv = c1/((xu(i)-xu(i-1))*(yu(jrow)-yu(jrow-1))*(zt(k)-zt(k-1)))
c
c-----------------------------------------------------------------------
c         construct velocity at position of particle by 3-d linear
c         interpolation. 
c-----------------------------------------------------------------------
c
          xeyszb = xe*ys*zb*dv
          xwyszb = xw*ys*zb*dv
          xeysza = xe*ys*za*dv
          xwysza = xw*ys*za*dv
          xeynzb = xe*yn*zb*dv
          xwynzb = xw*yn*zb*dv
          xeynza = xe*yn*za*dv
          xwynza = xw*yn*za*dv
c          
          uu = u(i-1,k-1,j,1,tau)*xeyszb   + u(i,k-1,j,1,tau)*xwyszb
     &        +u(i-1,k  ,j,1,tau)*xeysza   + u(i,k  ,j,1,tau)*xwysza
     &        +u(i-1,k-1,j-1,1,tau)*xeynzb + u(i,k-1,j-1,1,tau)*xwynzb
     &        +u(i-1,k  ,j-1,1,tau)*xeynza + u(i,k  ,j-1,1,tau)*xwynza
c
          vv = u(i-1,k-1,j,2,tau)*xeyszb   + u(i,k-1,j,2,tau)*xwyszb
     &        +u(i-1,k  ,j,2,tau)*xeysza   + u(i,k  ,j,2,tau)*xwysza
     &        +u(i-1,k-1,j-1,2,tau)*xeynzb + u(i,k-1,j-1,2,tau)*xwynzb
     &        +u(i-1,k  ,j-1,2,tau)*xeynza + u(i,k  ,j-1,2,tau)*xwynza
c
c         interpolate vertical velocities at the bases of
c         the "u" cells.
c
          if (pxyz(3,n) .gt. zw(k-1)) then
	    za = pxyz(3,n) - zw(k-1)
	    zb = zw(k) - pxyz(3,n)
	    dv = c1/((xu(i)-xu(i-1))*(yu(jrow)-yu(jrow-1))
     &              *(zw(k)-zw(k-1)))
            xeyszb = xe*ys*zb*dv
            xwyszb = xw*ys*zb*dv
            xeysza = xe*ys*za*dv
            xwysza = xw*ys*za*dv
            xeynzb = xe*yn*zb*dv
            xwynzb = xw*yn*zb*dv
            xeynza = xe*yn*za*dv
            xwynza = xw*yn*za*dv
            ww = adv_vbu(i-1,k-1,j)*xeyszb   + adv_vbu(i,k-1,j)*xwyszb
     &          +adv_vbu(i-1,k  ,j)*xeysza   + adv_vbu(i,k  ,j)*xwysza
     &          +adv_vbu(i-1,k-1,j-1)*xeynzb + adv_vbu(i,k-1,j-1)*xwynzb
     &          +adv_vbu(i-1,k  ,j-1)*xeynza + adv_vbu(i,k  ,j-1)*xwynza
	  else
	    if (k-2 .eq. 0) then
	      za = pxyz(3,n)
	      zb = zw(k-1) - pxyz(3,n)
	      dv = c1/((xu(i)-xu(i-1))*(yu(jrow)-yu(jrow-1))
     &              *(zw(k-1)-0.0))
	    else
	      za = pxyz(3,n) - zw(k-2)
	      zb = zw(k-1) - pxyz(3,n)
	      dv = c1/((xu(i)-xu(i-1))*(yu(jrow)-yu(jrow-1))
     &              *(zw(k-1)-zw(k-2)))
	    endif
            xeyszb = xe*ys*zb*dv
            xwyszb = xw*ys*zb*dv
            xeysza = xe*ys*za*dv
            xwysza = xw*ys*za*dv
            xeynzb = xe*yn*zb*dv
            xwynzb = xw*yn*zb*dv
            xeynza = xe*yn*za*dv
            xwynza = xw*yn*za*dv
            ww = adv_vbu(i-1,k-2,j)*xeyszb   + adv_vbu(i,k-2,j)*xwyszb
     &          +adv_vbu(i-1,k-1,j)*xeysza   + adv_vbu(i,k-1,j)*xwysza
     &          +adv_vbu(i-1,k-2,j-1)*xeynzb + adv_vbu(i,k-2,j-1)*xwynzb
     &          +adv_vbu(i-1,k-1,j-1)*xeynza + adv_vbu(i,k-1,j-1)*xwynza
	  endif
# ifdef lyapunov
c
c-----------------------------------------------------------------------
c     construct the shear and compression above and below the particle
c-----------------------------------------------------------------------
c
c     for du/dx and du/dy
c
c
      uxn  = (u(i,k-1,j,1,tau) - u(i-1,k-1,j,1,tau))*dxtr(i)*csu(jrow)
      uxs  = (u(i,k-1,j-1,1,tau) - u(i-1,k-1,j-1,1,tau))*dxtr(i)
     &                                                  *csu(jrow-1)
      uxn2 = (u(i,k,j,1,tau) - u(i-1,k,j,1,tau))*dxtr(i)*csu(jrow)
      uxs2 = (u(i,k,j-1,1,tau) - u(i-1,k,j-1,1,tau))*dxtr(i)*csu(jrow-1)
c
      uye   = (u(i,k-1,j,1,tau) - u(i,k-1,j-1,1,tau))*dytr(jrow)
      uyw   = (u(i-1,k-1,j,1,tau) - u(i-1,k-1,j-1,1,tau))*dytr(jrow)
      uye2  = (u(i,k,j,1,tau) - u(i,k,j-1,1,tau))*dytr(jrow)
      uyw2  = (u(i-1,k,j,1,tau) - u(i-1,k,j-1,1,tau))*dytr(jrow)
c
c     for dv/dx and dv/dy
c
      vxn  = (u(i,k-1,j,2,tau) - u(i-1,k-1,j,2,tau))*dxtr(i)*csu(jrow)
      vxs  = (u(i,k-1,j-1,2,tau) - u(i-1,k-1,j-1,2,tau))*dxtr(i)
     &                                                  *csu(jrow-1)
      vxn2 = (u(i,k,j,2,tau) - u(i-1,k,j,2,tau))*dxtr(i)*csu(jrow)
      vxs2 = (u(i,k,j-1,2,tau) - u(i-1,k,j-1,2,tau))*dxtr(i)*csu(jrow-1)
c
      vye   = (u(i,k-1,j,2,tau) - u(i,k-1,j-1,2,tau))*dytr(jrow)
      vyw   = (u(i-1,k-1,j,2,tau) - u(i-1,k-1,j-1,2,tau))*dytr(jrow)
      vye2  = (u(i,k,j,2,tau) - u(i,k,j-1,2,tau))*dytr(jrow)
      vyw2  = (u(i-1,k,j,2,tau) - u(i-1,k,j-1,2,tau))*dytr(jrow)

c
c     interpolate du/dx, du/dy, dv/dx, and dv/dy to particle position
c
      dxr   = c1/((xw+xe)*(zb+za))
      dyr   = c1/((ys+yn)*(zb+za))
c
      yszb = ys*zb*dyr
      ynzb = yn*zb*dyr
      ysza = ys*za*dyr
      ynza = yn*za*dyr
c      
      xwzb = xw*zb*dxr
      xezb = xe*zb*dxr
      xwza = xw*za*dxr
      xeza = xe*za*dxr
c
      dudx = uxn*yszb + uxs*ynzb + uxn2*ysza + uxs2*ynza
      dudy = uye*xwzb + uyw*xezb + uye2*xwza + uyw2*xeza
c
      dvdx = vxn*yszb + vxs*ynzb + vxn2*ysza + vxs2*ynza
      dvdy = vye*xwzb + vyw*xezb + vye2*xwza + vyw2*xeza
c
c-----------------------------------------------------------------------
c     integrate the deformation matrix. note that this is not quite
c     correct when a particle encounters a boundary and "slips" along
c     the boundary.
c
c     the Lyapunov exponents can be computed from the eigenvalues of
c     the deformation rate matrix "em". the exponents are given by:
c     log(abs(eigen(1..2)))/T where T is the integration time.
c     let c = (em(2,2)-em(1,1))**2 + 4*(em(1,2)*em(2,1)) 
c     if (c >= 0.0) then
c       abs(eigen(1..2)) = abs(((em(2,2)+em(1,1))**2 + or - sqrt(c))/2)
c     else
c       abs(eigen(1)) = abs(eigen(2)) = 
c       sqrt((em(2,2)+em(1,1))**2 + abs(c)))/2
c     endif
c-----------------------------------------------------------------------
c
      em11 = em(1,1,n) + (dudx*em(1,1,n) + dudy*em(2,1,n))*dtuv
      em12 = em(1,2,n) + (dudx*em(1,2,n) + dudy*em(2,2,n))*dtuv
      em21 = em(2,1,n) + (dvdx*em(1,1,n) + dvdy*em(2,1,n))*dtuv
      em22 = em(2,2,n) + (dvdx*em(1,2,n) + dvdy*em(2,2,n))*dtuv
c
      em(1,1,n) = em11
      em(1,2,n) = em12
      em(2,1,n) = em21
      em(2,2,n) = em22
# endif
c
c-----------------------------------------------------------------------
c         remember where the particle was
c-----------------------------------------------------------------------
c
          xold = pxyz(1,n)
          yold = pxyz(2,n)
          zold = pxyz(3,n)
c
c-----------------------------------------------------------------------
c         integrate the particle trajectory forward for one time step
c         taking convergence of meridians into account.
c-----------------------------------------------------------------------
c
          rcos      = cmdeg/cos(pxyz(2,n)*rrad)
          pxyz(1,n) = pxyz(1,n) + dtuv*uu*rcos
          pxyz(2,n) = pxyz(2,n) + dtuv*vv*cmdeg
          pxyz(3,n) = pxyz(3,n) - dtuv*ww
c
c-----------------------------------------------------------------------
c         update bottom most northeast index of bounding volume
c-----------------------------------------------------------------------
c
          if (pxyz(1,n) .ge. xu(i)) then
            pijk(1,n) = i + 1
          else if (pxyz(1,n) .lt. xu(i-1)) then
            pijk(1,n) = i - 1
          endif
c
          if (pxyz(2,n) .ge. yu(jrow)) then
            pijk(2,n) = jrow + 1
          else if (pxyz(2,n) .lt. yu(jrow-1)) then
            pijk(2,n) = jrow - 1
          endif
c
          if (pxyz(3,n) .ge. zt(k)) then
            pijk(3,n) = k + 1
          else if (pxyz(3,n) .lt. zt(k-1)) then
            pijk(3,n) = k - 1
          endif
c
c-----------------------------------------------------------------------
c         do not allow any component of the trajectory to enter
c         land. If it does, reset it to its previous value
c         thereby simulating free slip conditions. hey... not perfect, 
c         but beats loosing particles in land. Also if the grid has
c         isolated "T,S" cells (ones where all eight surrounding
c         velocities are on land), replace references to "tmask(i,k,j)"
c         by umask(i,k,j) + umask(i-1,k,j) + umask(i,k,j-1) + 
c         umask(i-1,k,j-1) to prevent stagnation of particles if this
c         is a problem.
c-----------------------------------------------------------------------
c       
c
# ifdef cyclic
          if (pijk(1,n) .gt. imt) then
            if (tmask(3,k,j) .ne. 0) then
              pijk(1,n) = pijk(1,n) - (imt-2)
              i         = pijk(1,n)
              pxyz(1,n) = xu(i-1) + (pxyz(1,n)-xu(imt))
            else
              pijk(1,n) = i
              pxyz(1,n) = xold
            endif
          endif
          if (pijk(1,n) .lt. 2) then
            if (tmask(imt-2,k,j) .ne. 0) then
              pijk(1,n) = pijk(1,n) + (imt-2)
              i         = pijk(1,n)
              pxyz(1,n) = xu(i-1) + (pxyz(1,n)-xu(1))
            else
              pijk(1,n) = i
              pxyz(1,n) = xold
            endif
          endif
# endif
c
c-----------------------------------------------------------------------
c         constrain particles vertically to lie within ocean
c-----------------------------------------------------------------------
c
          if (pijk(3,n) .ne. k) then
            if (pijk(3,n) .eq. 1 .or. pijk(3,n) .gt. km) then
              pxyz(3,n) = zold
              pijk(3,n) = k
            else if (pijk(3,n) .gt. k .and. tmask(i,k+1,j) .eq. c0) then
              pxyz(3,n) = zold
              pijk(3,n) = k
            endif
          endif
c
c-----------------------------------------------------------------------
c         constrain particles longitudinally to stay within ocean 
c-----------------------------------------------------------------------
c
          if (pijk(1,n) .ne. i) then
            if (pijk(1,n) .gt. i .and. tmask(i+1,k,j) .eq. c0) then
              pxyz(1,n) = xold
              pijk(1,n) = i
            else if (pijk(1,n) .lt. i .and. tmask(i-1,k,j) .eq. c0) then
              pxyz(1,n) = xold
              pijk(1,n) = i
            endif 
          endif 
c
c-----------------------------------------------------------------------
c         constrain particles latitudinally to stay within ocean 
c-----------------------------------------------------------------------
c
          if (pijk(2,n) .ne. jrow) then
            if (pijk(2,n) .gt. jrow .and. tmask(i,k,j+1) .eq. c0) then
              pxyz(2,n) = yold
              pijk(2,n) = jrow
            else if (pijk(2,n) .lt. jrow .and. tmask(i,k,j-1) .eq. c0)
     &        then
              pxyz(2,n) = yold
              pijk(2,n) = jrow
            endif 
          endif 
        endif
      enddo
c
c-----------------------------------------------------------------------
c     write particle positions
c     note: last positions are also written to "iorest" for restarting
c
c     coded August 1999 by SM Griffies smg@gfdl.gov      
c-----------------------------------------------------------------------
c
      if ((jrow .eq. jmtm1) .and. trajts) then

        reltim = relyr
        file_name = file_stamp ('particles', stop_stamp, '.dta')
        
        write (stdout,9000) nptraj, stamp, pxyz(1,1), pxyz(2,1)
     &,                       pxyz(3,1)*0.01

        if (num_traj .eq. 0) then

          call mpp_open(io, file_name, action=MPP_WRONLY, 
     &         form=OUTPUT_FORM, threading=MPP_MULTI, fileset=MPP_MULTI,
     &         iospec= '-F cachea')
        
          num_traj = num_traj + 1
          file_name = file_stamp ('particles', stop_stamp, '.dta')
          gvtitle   = 'Particle trajectories ' // momver

c-----------------------------------------------------------------------
c         Define axes           
c-----------------------------------------------------------------------          

          do n=1,nptraj
	    bufpart(n) = n
	  enddo
          call mpp_write_meta(
     &      io, axis_n, 'N', 'particle', 'trajectories',
     &      cartesian='X', data=bufpart(1:nptraj))
          call mpp_write_meta(
     &         io, axis_t, 'Time', timeunits, 
     &         'Time since initial condition', cartesian='T')

c-----------------------------------------------------------------------
c         attributes of data 
c-----------------------------------------------------------------------

          call mpp_write_meta(
     &      io, field_pos(1), (/axis_n,axis_t/), 
     &      'x_pos', 'degrees_E', 'longitude of particle',
     &      xt(1), xt(imt), -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &      io, field_pos(2), (/axis_n,axis_t/), 
     &      'y_pos', 'degrees_N', 'latitude of particle',
     &      yt(1), yt(jmt), -1.0E+34, +1.0E+34)
          call mpp_write_meta(
     &      io, field_pos(3), (/axis_n,axis_t/), 
     &      'z_pos', 'cm', 'depth of particle',
     &      0.0, zw(km), -1.0E+34, +1.0E+34)
#  ifdef lyapunov
          call mpp_write_meta(
     &      io, field_em(1,1), (/axis_n,axis_t/), 
     &      'em11', '? units', 'deformation component 11',
     &      -1.e25, 1.e25, -1.0E+34, +1.0E+34)
           call mpp_write_meta(
     &      io, field_em(2,1), (/axis_n,axis_t/), 
     &      'em21', '? units', 'deformation component 21',
     &      -1.e25, 1.e25, -1.0E+34, +1.0E+34)
           call mpp_write_meta(
     &      io, field_em(1,2), (/axis_n,axis_t/), 
     &      'em12', '? units', 'deformation component 12',
     &      -1.e25, 1.e25, -1.0E+34, +1.0E+34)
           call mpp_write_meta(
     &      io, field_em(2,2), (/axis_n,axis_t/), 
     &      'em22', '? units', 'deformation component 22',
     &      -1.e25, 1.e25, -1.0E+34, +1.0E+34)
#  endif
           
c-----------------------------------------------------------------------
c         averaging period 
c-----------------------------------------------------------------------
           
          call mpp_write_meta(
     &      io, field_period, (/axis_t/), 
     &      'period', 'years', 'averaging period',
     &      +0.0E+00, +1.0E+20, -1.0E+34, +1.0E+34)
          
c-----------------------------------------------------------------------
c         figure title 
c-----------------------------------------------------------------------
          
          call mpp_write_meta(io, 'title', cval=gvtitle )

        endif

c-----------------------------------------------------------------------
c       write axes now that all meta data has been written
c-----------------------------------------------------------------------
      
        call mpp_write (io,axis_n)

c-----------------------------------------------------------------------
c       write trajectory data 
c-----------------------------------------------------------------------
        
        if (rununits.eq.'days') then
          timrec = reltim*yrlen
        else
          timrec = reltim
        endif
        period = c0
        
        do m=1,3
          do n=1,nptraj
	    bufpart(n) = pxyz(m,n)
          enddo
          call mpp_write(io, field_pos(m), bufpart, timrec)
        enddo
#  ifdef lyapunov
        do j=1,2
          do i=1,2
            do n=1,nptraj
	      bufpart(n) = em(i,j,n)
            enddo
          call mpp_write(io, field_em(i,j), bufpart, timrec)
          enddo
        enddo
#  endif
        call mpp_write(io, field_period, period, timrec)

        write (stdout,'(/a,a,a,i10,a,a)') 
     &  ' => Particle trajectories written to file '     
     &, trim(file_name),' on ts = ',itt,' ', stamp

      endif

# ifdef timing
      call toc ('diagnostic', 'particle trajectories')
# endif

      return
9000  format (1x,'=> Following',i6, ' particles at ',a32
     &,       '. (lon,lat,depth) of 1st particle is ('
     &,f8.2,',',f8.2,',',f8.2,'m)')
      end

#else
      subroutine ptraj (j, jrow)
      return
      end
#endif
