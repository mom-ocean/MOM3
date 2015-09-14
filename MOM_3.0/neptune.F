       subroutine neptune (js, je)
c
c=======================================================================
c     calculate neptune (maximum entropy) velocities
c
c     the option "neptune" provides a subgridscale parameterization 
c     for the interaction of eddies and topography
c
c     reference:
c       Holloway, G., 1992: representing topographic stress for large
c     scale ocean models, J. Phys. Oceanogr., 22, 1033-1046
c
c     neptune is calculated as an equilibrium streamfunction given by 
c     pnep=-f*snep*snep*hnep and is applied through eddy viscosity
c
c     hnep = model streamfunction depth
c     snep = spnep + (senep-spnep)*(0.5 + 0.5*cos(2.0*latitude))
c
c     the neptune length scale snep has a value of senep at the
c     equator and smoothly changes to spnep at the poles
c
c     author: Greg Holloway e-mail zounds@ios.bc.ca (or check mail_list)
c             Micahel Eby   e-mail eby@uvic.ca (or check mail_list)
c    (mpp items)  v. balaji e-mail=> vb@gfdl.gov
c             
c=======================================================================
c
       use mom3_mpp_mod
#ifdef neptune
      use arrays_2d
# ifdef rot_grid
      use memory_window
# endif
# include "param.h"
# include "coord.h"
# include "cnep.h"
# include "emode.h"
# include "grdvar.h"
# include "mw.h"
# include "scalar.h"
# include "taskrows.h"
c
      dimension kmz(imt,js:je)
c
c     compute an array to indicate "interior" streamfunction grid boxes
c
      if (js .eq. 1) then
        do i=1,imt
          kmz(i,1) = 0
        enddo
      endif
c
      do jrow=js,je
        kmz(1,jrow) = 0
      enddo
c
      do jrow=js+1,je
        do i=2,imt
          kmz(i,jrow) = min(kmu(i-1,jrow-1), kmu(i,jrow-1), 
     &                      kmu(i-1,jrow), kmu(i,jrow))
        enddo
      enddo
c
c     calculate the topographic stress equilibrium streamfunction
c
c     snep = spnep + (senep-spnep)*(0.5 + 0.5*cos(2.0*latitude))
c     pnep  = -f*snep*snep*hnep 
c
# ifdef rot_grid
c
c     if using a rotated grid, read geographic latitude from a file
c
      write (stdout,'(/,(1x,a))')
     & '==> Warning: when using "rot_grid", rotation of datasets      '
     &,'             must be done by the researcher following the     '
     &,'             examples in rotation.F. After this has been done,'
     &,'             remove the call abort() in neptune.F to proceed.         '
      if (imt .eq. imt) call abort()
c
      write (stdout,'(//a,a/)')
     & ' =>Reading velocity grid latitude from file "lat.dta" '
      call getunit (iolat, 'lat.dta'
     &,               'unformatted sequential rewind ieee')
      read (iolat) imax, jmax
c     use array unep to temporarily hold latitudes to save memory
c     unep(i,jrow,1) is for velocity grid 
c     unep(i,jrow,2) is for tracer grid 
      read (iolat) ((unep(i,jrow,1),i=1,imt),jrow=1,jmt)
      read (iolat) ((unep(i,jrow,2),i=1,imt),jrow=1,jmt)
      call relunit (iolat)
c
      do jrow=jscomp,jecomp
        do i=2,imtm1
          tlat = unep(i,jrow,2)/radian
          f    = c2*omega*sin(tlat)
          snep = spnep + (senep - spnep)*
     &           (p5 + p5*cos(c2*tlat))
# else
      do jrow=jscomp,jecomp
        tlat = yt(jrow)/radian
        f    = c2*omega*sin(tlat)
        snep = spnep + (senep - spnep)*
     &         (p5 + p5*cos(c2*tlat))
        do i=2,imtm1
# endif
c
c         find depth on streamfunction grid
c
          hnep = 0
          if (kmz(i,jrow) .ne. 0) then
	    hnep = zw(kmz(i,jrow))
	  endif
c
          pnep(i,jrow) = -f*snep*snep*hnep
c
        enddo
# ifdef cyclic
        pnep(1,jrow) = pnep(imtm1,jrow)
        pnep(imt,jrow) = pnep(2,jrow)
# endif  
      enddo
# ifdef parallel_1d
      call mpp_update_domains (pnep, mom_domain)
# endif
c
c     calculate depth independent velocity components from pnep
c
      do jrow=jscomp,jecomp
        do i=2,imtm1
          diag1 = pnep(i+1,jrow+1) - pnep(i  ,jrow)
          diag0 = pnep(i  ,jrow+1) - pnep(i+1,jrow)
          unep(i,jrow,1) = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
          unep(i,jrow,2) =  (diag1-diag0)*dxu2r(jrow)*hr(i,jrow)
     &                      *csur(jrow)
        enddo
# ifdef cyclic
          unep(1,jrow,1) = unep(imtm1,jrow,1)
          unep(1,jrow,2) = unep(imtm1,jrow,2)
          unep(imt,jrow,1) = unep(2,jrow,1)
          unep(imt,jrow,2) = unep(2,jrow,2)
# endif  
      enddo
# ifdef parallel_1d
      call mpp_update_domains (unep, mom_domain)
# endif  
c
#endif
      return
      end

