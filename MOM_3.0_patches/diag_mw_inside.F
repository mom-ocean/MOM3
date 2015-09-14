      subroutine diag_mw_inside (joff, js, je, is, ie)
c
c=======================================================================
c     calculate diagnostics inside the mw
c
c     input:
c
c      joff   = offset between row j in the MW and latitude jrow on disk
c      js     = starting row for calculations
c      je     = ending row for calculations
c      is     = starting longitude index for calculations
c      ie     = ending longitude index for calculations
c=======================================================================
c
      use mpp_io_mod
#if defined bbl_ag || defined explicit_fresh_water_flux
      use arrays_2d
#endif
      use memory_window
      character*128 file_name, file_stamp
#include "param.h"
#include "taskrows.h"
#ifdef bbl_ag
# include "bbl.h"
#endif
#include "coord.h"
#ifdef matrix_sections
# include "cprnts.h"
#endif
#include "diag.h"
#include "docnam.h"
#include "grdvar.h"
#include "iounit.h"
#ifdef isoneutralmix
# include "isopyc.h"
#endif
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"
#include "vmixc.h"
c
#ifdef isoneutralmix
      dimension aibuf(imt,km)
#endif
#if defined meridional_overturning || defined meridional_overturning_basin
# include "state.h"
      dimension vbarx(0:km-1),sigm(imt,km)
# include "dens.h"
#endif
c-----------------------------------------------------------------------
c     bail out if starting row exceeds ending row
c-----------------------------------------------------------------------
c
      if (js .le. je) then
c
c-----------------------------------------------------------------------
c     limit longitudes
c-----------------------------------------------------------------------
c
      istrt  = max(2,is)
      iend   = min(imt-1,ie)
c
c
#ifdef time_averages
c
c-----------------------------------------------------------------------
c     diagnostic: accumulate "tau" data for time means 
c     author:  s.m. griffies e-mail  smg@gfdl.gov
c-----------------------------------------------------------------------
c
      if (timavgperts .and. .not. euler2) then
        call timeavg_out(joff, js, je, .true.)
      endif
#endif
c
#ifdef snapshots
c
c-----------------------------------------------------------------------
c     diagnostic: write instantaneous sample of MOM data
c     author:    r.c.pacanowski   e-mail  rcp@gfdl.gov
c                a. rosati        e-mail  ar@gfdl.gov
c-----------------------------------------------------------------------
c
      if (snapts .and. .not. euler2) then
	if (istrt .ne. 2 .and. iend .ne. imt-1) then
	  write (stdout,*) '=>Error: istrt = ',istrt,' and iend ='
     &,   iend,' are not allowed when calling "snapshots_out"'
          call abort()
	else
          call snapshots_out (joff, js, je, .true.)
	endif
      endif
#endif
#ifdef snapshots_regions
      do nreg=1, nregions	
        if (snapts_reg(nreg) .and. .not. euler2) then
          call region_out 
     &         (istrt, iend, joff, js, je, .true., nreg)
        endif
      enddo
#endif
c
#ifdef stability_tests
c
c-----------------------------------------------------------------------
c     diagnostic: compute stability diagnostics
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (stabts .and. eots) then
	if (istrt .ne. 2 .and. iend .ne. imt-1) then
	  write (stdout,*) '=>Error: istrt = ',istrt,' and iend ='
     &,   iend,' are not allowed when calling "stab"'
          call abort()
	else
          do j=js,je
            jrow = joff + j
	    call stab (j, jrow)
          enddo
	endif
      endif
#endif
c
#ifdef trajectories
c
c-----------------------------------------------------------------------
c     diagnostic: integrate particle trajectories
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (eots) then
	if (istrt .ne. 2 .and. iend .ne. imt-1) then
	  write (stdout,*) '=>Error: istrt = ',istrt,' and iend ='
     &,   iend,' are not allowed when calling "ptraj"'
          call abort()
	else
          do j=js,je
            jrow = joff + j
	    call ptraj (j, jrow)
          enddo
        endif
      endif
#endif
c
#ifdef meridional_overturning
c
c-----------------------------------------------------------------------
c     construct meridional overturning streamfunction
c     authors:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c               s. m. griffies   e-mail  smg@gfdl.gov
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = joff + j
        if (jrow .lt. jmtm1 .and. vmsfts .and. eots) then
# ifdef timing
          call tic ('diagnostic', 'meridional overturning')
# endif
c
c         Streamfunction with depth as vertical coordinate
c 
          rsverdrup = 1.0e-12
          do k=1,km
            vbarx(k-1) = c0
          enddo
	  vmsf(jrow,km) = c0

# ifdef bbl_ag
          do i=istrt,iend
	    kb=kmt(i,jrow)
	    kb1=kmt(i+1,jrow)
	    kbn=kmt(i,jrow+1)
	    kbn1=kmt(i+1,jrow+1)
c
	    kbmax = max(kb,kbn)
	    kbmax1 = max(kb1,kbn1)
	    kbmin = min(kb,kbn)
	    kbmin1 = min(kb1,kbn1)
c
c           Treat flow as occurring at maximum of kb and kb1, unless
c           it detrains, when it occurs at the minimum of kb and kb1
c
	    vbarx(kbmax)=vbarx(kbmax) - u(i,km,j,2,tau)
     &           *csu(jrow)*dxu(i)
     &           *dzt(km)*0.5
     &           *oppsigny(i,j)
	    vbarx(kbmin)=vbarx(kbmin) - u(i,km,j,2,tau)
     &           *csu(jrow)*dxu(i)
     &           *dzt(km)*0.5*(1.0-oppsigny(i,j))
	    vbarx(kbmax1)=vbarx(kbmax1) - u(i,km,j,2,tau)
     &           *csu(jrow)*dxu(i)
     &           *dzt(km)*0.5
     &           *oppsigny(i+1,j)
	    vbarx(kbmin1)=vbarx(kbmin1) - u(i,km,j,2,tau)
     &           *csu(jrow)*dxu(i)
     &           *dzt(km)*0.5*(1-oppsigny(i+1,j))
          enddo



          do k=km-1,1,-1
# else
          do k=km,1,-1
# endif
            do i=istrt,iend
              vtot = u(i,k,j,2,tau)
# ifdef gent_mcwilliams
     &         + umask(i,k,j)
     &         *(adv_vntiso(i,k,j)+adv_vntiso(i+1,k,j))/(2.0*csu(jrow))
# endif
              vbarx(k-1) = vbarx(k-1) - vtot*csu(jrow)*dxu(i)
# ifdef partial_cell
     &                                            *dhu(i,k,j)
# else
     &                                            *dzt(k)
# endif
            enddo
            vmsf(jrow,k-1) = vmsf(jrow,k) + vbarx(k-1)*rsverdrup
          enddo
c
c
c         when vmsf=0 at the ocean bottom, the plotting routines cannot
c         plot the contours correctly because of the "land" at the
c         bottom, so reset vmsf to non-zero at the bottom if the value
c         of vmsf at km-1 is non-zero
c
          if (vmsf(jrow,km-1) .ne. c0) vmsf(jrow,km) = 1.0e-12
c
# if defined gent_mcwilliams
          do k=1,km
            vbarx(k-1) = c0
          enddo
          vmsf_gm(jrow,km) = c0
#  ifdef bbl_ag
          vmsf_gm(jrow,km-1) = c0
          do k=km-1,1,-1
#  else
          do k=km,1,-1
#  endif
            do i=istrt,iend
              vstargm = umask(i,k,j)
     &         *(adv_vntiso(i,k,j)+adv_vntiso(i+1,k,j))/(2.0*csu(jrow))
              vbarx(k-1) = vbarx(k-1) - vstargm*csu(jrow)*dxu(i)
#  ifdef partial_cell
     &                                            *dhu(i,k,j)
#  else
     &                                            *dzt(k)
#  endif
            enddo
            vmsf_gm(jrow,k-1) = vmsf_gm(jrow,k) + vbarx(k-1)*rsverdrup
          enddo
          if (vmsf_gm(jrow,km-1) .ne. c0) vmsf_gm(jrow,km) = 1.0e-12
# endif
# if defined meridional_overturning_density
c
c-----------------------------------------------------------------------
c         Streamfunction with potential density 
c         as vertical coordinate: sig = (rho - 1)*1000.
c         Can employ various sigmas, referenced to different 
c         depth levels.  The trick is to define the vertical
c         sigma levels in a fashion to capture the essential 
c         aspects of the streamfunction.  Some tuning of sigAvert
c         will likely be necessary. 
c
c         NOTE: Effects from gent_mcwilliams have not yet been 
c         included in this streamfunction.  
c
c         I. First sigma vertical coordinate
c         USER INPUT:
c          (1) Define the reference level kref for use in defining
c              vertical coordinate.
c          (2) Define the values for the vertical coordinate sig(km)
c-----------------------------------------------------------------------
c
          depth_ref = 0.0e2
          kref = indp (depth_ref, zt, km)
c
          sigAvert(1) = 23.0
          do m=2,km
           sigAvert(m) = sigAvert(m-1) + .3
          enddo 
c
c         END OF USER INPUT for first sigma.
c
          csigma = 1000.0
c
c         Evaluate sigma for the ocean model grid points.
c
          do k=1,km
            do i=istrt,iend
              tprime = t(i,k,j,1,tau)-to(kref)
	      sprime = t(i,k,j,2,tau)-so(kref)
              sigm(i,k)=csigma*(ro0(kref)+dens(tprime,sprime,kref))
            enddo
          enddo          
c
c         Compute the overturning streamfunction by integrating 
c         zonally and vertically.  The vertical integration extends
c         only up to the depth of the particular density surface. 
c
          do m=1,km
            do i=istrt,iend
              do k=1,km 
                if (sigm(i,k) .ge. sigAvert(m)) then 
                  sigAmsf(jrow,m) = sigAmsf(jrow,m)
     &             - u(i,k,j,2,tau)*csu(jrow)*dxu(i)
#  ifdef partial_cell
     &                                              *dhu(i,k,j)
#  else
     &                                              *dzt(k)
#  endif
                endif
                if (sigm(i,k) .lt. sigAvert(1)) then 
                  sigAmsf(jrow,1) = sigAmsf(jrow,1)
     &             - u(i,k,j,2,tau)*csu(jrow)*dxu(i)
#  ifdef partial_cell
     &                                              *dhu(i,k,j)
#  else
     &                                              *dzt(k)
#  endif
                endif
              enddo
            enddo
          enddo
          do m=1,km
            sigAmsf(jrow,m) = sigAmsf(jrow,m)*rsverdrup
          enddo
c
c-----------------------------------------------------------------------
c         II. Second sigma
c         USER INPUT:
c           (1) Define the reference level kref for use in defining
c               vertical coordinate.
c           (2) Define the values for the vertical coordinate sig(km)
c-----------------------------------------------------------------------
c
          depth_ref = 2000.0e2
          kref = indp (depth_ref, zt, km)
c
          sigBvert(1) = 30.0
          do m=2,km
           sigBvert(m) = sigBvert(m-1) + .3
          enddo 
c
c         END OF USER INPUT for second sigma.
c
c         Evaluate sigma for the ocean model grid points.
c
          do k=1,km
            do i=istrt,iend
              tprime = t(i,k,j,1,tau)-to(kref)
	      sprime = t(i,k,j,2,tau)-so(kref)
              sigm(i,k)=csigma*(ro0(kref)+dens(tprime,sprime,kref))
            enddo
          enddo          
c
c         Compute the overturning streamfunction by integrating 
c         zonally and vertically.  The vertical integration extends
c         only up to the depth of the particular 
c         density surface. 
c
          do m=1,km
            do i=istrt,iend
              do k=1,km 
                if (sigm(i,k) .ge. sigBvert(m)) then 
                  sigBmsf(jrow,m) = sigBmsf(jrow,m)
     &             - u(i,k,j,2,tau)*csu(jrow)*dxu(i)
#  ifdef partial_cell
     &                                              *dhu(i,k,j)
#  else
     &                                              *dzt(k)
#  endif
                endif
                if (sigm(i,k) .lt. sigBvert(1)) then 
                  sigBmsf(jrow,1) = sigBmsf(jrow,1)
     &             - u(i,k,j,2,tau)*csu(jrow)*dxu(i)
#  ifdef partial_cell
     &                                              *dhu(i,k,j)
#  else
     &                                              *dzt(k)
#  endif
                endif
              enddo
            enddo
          enddo
          do m=1,km
            sigBmsf(jrow,m) = sigBmsf(jrow,m)*rsverdrup
          enddo
c
c-----------------------------------------------------------------------
c         Other sigmas can be used for defining the vertical coordinate.
c         Simply do the following: 
c         1. Emulate the above code block.
c         2. Add  common /cdiag/ sigCmsf(jmt,km),sigBvert(km) to diag.h.
c         3. Initialize sigCmsf(jrow,k) to zero inside diag_mw_inside.F
c         4. call vms_sigma_outC inside of diag_mw_outside.F
c         5. Create a subroutine vms_sigma_outC inside diag_mw_outside.F
c            like subroutines vms_sigma_outA, and vms_sigma_outB
c            which are already there. 
c-----------------------------------------------------------------------
c
# endif
# ifdef timing
          call toc ('diagnostic', 'meridional overturning')
# endif
        endif
      enddo
#endif
c
#ifdef meridional_overturning_basin
c
c-----------------------------------------------------------------------
c       meridional overturning of mass by basin number
c       0=global, 1=Atlantic, 2=Pacific
c       other basin number distributions can be setup by researcher.
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = joff + j
        if (jrow .lt. jmtm1 .and. vmsfts .and. eots) then
# ifdef timing
          call tic ('diagnostic', 'meridional overturning')
# endif
c
c         Streamfunction with depth as vertical coordinate
c
          rsverdrup = 1.0e-12
          do L=0,nbasin
c
            do k=1,km
              vbarx(k-1) = c0
            enddo
            vmsf(jrow,km,l) = c0
# ifdef bbl_ag
          do i=istrt,iend
	    kb=kmt(i,jrow)
	    kb1=kmt(i+1,jrow)
	    kbn=kmt(i,jrow+1)
	    kbn1=kmt(i+1,jrow+1)
c
	    kbmax = max(kb,kbn)
	    kbmax1 = max(kb1,kbn1)
	    kbmin = min(kb,kbn)
	    kbmin1 = min(kb1,kbn1)
c
c           Treat flow as occurring at maximum of kb and kb1, unless
c           it detrains, when it occurs at the minimum of kb and kb1
c
	    vbarx(kbmax)=vbarx(kbmax) - u(i,km,j,2,tau)
     &           *basinmask(i,jrow,l)
     &           *csu(jrow)*dxu(i)
     &           *dzt(km)*0.5
     &           *oppsigny(i,j)
	    vbarx(kbmin)=vbarx(kbmin) - u(i,km,j,2,tau)
     &           *basinmask(i,jrow,l)
     &           *csu(jrow)*dxu(i)
     &           *dzt(km)*0.5*(1.0-oppsigny(i,j))
	    vbarx(kbmax1)=vbarx(kbmax1) - u(i,km,j,2,tau)
     &           *basinmask(i,jrow,l)
     &           *csu(jrow)*dxu(i)
     &           *dzt(km)*0.5
     &           *oppsigny(i+1,j)
	    vbarx(kbmin1)=vbarx(kbmin1) - u(i,km,j,2,tau)
     &           *basinmask(i,jrow,l)
     &           *csu(jrow)*dxu(i)
     &           *dzt(km)*0.5*(1-oppsigny(i+1,j))
          enddo


            do k=km-1,1,-1
# else
            do k=km,1,-1
# endif
              do i=istrt,iend
              vtot = u(i,k,j,2,tau)
# ifdef gent_mcwilliams
     &         + umask(i,k,j)
     &         *(adv_vntiso(i,k,j)+adv_vntiso(i+1,k,j))/(2.0*csu(jrow))
# endif
               vbarx(k-1) = vbarx(k-1) - vtot*csu(jrow)*dxu(i)
# ifdef partial_cell
     &                                              *dhu(i,k,j)
# else
     &                                              *dzt(k)
# endif
     &                       *basinmask(i,jrow,L)
              enddo
              vmsf(jrow,k-1,L) = vmsf(jrow,k,L) + vbarx(k-1)*rsverdrup
            enddo
c
c           when vmsf=0 at the ocean bottom, the plotting routines can`t
c           plot the contours correctly because of the "land" at the
c           bottom, so reset vmsf to non-zero at the bottom if the value
c           of vmsf at km-1 is non-zero
c
            if (vmsf(jrow,km-1,L) .ne. c0) vmsf(jrow,km,L) = 1.0e-12
c
# if defined gent_mcwilliams
            do k=1,km
              vbarx(k-1) = c0
            enddo
            vmsf_gm(jrow,km,L) = c0
#  ifdef bbl_ag
            do k=km-1,1,-1
#  else
            do k=km,1,-1
#  endif
              do i=istrt,iend
                vstargm = umask(i,k,j)
     &          *(adv_vntiso(i,k,j)+adv_vntiso(i+1,k,j))/(2.0*csu(jrow))
                vbarx(k-1) = vbarx(k-1) - vstargm*csu(jrow)*dxu(i)
#  ifdef partial_cell
     &                                            *dhu(i,k,j)
#  else
     &                                            *dzt(k)
#  endif
     &                     *basinmask(i,jrow,L)
              enddo
              vmsf_gm(jrow,k-1,L) = vmsf_gm(jrow,k,L) + 
     &                              vbarx(k-1)*rsverdrup
            enddo
            if (vmsf_gm(jrow,km-1,L) .ne. c0) vmsf_gm(jrow,km,L)=1.0e-12
# endif

          enddo
# ifdef timing
          call toc ('diagnostic', 'meridional overturning')
# endif
        endif
      enddo
#endif
c
#ifdef show_zonal_mean_of_sbc
c
c-----------------------------------------------------------------------
c     construct zonal mean of surface b.c. and related items
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (zmbcts .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'show zonal mean of s.b.c.')
# endif
        do j=js,je
          jrow = joff + j
          zmau(jrow) = c0
          zmat(jrow) = c0
          do i=istrt,iend
            zma1            = umask(i,1,j)*csu(jrow)*dxu(i)*dyu(jrow)
            zmau(jrow)      = zmau(jrow)    + zma1
            zmsmf(jrow,1)   = zmsmf(jrow,1) + zma1*smf(i,j,1)
            zmsmf(jrow,2)   = zmsmf(jrow,2) + zma1*smf(i,j,2)
            zmsm(jrow,1)    = zmsm(jrow,1)  + zma1*u(i,1,j,1,tau)
            zmsm(jrow,2)    = zmsm(jrow,2)  + zma1*u(i,1,j,2,tau)
            zma2            = tmask(i,1,j)*cst(jrow)*dxt(i)*dyt(jrow)
            zmat(jrow)      = zmat(jrow) + zma2
            do n=1,nt
              zmstf(jrow,n) = zmstf(jrow,n) + zma2*stf(i,j,n)
              zmst(jrow,n)  = zmst(jrow,n)  + zma2*t(i,1,j,n,tau)
            enddo
# if defined explicit_fresh_water_flux
            zmsff(jrow)     = zmsff(jrow,n) + zma2*sfft(i,jrow)
# else
            zmsff(jrow)     = c0 
# endif
          enddo
        enddo
# ifdef timing
        call toc ('diagnostic', 'show zonal mean of s.b.c.')
# endif
      endif
#endif
c
#ifdef matrix_sections
c
c-----------------------------------------------------------------------
c     print "tau" (not "tau+1") variables  at specified latitudes
c
c     author:    r.c.pacanowski   e-mail  rcp@gfdl.gov
c                a. rosati        e-mail  ar@gfdl.gov
c     based on previous code by M. Cox
c-----------------------------------------------------------------------
c
      if (prxzts .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'matrix sections')
# endif
        reltim = relyr
        do j=js,je
          jrow = joff + j
          do jlat=1,nlatpr
            jj = indp (prlat(jlat), yt, jmt)
            if (jj .eq. jrow .and. prlat(jlat) .le. yt(jmt))
     &        then
              iss = indp (prslon(jlat), xt, imt)
              iee = indp (prelon(jlat), xt, imt)
              ks = indp (prsdpt(jlat), zt, km)
              ke = indp (predpt(jlat), zt, km)
              fx = 1.0e-2
              if (jlat .eq. 1) write(stdout,9000)
              do m=1,nt
                scl = c1
                if (m .eq. 2) scl=1.e-3
                if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
                  write (stdout,9100) trname(m), itt, jrow 
     &,           yt(jrow), xt(iss), xt(iee), fx*zt(ks), fx*zt(ke), scl
                  call matrix (t(1,1,j,m,tau), imt,iss,iee, ks, ke, scl)
                endif
                if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
c
                  file_name = file_stamp ('sections', xstamp,'.dta')
	          call getunit (io, trim(file_name)
     &,                        'unformatted sequential append ieee')
c
                  write (stdout,*) ' =>Zonal ',trname(m), ' slice: lat='
     &,           yt(jrow), ' written unformatted to ',trim(file_name)
     &,           ' on ts=', itt, stamp
                  write (stdout,'(///)')
                  iotext = ' read (ioprxz) imt, km, m, nt, reltim'
                  write (io) stamp, iotext, expnam
                  write (io) imt, km, m, nt, reltim
                  write(iotext,'(a10,i4,a4,i2)') ' for jrow=',jrow
     &,           '  m=',m
                  iotext(18:)=':read(ioprxz)((t(i,k,m),i=1,imt),k=1,km)'
                  write (io) stamp, iotext, expnam
                  call wrufio (io, t(1,1,j,m,tau), imt*km)
c
                  call relunit (io)
c
                endif
              enddo
              scl = 1.e-3
              if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
	        i1 = max(iss,2)
		i2 = min(iee,imtm1)
                write (stdout,9100)  'adv_vbu ', itt, jrow
     &,         yt(jrow), xt(i1), xt(i2), fx*zw(ks), fx*zw(ke), scl
                call matrix (adv_vbu(1,1,j), imt, i1, i2, ks, ke, scl)
                write (stdout,9100)  'adv_vbt ', itt, jrow
     &,         yt(jrow), xt(i1), xt(i2), fx*zw(ks), fx*zw(ke), scl
                call matrix (adv_vbt(1,1,j), imt, i1, i2, ks, ke, scl)
              endif
              if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
c
                file_name = file_stamp ('sections', xstamp,'.dta')
	        call getunit (io, trim(file_name)
     &,                      'unformatted sequential append ieee')
c
                write (stdout,*) ' => Zonal adv_vbu slice: lat='
     &,         yt(jrow), ' written unformatted to ',trim(file_name)
     &,         ' on ts=', itt, stamp
                write (stdout,'(///)')
                iotext = ' read (ioprxz) imt, km, reltim'
                write (io) stamp, iotext, expnam
                write (io) imt, km, reltim
                write(iotext,'(a10,i4)') ' for jrow=',jrow
                iotext(12:)=
     &          ': read(ioprxz)((adv_vbu(i,k),i=1,imt),k=0,km)'
                write (io) stamp, iotext, expnam
                call wrufio (io, adv_vbu(1,0,j), imt*(km+1))
c
                write (stdout,*) ' => Zonal adv_vbt slice: lat='
     &,         yt(jrow), ' written unformatted to ',trim(file_name)
     &,         ' on ts=', itt, stamp
                write (stdout,'(///)')
                iotext = ' read (ioprxz) imt, km, reltim'
                write (io) stamp, iotext, expnam
                write (io) imt, km, reltim
                write(iotext,'(a10,i4)') ' for jrow=',jrow
                iotext(12:)=
     &          ': read(ioprxz)((adv_vbt(i,k),i=1,imt),k=1,km)'
                write (io) stamp, iotext, expnam
                call wrufio (io, adv_vbt(1,0,j), imt*(km+1))
c
                call relunit (io)
c
              endif
c
              scl = c1
              if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
                write (stdout,9100) 'u velocity', itt
     &,        jrow, yt(jrow), xt(iss),xt(iee),fx*zt(ks), fx*zt(ke), scl
                call matrix (u(1,1,j,1,tau), imt, iss, iee, ks, ke, scl)
              endif
              if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
c
                file_name = file_stamp ('sections', xstamp, '.dta')
	        call getunit (io, trim(file_name)
     &,                      'unformatted sequential append ieee')
c
                write (stdout,*) ' => Zonal u velocity slice: lat='
     &,         yt(jrow), ' written unformatted to ',trim(file_name)
     &,         ' on ts=', itt, stamp
                write (stdout,'(///)')
                iotext = ' read (ioprxz) imt, km, reltim'
                write (io) stamp, iotext, expnam
                write (io) imt, km, reltim
                write(iotext,'(a10,i4)') ' for jrow=',jrow
                iotext(12:) = ': read (ioprxz)((u(i,k),i=1,imt),k=1,km)'
                write (io) stamp, iotext, expnam
                call wrufio (io, u(1,1,j,1,tau), imt*km)
c
                call relunit (io)
c
              endif
c
              scl = c1
              if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
                write (stdout,9100) 'v velocity', itt, jrow
     &,         yt(jrow), xt(iss), xt(iee), fx*zt(ks), fx*zt(ke), scl
                call matrix (u(1,1,j,2,tau), imt, iss, iee, ks, ke, scl)
              endif
              if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
c
                file_name = file_stamp ('sections', xstamp, '.dta')
	        call getunit (io, trim(file_name)
     &,                      'unformatted sequential append ieee')
c
                write (stdout,*) ' => Meridional v velocity slice: lat='
     &,         yt(jrow),' written unformatted to ',trim(file_name)
     &,          ' on ts=', itt, stamp
                write (stdout,'(///)')
                iotext = ' read (ioprxz) imt, km, reltim'
                write (io) stamp, iotext, expnam
                write (io) imt, km, reltim
                write(iotext,'(a10,i4)') ' for jrow=',jrow
                iotext(12:) = ': read (ioprxz)((v(i,k),i=1,imt),k=1,km)'
                write (io) stamp, iotext, expnam
                call wrufio (io, u(1,1,j,2,tau), imt*km)
c
                call relunit (io)
c
              endif
            endif
          enddo
        enddo
# ifdef timing
        call toc ('diagnostic', 'matrix sections')
# endif
      endif
9000  format(/' Zonal section printouts at specified latitudes:'/)
9100  format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     &,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     &,', scaling=',1pg10.3)
#endif
c
#if defined save_mixing_coeff
c
c-----------------------------------------------------------------------
c     diagnostic: save estimated mixing coefficients on east, north, and
c                 bottom face of T and U cells
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (cmixts .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'mixing coeff')
# endif
        call cmix_out (joff, js, je, .true.)
# ifdef timing
        call toc ('diagnostic', 'mixing coeff')
# endif
      endif
#endif
c
#if defined cross_flow && defined isoneutralmix
c
c-----------------------------------------------------------------------
c     diagnostic: compute diapycnal and isopycnal components of flow
c
c     author:    r.c.pacanowski   e-mail  rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (crossts .and. .not. euler2) then
	call cross_out (joff, js, je, .true.)
      endif
#endif
c
#ifdef save_pressure 
c
c-----------------------------------------------------------------------
c     diagnostic: write instantaneous pressure data
c     author:    s.m. griffies e-mail  smg@gfdl.gov
c-----------------------------------------------------------------------
c
      if (pressts .and. .not. euler2) then
	call pressure_out (joff, js, je, .true.)
      endif
#endif
c
#ifdef save_density_terms
c
c-----------------------------------------------------------------------
c     diagnostic: Save instantaneous terms in local 
c     potential density equation.
c     author:  S. M. Griffies  e-mail  smg@gfdl.gov
c-----------------------------------------------------------------------
c
      if (rhots .and. .not. euler2) then
	call density_out (joff, js, je, is, ie, .true.)
      endif
#endif
c
#if defined save_ncdiff
c
c-----------------------------------------------------------------------
c     diagnostic: Save fields related to nonconstant diffusivities.
c     Should only map fields on diffts time step, since otherwise,
c     certain of the fields are not known.  
c     author:  S. M. Griffies   e-mail  smg@gfdl.gov
c-----------------------------------------------------------------------
c
      if (diffsnapts .and. .not. euler2 .and. diffts) then
	call ncdiff_out (joff, js, je, .true.)
      endif
#endif
c
      endif
      return
      end

