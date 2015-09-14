#if defined velocity_horz_mix_var || defined velocity_horz_mix_smag
      subroutine smagnlc (joff, js, je, is, ie)
c
      use memory_window
# include "param.h"
# include "coord.h"
# include "grdvar.h"
# include "hmixc.h"
# include "mw.h"
# include "switch.h"
# if defined stability_tests
#  include "scalar.h"
# endif
      dimension f_smag(imt)
c
c=======================================================================
c     Compute Smagorinsky viscosity or 
c     specify a spatially variable viscosity.
c
c     Viscosities can be used for either 
c     laplacian or biharmonic operators.
c     
c     Also can compute horizontal diffusivites for tracers.
c
c     input:
c
c      joff   = offset between row j in the MW and latitude jrow on disk
c      js     = starting row for calculations
c      je     = ending row for calculations
c      is     = starting longitude index for calculations
c      ie     = ending longitude index for calculations
c
c     output:
c
c      am_varh  = viscosity at northeast corner of U cell
c      visc_ceu = viscosity at eastern face of U cell
c      visc_cnu = viscosity at northern face of U cell
c      diff_ceu = diffusivity at eastern face of T cell
c      diff_cnu = diffusivity at northern face of T cell
c
c      r.c. pacanowski   email: rcp@gfdl.gov
c      s.m. griffies     email: smg@gfdl.gov 
c=======================================================================
c
# ifdef timing
      call tic ('hmixc', 'smag')
# endif
c
c-----------------------------------------------------------------------
c     set bounds
c-----------------------------------------------------------------------
c
      istrt  = max(2,is)
      iend   = min(imt-1,ie)
      jstrt  = max(js-1,jbuf)
      jend   = je-1
      jsc    = max(js-1,jbuf+1)
c
# if defined velocity_horz_mix_smag
c-----------------------------------------------------------------------
c     compute rates of strain on northeast corner of U cell. The 
c     resulting viscosity am_varh(i,k,j) is defined on northeast 
c     corner of the U(i,k,j) cell. 
c
c     am_varh will be chosen as the max Smag value found in a column. 
c     Otherwise, the more quiescent interior can have unresolved 
c     boundary currents.   
c-----------------------------------------------------------------------
c
      pi  = 4.0*atan(1.0)
#  if defined velocity_horz_laplacian
      fx2    = (k_smag/pi)**2
      npower = 2
#  endif
#  if defined velocity_horz_biharmonic
      fx2    = p125*(k_smag/pi)**2
      npower = 4 
#  endif
c
      do j=jstrt,jend
        jrow = j + joff
        do i=istrt,iend
          f_smag(i) =  fx2*(min(cst(jrow+1)*dxt(i+1),dyt(jrow+1)))
     &                 **npower
        enddo
        f1   = cstr(jrow+1)
        f2   = dyt2r(jrow+1)*cst(jrow+1)
        do i=istrt,iend
          visc_vert_max = c0
          do k=1,km
            tension = ((u(i+1,k,j,1,taum1)+u(i+1,k,j+1,1,taum1))
     &               - (u(i  ,k,j,1,taum1)+u(i  ,k,j+1,1,taum1))
     &                )*f1*dxt2r(i+1)
     &      -(csur(jrow+1)*(u(i,k,j+1,2,taum1)+u(i+1,k,j+1,2,taum1)) 
     &       -csur(jrow  )*(u(i,k,j  ,2,taum1)+u(i+1,k,j  ,2,taum1))
     &       )*f2
            shear   = ((u(i+1,k,j,2,taum1)+u(i+1,k,j+1,2,taum1))
     &               - (u(i  ,k,j,2,taum1)+u(i  ,k,j+1,2,taum1))
     &                )*f1*dxt2r(i+1)
     &      +(csur(jrow+1)*(u(i,k,j+1,1,taum1)+u(i+1,k,j+1,1,taum1)) 
     &       -csur(jrow  )*(u(i,k,j  ,1,taum1)+u(i+1,k,j  ,1,taum1))
     &       )*f2
            am_varh(i,k,j) = f_smag(i)*sqrt(tension**2 + shear**2)
# if defined velocity_horz_mix_smag_vert_max
            if(am_varh(i,k,j) .gt. visc_vert_max)
     &          visc_vert_max = am_varh(i,k,j)
# endif
          enddo
# if defined velocity_horz_mix_smag_vert_max
          am_varh(i,:,j) = visc_vert_max
# endif
        enddo
        call setbcx (am_varh(1,1,j), imt, km)
      enddo
# endif
c
# if !defined velocity_horz_mix_smag
c
c     Variable horizontal viscosity
c     USER INPUT (specify variation in am_varh here)
c     Default values are constants from namelist.  
c
c     examples:
c
c     A. taper to zero as cosine
c     laplacian:      
c     am_varh(i,k,j) = am*csu(jrow)
c     biharmonic:
c     am_varh(i,k,j) = abs(ambi*csu(jrow)**3)
c
c     B. taper to am/2 as cosine      
c     laplacian:
c     am_varh(i,k,j) = am*(csu(jrow)+c1)
c     biharmonic:
c     am_varh(i,k,j) = abs(ambi*(csu(jrow)+c1)**3)
c
c     C. linearly taper to small value when poleward of critical lat.
c            amtaper_min=am/10.
c            amtaper_start_lat=60.
c            amtaper_end_lat=87.
c            amtaper = am
c            if (abs(yu(jrow)) .gt. amtaper_start_lat) then
c               amtaper = am -
c     &          (am-amtaper_min)*(amtaper_start_lat-abs(yu(jrow)))**2
c     &          /(amtaper_start_lat-amtaper_end_lat)**2
c            endif
c            if (amtaper .lt. amtaper_min) amtaper = amtaper_min
c            am_varh(i,k,j) = amtaper
      
      do j=jstrt,jend
        jrow = j + joff
        do k=1,km
          do i=istrt,iend
#  if defined velocity_horz_laplacian
            am_varh(i,k,j) = am
#   if defined am_cosine
            am_varh(i,k,j) = am*csu(jrow)
#   endif
#  endif
#  if defined velocity_horz_biharmonic
            am_varh(i,k,j) = abs(ambi)
#   if defined am_cosine
            am_varh(i,k,j) = abs(ambi*csu(jrow)**3)
#   endif
#  endif
          enddo
        enddo
        call setbcx (am_varh(1,1,j), imt, km)
      enddo
# endif
c
# if defined stability_tests
c
c     Check to see if the specified viscosity satisfies 
c     linear stability constraint for the diffusion equation.
c 
      if (stabts) then 
        ncount_stab_large_vis = 0
        do j=jstrt,jend
          jrow = j + joff
	  if (abs(yu(jrow)) .ge. 89.9) cycle           
          do k=1,km
            do i=istrt,iend
              dsmin = min(dxu(i)*csu(jrow),dyu(jrow))
              amtest = am_varh(i,k,j)
#  if defined velocity_horz_laplacian
              am_crit = p125*dsmin**2/(dtuv + epsln)
#  endif      
#  if defined velocity_horz_biharmonic
              am_crit = p125*p125*dsmin**4/(dtuv + epsln)
#  endif      
#  if defined velocity_horz_mix_smag
     &             +visc_c_back
#  endif      
              if (amtest .gt. am_crit
     &          .and. ncount_stab_large_vis .le. 50) then
                ncount_stab_large_vis = ncount_stab_large_vis+1
                write (stdout,'(a,3(i4,a),es14.7,a,es14.7)') 
     &          'Warning: Viscosity too large at (i,k,jrow) = ('
     &,          i,',',k,',',jrow,'), am_crit = ',am_crit
     &,         ', am = ',amtest
              endif
              if (ncount_stab_large_vis .eq. 50) then
                write(stdout,*) 'Stopped counting viscosity violations'
                write(stdout,*) 'after 50 points.'
              endif
            enddo
          enddo
        enddo  
      endif
# endif
c
c-----------------------------------------------------------------------
c     mixing coeffs on the east and north face of U and T cells
c     varh_cu is needed at the center of the U cell for the metric term
c-----------------------------------------------------------------------
c
# ifdef tracer_horz_mix_smag
      prandtlr = c1/prandtl
# endif
c
c     east face  
c
      do j=jsc,jend
        jm1 = max(j-1,jbuf)
	jrow = joff + j  
        do k=1,km
          do i=istrt,iend
            im1 = max(i-1,istrt)  
            visc_ceu(i,k,j) = p5*(am_varh(i,k,j)+am_varh(i,k,jm1))
# ifdef velocity_horz_mix_smag
     &                        + visc_c_back
# endif
# if defined tracer_horz_mix_smag
            diff_cet(i,k,j) = prandtlr*
     &                        p5*(am_varh(i,k,jm1)+am_varh(im1,k,jm1)) 
     &                        + diff_c_back
# endif
c
# if defined tracer_horz_mix_var && !defined tracer_horz_mix_smag
c
c           variable horizontal diffusion coefficient
c           USER INPUT (specify variation in diff_cet here)
c           Default values are constants from namelist
c           See examples for variable viscosity for suggestions.            
c
#  if defined tracer_horz_laplacian
            diff_cet(i,k,j) = ah
#  endif
#  if defined tracer_horz_biharmonic
            diff_cet(i,k,j) = abs(ahbi)
#  endif
# endif
c
            varh_cu(i,k,j)  = p25*(am_varh(i-1,k,j) + am_varh(i,k,j)
     &                            +am_varh(i,k,jm1)+am_varh(i-1,k,jm1)) 
          enddo
        enddo
	call setbcx (visc_ceu(1,1,j), imt, km)
	call setbcx (varh_cu(1,1,j),  imt, km)
# if defined tracer_horz_mix_var 
	call setbcx (diff_cet(1,1,j), imt, km)
# endif
      enddo
c
c     north face  
c
      do j=jstrt,jend
        jm1 = max(j-1,jbuf)
	jrow = joff + j  
        do k=1,km
          do i=istrt,iend
            im1 = max(i-1,istrt)  
            visc_cnu(i,k,j) = p5*(am_varh(im1,k,j) + am_varh(i,k,j))
# if defined velocity_horz_mix_smag
     &                        + visc_c_back
# endif
# if defined tracer_horz_mix_smag 
	    diff_cnt(i,k,j) = prandtlr*
     &                        p5*(am_varh(im1,k,j) + am_varh(im1,k,jm1))
     &                        + diff_c_back
# endif
c
# if defined tracer_horz_mix_var && !defined tracer_horz_mix_smag
c
c           variable horizontal diffusion coefficient
c           USER INPUT (specify variation in diff_cnt here)
c
#  if defined tracer_horz_laplacian
            diff_cnt(i,k,j) = ah
#  endif
#  if defined tracer_horz_biharmonic
            diff_cnt(i,k,j) = abs(ahbi)
#  endif
# endif
          enddo
        enddo
	call setbcx (visc_cnu(1,1,j), imt, km)
# if defined tracer_horz_mix_var 
	call setbcx (diff_cnt(1,1,j), imt, km)
# endif
      enddo
c
# ifdef matrix_sections
      if (prxzts .and. eots) then
        call diagnl (joff, jsc, jend)
      endif
# endif
# ifdef trace_indices
      write (stdout,'(4x,6(a,i4))')
     & "=> In smagnlc: js=",js," je=",je," joff=",joff," jstrt=",jstrt
     &," jend=",jend," jsc=",jsc
# endif
# ifdef timing
      call toc ('hmixc', 'smag')
# endif
      return
      end




# ifdef matrix_sections
      subroutine diagnl (joff, js, je)
      character*128 file_name, file_stamp
#  include "param.h"
#  include "coord.h"
#  include "cprnts.h"
#  include "hmixc.h"
#  include "iounit.h"
#  include "switch.h"
#  include "tmngr.h"
c
c=======================================================================
c     Show some diagnostics
c
c     input:
c
c      joff   = offset between row j in the MW and latitude jrow on disk
c      js     = starting row for calculations
c      je     = ending row for calculations
c=======================================================================
c
#  ifdef timing
      call tic ('diagnostic', 'matrix sections')
#  endif
      do j=js,je
	jrow = j + joff
        reltim = relyr
        do jlat=1,nlatpr
          jj = indp (prlat(jlat), yt, jmt)
          if (jj .eq. jrow .and. prlat(jlat) .le. yt(jmt)) then
            is = indp (prslon(jlat), xt, imt)
            ie = indp (prelon(jlat), xt, imt)
            ks = indp (prsdpt(jlat), zt, km)
            ke = indp (predpt(jlat), zt, km)
            fx = 1.0e-2
#  if defined tracer_horz_mix_var
c
c           write out the diffusion coeffs for tracers
c
            scl = 1.e7
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'diff_cet', itt, jrow 
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (diff_cet(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
              file_name = file_stamp ('sections', xstamp, '.dta')
	      call getunit (io, trim(file_name)
     &,                    'unformatted sequential append ieee')
              write (stdout,*) ' => diff_cet ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((diff_cet(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, diff_cet(1,1,j), imt*km)
              call relunit (io)
            endif
c
            scl = 1.e7
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'diff_cnt', itt, jrow 
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (diff_cnt(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
              file_name = file_stamp ('sections', xstamp, '.dta')
	      call getunit (io, trim(file_name)
     &,                    'unformatted sequential append ieee')
              write (stdout,*) ' => diff_cnt ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((diff_cnt(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, diff_cnt(1,1,j), imt*km)
              call relunit (io)
            endif
#  endif
#  if defined velocity_horz_mix_var
c
c           write out the viscosity coeffs for momentum
c
            scl = 1.e7
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'visc_ceu', itt, jrow 
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (visc_ceu(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
              file_name = file_stamp ('sections', xstamp, '.dta')
	      call getunit (io, trim(file_name)
     &,                    'unformatted sequential append ieee')
              write (stdout,*) ' => visc_ceu ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((visc_ceu(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, visc_ceu(1,1,j), imt*km)
              call relunit (io)
            endif
c
            scl = 1.e7
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'visc_cnu', itt, jrow 
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (visc_cnu(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
              file_name = file_stamp ('sections', xstamp, '.dta')
	      call getunit (io, trim(file_name)
     &,                    'unformatted sequential append ieee')
              write (stdout,*) ' => visc_cnu ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((visc_cnu(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, visc_cnu(1,1,j), imt*km)
              call relunit (io)
            endif
#  endif
          endif
        enddo
      enddo
#  ifdef timing
      call toc ('diagnostic', 'matrix sections')
#  endif
      return
9100    format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     &,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     &,', scaling=',1pg10.3)
      end
# endif

#else
      subroutine smagnl
      return
      end
#endif


      
      
