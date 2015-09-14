      subroutine ncdiffi (error)
c
      use arrays_2d
      use memory_window
      logical error
c
#include "param.h"
#ifdef nonconst_diffusivity
c
# include "coord.h"
# include "iounit.h"
# include "isopyc.h"
# include "tmngr.h"
# include "scalar.h"
# include "switch.h"
# include "taskrows.h"

c=======================================================================
c     Initialization for nonconstant diffusivity schemes.
c     There are two schemes available, one or both of which 
c     can be enabled.  If both are on, then will use the minimum
c     diffusivity between the two.
c
c     hl_diffusivity is based on 
c     Held and Larichev, 1996: A scaling theory for horizontally 
c     homogeneous baroclinically unstable flow on a beta plane,
c     Journal of Atmospheric Sciences, vol 53, pages 946-952.
c
c     vmhs_diffusivity is based on  
c     Visbeck, Marshall, Haine, and Spall, 1997:
c     Specification of eddy transfer coefficients in coarse resolution 
c     ocean circulation models. Journal of Physical Oceanography, 
c     vol 27, pages 381-402.
c
c     input:
c       error  = logical to signal problems
c
c     output:
c       depthbot         = Bottom depth for performing vertical average.
c                          Can be set in namelist.
c       depthtop         = Top depth for performing vertical average.  
c                          Can be set in namelist.
c       kncbot           = vertical level corresponding to depthbot.
c       knctop           = vertical level corresponding to depthtop
c       diffmin          = Minimum diffusivity (cm2/sec).
c                          Can be set in namelist.
c       diffmax          = Maximum diffusivity (cm2/sec).
c                          Can be set in namelist.
c       vmhs_alpha       = Dimensionless parameter for setting the 
c                          overall scale of the vmhs_diffusivity 
c                          diffusivity.
c                          Can be set in namelist.
c       vmhs_rate2_limit = minimum squared growth rate (units sec^-2)
c                          for defining a baroclinic zone in 
c                          vmhs_diffusivity.
c                          Can be set in namelist.
c       ijvmhs_maxlen    = (1/2) the maximum number of grid points 
c                          to be accepted for the width of the
c                          baroclinic zone in vmhs_diffusivity
c                          Can be set in namelist.
c       diffint          = Number of time steps (real) between which the 
c                          diffusivity is computed. 
c                          Can be set in namelist.
c       diffsnapint      = Period (days) between writing data file
c                          with various nonconstant diffusivity 
c                          related fields.  
c                          Can be set in namelist.
c       growth2max       = maximum squared growth rate allowed for use in 
c                          computing diffusivities.  
c                          Default is (1/day)^2
c       betaeff          = effective beta for use with hl_diffusivity
c                          This field is computed with an FIR filtered
c                          topography inside of topog.F.  Such filtering
c                          will smooth vertical walls.                                                      
c       numfltrgrth      = number of FIR filter passes applied to the 
c                          Eady growth rate field growth2a. 
c                          default=0. 
c       numfltrtopog     = number of FIR filter passes applied to the 
c                          bottom topog for computing betaeff. 
c                          This parameter is set inside topog.F.
c                          default=1. 
c
c     authors: SM Griffies:   smg@gfdl.gov
c              RC Pacanowski: rcp@gfdl.gov
c=======================================================================
c
c
      namelist /ncdiff/ depthbot,depthtop,diffmin,diffmax
     &,vmhs_alpha,vmhs_rate2_limit,ijvmhs_maxlen
     &,diffint, diffsnapint
c
      write (stdout,'(/,20x,a,/)')
     & 'N O N C O N S T A N T   D I F F U S I V I T Y   I N I T '
c
c-----------------------------------------------------------------------
c     initialize variables (all diffusivity units are cm**2/sec.)
c-----------------------------------------------------------------------
c
      depthbot         = 1000.e2
      depthtop         = 100.e2
      diffmin          = 1.e6
      diffmax          = 2.e7
      vmhs_alpha       = 0.015 
      vmhs_rate2_limit = (1.4e-6)**2
      ijvmhs_maxlen    = 10
      diffint          = 5.0
      idiffint         = nint(diffint)
      diffsnapint      = dtts*diffint/86400.0 
c
c-----------------------------------------------------------------------
c     provide for namelist over-ride of the above settings
c-----------------------------------------------------------------------
c
      call getunit (io, 'namelist.ncdiff'
     &     , 'formatted sequential rewind')
      read (io,ncdiff,end=100)
100   continue
c
c-----------------------------------------------------------------------
c     Non-namelist or derived quantities
c-----------------------------------------------------------------------
c
      numfltrgrth      = 0
      knctop           = indp(depthtop,zt,km)
      kncbot           = max (knctop,min(indp(depthbot,zt,km),km-1))
      ddepthtop        = zt(knctop) 
      ddepthbot        = zt(kncbot)
      growth2max       = (1/21600.0)**2
c
c-----------------------------------------------------------------------
c     write some details
c-----------------------------------------------------------------------
c
      write (stdout,ncdiff)
      call relunit (io)
      call getunit (iodoc, 'document.dta','formatted sequential append')
      write (iodoc, ncdiff)
      call relunit (iodoc)
c
      write(stdout,*)' '
      write(stdout,*)'------Parameters for nonconstant diffusivity-----'
      write(stdout,*)' Lower model depth to compute time scale(cm) = '
     &     ,ddepthbot
      write(stdout,*)' This is the model depth closest to(cm)      ='
     &     ,depthbot
      write(stdout,*)' It is at the model level                    ='
     &     ,kncbot
      write(stdout,*)' Upper model depth to compute time scale(cm) = '
     &     ,ddepthtop
      write(stdout,*)' This is the model depth closest to(cm)      ='
     &     ,depthtop
      write(stdout,*)' It is at the model level                    ='
     &     ,knctop
      write(stdout,*)' Max diffusivity allowed (cm2/sec)           = '
     &     ,diffmax
      write(stdout,*)' Min diffusivity allowed (cm2/sec)           = '
     &     ,diffmin
      write(stdout,*)' Max squared growth rate(sec^-2)             ='
     &     ,growth2max
# if defined vmhs_diffusivity
      write(stdout,*)' VMHS alpha which scales diffusivity         = '
     &     ,vmhs_alpha
      write(stdout,*)' VMHS min squared growth rate b/c zone (s^-2)= '
     &     ,vmhs_rate2_limit
      write(stdout,*)' VMHS max number grid points in 1/2 b/c zone = '
     &     ,ijvmhs_maxlen
# endif
      write(stdout,*)' Diffusivity updated every ',idiffint
     &,' time steps. Be aware of this frequency when comparing'
     &,' experiments using different tracer time steps dtts. '
# if defined save_ncdiff
      write(stdout,*)' Period on which ncdiff snapshot written     = '
     &     ,diffsnapint
# endif
      write(stdout,'(a)')
     & '--------------------------------------------------------------'
c
c-----------------------------------------------------------------------
c     check for problems
c-----------------------------------------------------------------------
c
# if !defined isoneutralmix
        write (stdout,'(/,(1x,a))')
     & '==> Error:"isoneutralmix" must be enabled since '
     &,'           "nonconst_diffusivity" predicts mixing'
     &,'            coeffs used by "isoneutralmix". '
        error = .true.
# endif
c
c-----------------------------------------------------------------------
c     Gather the depth dependent pieces needed for the vertical average
c     used to compute the growth rate of an Eady wave.
c     The p125 = p25*p5 factor arises from the following:
c       p25 is from the four triads.
c       p5  is from the two groups of Sx and Sy computed; one 
c       on east/north faces, and the other on the bottom face
c-----------------------------------------------------------------------
c
      do k=1,km
        fznoncst(k) = c0
        if(k .ge. knctop .and. k .le. kncbot) then
          fznoncst(k) = p125*(grav/(rho0*(ddepthbot-ddepthtop)))*dzt(k)
        endif
      enddo
c
c-----------------------------------------------------------------------
c     Initialize 
c-----------------------------------------------------------------------
c
      do j=1,jmw
        do i=1,imt
	  hldiff(i,j)     = diffmax
	  vmhsdiff(i,j)   = diffmax
          hllength(i,j)   = c0
          vmhslength(i,j) = c0
         enddo
      enddo
      if (init) then
        do j=jstask,jetask
          do i=1,imt
	    growth2a(i,j) = c0
          enddo
        enddo
      endif
# ifdef nonconst_diffusivity_initial
        do j=jstask,jetask
          do i=1,imt
            growth2a(i,j) = c0
            agm(i,j)      = athkdf
          enddo
        enddo
# endif
      do j=jstask,jetask
        do i=1,imt
          growth2b(i,j)   = c0
        enddo
      enddo
c
c
      return
      end


      subroutine hldiffusivity (joff, js, je, is, ie)
c
c=======================================================================
c     Compute mixing coefficients based on Held and Larichev (1996).
c
c     inputs:
c
c      joff   = offset between rows in the MW and latitude rows
c      js     = starting row for loading variables to calculate
c               coefficients. start at jstrt=max(js-1,1+jbuf)
c      je     = ending row for loading variables to calculate
c               coefficients. calculations end at je-1
c      is     = starting index for calculating coefficients in the 
c               longitude direction
c      ie     = ending index for calculating coefficients in the 
c               longitude direction
c
c     output:
c
c      hllength = held/larichev length scale (cm)
c      hldiff   = held/larichev diffusivity (cm2/sec)
c
c     author:  S. M. Griffies    e-mail=> smg@gfdl.gov
c              R. C. Pacanowski  e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use arrays_2d
      use memory_window
# include "param.h"
# include "coord.h"
# include "grdvar.h"
# include "isopyc.h"
# include "mw.h"
# include "scalar.h"
c
# ifdef timing
      call tic ('hldiff', 'hl_diffusivity')
# endif
c     
      pi   = c4*atan(c1)
      crad = pi/180.0
      do j=js,je
        jrow = j + joff
        do i = is,ie
          if(kmt(i,jrow) .ge. 1) then
            if(zt(kmt(i,jrow)) .gt. depthtop) then
              hllength(i,j) = tmask(i,1,j)*
     &             sqrt(growth2a(i,jrow))
     &             /(betaeff(i,jrow)+epsln)
              hldiff(i,j) = sqrt(growth2a(i,jrow))*hllength(i,j)**2
              if(hldiff(i,j) 
     &             .lt. diffmin) hldiff(i,j) = tmask(i,1,j)*diffmin
              if(hldiff(i,j)
     &             .gt. diffmax) hldiff(i,j) = tmask(i,1,j)*diffmax
            else
              hldiff(i,j)   = tmask(i,1,j)*athkdf
              hllength(i,j) = c0
            endif
          else
            hldiff(i,j)   = c0
            hllength(i,j) = c0
          endif
        enddo
        call setbcx(hllength(1,j),imt,1)
        call setbcx(hldiff(1,j),imt,1)
      enddo
c     
# ifdef timing
      call toc ('hldiffusivity', 'hl_diffusivity')
# endif
      return
      end


      subroutine vmhsdiffusivity (joff, js, je, is, ie)
c
c=======================================================================
c     Compute mixing coefficients based on Visbeck etal (1997).
c
c     inputs:
c
c      joff   = offset between rows in the MW and latitude rows
c      js     = starting row for loading variables to calculate
c               coefficients. start at jstrt=max(js-1,1+jbuf)
c      je     = ending row for loading variables to calculate
c               coefficients. calculations end at je-1
c      is     = starting index for calculating coefficients in the 
c               longitude direction
c      ie     = ending index for calculating coefficients in the 
c               longitude direction
c      tmask  = land/sea mask on "t" grid (land=0.0, sea=1.0)
c
c     output:
c
c      vmhslength = Visbeck etal length scale (cm)
c      vmhsdiff   = Visbeck etal diffusivity (cm2/sec)
c
c     author:  S. M. Griffies    e-mail=> smg@gfdl.gov
c              R. C. Pacanowski  e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use arrays_2d
      use memory_window
# include "param.h"
# include "coord.h"
# include "grdvar.h"
# include "isopyc.h"
# include "mw.h"
# include "scalar.h"
c
# ifdef timing
      call tic ('vmhsdiff', 'vmhs_diffusivity')
# endif
c
c     Determine the VMHS length scale according to the width
c     of the zone of significant baroclinicity.  
c
      do j = js,je
        jrow = j + joff 
        do i = is,ie
c         
          if(kmt(i,jrow) .ge. 1) then
c           
            if(zt(max(1,kmt(i,jrow))) .gt. depthtop) then
c             
              if(growth2a(i,jrow) .le. vmhs_rate2_limit) then
                vmhslength(i,j) = max(cst(jrow)*dxt(i),dyt(jrow))
              else
c               
c               scan north/south/east/west to get length scale
c               
                distn  = c0
                do jq = jrow,min(jmt-1,jrow+ijvmhs_maxlen),1
                  if(growth2a(i,jq) .ge. vmhs_rate2_limit) then
                    distn = distn + dyt(jq) 
                  else
                    exit
                  endif
                enddo
                dists  = c0
                do jq = jrow,max(2,jrow-ijvmhs_maxlen),-1
                  if(growth2a(i,jq) .ge. vmhs_rate2_limit) then
                    dists  = dists + dyt(jq) 
                  else
                    exit
                  endif 
                enddo
                diste  = c0
                do ip = i,min(imt-1,i+ijvmhs_maxlen),1
                  if(growth2a(ip,jrow) .ge. vmhs_rate2_limit) then
                    diste = diste + dxt(ip)*cst(jrow) 
                  else
                    exit
                  endif
                enddo
                distw  = c0
                do ip = i,max(2,i-ijvmhs_maxlen),-1
                  if(growth2a(ip,jrow) .ge. vmhs_rate2_limit) then
                    distw  = distw + dxt(ip)*cst(jrow) 
                  else
                    exit
                  endif 
                enddo
c               
c               Account for double counting the central t-cell.
c               
                distns = distn + dists - dyt(jrow)
                distew = diste + distw - dxt(i)*cst(jrow) 
c               
c               Define width of the baroclinic zone and get length scale
c               
                if(distns .le. distew) then
                  zone = distns
                  frac = min(distn,dists)/max(distn,dists)
                else
                  zone = distew
                  frac = min(diste,distw)/max(diste,distw)
                endif   
                vmhslength(i,j) = frac*zone
c               
              endif
c             
c             Compute the diffusivity 
c             
              vmhsdiff(i,j) = tmask(i,1,j)*vmhs_alpha*vmhslength(i,j)**2
     &             *sqrt(growth2a(i,jrow))
              if(vmhsdiff(i,j)
     &             .lt.diffmin) vmhsdiff(i,j)=tmask(i,1,j)*diffmin
              if(vmhsdiff(i,j)
     &             .gt.diffmax) vmhsdiff(i,j)=tmask(i,1,j)*diffmax
c             
c             Default in shallow waters to namelist parameter
c             
            else
              vmhslength(i,j) = max(cst(jrow)*dxt(i),dyt(jrow))
              vmhsdiff(i,j)   = tmask(i,1,j)*athkdf  
            endif
c           
c         Set length and diffusivity to zero over land
c
          else
            vmhslength(i,j) =  c0
            vmhsdiff(i,j)   =  c0
          endif
c         
        enddo
        call setbcx(vmhslength(1,j),imt,1)
        call setbcx(vmhsdiff(1,j),imt,1)
      enddo           
c
# ifdef timing
      call toc ('vmhsdiffusivity', 'vmhs_diffusivity')
# endif
      return
      end

      subroutine ncdiffusivity (joff, js, je, is, ie)
c
c=======================================================================
c     Accumulate mixing coefficients
c
c     inputs:
c
c      joff   = offset between rows in the MW and latitude rows
c      js     = starting row for loading variables to calculate
c               coefficients. start at jstrt=max(js-1,1+jbuf)
c      je     = ending row for loading variables to calculate
c               coefficients. calculations end at je-1
c      is     = starting index for calculating coefficients in the 
c               longitude direction
c      ie     = ending index for calculating coefficients in the 
c               longitude direction
c
c     output:
c
c      agm    = GM diffusivity (cm2/sec)
c
c     author:  S. M. Griffies    e-mail=> smg@gfdl.gov
c              R. C. Pacanowski  e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use memory_window
# include "size.h"
# include "isopyc.h"
c
# ifdef timing
      call tic ('ncdiffusivity', 'nc_diffusivity')
# endif
c
      do j=js,je
        jrow = j+joff
        do i=is,ie
          agm(i,jrow) = min(hldiff(i,j),vmhsdiff(i,j))
        enddo
        call setbcx(agm(1,jrow),imt,1)
      enddo
c
# ifdef timing
      call toc ('ncdiffusivity', 'nonconst_diffusivity')
# endif
c
      return
      end 


#else
      return
      end
#endif
