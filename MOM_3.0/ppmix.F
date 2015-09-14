#ifdef ppvmix
      subroutine ppmixi (error, cifdef, ifdmax, nifdef, vmixset)
c
      logical vmixset, error
      integer io_status
      character*(*) cifdef(ifdmax)
c
# include "size.h"
# include "accel.h"
# include "coord.h"
# include "iounit.h"
# include "rinum.h"
# include "scalar.h"
# include "stdunits.h"
# include "vmixc.h"
c
c=======================================================================
c     Initialization for the Pacanowski/Philander vertical mixing scheme
c     Pacanowski & Philander (JPO vol 11, #11, 1981). 
c    
c     input:
c       dzt    = thickness of vertical levels (cm)
c       km     = number of vertical levels
c       yt     = latitude of grid points (deg)
c       jmt    = number of latitudes
c       dtxcel = time step accelerator as a function of level
c       dtts   = density time step (sec)
c       dtuv   = internal mode time step (sec)
c       error  = logical to signal problems
c       cifdef = array of character strings for listing enabled "ifdefs"
c       ifdmax = size of "cifdef"
c       nifdef = current number of enabled "ifdefs" 
c       vmixset= logical to determine if a vertical mixing scheme was
c                chosen
c
c     output:
c       wndmix = min value for mixing at surface to simulate high freq
c                wind mixing (if absent in forcing). (cm**2/sec)
c       fricmx = maximum mixing (cm**2/sec)
c       diff_cbt_back = background "diff_cbt" (cm**2/sec)
c       visc_cbu_back = background "visc_cbu" (cm**2/sec)
c       diff_cbt_limit = largest "diff_cbt" (cm**2/sec)
c       visc_cbu_limit = largest "visc_cbu" (cm**2/sec)
c       cifdef = array of character strings for listing enabled "ifdefs"
c       nifdef = incremented by 1 if this routine is called 
c       error  = true if some inconsistency was found
c       vmixset= true
c       
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
c
      namelist /ppmix/ wndmix, fricmx, diff_cbt_back, visc_cbu_back
     &,                visc_cbu_limit, diff_cbt_limit
c
      write (stdout,'(/,20x,a,/)')
     & 'P P V M I X    I N I T I A L I Z A T I O N'
c
c-----------------------------------------------------------------------
c     initialize variables (all mixing units are cm**2/sec.)
c-----------------------------------------------------------------------
c
      wndmix    = 10.0
      fricmx    = 50.0
      diff_cbt_back =  0.1
      visc_cbu_back =  1.0
      dzmin     = 1.e10
      p25       = 0.25
      c0        = 0.0
c
# ifdef implicitvmix
c
c     simulate convective adjustment with large mixing coefficient
c     limits
c
      visc_cbu_limit = fricmx
      diff_cbt_limit = 1.0e6
# else
c
c     in regions of gravitational instability set mixing limits to the
c     maximum consistent with the "cfl" criterion. convective adjustment
c     will also act on the instability.
c
      visc_cbu_limit = fricmx
      diff_cbt_limit = fricmx
# endif
c
c-----------------------------------------------------------------------
c     provide for namelist over-ride of above settings + documentation
c-----------------------------------------------------------------------
c
      call getunit (io, 'namelist.ppmix'
     &,             'formatted sequential rewind')
      read (io,ppmix, IOSTAT=io_status)
      if (io_status .gt. 0) then
        write(stdout,*) 'IOSTAT reading ppmix namelist =',io_status
        call abort()
      endif
c
c-----------------------------------------------------------------------
c     set no-flux condition on density difference across bottom level
c-----------------------------------------------------------------------
c
      do j=1,jmw
        do i=1,imt
	  rhom1z(i,km,j) = c0
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     add character string to "ifdef option list" indicating that this
c     option is enabled
c-----------------------------------------------------------------------
c
      nifdef = nifdef + 1
      cifdef(nifdef) = 'ppvmix      '
c
c-----------------------------------------------------------------------
c     check for problems
c-----------------------------------------------------------------------
c
# if defined ppvmix && !defined implicitvmix && defined redi_diffusion
        write (stdout,'(/,(1x,a))')
     & '==> Error:  "ppvmix" must use "implicitvmix" if               '
     &,'            "redi_diffusion" is also enabled                  '
        error = .true.
# endif
      if (vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "ppvmix" cannot be enabled because another         '
     &,'            vertical mixing scheme has been enabled           '
        error = .true.
      else
        vmixset = .true.
      endif
c
      do k=1,km
        dzmin = min(dzmin,dzt(k))
      enddo
      if (dzmin .ge. 25.e2) then
        write (stdout,'(/,(1x,a))')
     & '==> Warning: "ppvmix" may not work well with coarse vertical  '
     &,'              resolution                                      '
      endif
c
      extlat = c0
      do jrow=1,jmt
        extlat = max(abs(yt(jrow)),extlat)
      enddo
      if (extlat .gt. 10.0) then
        write (stdout,'(/,(1x,a))')
     & '==> Warning: "ppvmix" may not work well outside the tropics   '
     &,'              where vertical shear is small unless solar      '
     &,'              shortwave penetration into the ocean is         '
     &,'              accounted for by enabeling  "shortwave"         '
      endif
c
# if !defined implicitvmix
      do k=1,km
        if ((dtts*dtxcel(k)*fricmx)/dzt(k)**2 .ge. p25) then
          write (stdout,'(/,(1x,a))')
     & '==> Error: vertical diffusive criteria exceeded for           '
     &,'           "fricmx".  use a smaller "dtts", "dtxcel", and/or  '
     &,'           "fricmx" .... or enable "implicitvmix"             '
          write (stdout,'(a48,i3)') ' at level =',k
          error = .true.
        endif
        if ((dtts*dtxcel(k)*diff_cbt_limit)/dzt(k)**2 .ge. p25) then
          write (stdout,'(/,(1x,a))')
     & '==> Error: vertical diffusive criteria exceeded for           '
     &,'            "diff_cbt_limit". use a smaller "dtts", "dtxcel"  '
     &,'         ,and/or  "diff_cbt_limit" ...or enable "implicitvmix"'
          write (stdout,'(a48,i3)') ' at level =',k
          error = .true.
        endif
      enddo
c
      if ((dtuv*fricmx)/dzmin**2 .ge. p25) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: vertical diffusive criteria exceeded for           '
     &,'            "fricmx". use a smaller "dtuv" and/or "fricmx"    '
     &,'             or enable "implicitvmix"                         '
        error = .true.
      endif
c
      if ((dtuv*visc_cbu_limit)/dzmin**2 .ge. p25) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: vertical diffusive criteria exceeded for           '
     &,'            "visc_cbu_limit". use a smaller "dtuv" or         '
     &,'            "visc_cbu_limit" or enable "implicitvmix"         '
        error = .true.
      endif
# else
        write (stdout,'(/,(1x,a))')
     & '==> Warning: enabeling "implicitvmix" with "ppvmix" uses      '
     &,'             variables defined at "tau" rather than at "tau-1"'
     &,'             as was done in MOM 1.x'
# endif
# ifdef bryan_lewis_vertical
        write (stdout,'(/,(1x,a/1x,a/1x,a/1x,a))')
     & '==> Warning: "bryan_lewis_vertical" tracer diffus coefficients'
     &,'              will  be added to "ppvmix" diffus coefficients  '
     &,'              Note that diff_cbt_back is being reset to zero  '
     &,'              while diff_cbu_back is unchanged                '
      diff_cbt_back =  0.0
# endif
c
c     write out namelist values
c
      write (stdout,ppmix)
      call relunit (io)
      call getunit (iodoc, 'document.dta'
     &,             'formatted sequential append')
      write (iodoc, ppmix)
      call relunit (iodoc)
      return
      end


      subroutine ppmix (joff, js, je, is, ie)
c
c=======================================================================
c     Compute vertical mixing coefficients based on...
c     Pacanowski & Philander (JPO vol 11, #11, 1981). 
c
c     Note: this parameterization was designed for equatorial models
c     and may not do a good job in mid or high latitudes. Simulations
c     in these regions (where vertical shear is small) are improved with
c     the addition of solar short wave penetration into the ocean which 
c     reduces buoyancy and enhances vertical mixing.
c
c     inputs:
c
c      joff   = offset between rows in the MW and latitude rows
c      js     = starting row for loading variables to calculate
c               coefficients. calcs start at jstrt=max(js-1,1+jbuf)
c      je     = ending row for loading variables to calculate
c               coefficients. calculations end at je-1
c      is     = starting index for calculating coefficients in the 
c               longitude direction
c      ie     = ending index for calculating coefficients in the 
c               longitude direction
c      km     = number of vertical levels
c      grav   = gravity (cm/sec**2)
c      umask  = land/sea mask on U grid (land=0.0, sea=1.0)
c      tmask  = land/sea mask on T grid (land=0.0, sea=1.0)
c      fricmx = max viscosity (cm**2/sec)
c      wndmix = min viscosity at bottom of 1st level to simulate 
c               missing high frequency windstress components (cm**2/sec)
c      visc_cbu_back  = background "visc_cbu" (cm**2/sec)
c      diff_cbt_back  = background "diff_cbt" (cm**2/sec)
c      visc_cbu_limit = largest "visc_cbu" in regions of gravitational
c                      instability (cm**2/sec)
c      diff_cbt_limit = largest "diff_cbt" in regions of gravitational 
c                       instability (cm**2/sec)
c      riu    = richardson number at bottom of U cells
c      rit    = richardson number at bottom of T cells
c
c     outputs:
c
c      visc_cbu = viscosity coefficient at bottom of U cells (cm**2/s)
c      diff_cbt = diffusion coefficient at bottom of T cells (cm**2/s)
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      use memory_window
# include "param.h"
# include "coord.h"
# include "grdvar.h"
# include "mw.h"
# include "rinum.h"
# include "scalar.h"
# include "switch.h"
# include "vmixc.h"
c
# ifdef timing
      call tic ('vmixc', 'ppmix')
# endif
c
c-----------------------------------------------------------------------
c     set local constants 
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
c-----------------------------------------------------------------------
c     compute gradient richardson number at base of T cells and U cells
c-----------------------------------------------------------------------
c
      call rinum (joff, js, je, is, ie)
c
c-----------------------------------------------------------------------
c     set starting and ending indices for computation of mixing
c     coefficients diff_cbt and visc_cbu
c-----------------------------------------------------------------------
c
      jstrt = max(js-1,1+jbuf)
      jend  = je-1
c
c-----------------------------------------------------------------------
c     viscosity and diffusivity coeffs are on bottom of T and U cells.
c-----------------------------------------------------------------------
c
c
      do j=jstrt,jend
        do k=1,km-1
          do i=istrt,iend
            t1                = c1/(c1 + c5*riu(i,k,j))
            visc_cbu(i,k,j)   = fricmx*t1**2 + visc_cbu_back
            t2                = c1/(c1 + c5*rit(i,k,j))
            diff_cbt(i,k,j,1) = fricmx*t2**3 + diff_cbt_back
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     limit coeffs on bottom of cells in unstable regions
c-----------------------------------------------------------------------
c
      do j=jstrt,jend
        do k=1,km-1
          do i=istrt,iend
            if (riu(i,k,j) .lt. c0) visc_cbu(i,k,j)   = visc_cbu_limit
            if (rit(i,k,j) .lt. c0) diff_cbt(i,k,j,1) = diff_cbt_limit
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     approximation for high freq wind mixing near the surface
c     set no flux through bottom of bottom level "km"
c-----------------------------------------------------------------------
c
      do j=jstrt,jend
        do i=istrt,iend
          if (diff_cbt(i,1,j,1) .lt. wndmix) diff_cbt(i,1,j,1) = wndmix
          if (visc_cbu(i,1,j)   .lt. wndmix) visc_cbu(i,1,j)   = wndmix
	  diff_cbt(i,km,j,1) = c0
	  visc_cbu(i,km,j) = c0
        enddo
      enddo
c
#ifdef bryan_lewis_vertical
c
c-----------------------------------------------------------------------
c     add Bryan-Lewis mixing if wanted
c-----------------------------------------------------------------------
c
      do j=jstrt,jend
        do k=1,km-1
          do i=istrt,iend
            diff_cbt(i,k,j,1) = diff_cbt(i,k,j,1) + Ahv(k)
          enddo
        enddo
      enddo
#endif
c
c-----------------------------------------------------------------------
c       set lateral bc
c-----------------------------------------------------------------------
c
      do j=jstrt,jend
        call setbcx (visc_cbu(1,1,j), imt, km)
        call setbcx (diff_cbt(1,1,j,1), imt, km)
      enddo
c
# ifdef matrix_sections
      if (prxzts .and. eots) then
        call diagvmix (joff, jstrt, jend)
      endif
# endif
# ifdef trace_indices
      write (stdout,'(2x,7(a,i4))')
     & "=> In ppmix: js=",js," je=",je," joff=",joff
     &," jstrt=",jstrt," jend=",jend," jrow=",jstrt+joff
     &," to ",jend+joff 
# endif
# ifdef timing
      call toc ('vmixc', 'ppmix')
# endif
c
      return
      end

#else
      subroutine ppmix
      return
      end
#endif
