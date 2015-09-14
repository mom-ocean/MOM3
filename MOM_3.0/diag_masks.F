      subroutine setglobalmasks (mask,maskvertical,whichmask)
#if defined tracer_averages || defined term_balances
# define traceravgbal
#endif
#if defined meridional_overturning_basin || defined traceravgbal || defined gyre_components
c=======================================================================
c
c     Subroutine setglobalmasks calculates the masks used in northward
c     meridional overturning calculation 
c     (ifdef meridional_overturning_basin) and
c     in volume weighted averages of tracers and surface tracer fluxes
c     (ifdefs tracer_averages, term_balances, gyre_components).
c
c     This subroutine is used for global models only.  it sets
c     northmasks to:
c     1 = arctic and atlantic
c     2 = indo-pacific
c     3 = inland seas, including Mediterranean.
c
c     and regionmasks to:
c     1 = southern ocean
c     2 = atlantic including Mediterranean and Hudson bay.
c     3 = pacific
c     4 = arctic
c     5 = indian
c
c     called by:
c     nor1st and reg1st
c
c     input: 
c     whichmask = a string: 'northmask' or 'regionmask'
c
c     output:
c     mask = a mask array
c
c     author: Eli Tziperman   e-mail => ett@gfdl.gov
c                                       eli@beach.weizmann.ac.il
c=======================================================================
c
      use arrays_2d
      use memory_window
#include "param.h"
#include "coord.h"
      integer mask(imt,jmt), maskvertical(km)
      character*(*) whichmask
      logical found_unset_mask_point, setglobalmasks_error

      write(*,*)" Entering setglobalmask, whichmask = ",whichmask
c
c-----------------------------------------------------------------------
c     use kmt to initialize land to zero, ocean to -999:
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do i=1,imt
          if (kmt(i,j).ne.0) then
            mask(i,j)=-999
          else
            mask(i,j)=0
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     set mask lines that separate between major ocean basins.  These
c     lines serve as "seeds" for setting the masks elesewhere:
c-----------------------------------------------------------------------
c
c     first, northmasks:
c
      if (whichmask .eq. 'northmask') then

c       Atlantic-Pacific divide:
        i1 = indp (360.0-69.0, xt, imt)
        j1 = indp (-54.0, yt, jmt)
        i2 = indp (360.0-69.0, xt, imt)
#if defined resolution_1x1
        j2 = indp (-67.5, yt, jmt)
#else
        j2 = indp (-85.0, yt, jmt)
#endif
        write(*,*)" Atlantic-Pacific divide, i1,j1,i2,j2=",i1,j1,i2,j2
        i = i1
        do j=j2,j1
          if (mask(i,j) .eq. -999) mask(i,j)=1
          if (mask(i-1,j) .eq. -999) mask(i-1,j)=2
        enddo

c       Atlantic-Indian divide:
        i1 = indp (20.0, xt, imt)
        j1 = indp (-34.0, yt, jmt)
        i2 = indp (20.0, xt, imt)
        j2 = indp (-85.0, yt, jmt)
        write(*,*)"Atlantic-Indian divide, i1,j1,i2,j2=",i1,j1,i2,j2
        i = i1
        do j=j2,j1
          if (mask(i,j) .eq. -999) mask(i,j)=1
          if (mask(i+1,j) .eq. -999) mask(i+1,j)=2
        enddo

c       Pacific-Arcric divide:
        i1 = indp (360.0-175.0, xt, imt)
        j1 = indp (66.0, yt, jmt)
        i2 = indp (360.0-164.0, xt, imt)
        j2 = indp (66.0, yt, jmt)
        write(*,*)"Pacific-Arcric divide, i1,j1,i2,j2=",i1,j1,i2,j2
        j = j1
        do i=i1,i2
          if (mask(i,j) .eq. -999) mask(i,j)=2
          if (mask(i,j+1) .eq. -999) mask(i,j+1)=1
        enddo

c       Mediterranean seed:
        i1 = indp (18.0, xt, imt)
        j1 = indp (36.0, yt, jmt)
        write(*,*)"Med seed, i1,j1=",i1,j1
        i = i1
        j = j1
        if (mask(i,j) .eq. -999) then
          mask(i,j) = 3
        else
          write(*,*) " setglobalmask error. Med seed over land?"
          call abort()
        endif

#if defined resolution_1x1
          write(*,*) "no Hudson seed"
#else
c       Hudson bay seed:
        i1 = indp (360.0-84.0, xt, imt)
        j1 = indp (60.0, yt, jmt)
        write(*,*)"Hudson bay seed, i1,j1=",i1,j1
        i = i1
        j = j1
        if (mask(i,j) .eq. -999) then
          mask(i,j)=3
        else
          write(*,*) " setglobalmask error. Hudson seed over land?"
          call abort()
        endif
#endif
      endif
c
c     second, region masks:
c
      if (whichmask .eq. 'regionmask') then

c       Southern Ocean-Atlantic divide #1:
        i1 = 1
        j1 = indp (-54.0, yt, jmt)
        i2 = indp (20.0, xt, imt)
        j2 = indp (-54.0, yt, jmt)
        write(*,*)"Southern Ocean-Atl divide, i1,j1,i2,j2=",i1,j1,i2,j2
        j = j1
        do i=i1,i2
          if (mask(i,j) .eq. -999) mask(i,j)=1
          if (mask(i,j+1) .eq. -999) mask(i,j+1)=2
        enddo
c       Southern Ocean-Atlantic divide #2:
        i1 = indp (360.0-69.0, xt, imt)
        j1 = indp (-54.0, yt, jmt)
        i2 = imt
        j2 = indp (-54.0, yt, jmt)
        write(*,*)"Southern Ocean-Atl divide, i1,j1,i2,j2=",i1,j1,i2,j2
        j = j1
        do i=i1,i2
          if (mask(i,j) .eq. -999) mask(i,j)=1
          if (mask(i,j+1) .eq. -999) mask(i,j+1)=2
        enddo

c       Southern Ocean-Indian divide:
        i1 = indp (20.0, xt, imt)+1
        j1 = indp (-54.0, yt, jmt)
        i2 = indp (146.0, xt, imt)
        j2 = indp (-54.0, yt, jmt)
        write(*,*)"Southern Ocn-Indian divide, i1,j1,i2,j2=",i1,j1,i2,j2
        j = j1
        do i=i1,i2
          if (mask(i,j) .eq. -999) mask(i,j)=1
          if (mask(i,j+1) .eq. -999) mask(i,j+1)=5
        enddo

c       Southern Ocean-Pacific divide:
        i1 = indp (146.0, xt, imt)+1
        j1 = indp (-54.0, yt, jmt)
        i2 = indp (360.0-69.0, xt, imt)
        j2 = indp (-54.0, yt, jmt)
        write(*,*)"Southern Ocn-Pacfc divide, i1,j1,i2,j2=",i1,j1,i2,j2
        j = j1
        do i=i1,i2
          if (mask(i,j) .eq. -999) mask(i,j)=1
          if (mask(i,j+1) .eq. -999) mask(i,j+1)=3
        enddo

c       Atlantic-Indian divide:
        i1 = indp (20.0, xt, imt)
        j1 = indp (-34.0, yt, jmt)
        i2 = indp (20.0, xt, imt)
        j2 = indp (-54.0, yt, jmt)+1
        write(*,*)"Atlantic-Indian divide, i1,j1,i2,j2=",i1,j1,i2,j2
        i = i1
        do j=j2,j1
          if (mask(i,j) .eq. -999) mask(i,j)=2
          if (mask(i+1,j) .eq. -999) mask(i+1,j)=5
        enddo
#if defined resolution_1x1
c       Pacific-Indian divide #1 (in southern ocean):
        i1 = indp (147.5, xt, imt)
        j1 = indp (-34.0, yt, jmt)
        i2 = indp (147.5, xt, imt)
        j2 = indp (-54.0, yt, jmt)+1
        write(*,*)"Pacific-Indian divide#1, i1,j1,i2,j2=",i1,j1,i2,j2
        i = i1
        do j=j2,j1
          if (mask(i,j) .eq. -999) mask(i,j)=5
          if (mask(i+1,j) .eq. -999) mask(i+1,j)=3
        enddo
c       Pacific-Indian divide #2 (between Australia and Asia):
        i1 = indp (101.5, xt, imt)
        j1 = indp (-0.8, yt, jmt)
        i2 = indp (130.5, xt, imt)
        j2 = indp (-0.8, yt, jmt)
        write(*,*)"Pacific-Indian divide#2, i1,j1,i2,j2=",i1,j1,i2,j2
        j = j1
        do i=i1,i2
          if (mask(i,j) .eq. -999) mask(i,j)=5
          if (mask(i,j+1) .eq. -999) mask(i,j+1)=3
        enddo

#else
c       Pacific-Indian divide #1 (in southern ocean):
        i1 = indp (146.0, xt, imt)
        j1 = indp (-34.0, yt, jmt)
        i2 = indp (146.0, xt, imt)
        j2 = indp (-54.0, yt, jmt)+1
        write(*,*)"Pacific-Indian divide#1, i1,j1,i2,j2=",i1,j1,i2,j2
        i = i1
        do j=j2,j1
          if (mask(i,j) .eq. -999) mask(i,j)=5
          if (mask(i+1,j) .eq. -999) mask(i+1,j)=3
        enddo
c       Pacific-Indian divide #2 (between Australia and Asia):
        i1 = indp (115.0, xt, imt)
        j1 = indp (-25.0, yt, jmt)
        i2 = indp (115.0, xt, imt)
        j2 = indp (0.0, yt, jmt)+1
        write(*,*)"Pacific-Indian divide#2, i1,j1,i2,j2=",i1,j1,i2,j2
        i = i1
        do j=j1,j2
          if (mask(i,j) .eq. -999) mask(i,j)=5
          if (mask(i+1,j) .eq. -999) mask(i+1,j)=3
        enddo
#endif

c       Pacific-Arcric divide:
        i1 = indp (360.0-175.0, xt, imt)
        j1 = indp (66.0, yt, jmt)
        i2 = indp (360.0-164.0, xt, imt)
        j2 = indp (66.0, yt, jmt)
        write(*,*)"Pacific-Arcric divide, i1,j1,i2,j2=",i1,j1,i2,j2
        j = j1
        do i=i1,i2
          if (mask(i,j) .eq. -999) mask(i,j)=3
          if (mask(i,j+1) .eq. -999) mask(i,j+1)=4
        enddo

c       Atlantic-Arcric divide #1:
        i1 = indp (360.0-45.0, xt, imt)
        j1 = indp (66.0, yt, jmt)
        i2 = imt
        j2 = indp (66.0, yt, jmt)
        write(*,*)"Atlantic-Arcric divide #1, i1,j1,i2,j2=",i1,j1,i2,j2
        j = j1
        do i=i1,i2
          if (mask(i,j) .eq. -999) mask(i,j)=2
          if (mask(i,j+1) .eq. -999) mask(i,j+1)=4
        enddo
c       Atlantic-Arcric divide #2:
        i1 = 1
        j1 = indp (66.0, yt, jmt)
        i2 = indp (16.0, xt, imt)
        j2 = indp (66.0, yt, jmt)
        write(*,*)"Atlantic-Arcric divide #2, i1,j1,i2,j2=",i1,j1,i2,j2
        j = j1
        do i=i1,i2
          if (mask(i,j) .eq. -999) mask(i,j)=2
          if (mask(i,j+1) .eq. -999) mask(i,j+1)=4
        enddo

c       Mediterranean seed:
        i1 = indp (18.0, xt, imt)
        j1 = indp (36.0, yt, jmt)
        write(*,*)"Med seed, i1,j1=",i1,j1
        i = i1
        j = j1
        if (mask(i,j) .eq. -999) then
          mask(i,j)=3
        else
          write(*,*) " setglobalmask error. Med seed over land?"
          call abort()
        endif

#if defined resolution_1x1
          write(*,*) "No Hudson seed"
#else
c       Hudson bay seed:
        i1 = indp (360.0-84.0, xt, imt)
        j1 = indp (60.0, yt, jmt)
        write(*,*)"Hudson bay seed, i1,j1=",i1,j1
        i = i1
        j = j1
        if (mask(i,j) .eq. -999) then
          mask(i,j)=3
        else
          write(*,*) " setglobalmask error. Hudson seed over land?"
          call abort()
        endif
#endif

      endif

      write(*,*) "northmask after setting the seeds:"
      call imatrx (mask, imt, 1, imt, -1, -jmt, 1)
c
c-----------------------------------------------------------------------
c     Iterate to complete mask values in all other locations by
c     propagating the values from the seeds:
c-----------------------------------------------------------------------
c
      found_unset_mask_point = .true.
      setglobalmasks_error = .false.
      iter = 0
      do while (found_unset_mask_point .and. iter.lt.imt*jmt
     &     .and. .not. setglobalmasks_error)
        iter=iter+1
        found_unset_mask_point=.false.
        do j=1,jmt
          do i=1,imt
            if (mask(i,j).eq.-999) then
              found_unset_mask_point=.true.
              do ip=-1,1
                do jp=-1,1
c                 go over neighboring points, but not diagonal neighboring points:
                  if (ip**2+jp**2.le.1) then
                    is=min(imt,i+ip)
                    is=max(is,1)
                    js=min(jmt,j+jp)
                    js=max(js,1)
c                   if (i,j) is still at -999, and if neighbor point
c                   (is,js) is set (.ne.-999), but not to land value
c                   (0), use it to set the mask at i,j:
                    if (mask(i,j).eq.-999) then
                      if (mask(is,js).ne.-999.and.mask(is,js).ne.0) then
                        mask(i,j)=mask(is,js)
                      endif
                    else if (mask(i,j).ne.mask(is,js)
     &                     .and.mask(is,js).ne.0
     &                     .and.mask(is,js).ne.-999) then
c                     if point (i,j) was already set during this
c                     iteration using some other (is,js), but the value
c                     that was set is not consistent with a non-land,
c                     neighbor with a set mask (.ne.-999) that is looked
c                     at now, complain:
                      write(*,*)"setglobalmasks: error at i,j=",i,j
                      setglobalmasks_error=.true.
                      go to 999
                    endif
                  endif
                enddo
              enddo
            endif
          enddo
        enddo
      enddo

 999  if (setglobalmasks_error) then
        write(*,*) "Stopping due to error in setglobalmasks: Specified "
        write(*,*) "divides and seeds do not lead to a unique mask."
        write(*,*) "Here is the mask when stopping:"
        call imatrx (mask, imt, 1, imt, -1, -jmt, 1)
        call abort()
      endif

c     check that all points are set:
      ifound_unset_mask = 0
      do j=1,jmt
        do i=1,imt
          if (mask(i,j).eq.-999) then
            ifound_unset_mask = 1
            if (whichmask .eq. 'northmask') then
              mask(i,j) = 3
            else if (whichmask .eq. 'regionmask') then
              mask(i,j) = 2
            endif
          endif
        enddo
      enddo
      if (ifound_unset_mask .eq. 1) then
        write(*,*) " setglobalmasks: Warning: mask is not set at all "
     &       ,"grid points.  These points were set arbitrarily to"
     &       ," atlantic (in regionmasks) & inland seas (in northmask)."
     &       ," Check for missing seed values in some of the basins."
      endif



c     set vertical region masks:
      if (whichmask .eq. 'regionmask') then
        k1 = indp (250.0, zt, km)
        k2 = kbot
        do k=1,k1
          maskvertical(k)=1
        enddo
        do k=k2,kbot
          maskvertical(k)=2
        enddo
      endif


#endif
      return
      end
