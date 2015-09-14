#ifdef explicit_fresh_water_flux
      subroutine get_ocean_bc 
!=======================================================================
!     This subroutine collects the sbc for the ocean model from the 
!     atmosphere and the ice model.
!     If no explicit fresh water is needed the surface fluxes are in 
!     sbcocn.
!     Note: sbcocn is updated prior an ocean segment and is kept 
!     constant during one ocean time step.
!     If explicit fresh water flux is required, sfft (and sffu) collects
!     the fresh water flux  from atmosphere, rivers and sea ice. strf
!     collects all tracer fluxes. If all fluxes are constant during one
!     ocean segment get_ocean_bc must be called befor an ocean segment
!     starts. With sea ice sfft and strf are updated prior each ocean 
!     time step and get_ocean_bc must be called within the segment loop. 
!
!     author:        m. schmidt e-mail=>martin.schmidt@io-warnemuende.de
!=======================================================================
# include "taskrows.h"
      call setsff   (jstask,jetask)
      call settracer(jstask,jetask)
      return
      end
      
      subroutine setsff (js, je)
!
!=======================================================================
!     collect all contributions of fresh water flux coming from
!     different component models:
!      - sbcocn comes from atmos and has just interpolated to the
!        ocean model grid in gosbc.
!      - rff is the river flow calculated in a hydrological module
!        (this may simply read a database)
!        for details see rivermod.F90, setriver.F, river.F, river.h
!      - other fluxes in connection with sea ice should be added here. 
!
!     js = global starting latitude index for the task
!     je = global ending latitude index for the task
!      - sfft collects the total fresh water flux on the t-grid
!      - sffu collects the total fresh water flux on the u-grid
!      - sbcocn(itflux) provides the total flux of tracer itflx
!        (useful only with coupled!)
!
!     This subroutine must be called before setvbc is called!  
!=======================================================================
!
      use arrays_2d
      use sbc_info
      use sbc_ocean
# ifdef river_inflow
      use river_mod
# endif
# include "derived_options.h"
# include "param.h"
# include "coord.h"
# include "grdvar.h"
# include "scalar.h"
!
# ifdef timing
      call tic ('setsff', 'setsff')
# endif
!
!----------------------------------------------------------------------
!     set freshwater flux onto U and T points 
!----------------------------------------------------------------------
!
      sfft(:,:) = 0.0
      sffu(:,:) = 0.0
#  ifdef need_sfft_expl
      sfft_expl(:,:) = 0.0
#  endif   
      
      do jrow=js,je
      
# ifdef simple_sbc
!       you may specify some fresh water flux here, for example
!       sfft(:,jrow) = sfft(:,jrow) + 3.e-6*min(kmt(i,jrow),1)
# else
!
!----------------------------------------------------------------------
!     add the freshwater flux from the atmosphere
!----------------------------------------------------------------------
!
        itflx = index_of_sbc ('fresh wtr')
        if (itflx .ne. 0) then
          sfft(:,jrow) = sfft(:,jrow)  
     &                         + sbcocn(:,jrow,itflx)*emsk(:,jrow)
#  ifdef atmos_tracer_explicit
          sfft_expl(:,jrow) = sfft_expl(:,jrow)  
     &                         + sbcocn(:,jrow,itflx)*emsk(:,jrow)
#  endif   
        endif
# endif   
      enddo
# ifdef river_inflow
!
!----------------------------------------------------------------------
!     add the river runoff to surface freshwater flux
!----------------------------------------------------------------------
!
      do nr=1, rivbox_in_pe
        ir = river_rec(nr)%ir
        jr = river_rec(nr)%jr
        sfft(ir,jr) = sfft(ir,jr) + river_rec(nr)%rff 
#  ifdef river_tracer_explicit
        sfft_expl(ir,jr) = sfft_expl(ir,jr) + river_rec(nr)%rff 
#  endif   
      enddo
# endif
# ifdef sea_ice
!
!----------------------------------------------------------------------
!     add the freshwater flux due to melting and freezing
!----------------------------------------------------------------------
!
# endif
!----------------------------------------------------------------------
!     calculate freshwater flux at u-points
!----------------------------------------------------------------------
!
      do jrow=js,je-1
	dyn  = dun(jrow)*cst(jrow+1)
	dys  = dus(jrow)*cst(jrow)
	dyr  = dyur(jrow)*csur(jrow)
        do i=2,imt-1
	  asw = duw(i)*dys
	  anw = duw(i)*dyn
	  ase = due(i)*dys
	  ane = due(i)*dyn
          sffu(i,jrow) = dyr*dxur(i)*(
     &                    sfft(i,jrow)*asw + sfft(i+1,jrow)*ase
     &                  + sfft(i,jrow+1)*anw + sfft(i+1,jrow+1)*ane)
        enddo
      enddo
      
!
# ifdef obc
!
!----------------------------------------------------------------------
!     set all surface fluxes on the open boundaries to zero
!----------------------------------------------------------------------
!
      do jrow=js,je
#  ifdef obc_south
        if (jrow .eq. 2) then
	  do i=1,imt
	    sfft(i,jrow) = c0
	  enddo
	endif
#  endif
#  ifdef obc_north
        if (jrow .eq. jmtm1) then
	  do i=1,imt
	    sfft(i,jrow) = c0
	  enddo
	endif
#  endif
#  ifdef obc_west
	sfft(2,jrow)     = c0
#  endif
#  ifdef obc_east
	sfft(imtm1,jrow) = c0
#  endif
      enddo
# endif
!
!----------------------------------------------------------------------
!     apply zonal boundary conditions
!----------------------------------------------------------------------
!
      call setbcx (sfft(1,js), imt, je-js+1)
!
# ifdef trace_indices
      write (stdout,'(2x,5(a,i4))') "=> In setsff: js=",js," je=",je
# endif
# ifdef timing
      call toc ('setsff', 'setsff')
# endif
!
      return
      end

      subroutine settracer (js, je)
!
!=======================================================================
!     collect all tracers coming with fresh water flux from
!     different model components:
!      - sbcocn comes from atmos and has just interpolated to the
!        ocean model grid in gosbc
!      - rtf is river tracer flow calculated in a hydrological module
!        (this may simply read a database)
!      - other fluxes in connection with sea ice will be added here. 
!
!     js = global starting latitude index for the task
!     je = global ending latitude index for the task
!      - sbcocn(itflux) collects the total flux of tracer itflx
!        (useful only with coupled!)
!
!     This subroutine must be called before setvbc is called!  
!
!=======================================================================
!
      use arrays_2d
      use sbc_info
      use sbc_ocean
# ifdef river_inflow
      use river_mod
# endif
# ifdef sea_ice
      use ice_mod
# endif
# include "derived_options.h"
# include "param.h"
# include "coord.h"
# include "grdvar.h"
# include "scalar.h"
# include "docnam.h"
!
# ifdef timing
      call tic ('get_ocean_bc', 'settracer')
# endif
!
# ifdef simple_sbc
!          The job is done in setvbc
# else
!     A general configuration is not possible.
!     Please check this code carefully. 
!
!     If sbcocn is the turbulent tracer flux (e.g., the sum of radiation
!     latent and sensible heat flux) disable atmos_tracer_explicit. 
!     In this case a "neutral" tracer flux with fresh water is added in
!     setvbc, so that a fresh water flux does not change tracer 
!     concentrationis but contributes to the total tracer. This is the
!     default. An exception is salt which is never neutral. If no saltflux
!     is specified, it is automatically zero.
!
!     If sbcocn is the total tracer flux (turbulent + tracer flux with 
!     fresh water) atmos_tracer_explicit must be enabled. 
!     
!     For two tracers (T,s) nothing has to be changed. If the salt flux is
!     zero, which is usually the case, the S.B.C. "saltflux" is not needed.
!     In this case strf(:,:,2) is set to zero below. Use namelist mbcin to
!     remove "saltflux" from the list of S.B.C.

!     For more than two tracers nloop must be modified!
!     If you remove the error stop, the additional tracer fluxes are
!     assumed to be zero (atmos_tracer_explicit enabled) or 
!     neutral (atmos_tracer_explicit disabled, this is default). 
!     If a flux must be prescribed, add lines following the scheme:
!        if (n.eq.1) itflx = index_of_sbc ('heatflux')
!        if (n.eq.2) itflx = index_of_sbc ('saltflux')
!        if (n.eq.3) itflx = index_of_sbc ('xxxxflux')
!        ...
!
      strf(:,:,:) = c0
      nloop: do n=1,nt
        itflx = 0
        if (n.eq.1) itflx = index_of_sbc ('heatflux')
        if (n.eq.2) itflx = index_of_sbc ('saltflux')
        if (n.gt.2) then
	   write(stdout,'(a)')
     &     '=> error in get_ocean_bc (atmosphere): '// 
     &     'Add definition of tracer fluxes for n>2!'
           call abort() 	    
	endif 
	if (itflx .ne. 0) then
#  ifdef trace_coupled_fluxes
#    ifdef atmos_tracer_explicit
          if (n.ne.2) write(stdout,'(a,i2.2,a,i2.2,a)') 
     &    '=> get_ocean_bc (atmosphere): '// 
     &    'The total '//trim(sbcname(itflx)) 
     &    //'['//trim(dunits(itflx))//'] from sbcocn(:,:,'
     &,   itflx,') goes to strf(:,:,',n,'). for '//trim(trname(n))
#    else
          if (n.ne.2) write(stdout,'(a,i2.2,a,i2.2,a/a)') 
     &    '=> get_ocean_bc (atmosphere): '// 
     &    'The turbulent '//trim(sbcname(itflx)) 
     &    //'['//trim(dunits(itflx))//'] from sbcocn(:,:,'
     &,   itflx,') goes to strf(:,:,',n,'). for '//trim(trname(n)) 
     &,   '   A neutral flux with fresh water is assumed.' 
#    endif
          if (n.eq.2) write(stdout,'(a,i2.2,a,i2.2,a)') 
     &    '=> get_ocean_bc (atmosphere): '// 
     &    'The total '//trim(sbcname(itflx)) 
     &    //'['//trim(dunits(itflx))//'] from sbcocn(:,:,'
     &,   itflx,') goes to strf(:,:,',n,'). for '//trim(trname(n))
#  endif
          do jrow=js,je
	    do i=1,imt
              strf(i,jrow,n) = sbcocn(i,jrow,itflx)
            enddo
	  enddo
#  ifdef trace_coupled_fluxes
	else
#    ifdef atmos_tracer_explicit
	  if (n.eq.1) write(stdout,'(a)') 
     &      '=> get_ocean_bc (atmosphere):'
     &      //' The total heatflux is set to zero.'
	  if (n.gt.2) write(stdout,'(a,i2.2,a)') 
     &      '=> get_ocean_bc (atmosphere): The total flux for tracer',n
     &,     ' is set to zero.'
#    else
	  if (n.eq.1) write(stdout,'(a)') 
     &      '=> get_ocean_bc (atmosphere):'
     &      //' The total heatflux is a neutral flux'
     &      //' with fresh water.'
	  if (n.gt.2) write(stdout,'(a,i2.2,a)') 
     &      '=> get_ocean_bc (atmosphere): The total flux for tracer',n
     &,     ' is a neutral flux with fresh water.'
#    endif
	  if (n.eq.2) write(stdout,'(a)') 
     &      '=> get_ocean_bc (atmosphere): '// 
     &      'The total saltflux is set to zero.'
#  endif
	endif
      end do nloop
#  ifdef river_inflow
!
!----------------------------------------------------------------------
!     add tracer fluxes due to river runoff
!----------------------------------------------------------------------
!
      do nr=1, rivbox_in_pe
        ir = river_rec(nr)%ir
        jr = river_rec(nr)%jr
        it = index_of_sbc ('heatflux')
        isst  = index_of_sbc ('sst')
        do n=1,nt
	  itflx = it + (n-1)
	  if (itflx .ne. 0.and.n.ne.2) then
#   if !defined river_tracer_explicit
#    ifdef trace_coupled_fluxes
	    if (nr.eq.1) then
	      write(stdout,'(a,i2.2,a)') 
     &        '=> get_ocean_bc (river): ' 
     &	      //'Assuming neutral river tracer flux '
     &        //'for strf(:,:,',n,'). '
     &        //trim(river_rec(nr)%trname(n))//' for '//trim(trname(n)) 
            endif
#    endif
#   else
#    ifdef trace_coupled_fluxes
	    if (nr.eq.1) then
	      write(stdout,'(a,a,i2.2,a,i2.2,a,a)') 
     &        '=> get_ocean_bc (river): ' 
     &	,     'Assigning rtf(',n,') to strf(:,:,',n,'). '
     &  ,     trim(river_rec(nr)%trname(n))//' for '//trim(trname(n))
            endif
#    endif
	    strf(ir,jr,n) = strf(ir,jr,n) + river_rec(nr)%rtf(n) 
#   endif			
	  endif
        enddo
      enddo
#  endif             
#  ifdef sea_ice
!
!----------------------------------------------------------------------
!     add the freshwater tracer flux due to melting and freezing
!----------------------------------------------------------------------
      do jrow=js,je
	do i=1,imt
          it = index_of_sbc ('heatflux')
          do n=1,nt
	    itflx = it + (n-1)
	    if (itflx .ne. 0) then
              strf(i,jrow,n) = strf(i,jrow,n) + 
!     &                                 something
	    endif
          enddo
	enddo
      enddo
!
#  endif
# endif /* ifdef simple_sbc */
!
# ifdef obc
!
!----------------------------------------------------------------------
!     set all surface fluxes on the open boundaries to zero
!----------------------------------------------------------------------
!
# endif
      return
      end
#else
      subroutine get_ocean_bc 
      return
      end
#endif
