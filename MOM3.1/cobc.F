#if defined obc_north || defined obc_south
      subroutine cobc (joff, j, is, ie)
      use arrays_2d
      use memory_window
      use sbc_info
      use sbc_ocean
c
c=======================================================================
c     compute tracers at "tau+1" for rows 2 and/or jmt-1
c     using open boundary conditions as in Stevens (1990).
c
c
c     input:
c       j    = boundary row in the MW
c       joff = offset relating "j" in the MW to latitude "jrow"
c            => jrow=j+joff  = boundary latitude (2 or jmt-1)
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c
c
c     Reference: Stevens, D.P., 1990: On open boundary condidtions for
c                three dimensional primitive equation ocean circulation 
c                models, Geophys. Astrophys. Fluid Dynamics, 51, 103-133
c
c     author:   Arne Biastoch   e-mail: abiastoch@ifm.uni-kiel.de
c               based on a version for MOM 1 by Rene Redler
c=======================================================================
c
# include "param.h"
# ifndef obc_north_anywhere
      parameter (istrt=2, iend=imt-1)
# endif
# include "accel.h"
# include "cobc.h"
# include "coord.h"
# include "emode.h"
# include "grdvar.h"
# include "hmixc.h"
# ifdef isoneutralmix
#  include "isopyc.h"
# endif
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# include "vmixc.h"
      dimension twodt(km)
      dimension ADV_Ty_obc(imt,km)
# include "fdift.h"
#ifdef partial_cell
c
c     additional statement functions for obc
c     (same as for interior but for tau)
c
      T_i_obc(i,k,j,n,ip) = t(i+ip,max(1,k-1),j,n,tau) - dhwe(i,k-1,j)
     &            *(t(i+ip,max(1,k-1),j,n,tau) - t(i+ip,k,j,n,tau))
     &               /dhwt(i+ip,k-1,j)
      T_j_obc(i,k,j,n,jp) = t(i,max(1,k-1),j+jp,n,tau) - dhwn(i,k-1,j)
     &            *(t(i,max(1,k-1),j+jp,n,tau) - t(i,k,j+jp,n,tau))
     &               /dhwt(i,k-1,j+jp)
#else
      T_i_obc(i,k,j,n,ip) = t(i+ip,k,j,n,tau)
      T_j_obc(i,k,j,n,jp) = t(i,k,j+jp,n,tau)
#endif
# ifdef timing
      call tic ('ocean', 'cobc')
# endif
c
c-----------------------------------------------------------------------
c     limit the longitude indices based on those from the argument list
c     Note: this is currently bypassed. istrt and iend are set as
c           parameters to optimize performance
c-----------------------------------------------------------------------
c
      jrow = j + joff
c
# ifdef obc_south
c
      if (jrow .eq. 2) then
#  ifdef obc_north_anywhere
         istrt = 2
         iend  = imt-1
#  endif
c
c-----------------------------------------------------------------------
c       1) compute the advective velocity "vad" 
c          at the north face of the T cells
c
c       2) compute phase velocity at the southern boundary: c1s
c-----------------------------------------------------------------------
c
        var = -dyu(jrow+1)/dtts 
c
        do k=1,km
          do i=2,imtm1
c	  
            vad(i,k) = (u(i,k,j,2,tau)*dxt(i)+u(i-1,k,j,2,tau)
     &                 *dxt(i-1))/(dxt(i)+dxt(i-1))*tmask(i,k,j)
            if (vad(i,k) .gt. c0) vad(i,k) = c0
c	     
	    do m=1,nt
	      var1(i,k,m) = t(i,k,j+2,m,taum1)-t(i,k,j+1,m,taum1)
	      if (var1(i,k,m) .eq. c0) then
		c1s(i,k,m) = var
	      else
		c1s(i,k,m) = var*(t(i,k,j+1,m,tau)-t(i,k,j+1,m,taum1))
     &                           /var1(i,k,m)*tmask(i,k,j)
                if (c1s(i,k,m) .gt. c0)  c1s(i,k,m) = c0
		if (c1s(i,k,m) .lt. var) c1s(i,k,m) = var
	      endif
	    enddo
c	     
          enddo
#  ifdef cyclic
          vad(  1,k) = vad(imtm1,k)
          vad(imt,k) = vad(    2,k)
	  do m=1,nt
              c1s(  1,k,m) = c1s(imtm1,k,m)
              c1s(imt,k,m) = c1s(    2,k,m)
	  enddo
#  else
          vad(imt,k) = c0
          vad(  1,k) = c0
#  endif
        enddo
      endif
# endif
c
# ifdef obc_north
c
#  ifdef obc_north_anywhere
      if (jrow .eq. jobcn) then
         istrt = iobcns
         iend  = iobcne
#  else
      if (jrow .eq. jmtm1) then
#  endif
c
c-----------------------------------------------------------------------
c       1) compute the advective velocity "vad" 
c          at the south face of the T cells
c
c       2) compute phase velocity at the northern boundary: c1n
c          (var1 is calculated in subroutine tracer because of 
c           running out of bounds)
c-----------------------------------------------------------------------
c
        var = dyu(jrow-1)/dtts
c
        do k=1,km
          do i=2,imtm1
c	  
            vad(i,k) = (u(i,k,j-1,2,tau)*dxt(i)+u(i-1,k,j-1,2,tau)
     &                 *dxt(i-1))/(dxt(i)+dxt(i-1))*tmask(i,k,j)
            if (vad(i,k) .lt. c0) vad(i,k) = c0
c
            do m=1,nt
	      if (var1(i,k,m) .eq. c0) then
 	        c1n(i,k,m) = var
	      else
 	        c1n(i,k,m) = -var*(t(i,k,j-1,m,tau)-t(i,k,j-1,m,taum1))
     &                           /var1(i,k,m)*tmask(i,k,j)
                if (c1n(i,k,m) .lt. c0)  c1n(i,k,m) = c0
	        if (c1n(i,k,m) .gt. var) c1n(i,k,m) = var
	      endif
	    enddo
c	   
          enddo
#  ifdef cyclic
          vad(  1,k) = vad(imtm1,k)
          vad(imt,k) = vad(    2,k)
	  do m=1,nt
	     c1n(  1,k,m) = c1n(imtm1,k,m)
	     c1n(imt,k,m) = c1n(    2,k,m)
	  enddo
#  else
          vad(imt,k) = c0
          vad(  1,k) = c0
#  endif
        enddo
c
      endif
# endif
c
c
      do n=1,nt
c
c-----------------------------------------------------------------------
c       diffusive flux across eastern face of T cells
c       (use constant horizontal diffusion)
c-----------------------------------------------------------------------
c
#ifdef tracer_horz_mix_const
        ah_csu_dyur =  diff_cnt*csu_dyur(jrow)
#endif
c
        do k=1,kbot
          do i=istrt-1,iend
c
c          only Laplacian diffusion 
c
           diff_fe(i,k,j) = (T_i_obc(i,k,j,n,1)
     &                      -T_i_obc(i,k,j,n,0))
#ifdef tracer_horz_mix_const
     &                      *ah_cstdxur(i,j)
#endif
#ifdef tracer_horz_mix_var
     &                      *diff_cet(i,k,j)*cstdxur(i,j)
#endif
#if defined partial_cell
     &                      *dhte(i,k,j)
#endif
c
            do jj=j-1,j
              jjrow = jj + joff
	      jp1 = min(jjrow+1,jmt-1)-joff
              diff_fn(i,k,jj) = (T_j_obc(i,k,jj,n,1)
     &                          -T_j_obc(i,k,jj,n,0))
#ifdef tracer_horz_mix_const
     &                         *ah_csu_dyur
#endif
#ifdef tracer_horz_mix_var
     &                         *diff_cnt(i,k,jj)*csu_dyur(jjrow)
#endif
#if defined partial_cell
     &                         *dhtn(i,k,jj)
#endif
            enddo
          enddo
        enddo
c
c-----------------------------------------------------------------------
c       diffusive flux across bottom face of T cells
c       use constant vertical diffusion
c-----------------------------------------------------------------------
c
        nt2 = min(n,nsmix)
        do k=1,kbot-1
          do i=istrt,iend
            diff_fb(i,k,j) = diff_cbt(i,k,j,nt2)*dz_wtr(i,k,j)*
     &                       (t(i,k,j,n,tau) - t(i,k+1,j,n,tau))
          enddo
        enddo
c
c       for vertical diffusion of tracers:
c       set k=0 elements to reflect surface tracer flux.
c       set bottom level "kz" to reflect insulation condition.
c       for vertical advection of tracers:
c       set k=0 elements of vertical advection array to the rigid
c       lid condition (w(0) = 0). also assume w(km)=0. If "kz" <
c       "km" and there is a bottom slope, w(kz) may be non zero.
c
        do i=istrt,iend
           kb              = kmt(i,jrow)
           diff_fb(i,0,j)  = stf(i,j,n)
           diff_fb(i,kb,j) = btf(i,j,n)
        enddo
c
c
c-----------------------------------------------------------------------
c       calculate meridional diffusion and advection of tracers
c-----------------------------------------------------------------------
c
# ifdef obc_south
c
        if (jrow .eq. 2) then
c
c         radiation condition at the southern wall
c         and meridional diffusion (del**2) of tracer
c
          do k=1,km
            do i=istrt,iend
              ADV_Ty_obc(i,k) = -(c1s(i,k,n)+vad(i,k))
     &                  *(t(i,k,j+1,n,tau)-t(i,k,j,n,tau))/dyu(jrow)
            enddo
          enddo
        endif
# endif
# ifdef obc_north
c
#  ifdef obc_north_anywhere
        if (jrow .eq. jobcn) then
#  else
        if (jrow .eq. jmtm1) then
#  endif
c
c         radiation condition at the northern wall
c         and meridional diffusion (del**2) of tracer
c
          do k=1,km
            do i=istrt,iend
              ADV_Ty_obc(i,k) = -(c1n(i,k,n)+vad(i,k))
     &                  *(t(i,k,j,n,tau)-t(i,k,j-1,n,tau))/dyu(jrow-1)
            enddo
          enddo
        endif
# endif
c
c-----------------------------------------------------------------------
c       construct tracer source terms here
c       force n/s wall to observed values
c-----------------------------------------------------------------------
c
        do k=1,km
          do i=istrt,iend
            source(i,k,j) = c0
          enddo
        enddo
c
# ifndef orlanski
        call obcsponge1 (j, jrow, istrt, iend, n, t(1,1,1,1,tau),
     &                   vad, source)
# endif
c
c-----------------------------------------------------------------------
c       calculate the new tracer quantities allowing for implicit
c       treatment of vertical diffusion
c-----------------------------------------------------------------------
c
        do k=1,km
          twodt(k) = dtts*dtxcel(k)
          do i=istrt,iend
            t(i,k,j,n,taup1) = t(i,k,j,n,tau) + twodt(k)*
     &                        (DIFF_Tx(i,k,j) + DIFF_Tz(i,k,j)
     &                       + DIFF_Ty(i,k,j)
     &                       + ADV_Ty_obc(i,k) 
     &                       + source(i,k,j)
     &                        )*tmask(i,k,j)
#ifdef trace_obc
# include "cobc_trace.inc"
#endif
          enddo
        enddo
c
#if defined implicitvmix || defined isoneutralmix
c
c not tested!
c-----------------------------------------------------------------------
c       add dT/dt component due to implicit vertical diffusion
c-----------------------------------------------------------------------
c
        call ivdift (joff, j, j, istrt, iend, n, twodt)
#endif
c
        call setbcx (t(1,1,j,n,taup1), imt, km)
c
c-----------------------------------------------------------------------
c       construct diagnostics associated with tracer "n"
c-----------------------------------------------------------------------
c
        call diagt1 (joff, j, j, istrt, iend, n, twodt)
c
c-----------------------------------------------------------------------
c       end of tracer component "n" loop
c-----------------------------------------------------------------------
c
      enddo
c
c-----------------------------------------------------------------------
c     construct diagnostics (for total dT/dt)
c-----------------------------------------------------------------------
c
      idiag = 1
      call diagt2 (joff, j, j, istrt, iend, idiag)
c
#if !defined simple_sbc && !defined smhi_sbc
c
c-----------------------------------------------------------------------
c     if needed, construct the Atmos S.B.C.(surface boundary conditions)
c     averaged over this segment
c     eg: SST and possibly SSS
c-----------------------------------------------------------------------
c
        isst  = index_of_sbc ('sst')
        isss  = index_of_sbc ('sss')
        if (isst .ne. 0 .or. isss .ne. 0) then
          tlev = tau
          if (euler2) tlev = taum1
        call asbct (joff, js, je, istrt, iend, isst, isss
     &,               t(1,1,1,1,tlev))
        endif
#endif
c
#ifdef trace_indices
      write (stdout,'(2x,3(a,i4))')
     & "=> In   cobc: j=",j," joff=",joff," jrow=",j+joff
# endif
# ifdef timing
      call toc ('ocean', 'cobc')
# endif
      return
      end

# ifndef orlanski
      subroutine obcsponge1 ( j, jrow, is, ie, n, tm1, vad, source)
c
c=======================================================================
c     newtonian damping variables for obc regions. damp to prescribed
c     values if there are inflow conditions. data must be prepared
c     using the "mkobc" routines included in the programs for working
c     with the MOM dataset.
c
c     input:
c       j    = open boundary row in the MW
c       jrow = (2,jmtm1) open boundary row
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       n    = (1,2) = (T,S) tracer component
c       tm1  = tracer at "tau"
c       vad  = advective vel. for tracers at open walls
c            ==> restore to boundary values for inflow conditions 
c
c     output:
c       source = newtonian damping term
c
c     author:   Arne Biastoch   e-mail: abiastoch@ifm.uni-kiel.de
c=======================================================================
c
      use xobc_mod
#  include "param.h"
#  include "iounit.h"
#  include "obc_data.h"
#  include "switch.h"
#  include "tmngr.h"
c
      dimension source(imt,km,2:jmw-1)
      dimension tm1(imt,km,jmw,nt)
      dimension vad(imt,km)
c
#  ifdef timing
        call tic ('cobc', 'obc sponge1')
#  endif
      if (n .eq. 1 .and. is .eq. 2) then
c
c-----------------------------------------------------------------------
c       decide whether to read sponge data or not
c-----------------------------------------------------------------------
c
        begtim = position_within_data (model_time, obc_start_time
     &,                              obc_end_time, .true.)
        methodobc = 3
        call timeinterp (begtim, indxob1, tobc1, obc1dpm, 12, .true.
     &,           methodobc, inextdobc1, iprevdobc1, wprevobc1,
     &            readob1, inextobc1, iprevobc1)
c
c       read in the next data record from disk when needed
c
        if (readob1) then
          call getunit (ionew4, 'obc1.ud', opt_obc1)
          read (ionew4, rec=inextdobc1) obctnext, spdpmn, im, kk, jm
#  ifdef obc_south	  
     &,  obcs,(((obbuf_south(i,k,nn,inextobc1),i=1,imt),k=1,km),nn=1,2)
#  endif     
#  ifdef obc_north	  
     &,  obcn,(((obbuf_north(i,k,nn,inextobc1),i=1,imt),k=1,km),nn=1,2)
#  endif     
          write (stdout,'(/a,i3,a,i2,a,i2,a,g14.7,1x,a/)')
     &    '=> read obc1 record =',inextdobc1,' into buffer =',inextobc1
     &,   ' method #',methodobc,' at day =', begtim, stamp
     
        call relunit (ionew4)
        endif
      endif
      if (n .le. 2) then
c
c-----------------------------------------------------------------------
c       construct newtonian damping term using obc data
c-----------------------------------------------------------------------
c
#  ifdef obc_south
	  if (jrow .eq. 2 .and. obcs .ne. c0) then
	    tnext = c1-wprevobc1
	    do k=1,km
	      do i=is,ie
                if (vad(i,k) .eq. 0.) then
	           data = tnext    *obbuf_south(i,k,n,inextobc1) 
     &                  + wprevobc1*obbuf_south(i,k,n,iprevobc1)
                   source(i,k,j) = source(i,k,j) - 
     &                               obcs*(tm1(i,k,j,n) - data)
                endif
              enddo
	    enddo
	  endif
#  endif
c
#  ifdef obc_north
#   ifdef obc_north_anywhere
 	  if (jrow .eq. jobcn .and. obcn .ne. c0) then
#   else
 	  if (jrow .eq. jmtm1 .and. obcn .ne. c0) then
#   endif
	    tnext = c1-wprevobc1
	    do k=1,km
	      do i=is,ie
                if (vad(i,k) .eq. 0.) then
	           data = tnext    *obbuf_north(i,k,n,inextobc1) 
     &                  + wprevobc1*obbuf_north(i,k,n,iprevobc1)
                   source(i,k,j) = source(i,k,j) - 
     &                               obcn*(tm1(i,k,j,n) - data)
                endif
              enddo
	    enddo
	  endif
#  endif
      endif
#  ifdef timing
        call toc ('cobc', 'obc sponge1')
#  endif
      return
      end

# endif
c
#else
      subroutine cobc (joff, j, is, ie)
      return
      end
#endif
