#ifdef time_varying_sbc_data
# if defined minimize_sbc_memory
      subroutine atmos (joff, js, je)
# else
      subroutine atmos
# endif
c
c-----------------------------------------------------------------------
c     this "atmos" model is an EXAMPLE of how to supply MOM with surface
c     boundary conditions which are a function of latitude, longitude,
c     and time (monthly mean data).
c     
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      use xctdbc_mod
      use sbc_info
      use sbc_ocean
      use timeinterp_mod
#include "param.h"
#include "ctdbc.h"
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"
c
#if defined minimize_sbc_memory
c
c-----------------------------------------------------------------------
c     Work only on rows js through je in the memory window
c-----------------------------------------------------------------------
c
#else
c
c-----------------------------------------------------------------------
c     Work on all rows assuming S.B.C. arrays are (imt,jmt)
c-----------------------------------------------------------------------
c
      joff = 0
      js   = 1
      je   = jmt
#endif
c
c-----------------------------------------------------------------------
c     determine the disk pointers, time weight interpolation factor,
c     and whether or not it is time to bring in new S.B.C. from disk
c     based on model time (days)
#if !defined minimize_sbc_memory
c     need to add "dt" to the model time because the call to 
c     atmos precedes the time stepping loop which calls mom, so the
c     model time has not yet been incremented when atmos executes.
#endif
c-----------------------------------------------------------------------
c
      if (joff .eq. 0) then
        do n=1,ntdbc
          a_time = model_time
#if !defined minimize_sbc_memory
            a_time = a_time + set_time(nint(dtts),0)
#endif
          daysbc(n) = position_within_data (a_time, data_start_time(n)
     &,                              data_end_time(n), period(n))
#ifdef debug_atmos
          write (stdout,*) ' atmos:n=',n,' daysbc(n)=',daysbc(n)
#endif
        enddo
      endif
c
c-----------------------------------------------------------------------
c     determine the disk pointers, time weight interpolation factor,
c     and whether or not it is time to bring in new S.B.C. from disk
c     based on the time (days) in MOM since dec 31, 1899 midnight.
c-----------------------------------------------------------------------
c
      if (joff .eq. 0) then
        do n=1,ntdbc
c
c         initialize the pointers (indices) to memory buffers and choose
c         interpolation method
c
          call timeinterp (daysbc(n), ntimeint(n)
     &,      tdrec(:,n), aprec(:,n)
     &,      ntdrec(n), period(n), method, inextd(n), iprevd(n)
     &,      wprev(n), rdtdbc(n), inextm(n), iprevm(n))
        enddo
      endif
c
c-----------------------------------------------------------------------
c     read in data for each S.B.C. when necessary
c     indices for "obc" which are used to relate how S.B.C. are
c     ordered within "sbcocn":
c
c     1 is for hellerman taux
c     2 is for hellerman tauy
c     3 is for levitus sst
c     4 is for levitus sss
c-----------------------------------------------------------------------
c
      do n=1,ntdbc
#if defined minimize_sbc_memory
c
c       read in the next data record and previous data record 
c       from disk for rows js through je
c
        call read_obc_data (iotdbc(n), iprevd(n), rstamp1
     &,                     obc(1,js,n,iprevm(n)), js, je, joff)
        call read_obc_data (iotdbc(n), inextd(n), rstamp2
     &,                     obc(1,js,n,inextm(n)), js, je, joff)
#else
        if (rdtdbc(n)) then
c
c         read in the next data record from disk
c
          call read_obc_data (iotdbc(n), inextd(n), rstamp2
     &,                       obc(1,js,n,inextm(n)), js, je, joff)
c
          write (stdout,'(/tr1,a,i2,a,i6,a,a,a,i3,a,a,a,a32)') 
     &      '==> Reading S.B.C. #',n, ' from rec ',inextd(n)
     &,     '[ending at ', rstamp2,'] into index ',inextm(n),'  for '
     &,     trim(nameout(n)), ' at MOM time ',stamp
          write (stdout,*)
     &   '  ...cksum =', checksum (obc(1,1,n,inextm(n)), imt*jmt, 1)
        endif
#endif
      enddo
c
c-----------------------------------------------------------------------
c     time dependent boundary condition section:
c
c     interpolate the monthly climatological values brought in
c     from disk to the time step level for MOM
c-----------------------------------------------------------------------
c
c
c     get indices related to S.B.C.
c     ordering of S.B.C. was set in driver.F 
c
      itaux = index_of_sbc ('taux')
      itauy = index_of_sbc ('tauy')
      ihflx = index_of_sbc ('heatflux')
      ipme  = index_of_sbc ('saltflux')
      isst  = index_of_sbc ('sst')
      isss  = index_of_sbc ('sss')
c
      n = 1
      wnext = c1-wprev(n)
      do j=js,je
        jrow = joff + j
        do i=1,imt
          sbcocn(i,jrow,itaux)  =  wprev(n)*obc(i,j,n,iprevm(n))
     &                             +  wnext*obc(i,j,n,inextm(n))
        enddo
      enddo
c
      n = 2
      wnext = c1-wprev(n)
      do j=js,je
        jrow = joff + j
        do i=1,imt
          sbcocn(i,jrow,itauy)  =  wprev(n)*obc(i,j,n,iprevm(n))
     &                            +   wnext*obc(i,j,n,inextm(n))
        enddo
      enddo
c
#ifdef restorst
c
c-----------------------------------------------------------------------
c     set flux for restoring surface tracers.
c     n=(3,4) restores model (SST,SSS) to levitus (SST,SSS) on
c     timescale "dampts" which is set via namelist.
c     note: set surface tracer flux to zero if not restoring. 
c-----------------------------------------------------------------------
c
      n = 3
      wnext = c1 - wprev(n)
      do j=js,je
        jrow = joff + j
        gamma = damp_factor(1)
        do i=1,imt
          datax           = wprev(n)*obc(i,j,n,iprevm(n))
     &                      +  wnext*obc(i,j,n,inextm(n))
          sbcocn(i,jrow,ihflx) = gamma*(datax - sbcocn(i,jrow,isst))
        enddo
      enddo
c
      n = 4
      wnext = c1 - wprev(n)
      do j=js,je
        jrow = joff + j
        gamma = damp_factor(2)
        do i=1,imt
          datax             = wprev(n)*obc(i,j,n,iprevm(n))
     &                        +  wnext*obc(i,j,n,inextm(n))
          sbcocn(i,jrow,ipme) = gamma*(datax - sbcocn(i,jrow,isss))
        enddo
      enddo
#else
      do j=js,je
        jrow = joff + j
        do i=1,imt
          sbcocn(i,jrow,ihflx)  = c0
          sbcocn(i,jrow,ipme)   = c0
        enddo
      enddo
#endif
c
c-----------------------------------------------------------------------
c     zero out SST and SSS accumulators for the next ocean segment
c-----------------------------------------------------------------------
c
      if (numasbc .ne. 0 .and. eots) then
        do n=1,numasbc
	  if (n .eq. 1) then
            m  = isst
	  elseif (n .eq. 2) then
            m  = isss
	  else
	    write (stdout,*) '=>Error: n=',n,' in atmos.F'
	    stop 'atmos.F'
	  endif
          do j=js,je
            jrow = joff + j
            do i=1,imt
              sbcocn(i,jrow,m) = c0
            enddo
          enddo
        enddo
      endif
c
      return
      end
#else
      subroutine nodataatmos
      return
      end
#endif
