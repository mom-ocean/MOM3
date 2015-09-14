#include "derived_options.h"
      subroutine addobcpsi
#if defined stream_function && defined obc
#ifndef orlanski
      use arrays_2d
      use xobc_mod
c
c-----------------------------------------------------------------------
c     add open boundary values to the solution (active open boundary)
c-----------------------------------------------------------------------
c
#  include "param.h"
#  include "emode.h"
#  include "iounit.h"
#  include "obc_data.h"
#  include "tmngr.h"
c
#if defined obc_south || defined obc_north
c
c-----------------------------------------------------------------------
c       decide whether to read obc data or not
c-----------------------------------------------------------------------
c
        begtim = position_within_data (model_time, obc_start_time
     &,                              obc_end_time, .true.)
        methodobc = 3
        call timeinterp (begtim, indxob1p, tobc1p, obc1pdpm, 12, .true.
     &,           methodobc, inextdobc1p, iprevdobc1p, wprevobc1p,
     &            readob1p, inextobc1p, iprevobc1p)
c
c       read in the next data record from disk when needed
c
        if (readob1p) then
          call getunit (ionew8, 'obcpsi1.ud', opt_obcpsi1)
          read (ionew8, rec=inextdobc1p) obctnext, spdpmn, im, kk, jm
# ifdef obc_south	  
     &,                 (psiwall_south(i,inextobc1p),i=1,imt)
# endif     
# ifdef obc_north	  
     &,                 (psiwall_north(i,inextobc1p),i=1,imt)
# endif     
          write (stdout,'(/a,i3,a,i2,a,i2,a,g14.7,1x,a/)')
     &    '=> read obcpsi1 record =',inextdobc1p,' into buffer ='
     &,   inextobc1p
     &,   ' method #',methodobc,' at day =', begtim, stamp
     
        call relunit (ionew8)
        endif
c	
        tnext = c1-wprevobc1p
# ifdef obc_south
      do i=1,imt
          data           = tnext     *psiwall_south(i,inextobc1p)
     &                   + wprevobc1p*psiwall_south(i,iprevobc1p)
          psi(i,1,1)     = data
          psi(i,2,1)     = data
      enddo
#  ifdef smooth_psi
c
c smooth in boundary data
c
      gama=0.15
      do i=2,imt-1
         do iawi=1,5
            psi(i,2+iawi,1) = (gama**iawi)*psi(i,2,1)+
     &            (1.-gama**iawi)*psi(i,2+iawi,1)
         enddo
      enddo
#  endif
# endif
# ifdef obc_north
#  ifdef obc_north_anywhere
      do i=iobcns,iobcne
          data             = tnext     *psiwall_north(i,inextobc1p)
     &                     + wprevobc1p*psiwall_north(i,iprevobc1p)
          psi(i,jobcn,1)   = data
          psi(i,jobcn+1,1) = data
      enddo
#  else
      do i=1,imt
          data           = tnext     *psiwall_north(i,inextobc1p)
     &                   + wprevobc1p*psiwall_north(i,iprevobc1p)
          psi(i,jmtm1,1) = data
          psi(i,jmt  ,1) = data
      enddo
#  endif
# endif
#endif
c
#if defined obc_west || defined obc_east
c
c-----------------------------------------------------------------------
c       decide whether to read obc data or not
c-----------------------------------------------------------------------
c
        begtim = position_within_data (model_time, obc_start_time
     &,                              obc_end_time, .true.)
        methodobc = 3
        call timeinterp (begtim, indxob2p, tobc2p, obc2pdpm, 12, .true.
     &,           methodobc, inextdobc2p, iprevdobc2p, wprevobc2p,
     &            readob2p, inextobc2p, iprevobc2p)
c
c       read in the next data record from disk when needed
c
        if (readob2p) then
          call getunit (ionew10, 'obcpsi2.ud', opt_obcpsi2)
          read (ionew10, rec=inextdobc2p) obctnext, spdpmn, im, kk, jm
# ifdef obc_west	  
     &,                 (psiwall_west(j,inextobc2p),j=1,jmt)
# endif     
# ifdef obc_east	  
     &,                 (psiwall_east(j,inextobc2p),j=1,jmt)
# endif     
          write (stdout,'(/a,i3,a,i2,a,i2,a,g14.7,1x,a/)')
     &    '=> read obcpsi2 record =',inextdobc2p,' into buffer ='
     &,   inextobc2p
     &,   ' method #',methodobc,' at day =', begtim, stamp
     
        call relunit (ionew10)
        endif
c	
      do j=1,jmt
	  tnext = c1-wprevobc2p
# ifdef obc_west
          data           = tnext     *psiwall_west(j,inextobc2p)
     &                   + wprevobc2p*psiwall_west(j,iprevobc2p)
          psi(1,j,1)     = data
          psi(2,j,1)     = data
# endif
# ifdef obc_east
          data           = tnext     *psiwall_east(j,inextobc2p)
     &                   + wprevobc2p*psiwall_east(j,iprevobc2p)
          psi(imtm1,j,1) = data
          psi(imt,  j,1) = data
# endif
      enddo
#endif
#if defined obc_west && defined obc_east
c
c-----------------------------------------------------------------------
c     In the following example there are two meridional open boundaries
c     with a netto tranport through the domain, e.g. the ACC flowing
c     through a South Atlantic model.
c
c     The two psi functions at i=2 and imt-1 are running from 0 in the 
c     north to the same in the south (different shapes allowed).
c     The southern land mass that connects these open boundaries has
c     to be set to the value of the netto throughflow. 
c     The index jpsimax in this case is the northernmost extension 
c     of the land mas.
c
c     For more infomation the user is referred to the MOM3 manual,
c     "Stevens Open Boundary Conditions"
c-----------------------------------------------------------------------
c
c      do jrow=1,jpsimax+1
c        do i=1,imt
c	  if (map(i,jrow) .ne. 0) then
c	    psi(i,jrow,1) = psi(1,1,1)
c	  endif
c	enddo
c      enddo
# if defined agape_p3xp3 || defined agape_1x1
      do jrow=1,jpsimax
        do i=3,imt-2
	  if (map(i,jrow) .ne. 0) then
	    psi(i,jrow,1) = psi(1,1,1)
	  endif
	enddo
      enddo
# endif
#endif
#endif
      return
      end
#else
      return
      end
#endif      


#if defined explicit_free_surface && defined obc
      subroutine addobceta(js,je,etarel)
#ifndef orlanski
      use arrays_2d
      use xobc_mod
c
c-----------------------------------------------------------------------
c     add open boundary values to the solution (active open boundary)
c-----------------------------------------------------------------------
c
#  include "param.h"
#  include "emode.h"
#  include "iounit.h"
#  include "obc_data.h"
#  include "tmngr.h"
      dimension etarel(imt,js:je)
c
#if defined obc_south || defined obc_north
c
c-----------------------------------------------------------------------
c       decide whether to read obc data or not
c-----------------------------------------------------------------------
c
        begtim = position_within_data (model_time, obc_start_time
     &,                              obc_end_time, .true.)
        methodobc = 3
        call timeinterp (begtim, indxob1p, tobc1p, obc1pdpm, 12, .true.
     &,           methodobc, inextdobc1p, iprevdobc1p, wprevobc1p,
     &            readob1p, inextobc1p, iprevobc1p)
c
c       read in the next data record from disk when needed
c
        if (readob1p) then
          call getunit (ionew8, 'obcpsi1.ud', opt_obcpsi1)
          read (ionew8, rec=inextdobc1p) obctnext, spdpmn, im, kk, jm
# ifdef obc_south	  
     &,                 (psiwall_south(i,inextobc1p),i=1,imt)
# endif     
# ifdef obc_north	  
     &,                 (psiwall_north(i,inextobc1p),i=1,imt)
# endif     
          write (stdout,'(/a,i3,a,i2,a,i2,a,g14.7,1x,a/)')
     &    '=> read obcpsi1 record =',inextdobc1p,' into buffer ='
     &,   inextobc1p
     &,   ' method #',methodobc,' at day =', begtim, stamp
     
        call relunit (ionew8)
        endif
      tnext = c1 - wprevobc1p
c	
# ifdef obc_south
      do i=2,imtm1
          etarel(i,2)    = tnext     *psiwall_south(i,inextobc1p)
     &                   + wprevobc1p*psiwall_south(i,iprevobc1p)
      enddo
# endif
# ifdef obc_north
#  ifdef obc_north_anywhere
      do i=iobcns,iobcne
          etarel(i,jobcn-1) = tnext     *psiwall_north(i,inextobc1p)
     &                   + wprevobc1p*psiwall_north(i,iprevobc1p)
      enddo
#  else
      do i=2,imtm1
          etarel(i,jmtm1) = tnext     *psiwall_north(i,inextobc1p)
     &                    + wprevobc1p*psiwall_north(i,iprevobc1p)
      enddo
#  endif
# endif
#endif
c
#if defined obc_west || defined obc_east
c
c-----------------------------------------------------------------------
c       decide whether to read obc data or not
c-----------------------------------------------------------------------
c
        begtim = position_within_data (model_time, obc_start_time
     &,                              obc_end_time, .true.)
        methodobc = 3
        call timeinterp (begtim, indxob2p, tobc2p, obc2pdpm, 12, .true.
     &,           methodobc, inextdobc2p, iprevdobc2p, wprevobc2p,
     &            readob2p, inextobc2p, iprevobc2p)
c
c       read in the next data record from disk when needed
c
        if (readob2p) then
          call getunit (ionew10, 'obcpsi2.ud', opt_obcpsi2)
          read (ionew10, rec=inextdobc2p) obctnext, spdpmn, im, kk, jm
# ifdef obc_west	  
     &,                 (psiwall_west(j,inextobc2p),j=1,jmt)
# endif     
# ifdef obc_east	  
     &,                 (psiwall_east(j,inextobc2p),j=1,jmt)
# endif     
          write (stdout,'(/a,i3,a,i2,a,i2,a,g14.7,1x,a/)')
     &    '=> read obcpsi2 record =',inextdobc2p,' into buffer ='
     &,   inextobc2p
     &,   ' method #',methodobc,' at day =', begtim, stamp
     
        call relunit (ionew10)
        endif
c	
      do j=2,jmtm1
	  tnext = c1-wprevobc2p
# ifdef obc_west
          etarel(2,j)    = tnext     *psiwall_west(j,inextobc2p)
     &                   + wprevobc2p*psiwall_west(j,iprevobc2p)
# endif
# ifdef obc_east
          etarel(imtm1,j) = tnext     *psiwall_east(j,inextobc2p)
     &                    + wprevobc2p*psiwall_east(j,iprevobc2p)
# endif
      enddo
#endif
#if defined obc_west && defined obc_east
c
c-----------------------------------------------------------------------
c     In the following example there are two meridional open boundaries
c     with a netto tranport through the domain, e.g. the ACC flowing
c     through a South Atlantic model.
c
c     The two psi functions at i=2 and imt-1 are running from 0 in the 
c     north to the same in the south (different shapes allowed).
c     The southern land mass that connects these open boundaries has
c     to be set to the value of the netto throughflow. 
c     The index jpsimax in this case is the northernmost extension 
c     of the land mas.
c
c     For more infomation the user is referred to the MOM3 manual,
c     "Stevens Open Boundary Conditions"
c-----------------------------------------------------------------------
c
c      do jrow=1,jpsimax+1
c        do i=1,imt
c	  if (map(i,jrow) .ne. 0) then
c	    psi(i,jrow,1) = psi(1,1,1)
c	  endif
c	enddo
c      enddo
#endif
#endif
      return
      end
#endif
