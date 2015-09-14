#ifdef river_inflow
      subroutine check_river(kmt, xu, yu)
!-----------------------------------------------------------------------
!     Checks whether the geographical setup of rivers fits in the model
!     topography.
!     author:      M. Schmidt  e-mail=> martin.schmidt@io-warnemuende.de
!     with a lot of code taken from util.F
!-----------------------------------------------------------------------
# include"size.h"
# include"stdunits.h"

!     The position of rivers is checked and a map with the rivers
!     is shown
      integer kmt(imt,jmt)
      real xu(imt), yu(jmt)
      integer, dimension(:,:), allocatable :: map_river
      real, dimension(:), allocatable :: dxt
      real, dimension(:), allocatable :: cst, dyt
c
      logical error, warning
      character*80 levels
      character*80 lev1
# include"river.h"
      real, dimension(:), allocatable :: areariv
      allocate (map_river(imt,jmt))
      allocate (dxt(imt))
      allocate (dyt(jmt), cst(jmt))
      allocate (areariv(nriv))

      do j=1,jmt
        do i=1,imt
	  if (kmt(i,j).eq.0) then
	    map_river(i,j) = 75
	  else
	    map_river(i,j) = 78
	  endif
	enddo
      enddo
      pi = 4.0 * atan(1.0)
      degrad = pi / 180.0
      radius = 0.6371e9
      do i=2,imt-1
        dxt(i) = radius * degrad * (xu(i) - xu(i-1))
      end do
      do jrow = 2,jmt-1
        cst(jrow) = cos(0.5*(yu(jrow)+yu(jrow-1))*degrad)
        dyt(jrow) = radius * degrad * (yu(jrow) - yu(jrow-1))
      end do
c
      write (stdout,*)  ' '
      write (stdout,'(a)') 
     &     'The following rivers have been defined: '	
      do n=1,nriv 
        error = .FALSE.
        areariv(n)=0.0
        do m=1,nboxriv(n)
           ir = iriv(n,m)
           jr = jriv(n,m)
           if (kmt(ir,jr).lt.2) then
             write (stdout,'(a,i2,a,i2,a)') 
     &         'error => River box',m,' of river ',n,' is on land!' 
             write (stdout,'(a,i4,a,i4,a)') 
     &         '  ir = ',ir,'  jr = ',jr
             error = .TRUE.
	     map_river(ir,jr) = 80
           else
             areariv(n) = areariv(n)+dxt(ir)*cst(jr)*dyt(jr)
	     map_river(ir,jr) = n
           endif
	   if (ir.eq.1.or.ir.eq.imt.or.jr.eq.1.or.jr.eq.jmt) then 
              write (stdout,'(a,i2,a,i2,a)') 
     &         'error => River box',m,' of river ',n
     &         ,' touches the model boundary!' 
	  
	   endif
# ifdef cyclic
	   if (ir.eq.2.or.ir.eq.imt-1) then 
              write (stdout,'(a,i2,a,i2,a)') 
     &         'error => River box',m,' of river ',n
     &         ,' touches a cyclic the model boundary!' 	  
	   endif
# endif
           if(kmt(ir-1,jr).ne.0.and.kmt(ir+1,jr).ne.0.
     &        and.kmt(ir,jr-1).ne.0.and.kmt(ir,jr+1).ne.0.
     &        and.kmt(ir-1,jr-1).ne.0.and.kmt(ir-1,jr+1).ne.0.
     &        and.kmt(ir+1,jr-1).ne.0.and.kmt(ir+1,jr+1).ne.0)then
              warning = .true.
	   else
              warning = .false.
	   endif
	enddo
        if (.not.error) then
          write (stdout,'(a,i2,a,a,a,i2,a,50(2i4,2x))') 
     &          'river ',n,': ',rivername(n)
     &        ,' defined on ',nboxriv(n),' boxes at ' 		
     &        ,(iriv(n,m),jriv(n,m),m=1,nboxriv(n))
          write (stdout,'(a,1pe12.5,a)')
     &	   '                 surface river area: ', areariv(n),' cm**2'
        else
          write (stdout,'(a,a,i2,a,50(2i4,2x))') rivername(n)
     &        ,' is defined in ',nboxriv(n),' boxes at' 		
     &        ,(iriv(n,m),jriv(n,m),m=1,nboxriv(n))
          write(stdout,'(a,a)') 'error => Check river ',rivername(n)	  
	endif
        if (warning) then
          write (stdout,'(a,a,a)') 'warning => river ',rivername(n)	  
     &         ,'does not come from land!'
          write (stdout,'(a,a)') 
     &         '          This works but check river ',rivername(n)
	endif
      enddo
c
c     set character markers
c
#  ifndef skip_river_map
      write (stdout,'(/a/)')
     & 'Use -Dskip_river_map to eliminate the following map.'
      write (stdout,*)  ' '
      write (stdout,'(/,38x,132a)') 'Map of model rivers'
      write (stdout,*)  ' '
      lev1(1:51) = '1234567890ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmno'
      levels = lev1(1:51)//'pqrstuvwxyz!#$%<>[]{}().+*-=@'
      linewidth = 125
      nfirst = 1
      nlast  = nriv
      inc    = 3
      do n=nfirst,nlast,inc
        nn = n-nfirst+1
        if (n+inc .le. nlast) then
          jinc = inc
        else
          jinc = nlast-n+1
        endif
        write (stdout,'(6(a,i6,a,a,3x))') 
     &  ('River',j+nfirst-1,' is printed as ',levels(j:j)
     &         ,j=nn,nn+jinc-1)
      enddo
      write (stdout,'(a,a,3x)') 
     &   'Land               : ',levels(75:75)
     &  ,'normal ocean points: ',levels(78:78)
     &  ,'Errors             : ',levels(80:80)
      inc=124
      do l=0,imt,inc
        incr = min(inc,imt-l)
        write (stdout,8800) (l+i,i=1,incr,4)
        do  jj=1,jmt
          j  = jmt+1-jj
          write (stdout,8900) j, 
     &    (levels(map_river(l+i,j):map_river(l+i,j))
     &           ,i=1,incr) 
        enddo
      enddo
      return
8800  format (/, 2x, 31i4)
8900  format (1x,i3,1x, 124a1)
      write (stdout,*)  ' '
!!      call iplot(map_river, imt, imt, jmt )
#  else
      write (stdout,'(/a/)') 
     &      'Remove -Dskip_river_map to show the rivers.'
#  endif
      deallocate (map_river, cst, dxt, dyt, areariv)
      return
      end
#else
      subroutine check_river
      return
      end
#endif
