#ifdef makesponge
# define driver_only
# include "iomngr.F"
# include "grids.F"
# include "size_check.F"
# include "util.F"

      program sponge
c
c-----------------------------------------------------------------------
c     prepare sponge for damping T and S back to data
c     at specified latitudes and on specified time scales
c
c     author:  r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
# include "size.h"
# include "stdunits.h"
      character*32 stamp
      character*80 iotext
c
      character*12 tempin, saltin
      dimension tempin(13), saltin(13)
      dimension iotin(13), iosin(13)
c
# include "coord.h"
c
      data tempin/
     &               'jan.temp.mom', 'feb.temp.mom', 'mar.temp.mom'
     &,              'apr.temp.mom', 'may.temp.mom', 'jun.temp.mom'
     &,              'jul.temp.mom', 'aug.temp.mom', 'sep.temp.mom'
     &,              'oct.temp.mom', 'nov.temp.mom', 'dec.temp.mom'
     &,              'ann.temp.mom' /
c
      data saltin/
     &               'jan.salt.mom', 'feb.salt.mom', 'mar.salt.mom'
     &,              'apr.salt.mom', 'may.salt.mom', 'jun.salt.mom'
     &,              'jul.salt.mom', 'aug.salt.mom', 'sep.salt.mom'
     &,              'oct.salt.mom', 'nov.salt.mom', 'dec.salt.mom'
     &,              'ann.salt.mom' /
      logical save_row
      dimension spngn(jmt), spngs(jmt), spng(jmt)
      integer damp_to_sponge_row
      dimension damp_to_sponge_row(jmt)
      dimension save_row(jmt)
      dimension temps(imt,km), salts(imt,km)
      dimension tempn(imt,km), saltn(imt,km) 
      real, dimension (:,:,:), allocatable :: t_sponge, s_sponge
c
      stdout = 6
c
c-----------------------------------------------------------------------
c     initialize the i/o units
c-----------------------------------------------------------------------
c
      call ioinit
      do n=1,13
        call getunit (iotin(n), tempin(n)
     &,               'unformatted sequential rewind ieee')
        call getunit (iosin(n), saltin(n)
     &,               'unformatted sequential rewind ieee')
      enddo
      call getunit (iosp, 'sponge.mom'
     &,             'unformatted sequential rewind ieee')
c
c-----------------------------------------------------------------------
c     construct the grid.
c-----------------------------------------------------------------------
c
      write (stdout,'(//1x,a,//)') 'Construct the MOM grid'
c
      maxlen = max(imt,jmt,km)
      call gcoord (maxlen, imt2, jmt2, km2, dxtdeg, dytdeg, dxudeg
     &,            dyudeg, dzt, dzw, xt, xu, yt, yu, zt, zw)
c
      call size_check (imt2, jmt2, km2, 'sponge', 'stop')
c
c-----------------------------------------------------------------------
c     construct sponge
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
        spngn(jrow) = 0.0
        spngs(jrow) = 0.0
        save_row(jrow) = .false.
	damp_to_sponge_row(jrow) = 1
      enddo
c
c     set sponge timescale (days) next to artificial walls
c
      sponge_timescale   = 1.0
c
      num_rows_south = 0
# if !defined no_southern_sponge
c
c     define southern sponge zone
c     USER INPUT: set the latitude rows (defined by latitudes y1
c                 through y2) on the T-grid
c                 and the latitude row (jrow_s) for damping
c                 the solution towards levitus data. All rows in the
c                 sponge will be damped to levitus data on
#  if defined sponges_use_all_rows
c                 each row
#  else
c                 row jrow_s
#  endif
c                 using a linearly interpolated e-folding time scale
c                 from "amp1" at "y1" to "amp2" at "y2".
c
      y1     = yt(2)
      y2     = -10.5
      y2     = yt(indp(y2, yt, jmt))
#  if defined sponges_old
      jrow_s = indp (y1, yt, jmt)
      js     = jrow_s
      je     = indp(y2, yt, jmt)
#  else
      jrow_s = indp (y2, yt, jmt)
      js     = 2
      je     = jrow_s
#  endif
c
      amp1  = sponge_timescale
      amp2  = 50.0*sponge_timescale
      call lintrp (y1, amp1, y2, amp2, yt, spngs, jmt)
c
      write (stdout,*) 'southern sponge is defined from latitude = ',y1
     &,                ' to ',y2,',   jrow=',js,' to ',je
#  if defined sponges_use_all_rows
      write (stdout,*) 'solution will be damped to data on each row'
#  else
      write (stdout,*) 'solution will be damped to data on row =',jrow_s
#  endif
c
      do jrow=js,je
        spngs(jrow)              = 1.0/(spngs(jrow)*86400.0)
#  if defined sponges_use_all_rows
        save_row(jrow)           = .true.      
        damp_to_sponge_row(jrow) = jrow-js+1
	num_rows_south           = num_rows_south + 1
#  else
        damp_to_sponge_row(jrow) = 1
        num_rows_south           = 1
        save_row(jrow_s)         = .true.      
#  endif
      enddo
# endif
# if !defined no_northern_sponge
c
c     define northern sponge zone
c     USER INPUT: set the latitude rows (defined by latitudes y1
c                 through y2) on the T-grid
c                 and the latitude row (jrow_n) for damping
c                 the solution towards levitus data. All rows in the
c                 sponge will be damped to levitus data on jrow_n
c                 using a linearly interpolated e-folding time scale
c                 from "amp1" to "amp2".
c
      y1     = yt(jmt-1)
      y2     = 10.5
      y2     = yt(indp(y2, yt, jmt))
#  if defined sponges_old
      jrow_n = indp (y1, yt, jmt)
      js     = indp(y2, yt, jmt)
      je     = jrow_n
#  else
      jrow_n = indp (y2, yt, jmt)
      js     = jrow_n
      je     = jmt-1
#  endif
c
      amp1  = sponge_timescale
      amp2  = 50.0*sponge_timescale
      call lintrp (y1, amp1, y2, amp2, yt, spngn, jmt)
c
      write (stdout,*) 'northern sponge is defined from latitude = ',y1
     &,                ' to ',y2,',   jrow=',js,' to ',je
#  if defined sponges_use_all_rows
      write (stdout,*) 'solution will be damped to data on each row'
#  else
      write (stdout,*) 'solution will be damped to data on row =',jrow_n
#  endif
      do jrow=js,je
        spngn(jrow)              =  1.0/(spngn(jrow)*86400.0)
#  if defined sponges_use_all_rows
        save_row(jrow)           = .true.      
        damp_to_sponge_row(jrow) = jrow-js+1 + num_rows_south
#  else
        damp_to_sponge_row(jrow) = num_rows_south + 1
        save_row(jrow_n)         = .true.      
#  endif
      enddo
# endif
c      
c
c     combine northern and southern damping regions into one
c
      do jrow=1,jmt
        if (spngn(jrow) .ne. 0.0 .and. spngs(jrow) .ne. 0.0) then
	  write (stdout,*) '=>Error: sponge conflict. jrow=',jrow
	  stop
	else
	  spng(jrow) = spngn(jrow) + spngs(jrow)
	endif
      enddo
c
c     show damping factors & which sponge row is associated with jrow
c
      write (stdout,'(/a/)') ' Sponge damping factors (1/sec):'
      max_rows = 0
      do jrow=1,jmt
        write (stdout,'(1x,a,i4,a,f10.5,a,e14.7,a, a, i4)')
     &  'jrow=',jrow,' lat=',yt(jrow),' spng=',spng(jrow),' (1/sec)'
     &, ', damp_to_sponge_row=', damp_to_sponge_row(jrow)
        if (spng(jrow) .ne. 0.0) then
	  write (stdout,'(10x,a,f10.2,a)')
     & '  1/spng = ',1.0/(86400.0*spng(jrow)), 'days.'
        endif
	if (save_row(jrow)) then
	  max_rows = max_rows + 1
	  write (stdout,'(a,i4,a,i4)')
     & ' Saving jrow =',jrow,' as damping row =',max_rows
	endif
      enddo
c
      allocate (t_sponge(imt,km,max_rows))
      allocate (s_sponge(imt,km,max_rows))
c
c     collect sponge rows
c
      do month=1,13
        num_row = 0
        do jrow=1,jmt
          if (save_row(jrow)) then
c
c           read data for jrow
c
            read (iotin(month)) 
            read (iotin(month)) stamp, dpm, im, kk, js, alt, m, xt, zt
     &,                         temps
c
            read (iosin(month))
            read (iosin(month)) stamp, dpm, im, kk, js, alt, m, xt, zt
     &,                         salts
c
            num_row = num_row + 1
	    do k=1,km
	      do i=1,imt
	        t_sponge(i,k,num_row) = temps(i,k)
# if defined salinity_psu                
	        s_sponge(i,k,num_row) = salts(i,k)*1000.0+35.0
# else
	        s_sponge(i,k,num_row) = salts(i,k)
# endif
	      enddo
	    enddo
	  else
c
c           skip jrow
c
            read (iotin(month))
            read (iotin(month))
c
	    read (iosin(month))        
	    read (iosin(month))        
	  endif
        enddo
        if (num_row .ne. max_rows) then
	  write (stdout,*) '=> Error: change parameter max_rows to '
     &,                        num_row
          stop
	endif
c
c       save the sponge data
c
        iotext =
     &'dim ss(j),sn(j),d(i,k,4);read(n)stmp,avg,i,k,j,js,jn,mon,ss,sn,d'
        write (iosp) iotext
        write (iosp) stamp, dpm, imt, km, jmt, max_rows, month, spng
     &,              damp_to_sponge_row,  t_sponge, s_sponge
        write (stdout,'(1x,a,i3,a,a,a,g14.7,a)')
     &   'Writing sponges for ',num_row,' latitude rows from ',stamp
     &,  ' averaged over previous ',dpm, ' days'
      enddo
      deallocate (t_sponge, s_sponge)
      stop
      end

      subroutine lintrp (y1, t1, y2, t2, y, a, len)
c
c     linearly interpolate between t1 at location y1 and t2 at 
c     location y2 for intermediate points within "a" 
c
c     y = array of coordinate locations 
c     a = interpolated values between t1 & t2 at indices between
c         coordinate locations corresponding to y1 and y2
c
# include "stdunits.h"
      dimension y(len), a(len)
c
      j1 = indp (y1, y, len)
      j2 = indp (y2, y, len)
      if (j1 .eq. j2) then
        write (stdout,*) 'Error in lintrp: j1=',j1,', j2=',j2
        stop '=>lintrp'
      endif 
      a(j1) = t1
      a(j2) = t2
      if (abs(j1-j2) .gt. 1) then
        if (j1 .lt. j2) then
          inc = 1
        else
          inc = -1
        endif
        do j=j1+inc,j2-inc,inc
          a(j) = ((y(j)-y(j1))*t2 + (y(j2)-y(j))*t1)/(y(j2) - y(j1))     
        enddo
      endif
      return
      end
#else
      subroutine sponge
      return
      end
#endif

