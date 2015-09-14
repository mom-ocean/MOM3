      subroutine setriver
#ifdef river_inflow
      use river_mod
      use time_manager_mod
      use timeinterp_mod 
# include"size.h"
# include"stdunits.h"
# include"grdvar.h"
# include"taskrows.h"
      real, dimension(:), allocatable :: area
# ifndef simple_sbc
      real, dimension(:), allocatable :: tracriv
# endif
      character*20 filename
      integer  ioriv, mrecriv
      real transport, aprec
      type(time_type) :: a_time
      character*32 rstamp1
      logical rdriv, perriv
# include"river.h"
# include"scalar.h"
# include"switch.h"
# include"tmngr.h"
!
!-----------------------------------------------------------------------
!     set up the time dependent data of river runoff
!
!     author:      M. Schmidt  e-mail=> martin.schmidt@io-warnemuende.de
!     inspired by the atmos code of r. c. pacanowski
!-----------------------------------------------------------------------
!     River data are expected to be in one data file for each river.
!     The fresh water and tracer flux for each river is added in one or 
!     more surface cells. For simple parallel computation the river flux
!     is distributed to the river cells in setriver and is stored in a  
!     data structure named river_data of type  river_data_type
!
!      type, public :: river_data_type     
!         integer               :: index    ! identifies river
!         character*16          :: name     ! name of the river 
!         integer               :: ir,jr    ! zonal and meridional grid index
!         real                  :: area     ! area of all boxes of the river
!         integer               :: mrecriv  ! number of data records
!         logical               :: perriv   ! true for periodic data
!         integer               :: ninterp  ! data index for timeinterp
!         real, pointer         :: aprec(:) ! period lengths
!         real, pointer         :: tdrec(:) ! times of data records
!         type(time_type)       :: start_time, end_time
!         character*32, pointer :: dstamp(:)! time stamps record end 
!         real, pointer         :: rff(:)   ! river fresh water flux(time)
!         real, pointer         :: rtf(:,:) ! river tracer flux(tracer index, time)
!         character*20, pointer :: trname(:)! Name of the tracers ( tracer index)
!      end type river_data_type
!
!     The data are interpolated with timeinterp in module tinterp.F90.
!     The time steps are not fixed. Since each river has its own
!     data file and the data sets of the river model boxes are stored
!     in separate data structures, average flows, monthly data or data on
!     even shorter time intervalls can be processed simultaneously. 
!      
!-----------------------------------------------------------------------
!
!     river.F is called outside the ocean segment, i.e. with respect to 
!     time interpolation rivers are constant during the ocean time steps
!     => river data are interpolated to model_time + river_time_offset 
!     
!     river_time_offset=0 if river.F runs within the ocean time stepping
!
!-----------------------------------------------------------------------
      if (rununits .eq. 'days') then
	segsec=segtim*86400.
      elseif (rununits .eq. 'years') then
	segsec=segtim*yrlen*86400.
      elseif (rununits .eq. 'months' .and. eqmon) then
	segsec=segtim*30.0*86400.
      else
        write (stdout,'(a/a)')
     & '=>Error: specified integration time units = ',trim(rununits)
     &,'units must be days, equal months or years'
     & //' to transform segtim into seconds in setriver.F'
        call abort()
      endif
      river_time_offset = set_time(nint(0.5*(segsec-dtts)),0)
!-----------------------------------------------------------------------
!
!     Determine the number of rivers in the domain of the PE 
!
!-----------------------------------------------------------------------
      rivbox_in_pe = 0
      allocate(area(nriv))
      do nr=1,nriv
        area(nr)=0.0
        do nb=1, nboxmax
          ir = iriv(nr,nb)
          jr = jriv(nr,nb)
          if (jr.ge.1) area(nr) = area(nr) + dxt(ir)*cst(jr)*dyt(jr)
          if (jr.ge.jstask.and.jr.le.jetask) then
             rivbox_in_pe = rivbox_in_pe + 1
          endif
        enddo
      enddo
!!      write(stdout,'(a,i2,a,i2)') 'River boxes on PN ',pn,': '
!!     &    , rivbox_in_pe
      write(stdout,*) 'River boxes on PN',pn,':',rivbox_in_pe
!     Do the following only if the PE is responsible for a river
      if(rivbox_in_pe.gt.0) then
!-----------------------------------------------------------------------
!
!     allocate and initialize the required data sets for the river boxes
!     on the PE 
!
!-----------------------------------------------------------------------
        allocate (river_data(rivbox_in_pe))
        nri = 0            ! 
        do nr=1,nriv       ! Loop over all rivers
!         Write out the geometry, name and index of the river        
	  do nb=1, nboxriv(nr)
            ir = iriv(nr,nb)
            jr = jriv(nr,nb)
            if (jr.ge.jstask.and.jr.le.jetask) then
              nri = nri + 1
              river_data(nri)%ir    = ir
              river_data(nri)%jr    = jr
	      river_data(nri)%area  = area(nr)
	      river_data(nri)%name  = rivername(nr)
	      river_data(nri)%index = riverindex(nr)
! Calculate the transformation factor from river discharge to fresh
! water velocity.
! q = Q/area 
!     q   : freshwater velocity in cm/s
!     Q   : river discharge in cm**3/s
!     area: surface of all river boxes of one river       
            endif
          enddo
        enddo              ! End of loop over all rivers
        write(stdout,'(a,i2)')'The following boxes have been defined: ' 
        do nr=1,rivbox_in_pe
          write(stdout,'(a,i3,a,i3,a,a,a,i2,a,1pe12.5)')  
     &       'ir = ',river_data(nr)%ir  
     &       ,' jr = ',river_data(nr)%jr 
     &       ,' for river ',river_data(nr)%name  
     &       ,' index = ',river_data(nr)%index  
     &       ,' area = ',river_data(nr)%area  
        enddo
        deallocate(area)
!       Now go through the river boxes and find the corresponding data
!       This is done by each PE. So the data sets for all river boxes are
!       prepared.  
!
!-----------------------------------------------------------------------
!
!       Analyze the structure of the data sets 
!
!-----------------------------------------------------------------------
!
# ifndef simple_sbc
        allocate(tracriv(ntriv))
# endif	
!       allocate fields needed for time interpolation
        allocate(iprevdriv(rivbox_in_pe),inextdriv(rivbox_in_pe)
     &          ,iprevmriv(rivbox_in_pe),inextmriv(rivbox_in_pe))
        allocate(wprevriv(rivbox_in_pe))
        do nb=1,rivbox_in_pe
!         Rivers may be distributed over two PE's. So prepare each river
!         box separately. 
!         Generate the file name  
          filename=trim(river_data(nb)%name)//'.dat'
          write (stdout,'(/)')
          write (stdout,*)
     &    '=> Preparing river box ir =',river_data(nb)%ir
     &,   ', jr =',river_data(nb)%jr
     &,   ', reading data from ',filename
!         Get the number of data records and tracers
          call getunit (ioriv, filename, 'formatted sequential rewind')
          read (ioriv,'(a)') iotext
          read (ioriv,*) mrecriv
          read (ioriv,*) ntriver
          if(mrecriv.ne.0) then
            river_data(nb)%mrecriv = mrecriv           
            write(stdout,'(a,a,a,i4,a)')
     &          'The file ',trim(filename),' consists of ',mrecriv,
     &          ' data records'
          else
            write(stdout,'(a,a)')
     &          '=>Error! Something is wrong in file ',filename 
	   stop '=>setriver.F'
          endif
          write(stdout,'(a,a,a,i3,a)')
     &          'The file ',trim(filename),' contains data for',ntriver
     &          ,' tracers.'
# ifndef simple_sbc
          if (ntriver.lt.ntriv) then
            write(stdout,'(a,a,a,i3,a)')
     &          '=> Error! The river file ',trim(filename)
     &          ,' contains less tracers than required by ntriv =',ntriv
     &          ,' (see river.h)'
            call abort()
          elseif (ntriver.gt.ntriv) then
            write(stdout,'(a,a,a,i3,a)')
     &          '=> Warning! The river file ',trim(filename)
     &          ,' contains more tracers than required by ntriv =',ntriv
     &          ,' (see river.h)'
            write(stdout,'(a)')
     &          'Continue anyway assuming everything is well organized!' 
          endif
# endif
!         assume non-periodic river data	
          perriv = .false.
          river_data(nb)%perriv = perriv
!
!         Allocate the data fields 
!
          allocate (river_data(nb)%rff(mrecriv))          
          allocate (river_data(nb)%aprec(mrecriv))          
          allocate (river_data(nb)%tdrec(mrecriv))          
          allocate (river_data(nb)%dstamp(mrecriv))          
!
# ifndef simple_sbc
!         Here a lot of comfortable things may be organized later using 
!         tracer names.       
	  allocate (river_data(nb)%rtf(ntriv,mrecriv))
	  allocate (river_data(nb)%trname(ntriv))
          do n=1,ntriv
	    river_data(nb)%trname(n) = rivertrname(n) 
          enddo
!
# endif
!         check that a calendar is defined
!
          if ( get_calendar_type().eq.no_calendar) then
	    write(stdout,'(a/a/a)') '=> Error: No calendar type found!'
     &,     ' Time interpolation is not working without calendar.'
     &,     ' Define a calendar in the program calling setriver!'
            call abort()
	  endif 	    
          if ( get_calendar_type().eq.thirty_day_months ) then
            write(stdout,'(a,a)') '=> Warning: Calendar type is ' 
     &,     ' "thirty_day_months".'
            if(monthly_data) then
              write(stdout,'(a/a)') '=> Warning: Assuming monthly data.' 
     &,     ' The data record period will be set to 30. for all month.'
	    endif 	    
	  endif 	    
          if ( get_calendar_type().eq.julian ) then
            write(stdout,'(a,a)') '=> Warning: Calendar type is ' 
     &,     ' "julian".'
            if(monthly_data.and.perriv) then
              write(stdout,'(a/a)') 
     & 	    '=> Warning: Assuming monthly periodic data.' 
     &,     ' The data record period for february will be set to 28.25.'
	    endif 	    
            if(monthly_data.and.(.not.perriv)) then
              write(stdout,'(a/a)') 
     & 	    '=> Warning: Assuming monthly periodic data.' 
     &,     ' The data record period for february will be set to 29'
     &,     ' for leap years.'
	    endif 	    
	  endif 	    
          if ( get_calendar_type().eq.no_leap ) then
            write(stdout,'(a,a)') '=> Warning: Calendar type is ' 
     &,     ' "no_leap".'
            if(monthly_data) then
              write(stdout,'(a/a)') '=> Warning: Assuming monthly data.' 
     &,     ' The data record period for february will be set to 28.'
	    endif 	    
	  endif 	    
          if(.not.monthly_data) then
            write(stdout,'(a/a)') '=> Warning: No monthly data. ' 
     &,     ' The data record period will accepted without checks.'
	  endif 	    
!        
!         assign the input data to river box data structures
!
          read(ioriv,'(a)') iotext
          do m=1,mrecriv
! In future, time records should have stamp format
!	  read (ioriv) dstamp(m), aprec(m), data_buff
            read(ioriv,*) iyr, imon, iday, ihr, imin, isec
     &                  , river_data(nb)%aprec(m), transport
# ifndef simple_sbc
     &                  , (tracriv(n),n=1,ntriv)
# endif
# ifdef debug_river
            write(stdout,'(i4,5i3,a8,f8.4,a6,f11.0,5(a,i2.2,a,f7.3))')
     &                  iyr, imon, iday, ihr, imin, isec
     &                  , ' period=', river_data(nb)%aprec(m)
     &                  ,' d, Q = ',transport
#  ifndef simple_sbc
     &                  , (' tracer_',n,'=',tracriv(n),n=1,ntriv)
#  endif
# endif
!
!           time manager cannot work with year 0 so add one to the year
!
!           call get_stamp (dstamp(m,n), iyr, imon, iday, ihr, imin, isec)
            if (iyr .eq. 0.and.river_data(nb)%perriv) then
	      iyr = iyr + 1
	    endif
!	    
!           Adjust date for monthly data to the defined calendar type
!
            call set_stamp(river_data(nb)%dstamp(m)
     &                    ,iyr, imon, iday, ihr, imin, isec)
!
            river_data(nb)%rff(m) = transport/river_data(nb)%area          
# ifdef debug_river
            write (stdout,*) 'fresh water flux = ',river_data(nb)%rff(m)
     &                      ,' cm/s'
# endif
# ifndef simple_sbc
            do n=1,ntriv
     	      river_data(nb)%rtf(n,m) = tracriv(n)*river_data(nb)%rff(m)
# ifdef debug_river
              write (stdout,*) 
     &           trim(river_data(nb)%trname(n))//'flux = '
     &,          river_data(nb)%rtf(n,m)           ,' cm/s'
# endif
            enddo
# endif
	  enddo ! All records for one river are read and stored in river_data(nb)
          call relunit (ioriv)
!-----------------------------------------------------------------------
!
!         construct time structures to mark boundaries of dataset
!
!-----------------------------------------------------------------------
          call get_stamp (river_data(nb)%dstamp(1) 
     &                   ,iyr, imon, iday, ihr, imin, isec)
          a_time = set_date (iyr, imon, iday, ihr, imin, isec)
          river_data(nb)%start_time = a_time 
     &          - set_time(0,nint(river_data(nb)%aprec(1)))
          call get_date (river_data(nb)%start_time,iyr, imon, iday
     &                   , ihr, imin, isec)
          call set_stamp (rstamp1, iyr, imon, iday, ihr, imin, isec)
	  write (stdout,*) 
     &      ' Start time for river dataset',nb
     &,     ':',rstamp1
!       
          k = river_data(nb)%mrecriv
          call get_stamp (river_data(nb)%dstamp(k)
     &                   ,iyr, imon, iday, ihr, imin, isec) 
          river_data(nb)%end_time = 
     &                    set_date (iyr, imon, iday, ihr, imin, isec)
	  write (stdout,*) 
     &      ' End time for river dataset  ',nb
     &,     ':',river_data(nb)%dstamp(k)
!
!         construct center of data records in days and show data bounds
!
          do m=1,river_data(nb)%mrecriv
!	    
!           Adjust date for monthly data to the defined calendar type
!
            if(monthly_data) then
              if (river_data(nb)%aprec(m).gt.31.01.or.
     &            river_data(nb)%aprec(m).lt.28.0) then
	        write(stdout,'(a,i4,a,2pE12.3/a,a)') 
     &	         '=> Error: Record length in data set: ',m,' aprec = '
     &,           river_data(nb)%aprec(m)
     &,          'This is not a monthly data record. ' 
     &,          'Check the data set or set monthly_data to .false.! ' 
                call abort()	      
	      endif
              if ( get_calendar_type().eq.thirty_day_months ) then
	        write(stdout,'(a,i4,a,2pE12.3,a)') 
     &	         'Record length in data set: ',m,' aprec = '
     &,           river_data(nb)%aprec(m),' changed to 30 days.' 
                river_data(nb)%aprec(m) = 30.
	      endif
              if ( get_calendar_type().eq.julian ) then
                call get_stamp (river_data(nb)%dstamp(m) 
     &                   ,iyr, imon, iday, ihr, imin, isec)
                if (imon.eq.3) then
	          if (abs(river_data(nb)%aprec(m)-28.).lt.1.001) then
                    if (river_data(nb)%perriv) then
	              write(stdout,'(a,i4,a,2pE12.3/a,a)') 
     &		      'Record length in data set: ',m,' aprec = '
     &,               river_data(nb)%aprec(m)
     &,             'A periodic data set has been defined, changed to ' 
     &,               '28.25 days.' 
                      river_data(nb)%aprec(m) = 28.25
	            else
                      a_time = 
     &                       set_date (iyr, imon, iday, ihr, imin, isec)
	              write(stdout,'(a,i4,a,2pE12.3)') 
     &		        'Record length in data set: ',m,' aprec = '
     &,                 river_data(nb)%aprec(m)
                      if (leap_year(a_time)) then
	                write(stdout,'(a,a)') 
     &                   'A leap year has been detected,' 
     &,                   ' changed to 29 days.' 
                        river_data(nb)%aprec(m) = 29.
		      else
	                write(stdout,'(a,a)') 
     &                   'A non leap year has been detected,' 
     &,                   ' changed to 28 days.' 
                        river_data(nb)%aprec(m) = 28.
		      endif !leap_year(a_time)
		    endif ! perriv
	          else
	            write(stdout,'(a,i4,a,2pE12.3/a)') 
     &		    '=>Error in data set: ',m,'! For february aprec = '
     &,             river_data(nb)%aprec(m)
     &,            ' has been found. Check the file'//trim(filename)
                    call abort()    
		  endif   ! abs(river_data(nb)%aprec(m)-28.).lt.1.01
	        endif     ! imon.eq.2
	      endif       ! get_calendar_type.eq.julian
              if ( get_calendar_type().eq.no_leap ) then
                call get_stamp (river_data(nb)%dstamp(m) 
     &                   ,iyr, imon, iday, ihr, imin, isec)
                if (imon.eq.3) then
	          if (abs(river_data(nb)%aprec(m)-28.).lt.1.001) then
	            write(stdout,'(a,i4,a,2pE12.3,a,a)') 
     &		      'Record length in data set: ',m,' aprec = '
     &,               river_data(nb)%aprec(m),' changed to 28 days.' 
                    river_data(nb)%aprec(m) = 28
	          else
	            write(stdout,'(a,i4,a,2pE12.3/a)') 
     &		    '=>Error in data set: ',m,'! For february aprec = '
     &,             river_data(nb)%aprec(m)
     &,            ' has been found. Check the file'//trim(filename)
                    call abort()    
		  endif   ! abs(river_data(nb)%aprec(m)-28.).lt.1.01
	        endif     ! imon.eq.2
	      endif       ! get_calendar_type.eq.no_leap
	    endif
	    if(m.eq.1) then
	      river_data(nb)%tdrec(1) = 0.5*river_data(nb)%aprec(1)
            else
              river_data(nb)%tdrec(m) = river_data(nb)%tdrec(m-1) 
     &                + 0.5*(river_data(nb)%aprec(m-1) 
     &                + river_data(nb)%aprec(m))
              rstamp1 = river_data(nb)%dstamp(m-1)
	    endif
# ifdef debug_river
	    write (stdout,'(a,i3,a,f10.4,4a)') ' data record # ',m
     &,     ' is centered at ',river_data(nb)%tdrec(m)
     &      , ' days, defined from '
     &,     rstamp1,' to ',river_data(nb)%dstamp(m)
# endif
	  enddo
!
!-----------------------------------------------------------------------
!       determine time interpolation factors and disk record indices for 
!       interpolating time dependent surface boundary conditions to the
!       model time. T data will be updated at the mid points of each 
!       record. For generality, each data set is treated separately.
!       This (for example) allows for a mixture of river data with 
!       different time structure.
!
!       note: to understand how the time interpolator works, use
!             run_timeinterp.
!-----------------------------------------------------------------------
!
!
!
!         express model time as number of days into dataset
!
          dayriv = position_within_data (model_time
     &           , river_data(nb)%start_time
     &           , river_data(nb)%end_time, river_data(nb)%perriv)
# ifdef debug_river
          write (stdout,*) ' dayriv =',dayriv
# endif
c
c         initialize the pointers (indices) to memory buffers and choose
c         interpolation method
c
          iprevmriv(nb) = 2
          inextmriv(nb) = 1
	  methodriv    = 3
c
c         find the disk indices and interpolation factor
c
          nset = - 1
	  call timeinterp (dayriv, nset, river_data(nb)%tdrec
     &        , river_data(nb)%aprec, river_data(nb)%mrecriv
     &        , river_data(nb)%perriv
     &        , methodriv, inextdriv(nb), iprevdriv(nb), wprevriv(nb)
     &        , rdriv, inextmriv(nb), iprevmriv(nb))
          river_data(nb)%ninterp = nset
c
# ifdef debug_river
          write (stdout,*) 'Time interpolation for river '
     &        ,river_data(nb)%name
          write (stdout,*) 'Number of records '
     &        ,river_data(nb)%mrecriv
          write (stdout,*) 'Contens of field tdrec '
     &        ,river_data(nb)%tdrec
          write (stdout,*) 'Contens of field aprec '
     &        ,river_data(nb)%aprec
          if(river_data(nb)%perriv) then
	    write (stdout,*) 
     &        'The data set is periodic'
          else
	    write (stdout,*) 
     &        'The data set is not periodic'
	  endif
          write (stdout,*) 'Pointer to data records: '
          write (stdout,*) 'iprevdriv ',iprevdriv(nb)
          write (stdout,*) 'iprevmriv ',iprevmriv(nb)
          write (stdout,*) 'inextdriv ',inextdriv(nb)
          write (stdout,*) 'inextmriv ',inextmriv(nb)
          write (stdout,*) 'wprevriv  ',wprevriv(nb)
# endif
        enddo    ! End of main loop over all boxes
# ifndef simple_sbc
        deallocate(tracriv)
# endif
!-----------------------------------------------------------------------
!
!     prepare the first river dataset
!
!-----------------------------------------------------------------------
     
        allocate (river_rec(rivbox_in_pe))
        do nb=1,rivbox_in_pe   ! Loop over all river boxes
          allocate (river_rec(nb)%rff)
!         Write out the geometry, name and index of the river        
          river_rec(nb)%ir    = river_data(nb)%ir
          river_rec(nb)%jr    = river_data(nb)%jr
	  river_rec(nb)%area  = river_data(nb)%area
          river_rec(nb)%name  = river_data(nb)%name
	  river_rec(nb)%index = river_data(nb)%index
	  in = inextdriv(nb)
	  ip = iprevdriv(nb)
          river_rec(nb)%rff    = wprevriv(nb)*river_data(nb)%rff(ip)
     &                       +(1-wprevriv(nb))*river_data(nb)%rff(in)
# ifndef simple_sbc
	  allocate (river_rec(nb)%rtf(ntriv))
	  allocate (river_rec(nb)%trname(ntriv))
          do n=1,ntriv
	    river_rec(nb)%trname(n) = river_data(nb)%trname(n) 
            river_rec(nb)%rtf(n) = wprevriv(nb)*river_data(nb)%rtf(n,ip)
     &                       +(1-wprevriv(nb))*river_data(nb)%rtf(n,in)

          enddo
# endif
        enddo
      endif
#endif
      return
      end
