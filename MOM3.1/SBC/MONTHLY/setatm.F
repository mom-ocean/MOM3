#ifdef time_varying_sbc_data
      subroutine setatm (dtatm)
c
c-----------------------------------------------------------------------
c     set up the monthly climatological data which has been prepared
c     using the run scripts in PREP_DATA
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      use mom3_mpp_mod
      use xctdbc_mod
      use sbc_info
      use sbc_ocean
      use timeinterp_mod
# include "param.h"
# include "ctdbc.h"
# include "tmngr.h"
# include "switch.h"
      character*32 seqrecl
      dimension xxu(imt), yyu(jmt), obc_data(imt,jmt)
c
c     names of sequential unformatted input data (S.B.C.)
c
      namein(1) = 'taux.mom'
      namein(2) = 'tauy.mom'
      namein(3) = 'sst.mom'
      namein(4) = 'sss.mom'
c
c     names of direct access data  (S.B.C.)
c
      nameout(1) = 'taux.ud'
      nameout(2) = 'tauy.ud'
      nameout(3) = 'sst.ud'
      nameout(4) = 'sss.ud'
c
       write (stdout,'(/a/)')
     &   ' ==> S.B.C. for MOM are comming from time varying S.B.C. data'
c
c-----------------------------------------------------------------------
c     transfer the monthly climatological surface boundary conditions
c     from the sequential unformatted files (prepared by PREP_DATA)
c     to direct access files.
c     disk resource factor of 13 is for 12 months + 1 annual mean
c-----------------------------------------------------------------------
c
c     define how dataset is to be used: (periodic or not) and number
c     of records for each dataset
c
      do n=1,ntdbc
        ntdrec(n) = 12
	period(n) = .true.
        if (ntdrec(n) .gt. maxbcr) then
          write (stdout,'(/,1x,a)') 'Error: "ntdrec" exceeds "maxbcr"'
          stop '=>setatm'
        endif
      enddo
c
      ntotrec = 13     
c               13 is for 12 monthly means + 1 annual mean
c
# if defined cray_ymp || defined cray_c90 || defined cray_t90 
c
c     lrec needs room for character*32 "stamp" + "period" + jrow
c     + latitude row of data = (4 + 1 + 1 + imt) words
c     8 characters per 64bit word = 4 words per character*32 "stamp" 
c
      lrec = 6 + imt
      write (opt_tdbc,'(a,i10,a,i10,a)')
     &' unformatted direct words =',lrec,' sds maxsize words='
     &,lrec*jmt*ntotrec,' novfl'
# else
c
c     lrec needs room for character*32 "stamp" + "period" + jrow
c     + latitude row of data = 2*(4 + 1 + 1 + imt) words
c     4 characters per 32bit word = 8 words per character*32 "stamp"
c     the factor of "2" is to account for 32bit platforms using
c     implicit double precision (real*8) 
c
      lrec = 2*(6 + imt)
      write (opt_tdbc,'(a,i10)') ' unformatted direct words =',lrec
# endif
c
      do n=1,ntdbc
        write (stdout,'(/a,a,a,1pg10.3,a)')
     & ' Sequential access disk resource for file ', namein(n)(1:13)
     &,' = ', lrec*jmt*ntotrec*1.e-6,' (MW)'
        write (seqrecl,'(a,i12)') ' words=',imt*jmt + imt + jmt + 4 + 8
        
        call getunit (nu, namein(n)
     &,               'unformatted sequential rewind ieee'//seqrecl)
# ifdef ssread_sswrite
#  ifdef minimize_sbc_memory
        call ostart (iotdbc(n), nameout(n), imt*jmt*ntotrec, imt, 1)
#  else
        call ostart (iotdbc(n), nameout(n), imt*jmt*ntotrec, imt*jmt, 1)
#  endif
# else
        write (stdout,'(/a,a,a,1pg10.3,a)')
     & ' Direct access disk resource for file ', nameout(n)(1:13),' = '
     &,  lrec*jmt*ntdrec(n)*1.e-6,' (MW)'
        call getunit (iotdbc(n), nameout(n), opt_tdbc)
# endif
c
c       transfer the sequential input data to direct access files
c
        write (stdout,'(/a,a)') 'Initializing ',nameout(n)
        cksum = 0.0
        do m=1,ntdrec(n)
          if (n .eq. 1 .or. n .eq. 2) then
c
c           read taux & tauy
c
	    read (nu)
	    read (nu) dstamp(m,n), aprec(m,n), iimt, jjmt, nn 
     &,            xxu, yyu, obc_data
          else
c
c           read sst & sss
c
	    read (nu)
	    read (nu) dstamp(m,n), aprec(m,n), iimt, jjmt, k, zz, nn
     &,               xxu, yyu, obc_data
#ifdef salinity_psu            
            if (n .eq. 4) then
              do i=1,imt
                do j=1,jmt
                  if (obc_data(i,j) .ne. 0.0) then
                    obc_data(i,j) = obc_data(i,j)*1.e3 + 35.0
                  endif
                enddo
              enddo
            endif
#endif
	  endif
c
c         time manager cannot work with year 0 so add one to the year
c         (add 2 years to the average for december)
c
          call get_stamp (dstamp(m,n), iyr, imon, iday, ihr, imin, isec)
          if (iyr .eq. 0 .and. imon .eq. 1 .and. iday .eq. 1) then
	    iyr = iyr + 2
	  else
	    iyr = iyr + 1
	  endif
	  call set_stamp (dstamp(m,n), iyr, imon, iday, ihr, imin, isec)
c
          cksum = cksum + checksum (obc_data, imt, jmt)
c
          call write_obc_data (iotdbc(n), m, dstamp(m,n), aprec(m,n)
     &,                        obc_data, n) 
	enddo
c
        write (stdout,*)
     & ' checksum for monthly climatology = ', cksum, ' on ', nameout(n)
        call relunit (nu)

        call mpp_sync()

c
c       construct time structures to mark boundaries of dataset
c
        call get_stamp (dstamp(1,n), iyr, imon, iday, ihr, imin, isec)
        a_time = set_date (iyr, imon-1, iday, ihr, imin, isec)
        ndim  = days_in_month(a_time)
        a_time = set_date (iyr, imon, iday, ihr, imin, isec)
	data_start_time(n) = decrement_date (a_time, 0, 0, ndim, 0,0,0)
        call get_date (data_start_time(n),iyr, imon, iday,ihr,imin,isec)
        call set_stamp (rstamp1, iyr, imon, iday, ihr, imin, isec)
c
        k = ntdrec(n)
        call get_stamp (dstamp(k,n), iyr, imon, iday, ihr, imin, isec) 
        data_end_time(n) = set_date (iyr, imon, iday, ihr, imin, isec)
c
c       construct center of data records in days and show data bounds
c
        do m=1,ntdrec(n)
	  if (m .eq. 1) then
	    tdrec(1,n) = 0.5*aprec(1,n)
	  else
	    tdrec(m,n) = tdrec(m-1,n) + 0.5*(aprec(m-1,n) + aprec(m,n))
	    rstamp1 = dstamp(m-1,n)
	  endif
	  write (stdout,'(a,i3,a,f10.4,4a)') ' data record # ',m
     &,   ' is centered at ',tdrec(m,n), 'days. defined from '
     &,   rstamp1,' to ',dstamp(m,n)
	enddo
      enddo
c
# if defined minimize_sbc_memory
c
      do n=1,ntdbc 
c
c       express model time as number of days into dataset
c
        daysbc(n) = position_within_data (model_time, data_start_time(n)
     &,                              data_end_time(n), period(n))
#ifdef debug_atmos
        write (stdout,*) ' n=',n,' daysbc(n)=',daysbc(n)
#endif
c
c       initialize the pointers (indices) to memory buffers and choose
c       interpolation method
c
        iprevm(n) = 1
        inextm(n) = 2
	method    = 3
        ntimeint(n) = -1
c
c       find the disk indices and interpolation factor
c
        call timeinterp (daysbc(n), ntimeint(n) 
     &,      tdrec(:,n), aprec(:,n), ntdrec(n)
     &,      period(n), method, inextd(n), iprevd(n), wprev(n)
     &,      rdtdbc(n), inextm(n), iprevm(n))
      enddo
# else
      write (stdout,'(/a,a,a/)') '==> Note: to save memory when using'
     &,' "time_varying_sbc_data", enable "minimize_sbc_memory"'
     &,' which will use disk instead'
c
c-----------------------------------------------------------------------
c     determine time interpolation factors and disk record indices for 
c     interpolating time dependent surface boundary conditions to the
c     model time. the S.B.C. will be updated at the mid points of each 
c     month. For generality, each S.B.C. is treated separately.
c     This (for example) allows for a mixture of monthly heatflux with
c     daily winds.
c     the test case assumes, however, that all S.B.C are monthly
c
c     note: to understand how the time interpolator works, use
c           run_timeinterp.
c-----------------------------------------------------------------------
c
      do n=1,ntdbc 
c
c
c       express model time as number of days into dataset
c
        daysbc(n) = position_within_data (model_time, data_start_time(n)
     &,                              data_end_time(n), period(n))
#ifdef debug_atmos
        write (stdout,*) ' n=',n,' daysbc(n)=',daysbc(n)
#endif
c
c       initialize the pointers (indices) to memory buffers and choose
c       interpolation method
c
        iprevm(n) = 1
        inextm(n) = 2
	method    = 3
        ntimeint(n) = -1
c
c       find the disk indices and interpolation factor
c
        call timeinterp (daysbc(n), ntimeint(n)
     &,      tdrec(:,n), aprec(:,n), ntdrec(n)
     &,      period(n), method, inextd(n), iprevd(n), wprev(n)
     &,      rdtdbc(n), inextm(n), iprevm(n))
c
c       read data just before and after current model time
c
        js   = 1
	je   = jmt
	joff = 0
        call read_obc_data (iotdbc(n), iprevd(n), rstamp1
     &,                     obc(1,1,n,iprevm(n)), js, je, joff)
        call read_obc_data (iotdbc(n), inextd(n), rstamp2
     &,                     obc(1,1,n,inextm(n)), js, je, joff)
c
        write (stdout,8800) n, iprevd(n), rstamp1, iprevm(n) 
     &,                        inextd(n), rstamp2, inextm(n)
     &,                        nameout(n), stamp
        write (stdout,'(/24x,a12,a,i2,a,i2/)')
     &  trim(nameout(n)),' is dataset index ',n
     &, ' for time interpolation using method #', method
        write (stdout,*)
     &   '  ...1st cksum =', checksum (obc(1,1,n,iprevm(n)), imt*jmt, 1)
     &,  '  ...2nd cksum =', checksum (obc(1,1,n,inextm(n)), imt*jmt, 1)
      enddo
8800  format (/1x,'==> Reading S.B.C. #',i2,' from rec ',i6, ' at ', a
     &,' into index ',i3,/24x,' and rec ',i6,' at ', a, ' into index '
     &,i3,/,24x, ' for ',a12, ' at MOM time ',a) 
# endif
c
c-----------------------------------------------------------------------
c     zero all fluxes in case a zero flux S.B.C. is desired
c     "m" relates how S.B.C. "n" is ordered within "sbcocn"
c-----------------------------------------------------------------------
c
      do n=1,numosbc
        do jrow=1,jmt
          do i=1,imt
            sbcocn(i,jrow,n) = 0.0
          enddo
        enddo
      enddo
c
# ifdef coupled
c
c     verify that "coupled" is not enabled
c
      write (stdout,*) '==> Error: "time_varying_sbc_data" and '
     &, '"coupled" are incompatible'
      stop '=>setatm'
# else
c
c-----------------------------------------------------------------------
c     since this is not an atmospheric model, the atmospheric 
c     time step "dtatm" does not have to be set. mom.F will set it.
c-----------------------------------------------------------------------
c
# endif
# ifdef simple_sbc
      write (stdout,'(/,a,a,/)')
     &'==> Error: "time_varying_sbc_data" and '
     &,' "simple_sbc" are incompatible'
      stop '=>setatm'
# endif
c
c-----------------------------------------------------------------------
c     initialize euler2 because it is tested within atmos but not set
c     until the first call to subroutine mom
c-----------------------------------------------------------------------
c
      euler2 = .false.
      return
      end

      subroutine write_obc_data (iou, m, stamp, period, obc, n)
# if defined sponges && !defined sponges_old
      use xsponge_mod
# endif
# include "size.h"
# include "stdunits.h"
      dimension obc(imt,jmt), buf(imt)
      character*(*) stamp
# if defined sponges && !defined sponges_old
      if (n .eq. 1 .or. n .eq. 2) then
        do jrow = 2,j_south_spng-2
	  do i=1,imt
	    obc(i,jrow) = obc(i,j_south_spng-1)
	  enddo
          if (n .eq. 1 .and. m .eq. 1) then
	    write (stdout,*) '=>Note: resetting stress on row'
     &,     jrow,' to stress on row', j_south_spng-1
          endif
	enddo
        do jrow = j_north_spng+1,jmt-1
	  do i=1,imt
	    obc(i,jrow) = obc(i,j_north_spng)
	  enddo
          if (n .eq. 1 .and. m .eq. 1) then
	    write (stdout,*) '=>Note: resetting stress on row'
     &,     jrow,' to stress on row', j_north_spng
          endif
	enddo
      elseif (n .eq. 3 .or. n .eq. 4) then
        do jrow = 2,j_south_spng-1
	  do i=1,imt
	    obc(i,jrow) = obc(i,j_south_spng)
	  enddo
          if (n .eq. 3 .and. m .eq. 1) then
	    write (stdout,*) '=>Note: resetting SST & SSS on row'
     &,     jrow,' to SST & SSS on row', j_south_spng
          endif
	enddo
        do jrow = j_north_spng+1,jmt-1
	  do i=1,imt
	    obc(i,jrow) = obc(i,j_north_spng)
	  enddo
          if (n .eq. 3 .and. m .eq. 1) then
	    write (stdout,*) '=>Note: resetting SST & SSS on row'
     &,     jrow,' to SST & SSS on row', j_north_spng
          endif
	enddo
      endif
# endif
# ifdef ssread_sswrite
      call oput (iou, imt*jmt, m, obc)
# else
      do j=1,jmt
        nrec = (m-1)*jmt + j
        do i=1,imt
	  buf(i) = obc(i,j)
        enddo
        write (iou, rec=nrec) stamp, period, j, buf
      enddo
# endif
      write (stdout,'(a,i2,1x,a,a,g14.7,a)') ' month=',m, stamp
     &,  '.  averaged over preceding ', period, 'days'
      return
      end


      subroutine read_obc_data (iou, m, stamp, obc, js, je, joff)
# include "size.h"
      dimension obc(imt,js:je), buf(imt)
      character*(*) stamp
# ifdef ssread_sswrite
#  ifdef minimize_sbc_memory
      do j=js,je
        jrow = j + joff
        nrec = (m-1)*jmt + jrow
        call oget (iou, imt, nrec, obc(1,j))
      enddo
#  else
      call oget (iou, imt*jmt, m, obc)
#  endif
      stamp = "(not avail with ssread_sswrite) "
# else
      do j=js,je
        jrow = j + joff
        nrec = (m-1)*jmt + jrow
        read (iou, rec=nrec) stamp, period, jj, buf
        do i=1,imt
	  obc(i,j) = buf(i)
        enddo
      enddo
# endif
      return
      end


#else
      subroutine nodatasetatm
      return
      end
#endif
