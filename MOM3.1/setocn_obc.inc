c====================== add-in "setocn_obc.inc" ========================
c
c initialize several variables
c
      do n=1,nt
        do k=1,km
#if defined obc_south || defined obc_north
          do i=1,imt
            var1(i,k,n) = c0
# ifdef obc_south
            c1s (i,k,n) = c0
# endif
# ifdef obc_north
            c1n (i,k,n) = c0
# endif
          enddo
#endif
#if defined obc_west || defined obc_east
          do j=1,jmw
            var2(k,j,n) = c0
# ifdef obc_west
            c1w (k,j,n) = c0
# endif
# ifdef obc_east
            c1e (k,j,n) = c0
# endif
          enddo
#endif
        enddo
      enddo
c      
      do k=1,km
#if defined obc_south || defined obc_north
        do i=1,imt
          vad(i,k) = c0
        enddo
#endif
#if defined obc_west || defined obc_east
        do j=1,jmw
          uad(k,j) = c0
        enddo
#endif
      enddo
c
#ifdef orlanski
      do i=1,imt
# ifdef obc_south
         c1ps(i)=0.
# endif	   
# ifdef obc_north
         c1pn(i)=0.
# endif	   
      enddo
c      
      do j=1,jmt
# ifdef obc_west
         c1pw(j)=0.
# endif	   
# ifdef obc_east
         c1pe(j)=0.
# endif	   
      enddo
#else
c
c-----------------------------------------------------------------------
c     read data for newtonian damping term at inflow regions of the
c     open boundary conditions and interpolate in time.
c     schematically: damping = - obc(j) * (t(i,k,j) - levitus(i,k,j))
c     all related data is assumed to have been prepared using the
c     "obc" routines in the MOM dataset.
c     disk resource factor of 13 is for 12 months + 1 annual mean
c-----------------------------------------------------------------------
c
#if defined obc_south || defined obc_north
      lrec = 4 + 4
# ifdef obc_south      
      lrec = lrec + 1 + 2*imt*km
# endif      
# ifdef obc_north      
      lrec = lrec + 1 + 2*imt*km
# endif      
      write (opt_obc1,'(a,1x,i8)') 'unformatted direct words =',lrec
      write (stdout,'(/a,1pg10.3,a)')
     & ' Sequential access disk resource for file "obc1.mom" = '
     &,lrec*13*1.e-6,' (MW)'
      call getunit (ionew3, 'obc1.mom'
     &,               'unformatted sequential rewind ieee')
      write (stdout,'(/a,1pg10.3,a)')
     & ' Direct access disk resource for file "obc1" = '
     &,(lrec/8)*13*1.e-6,' (MW)'
      call getunit (ionew4, 'obc1', opt_obc1) 
c
      sum = 0.0
      cksum = 0.0
      do m=1,13
        read (ionew3)
        read (ionew3) obcstamp, obdpm, im, kk, jm, j1, j2, mm
# ifdef obc_south	
     &,    obcs,(((obbuf_south(i,k,n,1),i=1,imt),k=1,km),n=1,2)
# endif
# ifdef obc_north
     &,    obcn,(((obbuf_north(i,k,n,1),i=1,imt),k=1,km),n=1,2)
# endif
        write (ionew4, rec=m) obcstamp, obdpm, im, kk, jm
# ifdef obc_south	
     &,      obcs, (((obbuf_south(i,k,n,1),i=1,imt),k=1,km),n=1,2)
# endif     
# ifdef obc_north	
     &,      obcn, (((obbuf_north(i,k,n,1),i=1,imt),k=1,km),n=1,2)
# endif     
        if (m .le. 12) then
          obc1dpm(m) = obdpm
	  sum = sum + obdpm
	  tobc1(m) = sum - 0.5*obdpm
	endif
	do n=1,2
# ifdef obc_south	
	  cksum = cksum + checksum (obbuf_south(1,1,n,1), imt, km)
# endif	  
# ifdef obc_north	
	  cksum = cksum + checksum (obbuf_north(1,1,n,1), imt, km)
# endif	  
	enddo
      enddo
      print *,' checksum for obc1 data = ',cksum
      if (annlevobc) then
        write (stdout,'(/a)')
     &  ' => Annual mean data will be used for obc1 '
        do m=1,12
	  write (ionew4, rec=m) obcstamp, obdpm, im, kk, jm
# ifdef obc_south	  
     &,         obcs, (((obbuf_south(i,k,n,1),i=1,imt),k=1,km),n=1,2)
# endif     
# ifdef obc_north	  
     &,         obcn, (((obbuf_north(i,k,n,1),i=1,imt),k=1,km),n=1,2)
# endif     
	enddo
      else
        write (stdout,'(/a)')
     &  ' => Monthly levitus data will be used for obc1'
      endif
      write (stdout,'(1a/)')  ' Newtonian damping setup:'
# ifdef obc_south      
        if (obcs .ne. c0) then
          write (stdout,'(a,f7.2,a,e10.3,a,e10.3)') ' jrow=2, 
     &,         , yt(2)=',yt(2), ', obcs =',obcs
     &,   ', Newtonian time scale (days)=',secday/obcs
        endif
# endif
# ifdef obc_north
        if (obcn .ne. c0) then
          write (stdout,'(a,f7.2,a,e10.3,a,e10.3)') ' jrow=jmt-1,
     &,         , yt(jmt-1)=',yt(jmt-1), ', obcn =',obcn
     &,   ', Newtonian time scale (days)=',secday/obcn
        endif
# endif	
      write (stdout,*) ' '
c
c     set the current model time in days and initialize interpolation
c     information
c
      begtim = (realdays(initial) - 1.0) + realdays(imodeltime)
c
      if (.not.eqyear) then
        write (stdout,*) '=>Warning leap year being used with '
     &, ' climatological sponges?'
      endif
      iprevobc1  = 1
      inextobc1  = 2
      indxob1    = 12
      methodobc  = 3
      call timeinterp(begtim, indxob1, tobc1, obc1dpm, 12, .true., 
     &                methodobc, inextdobc1, iprevdobc1, wprevobc1,
     &                readob1, inextobc1, iprevobc1)
c
      read (ionew4, rec=iprevdobc1) obctprev, spdpmp, im, kk, jm
# ifdef obc_south	  
     &,  obcs, (((obbuf_south(i,k,n,iprevobc1),i=1,imt),k=1,km),n=1,2)
# endif     
# ifdef obc_north	  
     &,  obcn, (((obbuf_north(i,k,n,iprevobc1),i=1,imt),k=1,km),n=1,2)
# endif     
      read (ionew4, rec=inextdobc1) obctnext, spdpmn, im, kk, jm
# ifdef obc_south	  
     &,  obcs, (((obbuf_south(i,k,n,inextobc1),i=1,imt),k=1,km),n=1,2)
# endif     
# ifdef obc_north	  
     &,  obcn, (((obbuf_north(i,k,n,inextobc1),i=1,imt),k=1,km),n=1,2)
# endif     
c
      write (stdout,'(2(/a,i3,1x,a,a,i2)/a,g14.7,1x,a,a,g14.7/)')
     & ' reading obc1 record ', iprevdobc1, obctprev, 
     & ' into buffer ', iprevobc1
     &,' reading obc1 record ', inextdobc1, obctnext, 
     & ' into buffer ', inextobc1
     &,' for day =', begtim, stamp, ' weight =',wprevobc1
      write (stdout,'(/a,i3,a,i2/)') 'obc1 is dataset index ',indxob1
     &,' for time interpolation using method #',methodobc
      call relunit (ionew3)
      call relunit (ionew4)
c
#endif
#if defined obc_west || defined obc_east
      lrec = 4 + 4
# ifdef obc_west      
      lrec = lrec + 1 + 2*jmt*km
# endif      
# ifdef obc_east      
      lrec = lrec + 1 + 2*jmt*km
# endif      
      write (opt_obc2,'(a,1x,i8)') 'unformatted direct words =',lrec
      write (stdout,'(/a,1pg10.3,a)')
     & ' Sequential access disk resource for file "obc2.mom" = '
     &,lrec*13*1.e-6,' (MW)'
      call getunit (ionew5, 'obc2.mom'
     &,               'unformatted sequential rewind ieee')
      write (stdout,'(/a,1pg10.3,a)')
     & ' Direct access disk resource for file "obc2" = '
     &,(lrec/8)*13*1.e-6,' (MW)'
      call getunit (ionew6, 'obc2', opt_obc2) 
c
      sum = 0.0
      cksum = 0.0
      do m=1,13
        read (ionew5)
        read (ionew5) obcstamp, obdpm, im, kk, jm, j1, j2, mm
# ifdef obc_west
     &,    obcw,(((obbuf_west (j,k,n,1),j=1,jmt),k=1,km),n=1,2)
# endif
# ifdef obc_east
     &,    obce,(((obbuf_east (j,k,n,1),j=1,jmt),k=1,km),n=1,2)
# endif
        write (ionew6, rec=m) obcstamp, obdpm, im, kk, jm
# ifdef obc_west	
     &,      obcw, (((obbuf_west (j,k,n,1),j=1,jmt),k=1,km),n=1,2)
# endif     
# ifdef obc_east	
     &,      obce, (((obbuf_east (j,k,n,1),j=1,jmt),k=1,km),n=1,2)
# endif     
        if (m .le. 12) then
          obc2dpm(m) = obdpm
	  sum = sum + obdpm
	  tobc2(m) = sum - 0.5*obdpm
	endif
	do n=1,2
# ifdef obc_west	
	  cksum = cksum + checksum (obbuf_west (1,1,n,1), jmt, km)
# endif	  
# ifdef obc_east	
	  cksum = cksum + checksum (obbuf_east (1,1,n,1), jmt, km)
# endif	  
	enddo
      enddo
      print *,' checksum for obc data2 = ',cksum
      if (annlevobc) then
        write (stdout,'(/a)')
     &  ' => Annual mean data will be used for obc2 '
        do m=1,12
	  write (ionew6, rec=m) obcstamp, obdpm, im, kk, jm
# ifdef obc_west	  
     &,         obcw, (((obbuf_west (j,k,n,1),j=1,jmt),k=1,km),n=1,2)
# endif     
# ifdef obc_east	  
     &,         obce, (((obbuf_east (j,k,n,1),j=1,jmt),k=1,km),n=1,2)
# endif     
	enddo
      else
        write (stdout,'(/a)')
     &  ' => Monthly levitus data will be used for obc2'
      endif
      write (stdout,'(1a/)')  ' Newtonian damping setup:'
# ifdef obc_west      
        if (obcw .ne. c0) then
          write (stdout,'(a,f7.2,a,e10.3,a,e10.3)') '    i=2, 
     &,         , xt(2)=',xt(2), ', obcw =',obcw
     &,   ', Newtonian time scale (days)=',secday/obcw
        endif
# endif
# ifdef obc_east
        if (obce .ne. c0) then
          write (stdout,'(a,f7.2,a,e10.3,a,e10.3)') '    i=imt-1,
     &,         , xt(imt-1)=',xt(imt-1), ', obce =',obce
     &,   ', Newtonian time scale (days)=',secday/obce
        endif
# endif	
      write (stdout,*) ' '
c
c     set the current model time in days and initialize interpolation
c     information
c
      begtim = (realdays(initial) - 1.0) + realdays(imodeltime)
c
      if (.not.eqyear) then
        write (stdout,*) '=>Warning leap year being used with '
     &, ' climatological sponges?'
      endif
      iprevobc2  = 1
      inextobc2  = 2
      indxob2    = 13
      methodobc  = 3
      call timeinterp(begtim, indxob2, tobc2, obc2dpm, 12, .true., 
     &                methodobc, inextdobc2, iprevdobc2, wprevobc2,
     &                readob2, inextobc2, iprevobc2)
c
      read (ionew6, rec=iprevdobc2) obctprev, spdpmp, im, kk, jm
# ifdef obc_west	  
     &,  obcw, (((obbuf_west (j,k,n,iprevobc2),j=1,jmt),k=1,km),n=1,2)
# endif     
# ifdef obc_east
     &,  obce, (((obbuf_east (j,k,n,iprevobc2),j=1,jmt),k=1,km),n=1,2)
# endif     
      read (ionew6, rec=inextdobc2) obctnext, spdpmn, im, kk, jm
# ifdef obc_west	  
     &,  obcw, (((obbuf_west (j,k,n,inextobc2),j=1,jmt),k=1,km),n=1,2)
# endif     
# ifdef obc_east
     &,  obce, (((obbuf_east (j,k,n,inextobc2),j=1,jmt),k=1,km),n=1,2)
# endif     
c
      write (stdout,'(2(/a,i3,1x,a,a,i2)/a,g14.7,1x,a,a,g14.7/)')
     & ' reading obc2 record ', iprevdobc2, obctprev, 
     & ' into buffer ', iprevobc2
     &,' reading obc2 record ', inextdobc2, obctnext, 
     & ' into buffer ', inextobc2
     &,' for day =', begtim, stamp, ' weight =',wprevobc2
      write (stdout,'(/a,i3,a,i2/)') 'obc2 is dataset index ',indxob2
     &,' for time interpolation using method #',methodobc
      call relunit (ionew5)
      call relunit (ionew6)
#endif
c-----------------------------------------------------------------------
c     read stream function data for at inflow regions of the
c     open boundary conditions.
c     all related data is assumed to have been prepared using the
c     "obcpsi" routines in the MOM dataset.
c-----------------------------------------------------------------------
c
#if defined obc_south || defined obc_north
      lrec = 2 + 3 
# ifdef obc_south      
      lrec = lrec + 2*imt
# endif      
# ifdef obc_north      
      lrec = lrec + 2*imt
# endif      
      write (opt_obcpsi1,'(a,1x,i8)') 'unformatted direct words =',lrec
      write (stdout,'(/a,1pg10.3,a)')
     & ' Sequential access disk resource for file "obcpsi1.mom" = '
     &,lrec*13*1.e-6,' (MW)'
      call getunit (ionew7, 'obcpsi1.mom'
     &,               'unformatted sequential rewind ieee')
      write (stdout,'(/a,1pg10.3,a)')
     & ' Direct access disk resource for file "obcpsi1" = '
     &,(lrec/8)*13*1.e-6,' (MW)'
      call getunit (ionew8, 'obcpsi1', opt_obcpsi1) 
c
      sum = 0.0
      cksum = 0.0
      do m=1,13
        read (ionew7)
        read (ionew7) obcstamp, obdpm, im, kk, jm
# ifdef obc_south	
     &,               (psiwall_south(i,1),i=1,imt)
# endif     
# ifdef obc_north	
     &,               (psiwall_north(i,1),i=1,imt)
#  ifdef cme
c
c multiply with 1e12 to get cm3/s
c
        do i=1,imt
	  psiwall_north(i,1) = psiwall_north(i,1)*1.e12
	enddo
#  endif
# endif     
        write (ionew8, rec=m) obcstamp, obdpm, im, kk, jm
# ifdef obc_south	
     &,               (psiwall_south(i,1),i=1,imt)
# endif     
# ifdef obc_north	
     &,               (psiwall_north(i,1),i=1,imt)
# endif     
        if (m .le. 12) then
          obc1pdpm(m) = obdpm
	  sum = sum + obdpm
	  tobc1p(m) = sum - 0.5*obdpm
	endif
# ifdef obc_south	
	  cksum = cksum + checksum (psiwall_south(1,1), imt, 1)
# endif	  
# ifdef obc_north	
	  cksum = cksum + checksum (psiwall_north(1,1), imt, 1)
# endif	  
      enddo
      print *,' checksum for obcpsi1 data = ',cksum
      if (annlevobc) then
        write (stdout,'(/a)')
     &  ' => Annual mean data will be used for obcpsi1 '
        do m=1,12
	  write (ionew8, rec=m) obcstamp, obdpm, im, kk, jm
# ifdef obc_south	  
     &,               (psiwall_south(i,1),i=1,imt)
# endif     
# ifdef obc_north	  
     &,               (psiwall_north(i,1),i=1,imt)
# endif     
	enddo
      else
        write (stdout,'(/a)')
     &  ' => Monthly data will be used for obcpsi1'
      endif
      write (stdout,*) ' '
c
c     set the current model time in days and initialize interpolation
c     information
c
      begtim = (realdays(initial) - 1.0) + realdays(imodeltime)
c
      if (.not.eqyear) then
        write (stdout,*) '=>Warning leap year being used with '
     &, ' climatological sponges?'
      endif
      iprevobc1p  = 1
      inextobc1p  = 2
      indxob1p    = 14
      methodobc   = 3
      call timeinterp(begtim, indxob1p, tobc1p, obc1pdpm, 12, .true., 
     &                methodobc, inextdobc1p, iprevdobc1p, wprevobc1p,
     &                readob1p, inextobc1p, iprevobc1p)
c
      read (ionew8, rec=iprevdobc1p) obctprev, spdpmp, im, kk, jm
# ifdef obc_south	  
     &,               (psiwall_south(i,iprevobc1p),i=1,imt)
# endif     
# ifdef obc_north	  
     &,               (psiwall_north(i,iprevobc1p),i=1,imt)
# endif     
      read (ionew8, rec=inextdobc1p) obctnext, spdpmp, im, kk, jm
# ifdef obc_south	  
     &,               (psiwall_south(i,inextobc1p),i=1,imt)
# endif     
# ifdef obc_north	  
     &,               (psiwall_north(i,inextobc1p),i=1,imt)
# endif     
c
      write (stdout,'(2(/a,i3,1x,a,a,i2)/a,g14.7,1x,a,a,g14.7/)')
     & ' reading obcpsi1 record ', iprevdobc1p, obctprev, 
     & ' into buffer ', iprevobc1p
     &,' reading obcpsi1 record ', inextdobc1p, obctnext, 
     & ' into buffer ', inextobc1p
     &,' for day =', begtim, stamp, ' weight =',wprevobc1p
      write (stdout,'(/a,i3,a,i2/)') 'obcpsi1 is dataset index'
     &, indxob1p,' for time interpolation using method #',methodobc
      call relunit (ionew7)
      call relunit (ionew8)
c
#endif
#if defined obc_west || defined obc_east
      lrec = 2 + 3
# ifdef obc_west      
      lrec = lrec + 2*jmt
# endif      
# ifdef obc_east      
      lrec = lrec + 2*jmt
# endif      
      write (opt_obcpsi2,'(a,1x,i8)') 'unformatted direct words =',lrec
      write (stdout,'(/a,1pg10.3,a)')
     & ' Sequential access disk resource for file "obcpsi2.mom" = '
     &,lrec*13*1.e-6,' (MW)'
      call getunit (ionew9, 'obcpsi2.mom'
     &,               'unformatted sequential rewind ieee')
      write (stdout,'(/a,1pg10.3,a)')
     & ' Direct access disk resource for file "obcpsi2" = '
     &,(lrec/8)*13*1.e-6,' (MW)'
      call getunit (ionew10, 'obcpsi2', opt_obcpsi2) 
c
      sum = 0.0
      cksum = 0.0
      do m=1,13
        read (ionew9)
        read (ionew9) obcstamp, obdpm, im, kk, jm, jpsimax, psimax
# ifdef obc_west
     &,               (psiwall_west(j,1),j=1,jmt)
# endif     
# ifdef obc_east	
     &,               (psiwall_east(j,1),j=1,jmt)
# endif     
        write (ionew10, rec=m) obcstamp, obdpm, im, kk, jm
# ifdef obc_west
     &,               (psiwall_west(j,1),j=1,jmt)
# endif     
# ifdef obc_east	
     &,               (psiwall_east(j,1),j=1,jmt)
# endif     
        if (m .le. 12) then
          obc2pdpm(m) = obdpm
	  sum = sum + obdpm
	  tobc2p(m) = sum - 0.5*obdpm
	endif
# ifdef obc_west	
	  cksum = cksum + checksum (psiwall_west (1,1), jmt, 1)
# endif	  
# ifdef obc_east	
	  cksum = cksum + checksum (psiwall_east (1,1), jmt, 1)
# endif	  
      enddo
      print *,' checksum for obc data2 = ',cksum
c
      if (annlevobc) then
        write (stdout,'(/a)')
     &  ' => Annual mean data will be used for obc2 '
        do m=1,12
	  write (ionew10, rec=m) obcstamp, obdpm, im, kk, jm
# ifdef obc_west	  
     &,               (psiwall_west(j,1),j=1,jmt)
# endif     
# ifdef obc_east	  
     &,               (psiwall_east(j,1),j=1,jmt)
# endif     
	enddo
      else
        write (stdout,'(/a)')
     &  ' => Monthly data will be used for obcpsi2'
      endif
c
      write (stdout,*) ' '
c
c     set the current model time in days and initialize interpolation
c     information
c
      begtim = (realdays(initial) - 1.0) + realdays(imodeltime)
c
      if (.not.eqyear) then
        write (stdout,*) '=>Warning leap year being used with '
     &, ' climatological sponges?'
      endif
      iprevobc2p  = 1
      inextobc2p  = 2
      indxob2p    = 15
      methodobc   = 3
      call timeinterp(begtim, indxob2p, tobc2p, obc2pdpm, 12, .true., 
     &                methodobc, inextdobc2p, iprevdobc2p, wprevobc2p,
     &                readob2p, inextobc2p, iprevobc2p)
c
      read (ionew10, rec=iprevdobc2p) obctprev, spdpmp, im, kk, jm
# ifdef obc_west	  
     &,               (psiwall_west(j,iprevobc2p),j=1,jmt)
# endif     
# ifdef obc_east	  
     &,               (psiwall_east(j,iprevobc2p),j=1,jmt)
# endif     
      read (ionew10, rec=inextdobc2p) obctnext, spdpmp, im, kk, jm
# ifdef obc_west	  
     &,               (psiwall_west(j,inextobc2p),j=1,jmt)
# endif     
# ifdef obc_east	  
     &,               (psiwall_east(j,inextobc2p),j=1,jmt)
# endif     
c
      write (stdout,'(2(/a,i3,1x,a,a,i2)/a,g14.7,1x,a,a,g14.7/)')
     & ' reading obcpsi2 record ', iprevdobc2p, obctprev, 
     & ' into buffer ', iprevobc2p
     &,' reading obcpsi2 record ', inextdobc2p, obctnext, 
     & ' into buffer ', inextobc2p
     &,' for day =', begtim, stamp, ' weight =',wprevobc2p
      write (stdout,'(/a,i3,a,i2/)') 'obcpsi2 is dataset index'
     &, indxob2p,' for time interpolation using method #',methodobc
      call relunit (ionew9)
      call relunit (ionew10)
c
#endif
#endif
