#ifdef makeimport
# include "iomngr.F"
      program import
c
c-----------------------------------------------------------------------
c     read the "formatted" MOM database and install it "unformatted"
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      call ioinit
      call top
      call stres
      call air
      call levit
      call showunits
      call showfiles
      stop
      end
      subroutine top
c
c-----------------------------------------------------------------------
c     read formatted Scripps 1 deg topography data & write unformatted
c-----------------------------------------------------------------------
c
      parameter (it=362, jt=180, jtp2=jt+2)
      dimension d(it,jt), cxd(it), cyd(jt), dd(it,jtp2), cydd(jtp2)
      character*80 iotext
c
      call getunit (nrdata, 'scripps.fmt'
     &,               'formatted sequential rewind')
      call getunit (nwdata, 'scripps.top'
     &,               'unformatted sequential rewind ieee')
      read (nrdata,9000) lons, lats, cxd, cyd, d
      print *,'  read Scripps topo:'
      print *,' lons=',lons, ' lats=',lats
      print *,' d(180,90) =',d(180,90)
      lats = jtp2
      do j=1,jt
        cydd(j+1) = cyd(j)
      enddo
      cydd(1) = cydd(2) - (cydd(3) - cydd(2))
      cydd(jtp2) = cydd(jt+1) + (cydd(jt+1) - cydd(jt))
      do j=1,jt
        do i=1,it
	  dd(i,j+1) = d(i,j)
	enddo
	dd(1,j+1) = dd(it-1,j+1)
	dd(it,j+1) = dd(2,j+1)
      enddo
      do i=1,it
	dd(i,1)    = dd(i,2)
	dd(i,jtp2) = dd(i,jtp2-1)
      enddo 
      print *,' expanded  topo'
      print *,' dd(180,91) =',dd(180,91)
c
      iotext =
     & 'dim x(i), y(j), d(i,j); read (nu) i, j, x, y, d'
      print *, iotext
      write (nwdata) iotext      
      write (nwdata)       lons, lats, cxd, cydd, dd
      write (*,'(/a/)') ' wrote Scripps 1 deg Topography '
c
      return
9000  format (2i3/(10e14.7))
      end


      subroutine stres
c
c-----------------------------------------------------------------------
c     read  formatted Sol Hellerman's monthly wind stress climatology 
c     write it unformatted
c-----------------------------------------------------------------------
c
      parameter (ib=182, jb=90, jbp2=jb+2)
      character*32 stamps, rstamp
      character*80 iotext
      dimension stamps(12), dpm(12)
c
      common /tau/ months, lons, lats, cxw(ib), cyw(jb)
      common /tau/ clmtau(ib,jb,12,2)
      dimension taux(ib,jb), tauy(ib,jb)
      dimension tx(ib,jbp2), ty(ib,jbp2), cywp2(jbp2)
c
      data stamps /
     &               'm/d/y= 2/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 3/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 4/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 5/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 6/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 7/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 8/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 9/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y=10/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y=11/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y=12/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 1/ 1/0000, h:m:s= 0: 0: 0'
     &            /
c
      data dpm /31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0
     &, 30.0, 31.0, 30.0, 31.0/
c
      call getunit (nrdata, 'sols.fmt'
     &,               'formatted sequential rewind')
      read (nrdata,9100) months, lons, lats, cxw, cyw, clmtau
      print *,'  read Hellerman stress:'
      print *,' months=',months,' lons=',lons, ' lats=',lats
c
      call getunit (nwdata, 'trash'
     &,               'unformatted sequential rewind ieee')
      call getunit (nw, 'hellerman.tau'
     &,               'unformatted sequential rewind ieee')
      do j=1,jb
        do i=1,ib
	  taux(i,j) = 0.0
	  tauy(i,j) = 0.0
	enddo
      enddo
      sum = 0.0
      do n=1,12
	sum = sum + dpm(n)
        write (nwdata) stamps(n), dpm(n), lons, lats
     &, cxw, cyw, (((clmtau(i,j,n,k),i=1,ib),j=1,jb),k=1,2)
        print *,' n=',n,' taux(90,45) =',clmtau(90,45,n,1)
        do j=1,jb
	  do i=1,ib
	    taux(i,j) = taux(i,j) + clmtau(i,j,n,1)*dpm(n)
	    tauy(i,j) = tauy(i,j) + clmtau(i,j,n,2)*dpm(n)
	  enddo
	enddo
        print *,' wrote hellerman.tau for = ',stamps(n), ' days ='
     &,  dpm(n)
      enddo
c
      rsum = 1.0/sum
      do j=1,jb
        do i=1,ib
	  taux(i,j) = taux(i,j)*rsum
	  tauy(i,j) = tauy(i,j)*rsum
	enddo
      enddo
c
      write (nwdata) stamps(12), sum, lons, lats, cxw, cyw, taux, tauy
      print *,' annual mean taux(90,45) =',taux(90,45)
      print *, ' 13th record is annual mean'
      print *,' wrote hellerman.tau for = ',stamps(12), ' days ='
     &,  sum
c
c     espand meridional dimension to cover poles
c
      rewind nwdata
      do n=1,13
        read (nwdata) rstamp, avg, lons, lats, cxw, cyw, taux, tauy
	do j=1,jb
	  cywp2(j+1) = cyw(j)
	  do i=1,ib
	    tx(i,j+1) = taux(i,j) 
	    ty(i,j+1) = tauy(i,j)
	  enddo
	  tx(1,j+1)  = tx(ib-1,j+1)
	  tx(ib,j+1) = tx(2,j+1)
	  ty(1,j+1)  = ty(ib-1,j+1)
	  ty(ib,j+1) = ty(2,j+1)
	enddo
	cywp2(1)    = cyw(1) - (cyw(2) - cyw(1))
	cywp2(jbp2) = cyw(jb) + (cyw(jb) - cyw(jb-1))
	lats = jbp2
	do i=1,ib
	  tx(i,1)    = tx(i,2)
	  tx(i,jbp2) = tx(i,jbp2-1)
	  ty(i,1)    = ty(i,2)
	  ty(i,jbp2) = ty(i,jbp2-1)
	enddo 
        print *,' expanded n=',n,' taux(90,46) =',tx(90,46)
        iotext =
     &'dim x(i),y(j),tx(i,j),ty(i,j);read(n)stamp,avg,i,j,mon,x,y,tx,ty'
        print *, iotext
        write (nw) iotext
        write (nw) rstamp, avg, lons, lats, n, cxw, cywp2, tx, ty
      enddo
c
      print *,' cyw='
      write (*,'(1x,10g12.5)') cyw            
      print *,' cywp2='
      write (*,'(1x,10g12.5)') cywp2            
      return
9100  format (3i3/(10e14.7))
      end



      subroutine air
c
c-----------------------------------------------------------------------
c     read in formatted Oort's monthly air temperature climatology and
c     write it unformatted
c-----------------------------------------------------------------------
c
      parameter (ia=73, ja=73, jap2=ja+2)
      character*32 stamps, rstamp
      character*80 iotext
      dimension stamps(12), dpm(12)
c
      common /airs/ months, lons, lats, cxa(ia), cya(ja)
     &,             clmair(ia,ja,12)
      dimension airt(ia,ja)
      dimension a2(ia,jap2),  cyap2(jap2)
c
      data stamps /
     &               'm/d/y= 2/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 3/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 4/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 5/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 6/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 7/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 8/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 9/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y=10/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y=11/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y=12/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 1/ 1/0000, h:m:s= 0: 0: 0'
     &            /
c
      data dpm /31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0
     &, 30.0, 31.0, 30.0, 31.0/
c
      call getunit (nrdata, 'oorts.fmt'
     &,               'formatted sequential rewind')
      read (nrdata,9100) months, lons, lats, cxa, cya, clmair 
      print *,'  Oort airt:'
      print *,' months=',months,' lons=',lons, ' lats=',lats
c
      call getunit (nwdata, 'trashair'
     &,               'unformatted sequential rewind ieee')
      call getunit (nw, 'oorts.air'
     &,               'unformatted sequential rewind ieee')
      do j=1,ja
        do i=1,ia
	  airt(i,j) = 0.0
	enddo
      enddo
      sum = 0.0
      do n=1,12
	sum = sum + dpm(n)
        write (nwdata) stamps(n), dpm(n), lons, lats
     &, cxa, cya, ((clmair(i,j,n),i=1,ia),j=1,ja)
        print *,' n=',n,' air(30,30) =',clmair(30,30,n)
        do j=1,ja
          do i=1,ia
	    airt(i,j) = airt(i,j) + clmair(i,j,n)*dpm(n)
	  enddo
        enddo
        print *,' wrote oorts.air for = ',stamps(n), ' days ='
     &,  dpm(n)
      enddo
c
      rsum = 1.0/sum
      do j=1,ja
        do i=1,ia
	  airt(i,j) = airt(i,j)*rsum
	enddo
      enddo
c
      write (nwdata) stamps(12), sum, lons, lats, cxa, cya, airt
      print *,' annual mean air(30,30) =',airt(30,30)
      print *, ' 13th record is annual mean'
      print *,' wrote oorts.air for = ',stamps(12), ' days ='
     &,  sum
c
c     espand meridional dimension to cover poles
c
      rewind nwdata
      do n=1,13
        read (nwdata) rstamp, avg, lons, lats, cxa, cya, airt
	do j=1,ja
	  cyap2(j+1) = cya(j)
	  do i=1,ia
	    a2(i,j+1) = airt(i,j) 
	  enddo
	  a2(1,j+1)  = a2(ia-1,j+1)
	  a2(ia,j+1) = a2(2,j+1)
	enddo
	cyap2(1)    = cya(1) - (cya(2) - cya(1))
	cyap2(jap2) = cya(ja) + (cya(ja) - cya(ja-1))
	lats = jap2
	do i=1,ia
	  a2(i,1)    = a2(i,2)
	  a2(i,jap2) = a2(i,jap2-1)
	enddo 
        print *,' expanded n=',n,' airt(30,31) =',a2(30,31)
        iotext =
     & 'dim x(i),y(j),airt(i,j);read(nu)stamp,avg,i,j,mon,x,y,airt'
        print *, iotext
        write (nw) iotext
        write (nw) rstamp, avg, lons, lats, n, cxa, cyap2, a2
      enddo
c
      print *,' cya='
      write (*,'(1x,10g12.5)') cya            
      print *,' cyap2='
      write (*,'(1x,10g12.5)') cyap2            
c
9100  format (3i3/(10e14.7))
      return
      end



      subroutine levit
      parameter (imd=362, jmd=180, kd=33)
      character*32 stamps
      character*80 iotext
      dimension stamps(12), dpm(12)
      common /bigone/ month, lons, lats, depths, cx(imd), cy(jmd)
     &,               dd(kd), data(imd,jmd)
c
      character*12 tempin, tempout, saltin, saltout
      dimension tempin(13), tempout(13), saltin(13), saltout(13)
      dimension iotin(13), iotout(13), iosin(13), iosout(13)
c
      data tempin /
     &               'jan.temp.fmt', 'feb.temp.fmt', 'mar.temp.fmt'
     &,              'apr.temp.fmt', 'may.temp.fmt', 'jun.temp.fmt'
     &,              'jul.temp.fmt', 'aug.temp.fmt', 'sep.temp.fmt'
     &,              'oct.temp.fmt', 'nov.temp.fmt', 'dec.temp.fmt'
     &,              'ann.temp.fmt' /
c
      data saltin /
     &               'jan.salt.fmt', 'feb.salt.fmt', 'mar.salt.fmt'
     &,              'apr.salt.fmt', 'may.salt.fmt', 'jun.salt.fmt'
     &,              'jul.salt.fmt', 'aug.salt.fmt', 'sep.salt.fmt'
     &,              'oct.salt.fmt', 'nov.salt.fmt', 'dec.salt.fmt'
     &,              'ann.salt.fmt' /
c
c
      data tempout /
     &               'jan.temp', 'feb.temp', 'mar.temp'
     &,              'apr.temp', 'may.temp', 'jun.temp'
     &,              'jul.temp', 'aug.temp', 'sep.temp'
     &,              'oct.temp', 'nov.temp', 'dec.temp'
     &,              'ann.temp' /
c
      data saltout /
     &               'jan.salt', 'feb.salt', 'mar.salt'
     &,              'apr.salt', 'may.salt', 'jun.salt'
     &,              'jul.salt', 'aug.salt', 'sep.salt'
     &,              'oct.salt', 'nov.salt', 'dec.salt'
     &,              'ann.salt' /
c
      data stamps /
     &               'm/d/y= 2/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 3/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 4/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 5/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 6/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 7/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 8/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 9/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y=10/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y=11/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y=12/ 1/0000, h:m:s= 0: 0: 0'
     &,              'm/d/y= 1/ 1/0000, h:m:s= 0: 0: 0'
     &            /
c
      data dpm /31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0
     &, 30.0, 31.0, 30.0, 31.0/
c
c
c-----------------------------------------------------------------------
c     read formatted levitus annual & monthly "t" and "s"  and
c     write unformatted files
c-----------------------------------------------------------------------
c
      do n=1,13
        call getunit (iotin(n), tempin(n)
     &,               'formatted sequential rewind')
        call getunit (iotout(n), tempout(n)
     &,               'unformatted sequential rewind ieee')
        call getunit (iosin(n), saltin(n)
     &,               'formatted sequential rewind')
        call getunit (iosout(n), saltout(n)
     &,               'unformatted sequential rewind ieee')
      enddo
c
c     monthly temperature (19 levels deep)
c
      sum = 0.0
      do n=1,12
        print *, ' '
	sum = sum + dpm(n)
        do k=1,19
          read (iotin(n),9200) month, lons, lats, depths, cx, cy, dd
     &,                        data
          call setcy (data, imd, jmd)
          iotext =
     &'dim x(i),y(j),z(k),t(i,j);read(nu)stamp,avg,i,j,k,kl,mon,x,y,z,t'
          print *, iotext
          write (iotout(n)) iotext
          write (iotout(n)) stamps(n), dpm(n), lons, lats, kd, k, n, cx
     &,                     cy, dd, data 
          print *, ' => converted seasonal T for level ', k, ' n=',n
     &,   ' month=',month,', lons=',lons,', lats=',lats
          if (k .eq. 3) print *,' monthly t(180,90,3)=',data(180,90)
        enddo
      enddo
      print *,' depths=',depths, ' lats=',lats, ' lons=',lons
c
c     annual mean temperature
c
      do k=1,kd
        read (iotin(13),9200) month, lons, lats, depths, cx, cy, dd
     &,                        data
        call setcy (data, imd, jmd)
        iotext =
     &'dim x(i),y(j),z(k),t(i,j);read(nu)stamp,avg,i,j,k,kl,mon,x,y,z,t'
        print *, iotext
        write (iotout(13)) iotext
	n = 13
        write (iotout(13)) stamps(12), sum, lons, lats, kd, k, n, cx
     &,                     cy, dd, data 
          if (k .eq. 3) print *,' annual t(180,90,3)=',data(180,90)
      enddo
      print *,' wrote ann.salt for = ',stamps(12), ' days =', sum
c
c     monthly salinity (24 levels deep)
c
      do n=1,12
        print *, ' '
        do k=1,24
          read (iosin(n),9200) month, lons, lats, depths, cx, cy, dd
     &,                        data
          call setcy (data, imd, jmd)
          iotext =
     &'dim x(i),y(j),z(k),s(i,j);read(nu)stamp,avg,i,j,k,kl,mon,x,y,z,s'
          print *, iotext
          write (iosout(n)) iotext
          write (iosout(n)) stamps(n), dpm(n), lons, lats, kd, k, n, cx
     &,                     cy, dd, data 
          print *, ' => converted seasonal S for level ', k, ' n=',n
     &,   ' month=',month,', lons=',lons,', lats=',lats
          if (k .eq. 3) print *,' monthly s(180,90,3)=',data(180,90)
        enddo
      enddo
c
c     annual mean salinity
c
      do k=1,kd
        read (iosin(13),9200) month, lons, lats, depths, cx, cy, dd
     &,                        data
        call setcy (data, imd, jmd)
        iotext =
     &'dim x(i),y(j),z(k),s(i,j);read(nu)stamp,avg,i,j,k,kl,mon,x,y,z,s'
        print *, iotext
        write (iosout(13)) iotext
	n = 13
        write (iosout(13)) stamps(12), sum, lons, lats, kd, k, n, cx
     &,                     cy, dd, data 
          if (k .eq. 3) print *,' annual s(180,90,3)=',data(180,90)
      enddo
      print *,' wrote ann.salt for = ',stamps(12), ' days =', sum
9200  format (3i3/(10e14.7))
      return
      end
      
      subroutine setcy (a, im, jm)
      dimension a(im,jm)
      do j=1,jm
        a(1,j)  = a(im-1,j)
	a(im,j) = a(2,j)
      enddo
      return
      end
#else
      subroutine import
      return
      end
#endif
