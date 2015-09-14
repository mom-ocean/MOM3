
      subroutine diagvmix (joff, js, je)
#ifdef matrix_sections
      character*128 file_name, file_stamp
# include "param.h"
# include "coord.h"
# include "cprnts.h"
# include "iounit.h"
# include "switch.h"
# include "tmngr.h"
# include "vmixc.h"
c
# ifdef timing
      call tic ('diagnostic', 'matrix sections')
# endif
      do j=js,je
	jrow = j + joff
        reltim = relyr
        do jlat=1,nlatpr
          jj = indp (prlat(jlat), yt, jmt)
          if (jj .eq. jrow .and. prlat(jlat) .le. yt(jmt)) then
            is = indp (prslon(jlat), xt, imt)
            ie = indp (prelon(jlat), xt, imt)
            ks = indp (prsdpt(jlat), zt, km)
            ke = indp (predpt(jlat), zt, km)
            fx = 1.0e-2
c
            scl = c1
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'diff_cbt', itt, jrow 
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (diff_cbt(1,1,j,1), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
              file_name = file_stamp ('sections', xstamp, '.dta')
	      call getunit (io, trim(file_name)
     &,                     'unformatted sequential append ieee')
              write (stdout,*) ' => diff_cbt ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((diff_cbt(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, diff_cbt(1,1,j,1), imt*km)
              call relunit (io)
            endif
c
            scl = c1
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'visc_cbu', itt, jrow 
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (visc_cbu(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
              file_name = file_stamp ('sections', xstamp, '.dta')
	      call getunit (io, trim(file_name)
     &,                     'unformatted sequential append ieee')
              write (stdout,*) ' => visc_cbu ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((visc_cbu(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, visc_cbu(1,1,j), imt*km)
              call relunit (io)
            endif
          endif
        enddo
      enddo
# ifdef timing
      call toc ('diagnostic', 'matrix sections')
# endif
9100    format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     &,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     &,', scaling=',1pg10.3)
#endif
      return
      end
