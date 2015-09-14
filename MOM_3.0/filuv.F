      subroutine filuv (joff, js, je)
      use arrays_2d
      use memory_window
c
c=====================================================================
c     filuv sets up input needed for fourier filtering
c     (when the "fourfil" option is defined) -or- symmetric finite
c     impulse response filtering (when ifdef "firfil" is defined) of
c     baroclinic velocities at the specified high latitude row "jrow".
c=====================================================================
c
# include "param.h"
# include "coord.h"
# include "cpolar.h"
# include "emode.h"
# include "grdvar.h"
# include "index.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
      dimension tempik(imt,km,2)
c
#ifdef timing
      call tic ('baroclinic', 'filtering')
#endif
      do n=1,2
        do j=js,je
          call setbcx (u(1,1,j,n,taup1), imt, km)
        enddo
      enddo
#if defined firfil
      rdtuv = 1.0/c2dtuv
      do j=js,je
        jrow = j+joff
        if (numflu(jrow) .ne. 0) then
	  fx = sign(1.0,phi(jrow))
	  do k=1,km
	    do i=1,imt
	      uonly = u(i,k,j,1,taup1)
	      vonly = u(i,k,j,2,taup1)
              tempik(i,k,1) = -fx*uonly*spsin(i) - vonly*spcos(i)
              tempik(i,k,2) =  fx*uonly*spcos(i) - vonly*spsin(i)
	    enddo
	  enddo
c
	  call fast_fir (tempik, imt, km, 2, umask(1,1,j), numflu(jrow))
c
	  do k=1,km
	    do i=1,imt
              uonly =fx*(-tempik(i,k,1)*spsin(i)+tempik(i,k,2)*spcos(i))
              vonly = -tempik(i,k,1)*spcos(i) - tempik(i,k,2)*spsin(i)
              u(i,k,j,1,taup1) = uonly
              u(i,k,j,2,taup1) = vonly
	    enddo
	  enddo
        endif
      enddo
#endif
#if defined fourfil
c
c---------------------------------------------------------------------
c     filter u and v at high latitudes
c---------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        if ((jrow.le.jfu1 .or. jrow.ge.jfu2) .and. jrow.ge.jfrst) then
          jj = jrow - jfrst + 1
c
          if (phi(jrow) .gt. c0) then
            jj = jj - jskpu + 1
            refcos = refcosn
            fx = c1
          else
            refcos = refcoss
            fx = -c1
          endif
c
          isave = 0
          ieave = 0
c
          do l=1,lsegf
            do k=1,km
              if (isuf(jj,l,k) .ne. 0) then
                is = isuf(jj,l,k)
                ie = ieuf(jj,l,k)
                iredo = 1
                if (is.ne.isave .or. ie.ne.ieave) then
                  iredo = 0
                  im = ie - is + 1
                  isave = is
                  ieave = ie
# ifdef cyclic
                  if (im .ne. imtm2) then
                    m = 2
                    n = nint(im*csu(jrow)/refcos)
                  else
                    m = 3
                    n = nint(im*csu(jrow)*p5/refcos)
                  endif
# else
                  m = 2
                  n = nint(im*csu(jrow)/refcos)
# endif
                endif
c
c               rotate vector to compensate for polar latitude
c
                ism1 = is - 1
                iea = ie
                if (ie .ge. imt) iea = imtm1
                do i=is,iea
                  tempik(i-ism1,k,1) = -fx*u(i,k,j,1,taup1)*spsin(i)
     &                             - u(i,k,j,2,taup1)*spcos(i)
                  tempik(i-ism1,k,2) =  fx*u(i,k,j,1,taup1)*spcos(i)
     &                             - u(i,k,j,2,taup1)*spsin(i)
                enddo
                if (ie .ge. imt) then
                  ieb = ie - imtm2
                  ii  = imtm1 - is
                  do i=2,ieb
                    tempik(i+ii,k,1) = -fx*u(i,k,j,1,taup1)*spsin(i)
     &                             - u(i,k,j,2,taup1)*spcos(i)
                    tempik(i+ii,k,2) =  fx*u(i,k,j,1,taup1)*spcos(i)
     &                            - u(i,k,j,2,taup1)*spsin(i)
                  enddo
                endif
# ifdef fourfil
                call filtr (tempik(1,k,1), im, m, n, iredo)
                call filtr (tempik(1,k,2), im, m, n, 1)
# endif
c
c               rotate vector back to compensate for polar latitude
c
                do i=is,iea
                  u(i,k,j,1,taup1) = fx*(-tempik(i-ism1,k,1)*spsin(i)
     &                   + tempik(i-ism1,k,2)*spcos(i))
                  u(i,k,j,2,taup1) = -tempik(i-ism1,k,1)*spcos(i)
     &                   - tempik(i-ism1,k,2)*spsin(i)
                enddo
                if (ie .ge. imt) then
                  do i=2,ieb
                    u(i,k,j,1,taup1) = fx*(-tempik(i+ii,k,1)*spsin(i)
     &                     + tempik(i+ii,k,2)*spcos(i))
                    u(i,k,j,2,taup1) = -tempik(i+ii,k,1)*spcos(i)
     &                     - tempik(i+ii,k,2)*spsin(i)
                  enddo
                endif
              endif
            enddo
          enddo
c
          if (isave .ne. 0 .and. ieave .ne. 0) then
            do i=1,imt
              tempik(i,1,1) = c0
              tempik(i,1,2) = c0
            enddo
c
            do k=1,km
              do i=1,imt
                tempik(i,1,1) = tempik(i,1,1) + u(i,k,j,1,taup1)*dzt(k)
                tempik(i,1,2) = tempik(i,1,2) + u(i,k,j,2,taup1)*dzt(k)
              enddo
            enddo
c
            do i=1,imt
              tempik(i,1,1) = tempik(i,1,1)*hr(i,jrow)
              tempik(i,1,2) = tempik(i,1,2)*hr(i,jrow)
            enddo
c
            do k=1,km
              do i=1,imt
                u(i,k,j,1,taup1) = u(i,k,j,1,taup1) - tempik(i,1,1)
                u(i,k,j,2,taup1) = u(i,k,j,2,taup1) - tempik(i,1,2)
              enddo
            enddo
c
            do k=1,km
              do i=1,imt
                u(i,k,j,1,taup1) = u(i,k,j,1,taup1)*umask(i,k,j)
                u(i,k,j,2,taup1) = u(i,k,j,2,taup1)*umask(i,k,j)
              enddo
            enddo
          endif
c
        endif
      enddo
#endif
#ifdef timing
      call toc ('baroclinic', 'filtering')
#endif
      return
      end
