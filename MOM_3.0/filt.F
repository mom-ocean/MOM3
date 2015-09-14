      subroutine filt (joff, js, je)
c
c=======================================================================
c     subroutine filt sets up input needed for fourier filtering
c     (when the "fourfil" option is defined) -or- symmetric finite
c     impulse response filtering (when the "firfil" option is defined)
c     of tracers at the specified high latitude row "jrow".
c=======================================================================
c
      use arrays_2d
      use memory_window
# include "param.h"
# include "grdvar.h"
# include "index.h"
# include "mw.h"
# include "scalar.h"
#if defined firfil
      dimension fik(imt,km,nt)
#endif
#if defined fourfil
      dimension tempik(imt,km)
#endif
#ifdef timing
      call tic ('tracer', 'filtering')
#endif
c
      do n=1,nt
        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo
c
#if defined firfil
      do j=js,je
        jrow = j+joff
        if (numflt(jrow) .ne. 0) then
          r2dt = c1/c2dtts
          do n=1,nt
	    do i=1,imt
	      do k=1,km
	        fik(i,k,n) = r2dt*(t(i,k,j,n,taup1)-t(i,k,j,n,taum1))
	      enddo
	    enddo
	  enddo
c
	  call fast_fir (fik, imt, km, nt, tmask(1,1,j), numflt(jrow))
c
          do n=1,nt
	    do i=1,imt
	      do k=1,km
	        t(i,k,j,n,taup1) = t(i,k,j,n,taum1) + c2dtts*fik(i,k,n)
	      enddo
	    enddo
	  enddo
        endif
      enddo
#endif
#if defined fourfil
      r2dt = c1/c2dtts
      do j=js,je
        jrow = j + joff
        if ((jrow.le.jft1.or.jrow.ge.jft2) .and. jrow.ge.jfrst) then
          jj = jrow-jfrst+1
c
          if (phi(jrow) .gt. c0) then
            jj = jj - jskpt + 1
            refcos = refcosn
          else
            refcos = refcoss
          endif
c
c        if previous strips were of same length, do not recompute
c        fourier coeffs
c
          isave = 0
          ieave = 0
          do l=1,lsegf
            do k=1,km
              if (istf(jj,l,k) .ne. 0) then
                is    = istf(jj,l,k)
                ie    = ietf(jj,l,k)
                iredo = 0
                if (is.ne.isave .or. ie.ne.ieave) then
                  iredo = -1
                  isave = is
                  ieave = ie
                  im = ie-is+1
# ifdef cyclic
                  if (im.ne.imtm2 .or. kmt(1,jrow).lt.k) then
                    m = 1
                    n = nint(im*cst(jrow)/refcos)
                  else
                    m = 3
                    n = nint(im*cst(jrow)*0.5/refcos)
                  endif
# else
                  m = 1
                  n = nint(im*cst(jrow)/refcos)
# endif
                endif
                do mm=1,nt
c
c                 filter tracer time tendency
c
                  idx  = iredo+mm
                  ism1 = is-1
                  iea  = ie
                  if (ie .ge. imt) iea = imtm1
                  do i=is,iea
                    tempik(i-ism1,k) =
     &                        (t(i,k,j,mm,taup1)-t(i,k,j,mm,taum1))*r2dt 
                  enddo
                  if (ie .ge. imt) then
                    ieb = ie-imtm2
                    ii  = imtm1-is
                    do i=2,ieb
                      tempik(i+ii,k) =
     &                        (t(i,k,j,mm,taup1)-t(i,k,j,mm,taum1))*r2dt 
                    enddo
                  endif
c
# ifdef fourfil
                  call filtr (tempik(1,k), im, m, n, idx)
# endif
c
c                 re-construct tracer at tau+1 from filtered tendency
c
                  do i=is,iea
                    t(i,k,j,mm,taup1) = 
     &                       t(i,k,j,mm,taum1) + c2dtts*tempik(i-ism1,k)
                  enddo
                  if (ie .ge. imt) then
                    do i=2,ieb
                      t(i,k,j,mm,taup1) =  
     &                       t(i,k,j,mm,taum1) + c2dtts*tempik(i+ii,k)
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
        endif
      enddo
c
#endif
#ifdef timing
      call toc ('tracer', 'filtering')
#endif
      return
      end
