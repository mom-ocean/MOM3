      subroutine fast_fir (f, im, km, nf, rmask, num)
c
c=======================================================================
c     simple finite impulse response filter with [.25, .5, .25] weights
c     using built in symmetric boundary conditions at land
c
c     input:
c             f     = functions to be filtered
c             im    = 1st dimension of "f" (longitudes)
c             km    = 2nd dimension of "f" (depth)
c             nf    = 3rd dimension of "f" (number of functions)
c             rmask = mask. must be (1.0,0.0) on (ocean,land) points
c             num   = number of filter applications
c
c     output:
c             f    = filtered quantities
c
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
      dimension f(im,km,nf), rmask(im,km), s(im,km)
      dimension cw(im,km), cc(im,km), ce(im,km)
c
      if (num .ne. 0) then
c
c     build weighting functions
c
      do k=1,km
        do i=2,im-1
	  cw(i,k) = 0.25*rmask(i-1,k)
	  ce(i,k) = 0.25*rmask(i+1,k)
	  cc(i,k) = 0.5*(1.0+rmask(i,k)) - cw(i,k) - ce(i,k)
	enddo
      enddo
c
c     filter "f"
c
      do npass=1,num
        do n=1,nf
	  do k=1,km
	    do i=2,im-1
	      s(i,k) = cw(i,k)*f(i-1,k,n) + cc(i,k)*f(i,k,n)
     &                +ce(i,k)*f(i+1,k,n)
	    enddo
	    s(1,k) = s(im-1,k)
	    s(im,k) = s(2,k)
	  enddo
	  do k=1,km
	    do i=1,im
	      f(i,k,n) = s(i,k)*rmask(i,k)
	    enddo
	  enddo
	enddo
      enddo
      endif
      return
      end
