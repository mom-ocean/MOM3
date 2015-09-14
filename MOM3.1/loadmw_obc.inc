c====================== add-in "loadmw_obc.inc" ========================
c
c       extrapolate tracer and velocity values onto boundary 
c       to prevent diffusion
c
#if defined obc_west || defined obc_east
c
      do j=js,je
        jrow = j + joff
        if (jrow .ge. 2 .and. jrow .le. jmtm2) then
	  do n=1,2
            do k=1,km
# ifdef obc_west
	      u(1,k,j,n,tau)   = u(2,k,j,n,tau)
	      u(1,k,j,n,taum1) = u(2,k,j,n,taum1)
# endif
# ifdef obc_east
	  u(imtm1,k,j,n,tau)   = u(imtm2,k,j,n,tau)
	  u(imtm1,k,j,n,taum1) = u(imtm2,k,j,n,taum1)
	  u(imt,  k,j,n,tau)   = u(imtm2,k,j,n,tau)
	  u(imt,  k,j,n,taum1) = u(imtm2,k,j,n,taum1)
# endif
	    enddo
	  enddo
        endif
        if (jrow .ge. 2 .and. jrow .le. jmtm1) then
	  do n=1,nt
            do k=1,km
# ifdef obc_west
	      t(1,k,j,n,tau)     = t(2,k,j,n,tau)
	      t(1,k,j,n,taum1)   = t(2,k,j,n,taum1)
# endif
# ifdef obc_east
	      t(imt,k,j,n,tau)   = t(imtm1,k,j,n,tau)
	      t(imt,k,j,n,taum1) = t(imtm1,k,j,n,taum1)
# endif
	    enddo
	  enddo
        endif
      enddo
#endif
#if defined obc_south || defined obc_north
c
      do j=js,je
        jrow = j + joff
# ifdef obc_south	
        if (jrow .eq. 2) then
	  do n=1,2
            do k=1,km
	      do i=1,imt
	        u(i,k,j-1,n,tau)   = u(i,k,j,n,tau)
	        u(i,k,j-1,n,taum1) = u(i,k,j,n,taum1)
	      enddo
	    enddo
	  enddo
	  do n=1,nt
            do k=1,km
	      do i=1,imt
	        t(i,k,j-1,n,tau)   = t(i,k,j,n,tau)
	        t(i,k,j-1,n,taum1) = t(i,k,j,n,taum1)
	      enddo
	    enddo
	  enddo
        endif
# endif
# ifdef obc_north
        if (jrow .eq. jmt) then
	  do n=1,2
            do k=1,km
	      do i=1,imt
	        u(i,k,j,  n,tau)   = u(i,k,j-2,n,tau)
	        u(i,k,j,  n,taum1) = u(i,k,j-2,n,taum1)
	        u(i,k,j-1,n,tau)   = u(i,k,j-2,n,tau)
	        u(i,k,j-1,n,taum1) = u(i,k,j-2,n,taum1)
	      enddo
	    enddo
	  enddo
	  do n=1,nt
            do k=1,km
	      do i=1,imt
	        t(i,k,j,n,tau)   = t(i,k,j-1,n,tau)
	        t(i,k,j,n,taum1) = t(i,k,j-1,n,taum1)
                rho(i,k,j)       = rho(i,k,j-1)
	      enddo
	    enddo
	  enddo
        endif
# endif
      enddo
#endif
