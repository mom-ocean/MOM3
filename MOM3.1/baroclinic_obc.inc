c====================== add-in "baroclinic_obc.inc" ====================
c
c       disable advective velocities at the open boundaries
c
#if defined obc_south || defined obc_north
# if defined obc_south && !defined obc_north
          if (jrow .eq. 2) then
# endif
# if !defined obc_south && defined obc_north
#  ifdef obc_north_anywhere
          if (jrow .eq. jobcn-1) then
#  else
          if (jrow .eq. jmtm2) then
#  endif
# endif
# if defined obc_south && defined obc_north
#  ifdef obc_north_anywhere
          if (jrow .eq. 2 .or. jrow.eq.jobcn-1) then
#  else
          if (jrow .eq. 2 .or. jrow.eq.jmtm2) then
#  endif
# endif
c
            do k=1,kbot
              do i=istrt,iend
                u(i,k,j,n,taup1) = (u(i,k,j,n,taup1) 
#  if !defined linearized_advection
     &             + ADV_Ux(i,k,j) + ADV_Uy(i,k,j) + ADV_Uz(i,k,j)
     &             - ADV_metric(i,k,j,n) 
#  endif
     &             )*umask(i,k,j)
	       enddo
	     enddo
	   endif
#endif
c
#if defined obc_west || defined obc_east
            do k=1,kbot
c
# ifdef obc_west	  
	      i = 2
	      u(i,k,j,n,taup1) = (u(i,k,j,n,taup1) 
#   if !defined linearized_advection
     &             + ADV_Ux(i,k,j) + ADV_Uy(i,k,j) + ADV_Uz(i,k,j)
     &             - ADV_metric(i,k,j,n) 
#   endif
     &             )*umask(i,k,j)
	      u(i-1,k,j,n,taup1) = u(i,k,j,n,taup1)

# endif
# ifdef obc_east	  
	      i = imtm2
	      u(i,k,j,n,taup1) = (u(i,k,j,n,taup1) 
#   if !defined linearized_advection
     &             + ADV_Ux(i,k,j) + ADV_Uy(i,k,j) + ADV_Uz(i,k,j)
     &             - ADV_metric(i,k,j,n) 
#   endif
     &             )*umask(i,k,j)
              u(i+1,k,j,n,taup1) = u(i,k,j,n,taup1)
              u(i+2,k,j,n,taup1) = u(i,k,j,n,taup1)
# endif
	     enddo
#endif
