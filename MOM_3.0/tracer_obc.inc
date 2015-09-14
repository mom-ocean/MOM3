c====================== add-in "tracer_obc.inc" ========================
#if defined obc_north
c
c-----------------------------------------------------------------------
c       calculate auxiliary field for cobc (only for biharmonic option)
c-----------------------------------------------------------------------
c
      if (jrow .eq. jmt-2) then
        do m=1,nt
          do k=1,km
	    do i=2,imtm1
              var1(i,k,m)=t(i,k,j,m,taum1)-t(i,k,j-1,m,taum1)
	    enddo
	  enddo
	enddo
      endif
#endif
c
c       choose start and end indices for "full" calculation
#ifdef obc_west
            iscalc = 3
#else
            iscalc = istrt
#endif
#ifdef obc_east
            iecalc = imtm2
#else
            iecalc = iend
#endif
c
#ifdef obc_west
        call cobc2 (joff, js, je, 2)
#endif
c	
        do j=js,je
          jrow   = j + joff
          do k=1,km
	    twodt(k) = c2dtts*dtxcel(k)
            do i=iscalc,iecalc
	      t(i,k,j,n,taup1) = t(i,k,j,n,taum1) + twodt(k)*(
     &                  DIFF_Tx(i,k,j) + DIFF_Ty(i,k,j) + DIFF_Tz(i,k,j)
     &                 - ADV_Tx(i,k,j) - ADV_Ty(i,k,j) - ADV_Tz(i,k,j)
#ifdef source_term
     &                 + source(i,k,j)
#endif
     &                           )*tmask(i,k,j)
            enddo
          enddo
        enddo
c
#ifdef trace_indices
      if (n.eq.1) write (stdout,'(2x,5(a,i4))')
     & "=> In tracer: js=",js," je=",je," joff=",joff
     &," jrows=",js+joff," to ",je+joff
#endif
c
#ifdef obc_east
        call cobc2 (joff, js, je, imtm1)
#endif
