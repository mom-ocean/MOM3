c====================== add-in "barotropic1_obc.inc" ===================
c
c-----------------------------------------------------------------------
c     calculate phase velocities for psi
c     Orlanski radiation condition (passive open boundary)   
c-----------------------------------------------------------------------
c
#ifdef obc_south
      var = -dyu(3)/dtts
      do i=1,imt
         var1(i,1,1) = psi(i,4,2) - psi(i,3,2)
         if (var1(i,1,1) .eq. c0) then
           c1ps(i) = var
         else
           c1ps(i) = var*(psi(i,3,1)-psi(i,3,2))/var1(i,1,1)
           if (c1ps(i) .ge. c0) c1ps(i) = c0
           if (c1ps(i) .lt. var) c1ps(i) = var
         endif
      enddo
#endif
#ifdef obc_north
# ifdef obc_north_anywhere
      var = dyu(jobcn-1)/dtts
      do i=iobcns,iobcne
         var1(i,1,1) = psi(i,jobcn-1,2) - psi(i,jobcn-2,2)
         if (var1(i,1,1) .eq. c0) then
           c1pn(i) = var
         else
           c1pn(i) =-var*(psi(i,jobcn-1,1)-psi(i,jobcn-1,2))/var1(i,1,1)
           if (c1pn(i) .le. c0) c1pn(i) = c0
           if (c1pn(i) .gt. var) c1pn(i) = var
         endif
      enddo
# else
      var = dyu(jmt-2)/dtts
      do i=1,imt
         var1(i,1,1) = psi(i,jmt-2,2) - psi(i,jmt-3,2)
         if (var1(i,1,1) .eq. c0) then
           c1pn(i) = var
         else
           c1pn(i) = -var*(psi(i,jmt-2,1)-psi(i,jmt-2,2))/var1(i,1,1)
           if (c1pn(i) .le. c0) c1pn(i) = c0
           if (c1pn(i) .gt. var) c1pn(i) = var
         endif
      enddo
# endif
#endif
c
#ifdef obc_west
c
      do j=1,jmt
         var = -dxu(3)*csu(j)/dtts
         var2(1,j,1) = psi(4,j,2)-psi(3,j,2)
         if (var2(1,j,1) .eq. c0) then
           c1pw(j) = var
         else
           c1pw(j) = var*(psi(3,j,1)-psi(3,j,2))/var2(1,j,1)
           if (c1pw(j) .ge. c0) c1pw(j) = c0
           if (c1pw(j) .lt. var) c1pw(j) = var
         endif
      enddo
c
#endif
#ifdef obc_east
c
      do j=1,jmt
         var = dxu(imt-2)*csu(j)/dtts
         var2(1,j,1) = psi(imt-2,j,2)-psi(imt-3,j,2)
         if (var2(1,j,1) .eq. c0) then
           c1pe(j) = var
         else
           c1pe(j) = -var*(psi(imt-2,j,1)-psi(imt-2,j,2))/var2(1,j,1)
           if (c1pe(j) .le. c0) c1pe(j) = c0
           if (c1pe(j) .gt. var) c1pe(j) = var
         endif
      enddo
c
#endif
