c====================== add-in "barotropic2_obc.inc" ===================
# ifdef orlanski
c
c-----------------------------------------------------------------------
c     passive Orlanski radiation condition (passive open boundary)   
c-----------------------------------------------------------------------
c
      do i=2,imtm1
#ifdef obc_south
         psi(i,2,1)    =psi(i,2,2)    - c1ps(i)*dtts*dyur(2)
     &                                 *(psi(i,3,2)    - psi(i,2,2))
         psi(i,1,1)    =psi(i,2,1)
#endif
#ifdef obc_north
# ifdef obc_north_anywhere
         psi(i,jobcn,1)=psi(i,jobcn,2)- c1pn(i)*dtts*dyur(jobcn)
     &                             *(psi(i,jobcn,2)- psi(i,jobcn-1,2))
         psi(i,jmt,1)  =psi(i,jobcn,1)
# else
         psi(i,jmt-1,1)=psi(i,jmt-1,2)- c1pn(i)*dtts*dyur(jmt-1)
     &                                 *(psi(i,jmt-1,2)- psi(i,jmt-2,2))
         psi(i,jmt,1)  =psi(i,jmt-1,1)
# endif
#endif
       enddo
c
      do j=2,jmtm1
#ifdef obc_west
         psi(2,j,1)    =psi(2,j,2)    - c1pw(j)*dtts*dxur(2)*csur(j)
     &                                 *(psi(3,j,2)    - psi(2,j,2))
         psi(1,j,1)    =psi(2,j,1)
#endif
#ifdef obc_east
         psi(imt-1,j,1)=psi(imt-1,j,2)- c1pe(j)*dtts*dxur(imt-1)*csur(j)
     &                                 *(psi(imt-1,j,2)- psi(imt-2,j,2))
         psi(imt,j,1)  =psi(imt-1,j,1)
#endif
       enddo
c
c
# else
      call addobcpsi
# endif
