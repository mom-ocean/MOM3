      subroutine unesco_density (tt, ss, p, rho, js, je)
!=======================================================================
!     compute unesco insitu density (rho in units of gm/cm^3) from
!     salinity (ss in psu), potential temperature (tt in deg C)
!     and pressure (p in bars)
!
!     Reference: Jackett and McDougall, Minimal Adjustment of 
!     Hydrographic Profiles to Achieve Static Stablilty, Journal of
!     Atmospehric and Oceanic Technology, Vol 12, 381-389,  April 1995
!     
!     R.C.Pacanowski rcp@gfdl.gov
!=======================================================================
#ifdef driver
      parameter (imt=122, km=15, jmw=3)
#else
# include "size.h"
#endif
      dimension tt(imt,km,jmw), ss(imt,km,jmw), p(imt,km,jmw)
      dimension rho(imt,km,jmw)
      

#if defined linearized_density
      do j=js,je
        rho(:,:,j) = 1.035 - 2.e-4*tt(:,:,j)
      enddo
#else
      dimension t(imt,km), t2(imt,km), t3(imt,km), t4(imt,km),t5(imt,km)
      dimension s(imt,km), s2(imt,km), s3(imt,km), s32(imt,km)
      dimension p1(imt,km), p2(imt,km)
      dimension rho0(imt,km), sk(imt,km)      
      
      c2cgs = 1.e-3 !convert kg/m^3 to gm/cm^3
      do j=js,je
        t(:,:)  = tt(:,:,j)
        s(:,:)  = ss(:,:,j)
        p1(:,:) = p(:,:,j)
        t2(:,:) = t(:,:)**2
        t3(:,:) = t(:,:)**3
        t4(:,:) = t(:,:)**4
        t5(:,:) = t(:,:)**5
        s2(:,:) = s(:,:)**2
        s3(:,:) = s(:,:)**3
        p2(:,:) = p1(:,:)**2
        s32(:,:) = sqrt(s3(:,:))

!       compute rho(s,theta,p=0) ! (same as rho(s,t_insitu,p=0) )

        rho0(:,:) = 999.842594 + 6.793952e-2*t(:,:)
     &              - 9.095290e-3*t2(:,:) + 1.001685e-4*t3(:,:)
     &              - 1.120083e-6*t4(:,:) + 6.536332e-9*t5(:,:)
     &              + s(:,:)*(0.824493 - 4.0899e-3*t(:,:)
     &              + 7.6438e-5*t2(:,:)
     &              - 8.2467e-7*t3(:,:) + 5.3875e-9*t4(:,:))
     &              + s32(:,:)*(-5.72466e-3 + 1.0227e-4*t(:,:)
     &              - 1.6546e-6*t2(:,:)) + 4.8314e-4*s2(:,:)

!       compute rho(s,theta,p)

        sk(:,:) = 1.965933e4 + 1.444304e2*t(:,:) - 1.706103*t2(:,:) 
     &            + 9.648704e-3*t3(:,:)  - 4.190253e-5*t4(:,:)
     &            + s(:,:)*(52.84855 - 3.101089e-1*t(:,:)
     &            + 6.283263e-3*t2(:,:) -5.084188e-5*t3(:,:))
     &            + s32(:,:)*(3.886640e-1 + 9.085835e-3*t(:,:)
     &            - 4.619924e-4*t2(:,:))
     &            + p1(:,:)*(3.186519 + 2.212276e-2*t(:,:)
     &            - 2.984642e-4*t2(:,:) + 1.956415e-6*t3(:,:)) 
     &            + p1(:,:)*s(:,:)*(6.704388e-3  -1.847318e-4*t(:,:)
     &            + 2.059331e-7*t2(:,:)) + 1.480266e-4*p1(:,:)*s32(:,:)
     &            + p2(:,:)*(2.102898e-4 - 1.202016e-5*t(:,:)
     &            + 1.394680e-7*t2(:,:)) +p2(:,:)*s(:,:)*(-2.040237e-6
     &            + 6.128773e-8*t(:,:) + 6.207323e-10*t2(:,:))

        rho(:,:,j) = c2cgs*rho0(:,:) / (1.0 - p1(:,:)/sk(:,:))
      enddo
#endif
      return
      end      

#ifdef driver
      program main
      parameter (imt=122, km=15, jmw=3)
      dimension t(imt,km,jmw), s(imt,km,jmw), p(imt,km,jmw)
      dimension rho(imt,km,jmw)
      
      t(:,:,:) = 3.0             ! potential temperature in deg C
      s(:,:,:) = 35.5            ! psu
      p(:,:,:) = 300.0           ! bars (1 bar = 10^6 dynes/cm^2)

!     test density: rho(t=3.0,s=35.5,p=300.0) = 1041.83267
!     for t=3.0,s=35.5,p=300.0, insitu temp = 3.25905 deg C

      do n=1,10000
        call unesco_density (t, s, p, rho, 1, jmw)
      enddo
      print *,' unesco_density rho=',rho(imt,km,1)
     &,rho(imt,km,2), rho(imt,km,3)
      stop
      end
#endif
