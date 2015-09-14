      subroutine tcmixi
#ifdef tcvmix
c
c=======================================================================
c     initialize mellor-yamada level 2.5 turbulence closure scheme
c
c     author:      a. rosati      e-mail=> ar@gfdl.gov
c=======================================================================
c
c
c     data initialize constants for turbulence closure vertical mixing
c     see "ctcmix.h"
c
      a1 = 0.92
      b1 = 16.6
      a2 = 0.74
      b2 = 10.1
      cc1 = 0.08
      e1 = 1.8
      e2 = 1.33
      e3 = 1.0
      aq = 2.0e7
# ifdef implicitvmix
      visc_cbu_limit = 1.e6
      diff_cbt_limit = 1.e6
      vdqlim = 1.e6
# else
      visc_cbu_limit = 1.e1
      diff_cbt_limit = 1.e1
      vdqlim = 1.e1
# endif
      alpha = 0.1
      vk    = 0.40
      sq    = 0.2
      small = 1.e-6
      big   = 1.e15
      deps  = 1.e-10
      c100  = 100.0
      p48   = 0.48
      p032  = 0.032
      c15   = 15.0
      p1    = 0.1
      p4    = 0.4
      p154  = 0.154
      do i=1,imtkm
        gh(i) = 0.154
	gmp1(i) = 0.154
	sm(i) = 0.39
	sh(i) = 0.49
      enddo
      rhoh2o = 1.025
c
      return
      end
      subroutine tcmix (j)
c
c=======================================================================
c
c     compute vertical mixing coefficients based on mellor-yamada level
c     2.5 turbulence closure scheme as implemented and described in
c     rosati & miyakoda (jpo,vol 18,#11,1988).
c
c     inputs:
c
c     j      = row j
c
c     outputs:
c
c     vertical mixing coefficient arrays "vdca,vvca or vdc,vvc"
c     in "tcslab.h"
c
c=======================================================================
c
      use arrays_2d
#include "param.h"
#include "coord.h"
#include "cprnts.h"
#include "chmix.h"
#include "ctask.h"
#include "ctcmix.h"
#include "ctmngr.h"
#include "cvbc.h"
#include "grdvar.h"
#include "iounit.h"
#include "scalar.h"
#include "slabs.h"
#include "switch.h"
#include "temp.h"
#include "vmixc.h"
      dimension drho(imt,km)
#ifdef lalg
      dimension sum0(imt), sum1(imt)
#endif
      dimension ekn(imt,km,nslabs)
      dimension bcl(imt), ustar(imt)
#ifndef implicitvmix
      dimension q2dif(imt,0:km)
#  ifdef leq
     $         ,q2ldif(imt,0:km)
#  endif
#endif
c
c----------------------------------------------------------------------
c     construct physical operators using statement functions
c----------------------------------------------------------------------
c
c     zonal,meridional & vertical advection of tke
c
      uq2x(i,k) = aux1(i+1,k) - aux1(i,k)
      vq2y(i,k) = p5*(fvn(i,k)+fvn(i,k+1))
     $              *(q2(i,k,jp,nc)+q2(i,k,jc,nc))
     $           -p5*(fvsq(i,k)+fvsq(i,k+1))
     $              *(q2(i,k,jc,nc)+q2(i,k,jm,nc))
      wq2z(i,k) = (fwb1(i,k-1) - fwb1(i,k))*dzw2r(k)*p5
c
c     zonal,meridional diffusion of tke
c
          q2xx(i,k)=bbtj*dxt4r(i)*
     $             (fm(i+1,k,jc)*aux3(i+1,k)-fm(i-1,k,jc)*aux3(i,k))
          q2yy(i,k) = cctj*fm(i,k,jp)*(q2(i,k,jp,nm)-q2(i,k,jc,nm))+
     $                 ddtj*fm(i,k,jm)*(q2(i,k,jm,nm)-q2(i,k,jc,nm))
#ifndef implicitvmix
c
c     vertical diffusion of tke
c
      q2zz(i,k) = eeq(i,k)*q2dif(i,k-1) - ffq(i,k)*q2dif(i,k)
#endif
#ifdef leq
c
c     zonal,meridional & vertical advection of tkel
c
      uq2lx(i,k) = aux2(i+1,k) - aux2(i,k)
      vq2ly(i,k) = p5*(fvn(i,k)+fvn(i,k+1))
     $              *(q2l(i,k,jp,nc)+q2l(i,k,jc,nc))
     $           -p5*(fvsq(i,k)+fvsq(i,k+1))
     $              *(q2l(i,k,jc,nc)+q2l(i,k,jm,nc))
      wq2lz(i,k) = (fwb2(i,k-1) - fwb2(i,k))*dzw2r(k)*p5
c
c     zonal,meridional diffusion of tke
c
          q2lxx(i,k)=bbtj*dxt4r(i)*
     $             (fm(i+1,k,jc)*aux4(i+1,k)-fm(i-1,k,jc)*aux4(i,k))
          q2lyy(i,k) = cctj*fm(i,k,jp)*(q2l(i,k,jp,nm)-q2l(i,k,jc,nm))+
     $                 ddtj*fm(i,k,jm)*(q2l(i,k,jm,nm)-q2l(i,k,jc,nm))
#  ifndef implicitvmix
c
c     vertical diffusion of tke
c
      q2lzz(i,k) = eeq(i,k)*q2ldif(i,k-1) - ffq(i,k)*q2ldif(i,k)
#  endif
#endif
c
c     initialize mixing coefficients to constant values
c
      if (init .and. first) then
        do  100 k=1,km
          do  90 i=1,imt
            vvc (i,k,jc,nc) = kappa_m
            vdc (i,k,jc,nc) = kappa_h
            vdqc(i,k,jc,nc) = kappa_h
            vvc (i,k,jc,nm) = kappa_m
            vdc (i,k,jc,nm) = kappa_h
            vdqc(i,k,jc,nm) = kappa_h
90        continue
100     continue
      endif
      if (first) then
c
c     set up some print locations for diagnostic printout
c
        jpr1 = jmt/2
        ipr1 = imt/2
        jpr2 = (3*jmt)/4
        ipr2 = imt/2
      endif
      if (j .eq. jsrow) then
c
c-----------------------------------------------------------------------
c       on the first row of a task:
c
c       1) compute the advective coefficients at the south face of the
c          "t" grid box for the next row (jsrow+1)
c
c       2) return
c-----------------------------------------------------------------------
c
        fxb = cstr(j+1)*dytr(j+1)*csu(j)
        do 200 k=1,km
          do 190 i=2,imt
            fvsq(i,k)  = (v(i,k,jc,nc)*dxu(i) + v(i-1,k,jc,nc)*dxu(i-1))
     $               *fxb*dxt4r(i)
190       continue
200     continue
c
        return
      else
c
c-----------------------------------------------------------------------
c       on all but the first row of a task:
c
c       1) compute the advective coefficients "fuw" at the west face of
c         the "t" grid box & "fvn" at the north face of the "t" grid box
c
c-----------------------------------------------------------------------
c
        fxa = cstr(j)*dytr(j)
        fxb = fxa*csu(j)
        do 300 k=1,km
          do 290 i=2,imt
            fuw(i,k) = (u(i-1,k,jc,nc)*dyu(j) + u(i-1,k,jm,nc)*dyu(j-1))
     $               *fxa*dxt4r(i)
            fvn(i,k) = (v(i  ,k,jc,nc)*dxu(i) + v(i-1,k,jc,nc)*dxu(i-1))
     $               *fxb*dxt4r(i)
290       continue
300     continue
c
      endif
c
      tdt=c2dtts
      do  400 k=1,km
        do  390 i=1,imt
          sprod(i,k) = c0
          bprod(i,k) = c0
          dissp(i,k) = c0
          tke(i,k)   = c0
          hdq2(i,k)  = c0
#ifdef leq
          tkel(i,k)  = c0
          hdq2l(i,k) = c0
#endif
390     continue
400   continue
c
c-----------------------------------------------------------------------
c    set time level ln depending on implicit or explicit treatment
c-----------------------------------------------------------------------
#ifdef implicitvmix
      ln = nc
#else
      ln = nm
#endif
c
c-----------------------------------------------------------------------
c     compute vertical difference of density on "t" grid row j+1 at
c     time level tau for implicit time differencing.
c     time level tau-1 for explicit  time differencing.
c     set no-flux condition on density at bottom. bottom boundary
c     condition on velocity will determine if there is drag.
c     bottom boundary condition on tracer is no flux.
c-----------------------------------------------------------------------
c
      do 500 ks=1,2
        call statec (t(1,1,jc,ln,1), t(1,1,jc,ln,2), tempik(1,1,1), ks)
        do 490 k=ks,kmm1,2
          do 480 i=1,imt
            drho(i,k) = tempik(i,k,1) - tempik(i,k+1,1)
480       continue
490     continue
500   continue
c
      do 600 i=1,imt
        if (kmt(i,j) .ne. 0) then
          drho(i,kmt(i,j)) = c0
        endif
600   continue
c
#ifdef dataassim
c
c     increment rho to increase stability for data assimilation
c
      do 700 k=1,kmm1
        do 690 i=1,imt
          drho(i,k)=drho(i,k)+2.e-6*dzwr(0)/dzwr(k-1)
690     continue
700   continue
#endif
c
c      compute vertical diffusion coefficients for q2
c
      do 800 k=2,km
        do 790 i=1,imt
          eeq(i,k) = p5*(vdqc(i,k-1,jc,ln)+vdqc(i,k,jc,ln))*dzwur(k)
790     continue
800   continue
      do 900 k=1,kmm1
        do 890 i=1,imt
          ffq(i,k) = p5*(vdqc(i,k,jc,ln)+vdqc(i,k+1,jc,ln))*dzwlr(k)
890     continue
900   continue
      do 1000 i=1,imt
          eeq(i,1)   = p5*vdqc(i,1,jc,ln)*dzwur(1)
          ffq(i,km)=c0
1000  continue
c
#ifdef leq
c     set top boundary condition for q2l
      do  1100 i=1,imt
        bcl(i) = c0
1100  continue
#endif
c
c-----------------------------------------------------------------------
c     solve for turbulant kinetic energy
c-----------------------------------------------------------------------
c
c
c     compute friction velocity ustar
c
      do 1200 i = 1,imt
        ustar(i) = smf(i,1)**2 + smf(i,2)**2
1200  continue
      b1tt = (b1**(2./3.))/rhoh2o
c    note that ustar is really b1**2/3*ustar**2
      do 1300 i=1,imt
        ustar(i) = b1tt*sqrt(ustar(i))
1300  continue
#ifdef lalg
c
c     compute length scale from algebraic eqn (lalg)
c
      do 1400 i = 1,imt
        sum0(i) = c0
        sum1(i) = c0
1400  continue
      do 1500 k=1,km
        do 1490 i=1,imt
          sum0(i) = sqrt(q2(i,k,jc,nc))*dzt(k) + sum0(i)
          sum1(i) = sqrt(q2(i,k,jc,nc))*dzt(k)*zw(k) + sum1(i)
1490    continue
1500  continue
      alpha=0.2
      do 1600 k=1,km
        do 1590 i=1,imt
          el0        = alpha*abs(sum1(i)/(sum0(i)+epsln))
          el(i,k,jc) = (el0*vk*abs(zw(k)))/(vk*abs(zw(k))+el0)
1590    continue
1600  continue
#endif
#ifdef leq
      do 1700 k=1,km
        do 1690 i=1,imt
          el(i,k,jc) = q2l(i,k,jc,nc)/(q2(i,k,jc,nc)+epsln)
          el(i,k,jp) = q2l(i,k,jp,nc)/(q2(i,k,jp,nc)+epsln)
1690    continue
1700  continue
#endif
c
c-------------------------------------------------------------------
c   compute vertical velocity in q2 columns
c-------------------------------------------------------------------
c
      do 1800 k=1,kmm1
        do 1790 i=2,imtm1
          wt(i,k) = c2dzt(k)*((fuw(i+1,k)-fuw(i,k))
     $                     +fvn(i  ,k)-fvsq(i,k))
1790    continue
1800  continue
      do 1900 k=2,kmm1
        do 1890 i=2,imtm1
          wt(i,k) = wt(i,k)+wt(i,k-1)
1890    continue
1900  continue
c
c     set vertical velocity at the surface to zero (rigidlid)
c     set vertical velocity at the bottom to zero
c
      do 2000 i=2,imtm1
        wt(i,0)        = c0
        wt(i,kmt(i,j)) = c0
2000  continue
c
c------------------------------------------------------------------
c     calculate quantities for the computation of
c     vertical advevtion of tke
c-----------------------------------------------------------------
c
      do 2100 k=1,kmm1
        do 2090 i=2,imtm1
          fwb1(i,k) = (wt(i,k)+wt(i,k+1))
     $               *(q2(i,k,jc,nc)+q2(i,k+1,jc,nc))
#ifdef leq
          fwb2(i,k) = (wt(i,k)+wt(i,k+1))
     $               *(q2l(i,k,jc,nc)+q2l(i,k+1,jc,nc))
#endif
2090    continue
2100  continue
c
c     for vertical advection of tke:
c     set k=0 elements of vertical advection array to
c     (wt(0)+wt(1))*(u*+q2(1))
#ifdef leq
c     for tkel:
c     (wt(0)+wt(1))*(bcl+q2l(1))
#endif
c     also assume wt(km)=0.if "kz" < "km"
c     and there is a bottom slope, w(kz) may be non zero.
c
      do 2200 i=2,imtm1
        fwb1(i,0) = (wt(i,0)+wt(i,1))
     $              *(ustar(i)+q2(i,1,jc,nc))
        fwb1(i,km) = c0
#ifdef leq
        fwb2(i,0) = (wt(i,0)+wt(i,1))
     $              *(bcl(i)+q2l(i,1,jc,nc))
        fwb2(i,km) = c0
#endif
2200  continue
c
c------------------------------------------------------------------
c     construct auxiliary arrays for statement functions to
c     prevent unnecessary computation
c------------------------------------------------------------------
c
      do 2300 k=1,kmm1
        do 2290 i=2,imt
          aux1(i,k) = p5*(fuw(i,k)+fuw(i,k+1))
     $                  *(q2(i,k,jc,nc)+q2(i-1,k,jc,nc))
          aux3(i,k) = dxu2r(i-1)*(q2 (i,k,jc,nm)-q2 (i-1,k,jc,nm))
#ifdef leq
          aux2(i,k) = p5*(fuw(i,k)+fuw(i,k+1))
     $                  *(q2l(i,k,jc,nc)+q2l(i-1,k,jc,nc))
          aux4(i,k) = dxu2r(i-1)*(q2l(i,k,jc,nm)-q2l(i-1,k,jc,nm))
#endif
2290    continue
2300  continue
c
c---------------------------------------------------------------------
c     set horizontal diffusion coefficients
c---------------------------------------------------------------------
c
      bbtj = 8.0*cstr(j)*cstr(j)
      cctj = csu(j  )*dyur(j  )*dytr(j)*cstr(j)
      ddtj = csu(j-1)*dyur(j-1)*dytr(j)*cstr(j)
      bbtj = bbtj*aq
      cctj = cctj*aq
      ddtj = ddtj*aq
c
c   store contribution from  advection of tke into advq2
c   store contribution from hor. diff. of tke into hdq2
c
      do 2400 k=1,kmm1
        do 2390 i=2,imtm1
          advq2(i,k) = uq2x(i,k)  + vq2y(i,k)  + wq2z(i,k)
          hdq2(i,k)  = q2xx(i,k)  + q2yy(i,k)
#ifdef leq
          advq2l(i,k) = uq2lx(i,k) + vq2ly(i,k) + wq2lz(i,k)
          hdq2l(i,k)  = q2lxx(i,k) + q2lyy(i,k)
#endif
2390    continue
2400  continue
c
      do 2500 k=1,kmm1
        do 2490 i=2,imtm1
c
c--------------------------------------------------------------------
c       calculate shear and buoyancy production
c--------------------------------------------------------------------
c
        sprod(i,k) = (p25*(u(i-1,k+1,jm,nc)+u(i-1,k+1,jc,nc)
     $                  +u(i,k+1,jm,nc)+u(i,k+1,jc,nc))
     $             -p25*(u(i-1,k,jm,nc)+u(i-1,k,jc,nc)
     $                  +u(i,k,jm,nc)+u(i,k,jc,nc)))**2
        sprod(i,k) = sprod(i,k)
     $        +(p25*(v(i-1,k+1,jm,nc)+v(i-1,k+1,jc,nc)
     $              +v(i,k+1,jm,nc)+v(i,k+1,jc,nc))
     $        -p25*(v(i-1,k,jm,nc)+v(i-1,k,jc,nc)
     $             +v(i,k,jm,nc)+v(i,k,jc,nc)))**2
        sprod(i,k) = sprod(i,k)*vvc(i,k,jc,nc)*dzwr(k)*dzwr(k)
c       betag      = grav*(8.75e-6*(t(i,k,1)+9.))
c       bprod(i,k) = vdc(i,k,jc,nc)*betag*(t(i,k,jc,nc,1)
c    $                                       -t(i,k+1,jc,nc,1))*dzwr(k)
        bprod(i,k) = vdc(i,k,jc,nc)*grav*drho(i,k)*dzwr(k)
c
c-----------------------------------------------------------------------
c       calculate dissipation
c-----------------------------------------------------------------------
c
        dissp(i,k) = c2*(sqrt(q2(i,k,jc,ln))*q2(i,k,jc,ln))
c       if explicit calc of dissp no need to divide by q2 and use q2b
#ifdef implicitvmix
#  ifdef leq
        dissp(i,k) = dissp(i,k)/(b1*q2l(i,k,jc,ln)+epsln)
#  endif
#  ifdef lalg
        dissp(i,k) = dissp(i,k)/(b1*el(i,k,jc)*q2(i,k,jc,ln)+epsln)
#  endif
#else
        dissp(i,k) = dissp(i,k)/(b1*el(i,k,jc)+epsln)
#endif
c
c-----------------------------------------------------------------------
c       compute partial tendency
c-----------------------------------------------------------------------
c
          tke(i,k) = q2(i,k,jc,nm)
     $                 +tdt*c2*(sprod(i,k)+bprod(i,k))
     $                 -tdt*advq2(i,k)+tdt*hdq2(i,k)
2490    continue
2500  continue
#ifdef implicitvmix
c
c     include vertical diffusion and solve implicitly
c
      call implq (q2a, tke, tdt, ustar, j)
#else
c
c     calculate tke explictly add on dissipation
c     and vertical diffusion
c
      do 2600 k=1,kmm1
        do 2590 i=2,imtm1
          q2dif(i,k) = q2(i,k,jc,nm) - q2(i,k+1,jc,nm)
2590    continue
2600  continue
c
c     for vertical diffusion of q2:
c     set k=0 to b1**2/3*ustar**2
c     set bottom level "kz" =0
c
      do 2700 i=2,imtm1
        q2dif(i,0)        = ustar(i) - q2(i,1,jc,nm)
        q2dif(i,kmt(i,j)) = c0
2700  continue
      do 2800 k=1,km
        do 2790 i=2,imtm1
          q2a(i,k) = tke(i,k)-tdt*dissp(i,k)
     $                     +tdt*q2zz(i,k)
2790    continue
2800  continue
#endif
#ifdef lalg
      do 2900 k=1,km
        do 2890 i=1,imt
          if (q2a(i,k) .lt. small) q2a(i,k) = small
2890    continue
2900  continue
#endif
c
c     mask out land values   q2
c
      do 3000 k=1,km
        do 2990 i=1,imt
          q2a(i,k) = q2a(i,k)*fm(i,k,jc)
2990    continue
3000  continue
#ifdef leq
c
c----------------------------------------------------------------------
c      compute master length scale  (leq)
c----------------------------------------------------------------------
c
      do 3100 k=1,kmm1
        do 3090 i=1,imt
c
c     compute wall proximity function w = 1 + e2(l/L)**2
c
          wallc = (c1/abs(zw(k))+c1/(abs(zw(k)-zw(kmt(i,j)))+deps))/vk
          wallf = c1+e2*(wallc*el(i,k,jc))**2
c
c-----------------------------------------------------------------------
c       calculate dissipation
c-----------------------------------------------------------------------
c
#ifdef implicitvmix
          dissp(i,k) = p5*dissp(i,k)*wallf
#else
          dissp(i,k) = (sqrt(q2(i,k,jc,ln))*q2(i,k,jc,ln)*wallf)/b1
#endif
c
c
c     compute partial tendency of tkel
c
          tkel(i,k) = q2l(i,k,jc,nm)
     $                +tdt*(el(i,k,jc)*e1*(sprod(i,k)+e3*bprod(i,k)))
     $                -tdt*advq2l(i,k)+tdt*hdq2l(i,k)
3090    continue
3100  continue
#ifdef implicitvmix
c
c     include vertical diffusion and solve implicitly
c
      call implq (q2la, tkel, tdt, bcl, j)
#else
c
c     calculate tkel explictly add on dissipation
c     and vertical diffusion
c
      do 3200 k=1,kmm1
        do 3190 i=2,imtm1
          q2ldif(i,k) = q2l(i,k,jc,nm) - q2l(i,k+1,jc,nm)
3190    continue
3200  continue
c
c     for vertical diffusion of q2l:
c     set k=0 to bcl=0
c     set bottom level "kz" =0
c
      do 3300 i=2,imtm1
        q2ldif(i,0)        = bcl(i)
        q2ldif(i,kmt(i,j)) = c0
3300  continue
      do 3400 k=1,km
        do 3390 i=2,imtm1
          q2la(i,k) = tkel(i,k)-tdt*dissp(i,k)
     $                       +tdt*q2lzz(i,k)
3390    continue
3400  continue
#endif
      do 3500 k=1,km
        do 3490 i=1,imt
          if (q2a(i,k) .lt. small .or. q2la(i,k) .lt. small) then
            q2a(i,k)  = small
            q2la(i,k) = small*c100
          endif
3490    continue
3500  continue
c
c     compute max l  (could modify e3 instead) for stable case
c
      do 3600 k=1,km
        do 3590 i=1,imt
          elmax(i,k) = big
3590    continue
3600  continue
c
      fxb=.1
      do 3700 k=1,kmm1
        do 3690 i=1,imt
          if(drho(i,k).le.c0) then
          elmax(i,k) = sqrt(q2a(i,k))
          bprol(i,k) = abs(grav*drho(i,k)*dzwr(k))
          bprol(i,k) = sqrt(bprol(i,k))
          elmax(i,k) = q2a(i,k)*fxb*elmax(i,k)
     $                 /(bprol(i,k)+epsln)
          endif
3690    continue
3700  continue
      do 3800 k=1,km
        do 3790 i=1,imt
          if (q2la(i,k) .gt. elmax(i,k) .and. q2a(i,k) .ne. small)
     $              q2la(i,k) = elmax(i,k)
c
c       mask out land values  q2l
c
          q2la(i,k) = q2la(i,k)*fm(i,k,jc)
3790    continue
3800  continue
#endif
c
c-----------------------------------------------------------------
c      solve for km , kh & kq
c-----------------------------------------------------------------
c
      d1     = a1*(c1-3.*cc1)
      d2     = a2
      a1a16  = a1*a1*6.
      a1a112 = a1*a1*12.
      a1a29  = a1*a2*9.
      a1a26  = a1*a2*6.
      a2b23  = a2*b2*3.
      a1a212 = a1*a2*12.
c     el     = q2la/(q2a+epsln)
      do 3900 k=1,kmm1
        do 3890 i=1,imt
          gh(i,k) = sprod(i,k)*el(i,k,jc)*el(i,k,jc)/
     $            (q2(i,k,jc,nc)*vvc(i,k,jc,nc)+epsln)
3890    continue
3900  continue
c
c      the following filter on gmp1 does not change the calculated
c      results significantly but does eliminate occasional spurious
c      vertical fluctuations in km and kh
c
      do 62 i=1,imt
        gmp1(i,1) = p25*p154+p5*gh(i,1)+p25*gh(i,2)
  62  continue
      do 4000 k=2,kmm1
        do 3990 i=1,imt
          gmp1(i,k) = p25*gh(i,k-1)+p5*gh(i,k)+p25*gh(i,k+1)
3990    continue
4000  continue
      do 4100 k=1,kmm1
        do 4090 i=1,imt
          gh(i,k)   = bprod(i,k)*el(i,k,jc)*el(i,k,jc)/
     $                (q2(i,k,jc,nc)*vdc(i,k,jc,nc)+epsln)
          gh(i,k)   = amin1(gh(i,k),p032)
          gh(i,k)   = amax1(gh(i,k),-.1e9)
          gmp1(i,k) = amin1(gmp1(i,k),p48-c15*gh(i,k))
4090    continue
4100  continue
      do 4200 k=1,kmm1
        do 4190 i=1,imt
          c11     = c1+a1a16*gmp1(i,k)-a1a29*gh(i,k)
          c12     = (-a1a112-a1a29)*gh(i,k)
          c21     = a1a26*gmp1(i,k)
          c22     = c1-(a2b23+a1a212)*gh(i,k)
          denom   = c1/(c11*c22-c12*c21+epsln)
          sm(i,k) = (c22*d1-c12*d2)*denom
          sh(i,k) = (c11*d2-c21*d1)*denom
4190    continue
4200  continue
      do  4300 k=1,kmm1
        do 4290 i=1,imt
          sm(i,k) = (sm(i,k)+abs(sm(i,k)))*p5
          sh(i,k) = (sh(i,k)+abs(sh(i,k)))*p5
4290    continue
4300  continue
      do 4400 k=1,km
        do 4390 i=1,imt
          ekn(i,k,jc) = el(i,k,jc)*sqrt(q2(i,k,jc,nc))
4390    continue
4400  continue
      fx1 = sq/0.392
      do 4500 k=1,kmm1
        do 4490 i=1,imt
          vdqca(i,k) = p25*vdqc(i,k,jc,nc)
     $                +p75*(ekn(i,k,jc)*fx1*sm(i,k))
          vdca(i,k)  = p25*vdc(i,k,jc,nc)
     $                +p75*(ekn(i,k,jc)*sh(i,k))
4490    continue
4500  continue
c
c     average l*q*sm to obtain a value over u v pts
c        (at eknp should use sm at j+1-but not available)
c
      do 4600 k=1,km
        do 4590 i=1,imt
          if (j .eq. jmtm1) then
             ekn(i,k,jp) = c0
          else
            ekn(i,k,jp) = el(i,k,jp)*sqrt(q2(i,k,jp,nc))
          endif
4590    continue
4600  continue
      do 4700 k=1,kmm1
        do 4690 i=1,imtm1
          vvca(i,k) = p25*vvc(i,k,jc,nc)
     $        +p75*(p25*(ekn(i,k,jp)*sm(i,k)+ekn(i+1,k,jp)*sm(i+1,k)
     $             +ekn(i,k,jc)*sm(i,k)+ekn(i+1,k,jc)*sm(i+1,k)))
4690    continue
4700  continue
c
c     set q2=q2l=0. at kz (bottom b.c.) also vvc,vdc=0.
c
      do 4800 i=1,imt
        kz = kmt(i,j)
        if (kz .eq. 0) then
          q2a(i,kz)   = c0
#ifdef leq
          q2la(i,kz)  = c0
#endif
          vdqca(i,kz) = c0
          vdca(i,kz)  = c0
        endif
4800  continue
c
c     limit mixing coeff
c
      do 4900 k=1,km
        do 4890 i=1,imt
          if (vdqca(i,k) .ge. vdqlim) vdqca(i,k) = vdqlim
          if (vvca(i,k) .ge. visc_cbu_limit) vvca(i,k)  = visc_cbu_limit
          if (vdca(i,k) .ge. diff_cbt_limit) vdca(i,k) = diff_cbt_limit
          if (vdqca(i,k) .lt. diff_cbt_back) vdqca(i,k) = diff_cbt_back
          if (vvca(i,k) .lt. visc_cbu_back) vvca(i,k)   = visc_cbu_back
          if (vdca(i,k) .lt. diff_cbt_back) vdca(i,k)   = diff_cbt_back
4890    continue
4900  continue
c
c
#ifdef cyclic
      do 5600 k=1,km
        q2a(1,k)    = q2a(imtm1,k)
        q2a(imt,k)  = q2a(2,k)
#ifdef leq
        q2la(1,k)   = q2la(imtm1,k)
        q2la(imt,k) = q2la(2,k)
#endif
#endif
5600  continue
c
c     set up advection coeff for next call
c
      factor = cst(j)*dyt(j)*cstr(j+1)*dytr(j+1)
      do 5700 k=1,km
        do 5690 i=1,imt
          fvsq(i,k) = fvn(i,k)*factor
5690    continue
5700  continue
c
c     printout some diagnostics
c
      if (tsits) then
c
      kmd2=km/2
      if (j .eq. jpr1 .or. j .eq. jpr2) then
      ii = ipr1
      if (j .eq. jpr2)   ii = ipr2
      if (jpr1 .eq. jpr2)ii = ipr1
      write(stdout,755)
      do 5800 k=1,kmd2
        tendq = (q2a(ii,k)-q2(ii,k,jc,nm))/tdt
        dissp(ii,k) = sqrt(q2(ii,k,jc,nm))*q2(ii,k,jc,nm)
        dissp(ii,k) = dissp(ii,k)/(b1*el(ii,k,jc)+epsln)
#ifdef implicitvmix
        diff = advq2(ii,k)-c2*sprod(ii,k)-c2*bprod(ii,k)+tendq
     1     +c2*dissp(ii,k)-hdq2(ii,k)
#else
        diff = q2zz(ii,k)
#endif
        rf = -bprod(ii,k)/(sprod(ii,k)+epsln)
        write(stdout,750)tendq,advq2(ii,k),diff,sprod(ii,k),bprod(ii,k)
     1         ,dissp(ii,k),vdca(ii,k),vvca(ii,k),rf,q2a(ii,k)
5800  continue
#ifndef leq
      write(stdout,777)
      do 5900 k=1,kmd2
        write(stdout,778)hdq2(ii,k)
5900  continue
#endif
#ifdef leq
      write(stdout,777)
      do 6000 k=1,kmd2
        pod = b1*(sm(ii,k)*gmp1(ii,k)+sh(ii,k)*gh(ii,k))
        rii = -gh(ii,k)/(gmp1(ii,k)+epsln)
        rif = -sh(ii,k)*gh(ii,k)/(sm(ii,k)*gmp1(ii,k)+epsln)
        write(stdout,778)hdq2(ii,k),hdq2l(ii,k),el(ii,k,jc),elmax(ii,k)
     1         ,rii,rif,pod
6000  continue
      write(stdout,767)
      do 6100 k=1,kmd2
        tendq = (q2la(ii,k)-q2l(ii,k,jc,nm))/tdt
        wallf = c1+e2*((c1/abs(zw(k))+c1/(abs(zw(k)-zw(kmt(i,j)))
     $                                    +deps)) *el(ii,k,jc)/vk)**2
        dissp(ii,k) = sqrt(q2(ii,k,jc,ln))*q2(ii,k,jc,ln)
        dissp(ii,k) = (dissp(ii,k)/b1)*wallf
        diff1 = (q2la(ii,k)-tkel(ii,k))/tdt+dissp(ii,k)
#ifdef implicitvmix
        diff2 = tendq+dissp(ii,k)-el(ii,k,jc)*e1*(sprod(ii,k)+
     $       bprod(ii,k))+advq2l(ii,k)-hdq2l(ii,k)
#else
        diff2 = q2lzz(ii,k)
#endif
        write(stdout,765)tendq,dissp(ii,k),diff1,diff2,tke(ii,k)
     $                ,vdqca(ii,k),advq2l(ii,k),q2la(ii,k),wallf
  765 format(5x,9e12.4)
  767 format(10x,'tendql=',7x,'dissp',7x,'diff1',7x,'diff2',7x,'tke '
     $      ,7x,' vdqc ',7x,'advq2l',6x,'q2la',7x,'wallf')
6100  continue
#endif
      endif
      endif
      reltim = relyr - deltat
      if (prxzts .and. eots) then
c
c---------------------------------------------------------------------
c       print the solution at specified latitudes on diagnostic tsteps
c---------------------------------------------------------------------
c
        do 6200 jlat=1,nlatpr
          jj = indp (prlat(jlat), yt, jmt)
          if (jj .eq. j .and. prlat(jlat) .le. yt(jmt)) then
            is  = indp (prslon(jlat), xt, imt)
            ie  = indp (prelon(jlat), xt, imt)
            ks  = indp (prsdpt(jlat), zt, km)
            ke  = indp (predpt(jlat), zt, km)
            fx  = 1.0e-2
            scl = c1
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'q2', itt
     $,        j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (q2(1,1,jc,nc), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
              io = iabs(ioprxz)
              write (stdout,*) ' => Zonal q2 slice: lat=',yt(j)
     $,       ' written unformatted to unit', io, ' on ts=',itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (ioprxz) imt, km, reltim
              write(iotext,'(a7,i4)') ' for j=',j
              iotext(12:) = ': read (ioprxz) ((q2(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, q2(1,1,jc,nc), imt*km)
            endif
c
#ifdef leq
            scl = c1
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'q2l ', itt
     $,        j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (q2l(1,1,jc,nc), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
              io = iabs(ioprxz)
              write (stdout,*) ' => Zonal q2l slice: lat=',yt(j)
     $,       ' written unformatted to unit', io, ' on ts=',itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (ioprxz) imt, km, reltim
              write(iotext,'(a7,i4)') ' for j=',j
              iotext(12:) =': read (ioprxz) ((q2l(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, q2l(1,1,jc,nc), imt*km)
            endif
c
#endif
            scl = c1
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'kh ', itt
     $,        j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (vdc(1,1,jc,nc), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
              io = iabs(ioprxz)
              write (stdout,*) ' => Zonal kh slice: lat=',yt(j)
     $,       ' written unformatted to unit', io, ' on ts=',itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (ioprxz) imt, km, reltim
              write(iotext,'(a7,i4)') ' for j=',j
              iotext(12:) =': read (ioprxz) ((vdc(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, vdc(1,1,jc,nc), imt*km)
            endif
c
          endif
6200    continue
      endif
   72 format(5x,'sm,sh,gmp1,gh,km,kh= ',6e12.4)
  750 format(5x,10e12.4)
  755 format(10x,'tendq=',7x,'advq',7x,'+diff',7x,'+prods',7x,'+prodb'
     1      ,7x,'-dissp',9x,'kh',9x,'km',9x,'rf',9x,'q2')
  777 format(10x,'hdiff ',7x,'hdiffl',7x,'   l',7x,'  lmax',7x,
     1 '  rii',7x,'   rif',7x,'   pod')
  778 format(5x,7e12.4)
9100  format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     $,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     $,', scaling=',1pg10.3)
#endif
      return
      end
