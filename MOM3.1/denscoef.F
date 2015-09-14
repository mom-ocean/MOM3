#ifdef drive_dencoef

# define driver_only
# include "grids.F"
# include "iomngr.F"
# include "util.F"
# include "size_check.F"

      program driver
#include "size.h"
#include "stdunits.h"
#include "coord.h"
#include "state.h"
c
      stdout = 6
c
c-----------------------------------------------------------------------
c     construct the grid.
c-----------------------------------------------------------------------
c
      maxlen = max(imt, jmt, km)
      call gcoord (maxlen, imt2, jmt2, km2
     &,            dxtdeg, dytdeg, dxudeg, dyudeg
     &,            dzt, dzw, xt, xu, yt, yu, zt, zw)
c
c     check that returned grid sizes match those in file "size.h"
c
      call size_check (imt2, jmt2, km2, 'eqstate', 'stop')
c
c-----------------------------------------------------------------------
c     construct coefficients for equation of state
c-----------------------------------------------------------------------
c
      call eqstate (zt, km, ro0, to, so, c, tmink, tmaxk, smink, smaxk)
c
      cksum = checksum(c, km, 9)
      write (stdout
     &,'(/6x,"Checksum for density coefficients =",e14.7)') cksum
c
      write (stdout,'(/a/a/a)') 
     & 'Note: it is no longer necessary to create density coeffs'
     &,'      before running the model because the model calculates'
     &,'      them.'
      stop
      end
#endif

      subroutine eqstate (zt, km, ro0_den, to_den, so_den, c_den
     &,                   tmink, tmaxk, smink, smaxk)
c
c=======================================================================
c
c             E Q U A T I O N   O F   S T A T E   M O D U L E
c
c     Calculate polynomial coefficients for density computations in MOM.
c
c     To generate the coefficients: 
c     
c        1) set up the grid in "grids.F" module
c
c        2) compile and run this module by setting the options in and
c           executing the "run_denscoef" script
c
c     To install the coefficients in MOM: 
c
c        3) follow the directions at the end of the output from 2)
c
c
c     This program calculates the 9 coefficients of a third order
c     polynomial approximation to the equation of state for sea water.
c       The program yields coefficients that will compute density as a
c     function of temperature, and salinity, at predetermined depths,
c     as used in the MOM subroutine "state".
c     More specifically, the densities calculated from the polynomial
c     formula are in the form of sigma anomalies.  The method is
c     taken from that described by Bryan & Cox (1972).
c       By default, the program uses the equation of state set by the
c     Joint Panel on Oceanographic Tables & Standards (UNESCO, 1981)
c     an described by Gill (1982).  An option exists to use the older
c     Knudsen-Ekman equation of state, as described by Fofonoff (1962),
c     if the user prefers.
c       Subroutine "lsqsl2" performs the iterative least-squares
c     polynomial fitting for the overdetermined system.  The algorithm
c     is outlined by Hanson and Lawson (1969), and the code looks as if
c     it has not be touched since that time.
c
c     references:
c        Bryan, K. & M. Cox, An approximate equation of state
c          for numerical models of ocean circulation, J. Phys.
c          Oceanogr., 2, 510-514, 1972.
c        Fofonoff, N., The Sea: Vol 1, (ed. M. Hill). Interscience,
c          New York, 1962, pp 3-30.
c        Gill, A., Atmosphere-Ocean Dynamics: International Geophysical
c          Series No. 30. Academic Press, London, 1982, pp 599-600.
c        Hanson, R., & C. Lawson, Extensions and applications of the
c          Householder algorithm for solving linear least squares
c          problems. Math. Comput., 23, 787-812, 1969.
c        UNESCO, 10th report of the joint panel on oceanographic tables
c          and standards. UNESCO Tech. Papers in Marine Sci. No. 36,
c          Paris, 1981.
c
c    ifdef options:
c       "knudsen"
c       To over-ride the default of using the UNESCO equation of state
c     and to instead employ the Knudsen-Ekman formula.
c       "insitu"
c       If the user desires the polynomial approximations to calculate
c     density as a function of in situ temperature, salinity, and depth,
c     then the ifdef option "insitu" must be defined.  Otherwise, the
c     default assumption is that potential temperatures will be used (as
c     in the ocean model code).
c       "extras"
c       If the user wishes to have a detailed report of the inputs and
c     results of the curve fitting processes written to the standard
c     output unit (stdout), then the ifdef option "extras" should be
c     defined.  The default is for a rather short summary to be written.
c
c
c     author:          k. dixon              e-mail=> kd@gfdl.gov
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
c
c-----------------------------------------------------------------------
c
#include "stdunits.h"
      parameter (kmax=200)
      character*10 fname
c
      parameter (kx = 5, kxx = 2*kx, kk = kx*kxx )
      parameter (ksdim = kk+72 , krdim = kk+36 )
c
      dimension a(kk,9), sigma(kk), sigman(kk), c(kk,9), x(9),
     &          sb(ksdim), r(krdim)
      dimension tmin(kmax), smin(kmax), tmax(kmax), smax(kmax),
     &          z(kmax), dd(kmax), ss(kmax), ab(13,kmax), ts(33,4),
     &          ta(kxx), sa(kxx), tp(kk), sp(kk), th(kk)
c
      real ro0_den, to_den, so_den, c_den, zt
      real tmink, tmaxk, smink, smaxk
      dimension ro0_den(km), to_den(km), so_den(km), c_den(km,9)
      dimension tmink(km), tmaxk(km), smink(km), smaxk(km)
      dimension zt(km)
c
      real realz
c
      data fname /'dncoef.new'/
c
c
c  enter bounds for polynomial fit: at 33 levels from sfc to 8000 m.
c           ts(k,1)=lower bnd of t at z=(k-1)*250 meters
c           ts(k,2)=upper bnd of t          
c           ts(k,3)=lower bnd of s          
c           ts(k,4)=upper bnd of s          
c  The user should review the appropriateness of the "ts" values set
c  below, and modify them if the intended modelling application could
c  be expected to yield temperature and salinity values outside of the
c  "ts" ranges set by default.
c
#ifdef bbl_uniform_ts_range_in_state
      data (ts(k,1),k=1,33) / 33*-2.0/
      data (ts(k,2),k=1,33) / 33*30.0/
      data (ts(k,3),k=1,33) / 33*28.0/
      data (ts(k,4),k=1,33) / 33*37.0/
#else
      data (ts(k,1),k=1,33) / 4*-2.0, 15*-1.0, 14*0.0 /
      data (ts(k,2),k=1,33) / 29.0, 19.0, 14.0, 11.0, 9.0, 28*7.0 /
      data (ts(k,3),k=1,33) / 28.5, 33.7, 34.0, 34.1, 34.2, 34.4,
     &                        2*34.5, 15*34.6, 10*34.7 /
      data (ts(k,4),k=1,33) / 37.0, 36.6, 35.8, 35.7, 35.3, 2*35.1,
     &                        26*35.0 /
#endif
c
c     z       = model levels (midpoint of model layers)
c     tmin, tmax, smin, smax = minimum and maximum in situ temperature
c               and salinity values which define the ranges to be used
c               when computing the polynomials at each model level
c     dd, ds  = increment between temperature and salinity values at
c               each model level to be used in constructing array of
c               temperature, salinity and density for curve fitting
c     ta, sa  = in situ temperature and salinity values available for
c               constructing array of data for curve fitting at each
c               model level
c     tp, sp  = in situ temperature and salinity values constructed from
c               all combinations of ta & sa
c     th      = potential temperature values associated with "tp" at a
c               given level and salinity
c     t1, s1, tot1, th1 = level mean insitu temp., salinity, density,
c               and potential temp. used in polynomial fitting
c     tot     = density (in sigma units) calculate from t1 and s1 at a
c               given model level
c     sigma   = insitu densities (in sigma units) calculated from "tp"
c               and "sp" values
c     sigman  = insitu density anomalies at a given level (formed by
c               subtracting "tot" from sigma)
c     tanom, sanom = temperature and salinity anomalies used in loading
c               array "a" for use in lsqsl2 curve fitting
c     x       = the 9 polynomial coefficients
c     r, sb   = used only in lsqsl2
c
c=======================================================================
c
#if defined read_denscoef || defined save_denscoef
      call ioinit
#endif
#ifdef read_denscoef
      write (stdout,'(/a/)')
     & ' =>Reading denscoef fields from file "denscoef.dta" '
      call getunit (io, 'denscoef.dta'
     &,             'unformatted sequential rewind')
c
      read (io) ro0_den
      read (io) to_den
      read (io) so_den
      read (io) c_den
      read (io) tmink
      read (io) tmaxk
      read (io) smink
      read (io) smaxk
      call relunit (io)
#else
      if (km .gt. kmax) then
        write (stdout,*) '=>Error: increase "kmax" > ',km,' in eqstate'
        call abort()
      endif
c
c     set some constants
c
      c0 = 0.0
      c1 = 1.0
      c2 = 2.0
c
c     construct depths (meters) from surface to midpoint of levels
c
      cmtocm = 1.0d-2
      do k=1,km
        z(k) = zt(k) * cmtocm
        if (z(k) .gt. 8000.0) then
          write (stdout,*) '=>Error:depth can`t exceed 8000m in eqstate'
          call abort()
        endif
      enddo
c
c     set the temperature and salinity ranges to be used for each
c     model level when performing the polynomial fitting
c
      do k=1,km
        realz = z(k)/250.0
        i = ifix (realz) + 1
        tmin(k) = ts(i,1)
        tmax(k) = ts(i,2)
        smin(k) = ts(i,3)
        smax(k) = ts(i,4)
      enddo
c
c  set temperature and salinity increments to be used in creating
c  curve fitting array at each level (twice as many temperature values
c  than salinity values)
c
      fkx = kx
      do k=1,km
        dd(k) = (tmax(k)-tmin(k)) / (c2*fkx-c1)
        ss(k) = (smax(k)-smin(k)) / (fkx-c1)
      enddo
c
c     loop over all model levels
c
      do k=1,km
c
        do i=1,kxx
          fi = i
          ta(i) = tmin(k) + (fi-c1)*dd(k)
          sa(i) = smin(k) + (fi-c1)*ss(k)
        enddo
c
c       load the "kxx" combinations of the 2*"kx" insitu temp. and "kx"
c       salinity values into "tp" and "sp"
c
        do i=1,kxx
          do j=1,kx
            ka = kx*i + j - kx
            tp(ka) = ta(i)
            sp(ka) = sa(j)
          enddo
        enddo
c
        t1  = c0
        s1  = c0
        tot = c0
        th1 = c0
        fkk = kk
c
c       calculate insitu density "sigma" for each t,s combination at
c       this depth "d"
c
        do ka=1,kk
          d = z(k)
          s = sp(ka)
          t = tp(ka)
c
#ifdef knudsen
c         "knuekm" returns density (in sigma units) from insitu
c         temperature, salinity, & depth (pressure) using the 
c         Knudsen-Ekman formula
c
          call knuekm(t,s,d,densit)
c
          sigma(ka) = densit
#else
c         "unesco" returns density (kg per m**3) from insitu
c         temperature, salinity, & depth (pressure) using the UNESCO
c         equation of state
c
          call unesco(t,s,d,densit)
c
          sigma(ka) = densit - 1.0d3 + 2.5d-2
#endif
c
c         "potem" returns potential temp. from from insitu temperature,
c         salinity, & depth (pressure)
c
          call potem(t,s,d,theta)
c
          th(ka) = theta
          t1 = t1 + tp(ka)
          s1 = s1 + sp(ka)
          tot = tot + sigma(ka)
          th1 = th1 + th(ka)
        enddo
c
c       form layer averages "t1", "s1", "th1", and "tot1", and compute
c       reference density "tot" from "t1" and "s1" at this depth "d"
c
        t1 = t1/fkk
        s1 = s1/fkk
        th1 = th1/fkk
        tot1 = tot/fkk
#ifdef knudsen
c
c       "knuekm" returns density from insitu temp., salinity, & depth
c       (pressure) using the Knudsen-Ekman formula
c
        call knuekm (t1, s1, d, densit)
c
        tot = densit
#else
c
c       "unesco" returns density from insitu temp., salinity, & depth
c       (pressure) using the UNESCO equation of state
c
        call unesco (t1, s1, d, densit)
        tot = densit - 1.0d3 + 2.5d-2

#endif
c
#ifdef extras
c
c       define "extras" for voluminous printout of calculation info.
c
        write (stdout,'(a49)')
     &  ' insitu temperatures used in polynomial fit & avg'
        write (stdout, 9071) kk, (tp(ka),ka=1,kk)
        write (stdout, 9072) t1, k
        write (stdout,'(a40)')
     &  ' salinities used in polynomial fit & avg'
        write (stdout, 9071) kk, (sp(ka),ka=1,kk)
        write (stdout, 9072) s1, k
        write (stdout,'(a53)')
     &  ' densities (sigma units) used in polynomial fit & avg'
        write (stdout, 9071) kk, (sigma(ka),ka=1,kk)
        write (stdout, 9072) tot1, k
        write (stdout,'(a54)')
     &  ' density calculated from level avg insitu T & salinity'
        write (stdout, 9072) tot, k
        write (stdout,'(a52)')
     &  ' potential temperatures used in polynomial fit & avg'
        write (stdout, 9071) kk, (th(ka),ka=1,kk)
#ifndef bbl_uniform_ts_range_in_state
        write (stdout, 9072) th1, k
#else
        write (stdout, 9072) t1, k
#endif
#ifdef insitu
        write (stdout,'(a47)')
     &  ' >> insitu temps were used in polynomial fit <<'
#else
        write (stdout,'(a50)')
     &  ' >> potential temps were used in polynomial fit <<'
#endif
c
#endif
#ifndef insitu
c
c       define insitu if using insitu temperatures (removes this line)
c
#ifndef bbl_uniform_ts_range_in_state
        t1 = th1
#endif
c
#endif
c
c       begin loading "ab" array with level averages
c
        ab(1,k) = z(k)
        ab(2,k) = tot
        ab(3,k) = t1
        ab(4,k) = s1
c
        do ka=1,kk
#ifndef insitu
c
c         define insitu (removes this line) if using insitu temperatures
c
          tp(ka) = th(ka)
#endif
c
c         create anomalies for temperature, salinity & density and
c         load work array "a" with the anomalies and their products
c
          tanom = tp(ka) - t1
          sanom = sp(ka) - s1
          sigman(ka) = sigma(ka) - tot
          a(ka,1) = tanom
          a(ka,2) = sanom
          a(ka,3) = tanom * tanom
          a(ka,4) = tanom * sanom
          a(ka,5) = sanom * sanom
          a(ka,6) = a(ka,3) * tanom
          a(ka,7) = a(ka,5) * tanom
          a(ka,8) = a(ka,3) * sanom
          a(ka,9) = a(ka,5) * sanom
        enddo
c
c       set the arguments used in call to "lsqsl2"
c       ndim = first dimension of array a
c       nrow =number of rows of array a
c       ncol = number of columns of array a
c       in = option number of lsqsl2
c       itmax = number of iterations
c
        ndim = 50
        nrow = kk
        ncol = 9
        in = 1
        itmax = 4
c
        it = 0
        ieq = 2
        irank = 0
        eps = 1.0e-7
        nhdim = 9
c
c       LSQL2 is  a Jet Propulsion Laboratory subroutine that
c       computes the least squares fit in an iterative manner for
c       overdetermined systems.
c
        call lsqsl2 (ndim, a, nrow, ncol, sigman, x, irank, in, itmax,
     &               it, ieq, enorm, eps, nhdim, h, c, r, sb)
c
#ifdef extras
        write (stdout, 9081) k, (x(i),i=1,9)
        write (stdout, 9082) tot
        write (stdout, '(///)')
c
#endif
        do i=1,ncol
          ab(i+4,k) = x(i)
        enddo
c
c       end loop on model levels
c
      enddo
c
      nn = ncol + 4
#ifdef extras
      write (stdout, 9091)
      write (stdout, 9092) ((ab(i,j),i=1,nn),j=1,km)
#endif
c
#if !defined salinity_psu
      do k=1,km
        ab(2,k)  = 1.e-3 * ab(2,k)
        ab(4,k)  = 1.e-3 * ab(4,k) - 0.035
        ab(5,k)  = 1.e-3 * ab(5,k)
        ab(7,k)  = 1.e-3 * ab(7,k)
        ab(10,k) = 1.e-3 * ab(10,k)
        ab( 9,k) = 1.e+3 * ab( 9,k)
        ab(11,k) = 1.e+3 * ab(11,k)
        ab(13,k) = 1.e+6 * ab(13,k)
      enddo
#else 
      do k=1,km
        ab(2,k) = 1.e-3 * ab(2,k)
        ab(5,k) = 1.e-3 * ab(5,k)
        ab(6,k) = 1.e-3 * ab(6,k)
        ab(7,k) = 1.e-3 * ab(7,k)
        ab(8,k) = 1.e-3 * ab(8,k)
        ab(9,k) = 1.e-3 * ab(9,k)
        ab(10,k)= 1.e-3 * ab(10,k)
        ab(11,k)= 1.e-3 * ab(11,k)
        ab(12,k)= 1.e-3 * ab(12,k)
        ab(13,k)= 1.e-3 * ab(13,k)
      enddo
#endif
c
c     save this data for a model execution converting insitu
c     temperatures to potential temperatures
c
      do k=1,km
        ro0_den(k) = ab(2,k)
        to_den(k)  = ab(3,k)
        so_den(k)  = ab(4,k)
        d = z(k)
        s = smin(k)
        t = tmin(k)
        call potem(t,s,d,theta)
        tmink(k)   = theta
        s = smax(k)
        t = tmax(k)
        call potem(t,s,d,theta)
        tmaxk(k)   = theta
        smink(k)   = smin(k)
        smaxk(k)   = smax(k)
        do i=5,13
          c_den(k,i-4) = ab(i,k)
        enddo
      enddo

#ifdef linearized_density
c
c     Density a linear function of temp. No pressure effects.
c     Take following numbers from Appendix 3 of Gill (1982)
c     with T=10C, S=35psu, p=0atm. 
c     Recall that unless salinity_psu is enabled, 
c     model salinity is S(model) = (S(psu)-35)/1000
c     
      do k=1,km
        c_den(k,1) = -2.55e-4
        ab(1+4,k)   = c_den(k,1)
        ro0_den(k) = 0.025022
        ab(2,k)    = ro0_den(k)
        to_den(k)  = 19.0
        ab(3,k)    = to_den(k) 
# if defined salinity_psu
        so_den(k)  = 35.0
# else
        so_den(k)  = 0.0
# endif
        ab(4,k)    = so_den(k)
        do i=2,9
           c_den(k,i) = c0
           ab(i+4,k) = c_den(k,i)
        enddo
      enddo
#endif
c
      write (stdout,9531)
      write (stdout,9532) 
     & (i,z(i),tmin(i),tmax(i),smin(i),smax(i),tmink(i),tmaxk(i),i=1,km)
c
#ifdef drive_dencoef
      write (stdout,9533)
c
c     show "ro0", "to" & "so" data statements
c
      do nx=2,4
        if (nx .eq. 2) write(stdout,9502) 'data ro0 /'
        if (nx .eq. 3) write(stdout,9502) 'data to /'
        if (nx .eq. 4) write(stdout,9502) 'data so /'
        n = 0
        do ii=1,1000
          is = n+1
          ie = n+5
          if (ie .lt. km) then
            write(stdout,9510) (ab(nx,i),i=is,ie)
            n = ie
          else
            ie = km
            n = ie-is+1
            if (n .eq. 1) write(stdout,9511) (ab(nx,i),i=is,ie)
            if (n .eq. 2) write(stdout,9512) (ab(nx,i),i=is,ie)
            if (n .eq. 3) write(stdout,9513) (ab(nx,i),i=is,ie)
            if (n .eq. 4) write(stdout,9514) (ab(nx,i),i=is,ie)
            if (n .eq. 5) write(stdout,9515) (ab(nx,i),i=is,ie)
            goto 600
          endif
        enddo
600     continue
      enddo
      do k=1,km
        write(stdout,9521) k
        write(stdout,9522) (ab(i,k),i=5,8)
        write(stdout,9522) (ab(i,k),i=9,12)
        write(stdout,9523) ab(13,k)
      enddo
#endif
c
#ifdef knudsen
      write (stdout,'(/a/a/)')
     & 'Note: density coefficients were calculated using the "knudsen"'
     &,'      formulation.' 
#else
      write (stdout,'(/a/a/)')
     & 'Note: density coefficients were calculated using the UNESCO'
     &,'      formulation.' 
#endif
#ifdef save_denscoef
      write (stdout,'(/a/)')
     & ' =>Writing denscoef fields to file "denscoef.dta" '
      call getunit (io, 'denscoef.dta'
     &,             'unformatted sequential rewind')
c
      write (io) ro0_den
      write (io) to_den
      write (io) so_den
      write (io) c_den
      write (io) tmink
      write (io) tmaxk
      write (io) smink
      write (io) smaxk
      call relunit (io)
#endif
c
c =====================================================================
c
#endif
      return
c
#ifdef extras
 9091 format(//,
     &' calculating coefficients for "MOM" density computations'/
     &'   z    sig0    t    s       x1          x2         ',
     &'x3          x4          x5          x6          x7          x8',
     &'          x9',/)
 9092 format(//,f5.0,f8.4,f5.1,f6.2,9e12.5)

 9071 format(/' kk = # of pts going into interpltn =',i4,/
     &       (1x,5e14.7))
 9072 format(5x,' avg =',e14.7,' for level ',i4,/)
 9081 format(' model level ',i3,': before scaling (x(i),i=1,9)='/
     &        1x,5e14.7,/,1x,4e14.7)
 9082 format(' reference sigma, about which density anomalies are ',
     &       'computed'/1x,e14.7)
#endif
c
#ifdef drive_dencoef
 9502 format(6x,a)
 9510 format(5x,'&',8x,5(f10.7,','))
 9511 format(5x,'&',8x,f10.7,'/',/'c')
 9512 format(5x,'&',8x,f10.7,',',f10.7,'/',/'c')
 9513 format(5x,'&',8x,2(f10.7,','),f10.7,'/',/'c')
 9514 format(5x,'&',8x,3(f10.7,','),f10.7,'/',/'c')
 9515 format(5x,'&',8x,4(f10.7,','),f10.7,'/',/'c')
 9521 format(6x,'data (c(',i3,',n),n=1,9)/')
 9522 format(5x,'&',9x,4(e13.7,','))
 9523 format(5x,'&',9x,e13.7,'/',/,'c')
#endif
 9531 format('    density coefficients were calculated  ',
#ifdef knudsen
     &       /'     employing the Knudsen-Ekman equation of state ',
#else
     &       /'    (employing the UNESCO equation of state)',
#endif
     &       /'    and are valid for the following depths and',
     &       ' both INSITU and POTENTIAL T and S ranges'
     &       /' ',t7,'k',t14,'depth',t27,'tmin',t37,
     &       'tmax',t52,'smin',t62,'smax',t70,'Pot tmin',t80,'Pot tmax')
 9532 format(' ',t5,i3,t12,f7.2,'e2',t25,f7.3,t35,f7.3,t50,f7.4,
     &       t60,f7.4,t70,f7.3,t80,f7.3)
 9533 format(' '/'   the 3rd order polynomial will return density ',
     &       'departures [gm/cm**3] as',/,'   a function of',
#ifdef insitu
     &       ' insitu ',
#else
     &       ' potential ',
#endif
     &       'temperature [deg C] & salinity [model units]'
     &/'   density coefficients are:'/)
 9534 format('   ',8f8.4)
c
      end




      subroutine knuekm (t, s, d, rho)
c=======================================================================
c     this subroutine calculates the density of seawater using the
c     Knudsen-Ekman equation of state.
c
c     input [units]:
c       in-situ temperature (t): [degrees centigrade]
c       salinity (s): [per mil]
c       depth (d): [meters of depth, to approximate pressure]
c     output [units]:
c       density (rho): sigma units
c
c     reference:
c        Fofonoff, N., The Sea: Vol 1, (ed. M. Hill). Interscience,
c          New York, 1962, pp 3-30.
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c
      t2 = t*t
      t3 = t2*t
      s2 = s*s
      s3 = s2*s
      f1 = -1.0d0 * (t - 3.98d0)**2 * (t + 2.83d2) /
     &     (5.0357d2*(t + 6.726d1))
      f2 = t3*1.0843d-6 - t2*9.8185d-5 + t*4.786d-3
      f3 = t3*1.6670d-8 - t2*8.1640d-7 + t*1.803d-5
      fs = s3*6.76786136d-6 - s2*4.8249614d-4 + s*8.14876577d-1
c
      sigma= f1 + (fs + 3.895414d-2)*
     &      (1.0d0 - f2 + f3*(fs - 2.2584586d-1))
c
      a= d*1.0d-4*(1.055d2 + t*9.50d0 - t2*1.58d-1 - d*t*1.5d-4)  -
     &   (2.27d2 + t*2.833d1 - t2*5.51d-1 + t3*4.0d-3)
      b1 = (fs - 2.81324d1)*1.d-1
      b2 = b1 * b1
      b  = -b1* (1.473d2 - t*2.72d0 + t2*4.0d-2 - d*1.0d-4*
     &     (3.24d1 - 0.87d0*t + 2.0d-2*t2))
      b  = b + b2*(4.5d0 - 1.0d-1*t - d*1.0d-4*(1.8d0 - 6.0d-2*t))
      co = 4.886d3/(1.0d0 + 1.83d-5*d)
c
      alpha = d*1.0d-6*(co + a + b)
c
      rho = (sigma + alpha)/(1.d0 - 1.0d-3*alpha)
c
      return
      end


      subroutine lsqsl2 (ndim,a,d,w,b,x,irank,in,itmax,it,ieq,enorm,eps1
     &,nhdim,h,aa,r,s)
c
c     this routine is a modification of lsqsol. march,1968. r. hanson.
c     linear least squares solution
c
c     this routine finds x such that the euclidean length of
c     (*) ax-b is a minimum.
c
c     here a has k rows and n columns, while b is a column vector with
c     k components.
c
c     an orthogonal matrix q is found so that qa is zero below
c     the main diagonal.
c     suppose that rank (a)=r
c     an orthogonal matrix s is found such that
c     qas=t is an r x n upper triangular matrix whose last n-r columns
c     are zero.
c     the system tz=c (c the first r components of qb) is then
c     solved. with w=sz, the solution may be expressed
c     as x = w + sy, where w is the solution of (*) of minimum euclid-
c     ean length and y is any solution to (qas)y=ty=0.
c
c     iterative improvements are calculated using residuals and
c     the above procedures with b replaced by b-ax, where x is an
c     approximate solution.
c
      implicit real*8 (a-h,o-z)
c
      real*8 sj,dp,up,bp,aj
      logical erm
      integer d,w
#include "stdunits.h"
c
c     in=1 for first entry.
c                   a is decomposed and saved. ax-b is solved.
c     in = 2 for subsequent entries with a new vector b.
c     in=3 to restore a from the previous entry.
c     in=4 to continue the iterative improvement for this system.
c     in = 5 to calculate solutions to ax=0, then store in the array h.
c     in  =  6   do not store a  in aa.  obtain  t = qas, where t is
c     min(k,n) x min(k,n) and upper triangular. now return.do not obtain
c     a solution.
c     no scaling or column interchanges are performed.
c     in  =  7   same as with  in = 6  except that soln. of min. length
c                is placed into x. no iterative refinement.  now return.
c     column interchanges are performed. no scaling is performed.
c     in  = 8    set addresses. now return.
c
c     options for computing  a matrix product   y*h  or  h*y are
c     available with the use of the entry points  myh and mhy.
c     use of these options in these entry points allow a great saving in
c     storage required.
c
c
c     dimension a(ndim,ndim),b(1),aa(d,w),s(1), x(1),h(nhdim,nhdim),r(1)
      dimension a(d,w),b(d),aa(d,w),s(d+72),x(w),h(nhdim,nhdim),r(d+36)
#if !defined read_denscoef
c     d = depth of matrix.
c     w = width of matrix.
      k=d
      n=w
      erm=.true.
c
c     if it=0 on entry, the possible error message will be suppressed.
c
      if (it.eq.0) erm=.false.
c
c     ieq = 2      if column scaling by least max. column length is
c     to be performed.
c
c     ieq = 1       if scaling of all components is to be done with
c     the scalar max(abs(aij))/k*n.
c
c     ieq = 3 if column scaling as with in =2 will be retained in
c     rank deficient cases.
c
c     the array s must contain at least max(k,n) + 4n + 4min(k,n) cells
c        the   array r must contain k+4n s.p. cells.
c
      data eps2/1.d-16/
c     the last card controls desired relative accuracy.
c     eps1  controls  (eps) rank.
c
      isw=1
      l=min0(k,n)
      m=max0(k,n)
      j1=m
      j2=n+j1
      j3=j2+n
      j4=j3+l
      j5=j4+l
      j6=j5+l
      j7=j6+l
      j8=j7+n
      j9=j8+n
      lm=l
      if (irank.ge.1.and.irank.le.l) lm=irank
      if (in.eq.6) lm=l
      if (in.eq.8) return
c
c     return after setting addresses when in=8.
c
      go to (10,360,810,390,830,10,10), in
c
c     equilibrate columns of a (1)-(2).
c
c     (1)
c
   10 continue
c
c     save data when in = 1.
c
      if (in.gt.5) go to 30
        do j=1,n
          do i=1,k
            aa(i,j)=a(i,j)
          enddo
        enddo
   30 continue
      if (ieq.eq.1) go to 60
        do j=1,n
!          am=0.e0
          am=1.e-30
          do i=1,k
            am = max(am,abs(a(i,j)))
          enddo
c
c         s(m+n+1)-s(m+2n) contains scaling for output variables.
c
          n2=j2+j
          if (in.eq.6) am=1.d0
          s(n2)=1.d0/am
          do i=1,k
            a(i,j)=a(i,j)*s(n2)
          enddo
        enddo
        go to 100
60      continue
!        am=0.d0
          am=1.e-30
        do j=1,n
          do i=1,k
            am= max(am,abs(a(i,j)))
          enddo
        enddo
        am=am/float(k*n)
        if (in.eq.6) am=1.d0
        do j=1,n
          n2=j2+j
          s(n2)=1.d0/am
        enddo
        do j=1,n
          n2=j2+j
          do i=1,k
            a(i,j)=a(i,j)*s(n2)
          enddo
        enddo
c       compute column lengths with d.p. sums finally rounded to s.p.
c
c     (2)
c
100     continue
        do j=1,n
          n7=j7+j
          n2=j2+j
          s(n7)=s(n2)
        enddo
c
c       s(m+1)-s(m+ n) contains variable permutations.
c
c       set permutation to identity.
c
        do j=1,n
          n1=j1+j
          s(n1)=j
        enddo
c
c       begin elimination on the matrix a with orthogonal matrices .
c
c       ip=pivot row
c
        do ip=1,lm
          dp=0.d0
          km=ip
          do j=ip,n
            sj=0.d0
            do i=ip,k
              sj=sj+a(i,j)**2
            enddo
            if (dp.gt.sj) go to 140
              dp=sj
              km=j
              if (in.eq.6) go to 160
140         continue
          enddo
c
c         maximize (sigma)**2 by column interchange.
c
c         suppress column interchanges when in=6.
c
c
c         exchange columns if necessary.
c
          if (km.eq.ip) go to 160
            do i=1,k
              a1=a(i,ip)
              a(i,ip)=a(i,km)
              a(i,km)=a1
            enddo
c
c           record permutation and exchange squares of column lengths.
c
            n1=j1+km
            a1=s(n1)
            n2=j1+ip
            s(n1)=s(n2)
            s(n2)=a1
            n7=j7+km
            n8=j7+ip
            a1=s(n7)
            s(n7)=s(n8)
            s(n8)=a1
160       continue
          if (ip.eq.1) go to 180
            a1=0.d0
            ipm1=ip-1
            do i=1,ipm1
              a1=a1+a(i,ip)**2
            enddo
            if (a1.gt.0.d0) go to 190
180       continue
          if (dp.gt.0.d0) go to 200
c
c           test for rank deficiency.
c
190         continue
            if (sqrt(dp/a1).gt.eps1) go to 200
            if (in.eq.6) go to 200
            ii=ip-1
            if (erm) write (stdout,1140) irank,eps1,ii,ii
            irank=ip-1
            erm=.false.
            go to 260
c
c           (eps1) rank is deficient.
c
200       continue
          sp=sqrt(dp)
c
c         begin front elimination on column ip.
c
c         sp=sqroot(sigma**2).
c
          bp=1.d0/(dp+sp*abs(a(ip,ip)))
c
c         store beta in s(3n+1)-s(3n+l).
c
          if (ip.eq.k) bp=0.d0
          n3=k+2*n+ip
          r(n3)=bp
          up=sign(dble(sp)+abs(a(ip,ip)),dble(a(ip,ip)))
          if (ip.ge.k) go to 250
            ipp1=ip+1
            if (ip.ge.n) go to 240
              do j=ipp1,n
                sj=0.d0
                do i=ipp1,k
                  sj=sj+a(i,j)*a(i,ip)
                enddo
                sj=sj+up*a(ip,j)
                sj=bp*sj
c  
c               sj=yj now
c
                do i=ipp1,k
                  a(i,j)=a(i,j)-a(i,ip)*sj
                enddo
                a(ip,j)=a(ip,j)-sj*up
              enddo
240         continue
            a(ip,ip)=-sign(sp,a(ip,ip))
c
            n4=k+3*n+ip
            r(n4)=up
250       continue
        enddo
        irank=lm
260     continue
        irp1=irank+1
        irm1=irank-1
        if (irank.eq.0.or.irank.eq.n) go to 360
        if (ieq.eq.3) go to 290
c
c       begin back processing for rank deficiency case
c       if irank is less than n.
c
        do j=1,n
          n2=j2+j
          n7=j7+j
          l=min0(j,irank)
c
c         unscale columns for rank deficient matrices when ieq.ne.3.
c
          do i=1,l
            a(i,j)=a(i,j)/s(n7)
          enddo
          s(n7)=1.d0
          s(n2)=1.d0
        enddo
290     continue
        ip=irank
300     continue
        sj=0.d0
        do j=irp1,n
          sj=sj+a(ip,j)**2
        enddo
        sj=sj+a(ip,ip)**2
        aj=sqrt(sj)
        up=sign(dble(aj)+abs(a(ip,ip)),dble(a(ip,ip)))
c
c       ip th element of u vector calculated.
c
        bp=1.d0/(sj+abs(a(ip,ip))*aj)
c
c       bp = 2/length of u squared.
c
        ipm1=ip-1
        if (ipm1.le.0) go to 340
        do i=1,ipm1
          dp=a(i,ip)*up
          do j=irp1,n
            dp=dp+a(i,j)*a(ip,j)
          enddo
          dp=dp/(sj+abs(a(ip,ip))*aj)
c
c         calc. (aj,u), where aj=jth row of a
c
          a(i,ip)=a(i,ip)-up*dp
c
c         modify array a.
c
          do j=irp1,n
            a(i,j)=a(i,j)-a(ip,j)*dp
          enddo
        enddo
340     continue
        a(ip,ip)=-sign(dble(aj),dble(a(ip,ip)))
c
c       calc. modified pivot.
c
c
c       save beta and ip th element of u vector in r array.
c
        n6=k+ip
        n7=k+n+ip
        r(n6)=bp
        r(n7)=up
c
c       test for end of back processing.
c
        if (ip-1) 360,360,350
350     continue
        ip=ip-1
        go to 300
360     continue
        if (in.eq.6) return
        do j=1,k
          r(j)=b(j)
        enddo
        it=0
c
c       set initial x vector to zero.
c
        do j=1,n
          x(j)=0.d0
        enddo
        if (irank.eq.0) go to 690
c
c       apply q to rt. hand side.
c
390     continue
        do ip=1,irank
          n4=k+3*n+ip
          sj=r(n4)*r(ip)
          ipp1=ip+1
          if (ipp1.gt.k) go to 410
            do i=ipp1,k
              sj=sj+a(i,ip)*r(i)
            enddo
410       continue
          n3=k+2*n+ip
          bp=r(n3)
          if (ipp1.gt.k) go to 430
            do i=ipp1,k
              r(i)=r(i)-bp*a(i,ip)*sj
            enddo
430       continue
          r(ip)=r(ip)-bp*r(n4)*sj
        enddo
        do j=1,irank
          s(j)=r(j)
        enddo
        enorm=0.d0
        if (irp1.gt.k) go to 510
          do j=irp1,k
            enorm=enorm+r(j)**2
          enddo
          enorm=sqrt(enorm)
          go to 510
460       continue
          do j=1,n
            sj=0.d0
            n1=j1+j
            ip=s(n1)
            do i=1,k
              sj=sj+r(i)*aa(i,ip)
            enddo
c
c           apply at to rt. hand side.
c           apply scaling.
c
            n7=j2+ip
            n1=k+n+j
            r(n1)=sj*s(n7)
          enddo
          n1=k+n
          s(1)=r(n1+1)/a(1,1)
          if (n.eq.1) go to 510
          do j=2,n
            n1=j-1
            sj=0.d0
            do i=1,n1
              sj=sj+a(i,j)*s(i)
            enddo
            n2=k+j+n
            s(j)=(r(n2)-sj)/a(j,j)
          enddo
c
c         entry to continue iterating.  solves tz = c = 1st irank
c         components of qb .
c
510     continue
        s(irank)=s(irank)/a(irank,irank)
        if (irm1.eq.0) go to 540
          do j=1,irm1
            n1=irank-j
            n2=n1+1
            sj=0.
            do i=n2,irank
              sj=sj+a(n1,i)*s(i)
            enddo
            s(n1)=(s(n1)-sj)/a(n1,n1)
          enddo
c
c         z calculated.  compute x = sz.
c
540     continue
        if (irank.eq.n) go to 590
          do j=irp1,n
            s(j)=0.d0
          enddo
          do i=1,irank
            n7=k+n+i
            sj=r(n7)*s(i)
            do j=irp1,n
              sj=sj+a(i,j)*s(j)
            enddo
            n6=k+i
            do j=irp1,n
              s(j)=s(j)-a(i,j)*r(n6)*sj
            enddo
            s(i)=s(i)-r(n6)*r(n7)*sj
          enddo
c
c         increment for x of minimal length calculated.
c
590     continue
        do i=1,n
          x(i)=x(i)+s(i)
        enddo
        if (in.eq.7) go to 750
c
c         calc. sup norm of increment and residuals
c
          top1=0.d0
          do j=1,n
            n2=j7+j
            top1= max(top1,abs(s(j))*s(n2))
          enddo
          do i=1,k
            sj=0.d0
            do j=1,n
              n1=j1+j
              ip=s(n1)
              n7=j2+ip
              sj=sj+aa(i,ip)*x(j)*s(n7)
            enddo
            r(i)=b(i)-sj
          enddo
        if (itmax.le.0) go to 750
c
c         calc. sup norm of x.
c
          top=0.d0
          do j=1,n
            n2=j7+j
            top= max(top,abs(x(j))*s(n2))
          enddo
c
c         compare relative change in x with tolerance eps .
c
          if (top1-top*eps2) 690,650,650
650   continue
      if (it-itmax) 660,680,680
660   continue
      it=it+1
      if (it.eq.1) go to 670
      if (top1.gt..25*top2) go to 690
670   continue
      top2=top1
      go to (390,460), isw
680   continue
      it=0
690   continue
      sj=0.d0
      do j=1,k
        sj=sj+r(j)**2
      enddo
      enorm=sqrt(sj)
      if (irank.eq.n.and.isw.eq.1) go to 710
      go to 730
710   continue
      enm1=enorm
c
c     save x array.
c
      do j=1,n
        n1=k+j
        r(n1)=x(j)
      enddo
      isw=2
      it=0
      go to 460
c
c     choose best solution
c
730   continue
      if (irank.lt.n) go to 750
      if (enorm.le.enm1) go to 750
      do j=1,n
        n1=k+j
        x(j)=r(n1)
      enddo
      enorm=enm1
c
c     norm of ax - b located in the cell enorm .
c
c
c     rearrange variables.
c
750   continue
      do j=1,n
        n1=j1+j
        s(j)=s(n1)
      enddo
      do j=1,n
        do i=j,n
          ip=s(i)
          if (j.eq.ip) go to 780
        enddo
780     continue
        s(i)=s(j)
        s(j)=j
        sj=x(j)
        x(j)=x(i)
        x(i)=sj
      enddo
c
c     scale variables.
c
      do j=1,n
        n2=j2+j
        x(j)=x(j)*s(n2)
      enddo
      return
c
c     restore a.
c
810   continue
      do j=1,n
        n2=j2+j
        do i=1,k
          a(i,j)=aa(i,j)
        enddo
      enddo
      return
c
c     generate solutions to the homogeneous equation ax = 0.
c
  830 if (irank.eq.n) return
      ns=n-irank
      do i=1,n
        do j=1,ns
          h(i,j)=0.d0
        enddo
      enddo
      do j=1,ns
        n2=irank+j
        h(n2,j)=1.d0
      enddo
      if (irank.eq.0) return
      do j=1,irank
        do i=1,ns
          n7=k+n+j
          sj=r(n7)*h(j,i)
c
c         this part of the code should not be executed. However some
c         compilers generate warning msgs that "irp1" is not defined
c         so it is arbitrarily set here to fool the compilers
c
          irp1 = 1
          if (i .gt. 0) then
            write (stdout,*) ' Error in lsqsl2: search for this msg'
            call abort()
          endif
          do k1=irp1,n
            sj=sj+h(k1,i)*a(j,k1)
          enddo
          n6=k+j
          bp=r(n6)
          dp=bp*r(n7)*sj
          a1=dp
          a2=dp-a1
          h(j,i)=h(j,i)-(a1+2.*a2)
          do k1=irp1,n
            dp=bp*a(j,k1)*sj
            a1=dp
            a2=dp-a1
            h(k1,i)=h(k1,i)-(a1+2.*a2)
          enddo
        enddo
      enddo
c
c     rearrange rows of solution matrix.
c
      do j=1,n
        n1=j1+j
        s(j)=s(n1)
      enddo
      do j=1,n
        do i=j,n
          ip=s(i)
          if (j.eq.ip) go to 900
        enddo
900     continue
        s(i)=s(j)
        s(j)=j
        do k1=1,ns
          a1=h(j,k1)
          h(j,k1)=h(i,k1)
          h(i,k1)=a1
        enddo
      enddo
#endif      
      return
c
 1140 format (/'warning. irank has been set to',i4,'  but(',1pe10.3,
     1 ') rank is',i4,'.  irank is now taken as ',i4)
      end



      subroutine potem (t, s, p, theta)
c
c=======================================================================
c     this subroutine calculates potential temperature as a function
c     of in-situ temperature, salinity, and pressure.
c
c     input [units]:
c       in-situ temperature (t): [degrees centigrade]
c       salinity (s): [per mil]
c       pressure (p): [decibars, approx. as meters of depth]
c     output [units]:
c       potential temperature (theta): [degrees centigrade]
c
c     references:
c        based on Fofonoff and Froese (1958) as shown in ...
c        Fofonoff, N., The Sea: Vol 1, (ed. M. Hill). Interscience,
c          New York, 1962, page 17, table iv.
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c
      b1    = -1.60d-5*p
      b2    = 1.014d-5*p*t
      t2    = t*t
      t3    = t2*t
      b3    = -1.27d-7*p*t2
      b4    = 2.7d-9*p*t3
      b5    = 1.322d-6*p*s
      b6    = -2.62d-8*p*s*t
      s2    = s*s
      p2    = p*p
      b7    = 4.1d-9*p*s2
      b8    = 9.14d-9*p2
      b9    = -2.77d-10*p2*t
      b10   = 9.5d-13*p2*t2
      b11   = -1.557d-13*p2*p
      potmp = b1+b2+b3+b4+b5+b6+b7+b8+b9+b10+b11
      theta = t-potmp
c
      return
      end



      subroutine unesco (t, s, pin, rho)
c
c=======================================================================
c     this subroutine calculates the density of seawater using the
c     standard equation of state recommended by unesco(1981).
c
c     input [units]:
c       in-situ temperature (t): [degrees centigrade]
c       salinity (s): [practical salinity units]
c       pressure (pin): [decibars, approx. as meters of depth]
c     output [units]:
c       density(rho): kilograms per cubic meter
c
c     references:
c        Gill, A., Atmosphere-Ocean Dynamics: International Geophysical
c         Series No. 30. Academic Press, London, 1982, pp 599-600.
c        UNESCO, 10th report of the joint panel on oceanographic tables
c          and standards. UNESCO Tech. Papers in Marine Sci. No. 36,
c          Paris, 1981.
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c
      c1p5 = 1.5d0
c
c  convert from depth [m] (decibars) to bars
      p = pin * 1.0d-1
c
      rw =     9.99842594d2 + 6.793952d-2*t - 9.095290d-3*t**2
     &        + 1.001685d-4*t**3 - 1.120083d-6*t**4 + 6.536332d-9*t**5
c
      rsto =   rw + (8.24493d-1 - 4.0899d-3*t + 7.6438d-5*t**2
     &        - 8.2467d-7*t**3 + 5.3875d-9*t**4) * s
     &       + (-5.72466d-3 + 1.0227d-4*t - 1.6546d-6*t**2) * s**c1p5
     &       + 4.8314d-4 * s**2
c
      xkw =     1.965221d4 + 1.484206d2*t - 2.327105d0*t**2 +
     &         1.360477d-2*t**3 - 5.155288d-5*t**4
c
      xksto =   xkw + (5.46746d1 - 6.03459d-1*t + 1.09987d-2*t**2
     &        - 6.1670d-5*t**3) * s
     &       + (7.944d-2 + 1.6483d-2*t - 5.3009d-4*t**2) * s**c1p5
c
      xkstp =   xksto + (3.239908d0 + 1.43713d-3*t + 1.16092d-4*t**2
     &        - 5.77905d-7*t**3) * p
     &       + (2.2838d-3 - 1.0981d-5*t - 1.6078d-6*t**2) * p * s
     &       + 1.91075d-4 * p * s**c1p5
     &       + (8.50935d-5 - 6.12293d-6*t + 5.2787d-8*t**2) * p**2
     &       + (-9.9348d-7 + 2.0816d-8*t + 9.1697d-10*t**2) * p**2 * s
c
      rho =    rsto / (1.0d0 - p/xkstp)
c
      return
      end
