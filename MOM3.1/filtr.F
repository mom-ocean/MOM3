      subroutine filtr (s, im, mm, n, iss)
#ifdef fourfil
c
c=======================================================================
c                                                                    ===
c  filter fourier analyses the arrays of various                     ===
c         physical quantities, then truncates the series and         ===
c         resynthesizes the filtered quantities where:               ===
c             s  =the string to be filtered                          ===
c             im =the length of s                                    ===
c             mm =1 (cosine series, deriv at bndry pts=0)            ===
c                =2 (  sine series,          bndry pts=0)            ===
c                =3 (full series, cyclic)                            ===
c             n  =number of waves to keep                            ===
c             iss=0 (cant use fourier coefs from previous call)      ===
c             iss>0 (can  use fourier coefs from previous call)      ===
c
c  author: descendent from Mike Cox                                  ===
c                                                                    ===
c=======================================================================
c
c---------------------------------------------------------------------
c     define global data
c---------------------------------------------------------------------
c
# include "param.h"
# include "filtr.h"
# include "ndcon.h"
# include "switch.h"
c
c---------------------------------------------------------------------
c     define local data and dimension argument arrays
c---------------------------------------------------------------------
c
      dimension temp(imtx4)
      dimension indx(imtx8),cof(imtx8)
      dimension cosine(imtx8)
      dimension denom(imtx4)
      dimension s(imt),sprime(imt)
c
c---------------------------------------------------------------------
c  begin executable code
c---------------------------------------------------------------------
c
# ifdef timing
      call tic ('filtering', 'filtr (fourier) ')
# endif
c
      if (im.lt.1 .or. mm.lt.1 .or. mm.gt.3 .or. n.lt.0 .or. iss.lt.0)
     $  then
        write (stdout,99) im, mm, n, iss
        write (stderr,99) im, mm, n, iss
        call abort()
      endif
c
c     calculate some useful constants
c
      if(mm.eq.2 .and. n.eq.0) then
c
        do 400 i=1,im
          s(i) = c0
400     continue
c
        goto 3201
      endif
c
      if (mm .eq. 1) then
        nmax = n - 1
      else
        nmax = n
      endif
c
      nmaxp1 = nmax + 1
      cc1 = p5*float(nmax) + p25
      cc2 = float(nmax) + p5
c
      if (mm .eq. 2) then
        lcy   = 2*(im + 1)
        fnorm = c2/float(im + 1)
      else
        lcy   = 2*im
        fnorm = c2/float(im)
      endif
c
      lh    = lcy/2
      lhm1  = lh - 1
      lqm   = (lh - 1)/2
      l2cy  = 2*lcy
      lcym1 = lcy - 1
      lcyp1 = lcy + 1
      imx2  = im*2
      imx4  = im*4
      imx8  = im*8
c
c     average incoming array
c
      ssum = c0
c
      do 500 i=1,im
        ssum = ssum + s(i)
500   continue
c
c     mm = 1  derivative must be zero at boundaries (cosine)
c     mm = 2  value must be zero at boundaries (sine)
c     mm = 3  cyclic boundary conditions (general fourier series)
c
      fim   = float(im)
      fimr  = c1/fim
      stemp = ssum*fimr
c
      if (n.gt.1 .or. mm.ne.1) goto 601
c
      do 600 i=1,im
        s(i)=stemp
600   continue
c
      go to 3201
601   continue
c
      if (mm .ne. 2) then
c
        do 700 i=1,im
          s(i) = s(i) - stemp
700     continue
c
      endif
c
      if (iss .gt. 0) goto 2501
c
c     assemble appropriate 1-cycle (2*pi) cosine array
c
c     use stored 1/4 cycle to calculate first 1/2 cycle
c
      jbase = icbase(lh)
c
      do 800 i=1,lqm
        cosine(i) = cossav(jbase+i)
800   continue
c
      do 900 i=1,lqm
        cosine(lh-i) = -cossav(jbase+i)
900   continue
c
c     fill in cos(pi/2) if lh is even
c
      if (2*(lqm+1) .eq. lh) cosine(lqm+1) = c0
c
c     fill in cos(pi) in any case
c
      cosine(lh) = -c1
c
c     fill in rest of cycle
c
      do 1000 i=1,lh
        cosine(lh+i) = -cosine(i)
1000  continue
c
c     assemble denominator array
c
      ibase = idbase(lh)
c
      do 1100 i=1,lhm1
        denom(i) = p25*denmsv(ibase+i)
1100  continue
c
      denom(lh) = 0.125
c
      do 1200 i=1,lhm1
        temp(i) = denom(lh-i)
1200  continue
c
      do 1300 i=1,lhm1
        denom(lh+i) = temp(i)
1300  continue
c
      nprint = 0
      denom(lcy) = c0
c
      do 1400 i=lcyp1,imx4
        denom(i) = denom(i-lcy)
1400  continue
c
c     assemble appropriate subscript arrays
c
c     calculate needed indices
c
      if (mm.eq.3) then
        fact1 = 2*nmax
        fact2 = 2*nmaxp1
      else
        fact1 = nmax
        fact2 = nmaxp1
      endif
c
      do 1500 i=1,imx4
        indx(i) = ind(i)*fact1
1500  continue
c
      do 1600 i=1,imx4
        indx(imx4+i) = ind(i)*fact2
1600  continue
c
c     calculate parameters for reducing indices
c
      maxind = imx4*fact2
      ncyc   = (maxind-1)/lcy + 1
      maxndx = lcy
      if (maxndx .ge. maxind) goto 1801
c
      do 1700 npwr=1,ncyc+2
        maxndx = 2*maxndx
        if (maxndx .ge. maxind) goto 1701
1700  continue
c
      write (stdout,999)
      write (stderr,999)
      call abort()
c
1701  continue
c
      do 1800 np=1,npwr
        maxndx = maxndx/2
        do 1790 i=1,imx8
          if (indx(i) .gt. maxndx) indx(i) = indx(i) - maxndx
1790    continue
1800  continue
c
1801  continue
c
c     gather coefficients
c
      do 1900 j=1,imx8
        cof(j) = cosine(indx(j))
1900  continue
c
c     assemble transformation array which will filter s
c
      if(mm.eq.1) then
c
c       cosine transform
c
        ioff1 = lcy
        ioff2 = lcy + imx4
c
        do 2000 j=1,im
          joffs = (j-1)*imt
          do 1990 i=1,im
            ftarr(joffs+i) =
     $         (cof(i-j+ioff1) - cof(i-j+ioff2)) *denom(i-j+ioff1) +
     $         (cof(i + j - 1) - cof(imx4+i+j-1))*denom(i+j-1) - p5
1990      continue
2000    continue
c
        do 2100 j=1,im
          ftarr(j*(imt+1)-imt) = ftarr(j*(imt+1)-imt) + cc1
2100    continue
c
      elseif (mm .eq. 2) then
c
c       sine transform
c
        ioff1 = lcy
        ioff2 = lcy + imx4
c
        do 2200 j=1,im
          joffs = (j-1)*imt
          do 2190 i=1,im
            ftarr(joffs+i) =
     $         (cof(i-j+ioff1) - cof(i-j+ioff2))*denom(i-j+ioff1) -
     $         (cof(i + j)     - cof(imx4+i+j)) *denom(i+j)
2190      continue
2200    continue
c
        do 2300 j=1,im
          ftarr(j*(imt+1)-imt) = ftarr(j*(imt+1)-imt) + cc1
2300   continue
c
      else if(mm.eq.3) then
c
c       general fourier transform
c
        if (2*n .eq. im) then
          genadj = p5
        else
          genadj = c0
        endif
c
        ioff1 = lcy
        ioff2 = lcy + imx4
c
        do 2400 j=1,im
          joffs = (j-1)*imt
          do 2390 i=1,im
            ftarr(joffs+i) = (c2*(cof(i-j+ioff1) - cof(i-j+ioff2)))
     $          *denom(2*i-2*j+ioff1) - p5 - genadj*cosnpi(i)*cosnpi(j)
2390      continue
2400    continue
c
        do 2500 j=1,im
          ftarr(j*(imt+1)-imt) = ftarr(j*(imt+1)-imt) + cc2
2500    continue
c
      endif
c
c     filter s
c
2501  continue
c
      do 2600 i=1,im
        sprime(i) = c0
2600  continue
c
c     note that ftarr(j,i)=ftarr(i,j), so following is legal
c
      do 2700 i=1,im
        ioff = (i-1)*imt
        do 2690 j=1,im
          sprime(j) = sprime(j) + s(i)*ftarr(ioff+j)
2690    continue
2700  continue
c
      do 2800 i=1,im
        sprime(i) = fnorm*sprime(i)
2800  continue
c
      if(mm.eq.2) then
c
        do 2900 i=1,im
          s(i) = sprime(i)
2900    continue
c
        goto 3201
      endif
c
3000  continue
      ssm = c0
c
      do 3100 i=1,im
        ssm = ssm + sprime(i)
3100  continue
c
      ssm = (ssum-ssm)*fimr
c
      do 3200 i=1,im
        s(i) = ssm+sprime(i)
3200  continue
c
3201  continue
c
   99 format (/' error => bad argument(s) in call to filtr'
     $       /' im,mm,n,iss = ',4i10)
  999 format (/' error => can not calculate parameters for reducing',
     $        ' indices in filtr')
# ifdef timing
      call toc ('filtering', 'filtr (fourier) ')
# endif
#endif
      return
      end

