c======================== include file "filtr.h" =======================
c
c     variables used in Fourier filtering
c
      parameter (imtx2=imt*2,ni=imt)
      parameter (imtd2=imt/2,lqmsum=imtd2*(imt-imtd2))
      parameter (lhsum=imt*(imt+1)/2)
      parameter (imtx4=imt*4,imtx8=imt*8,imtimt=imt*imt)
      parameter (imp1x2=(imt+1)*2)
c
c     cossav must remain full precision if most of filter is made half-p
      real cossav
c
      dimension icbase((imt+1)),idbase((imt+1)),ind(imtx8)
      dimension cossav(lqmsum),denmsv(lhsum),cosnpi(imt)
      dimension ftarr(imtimt)
      dimension circle(4)
      common /cfiltr/ ind,  denmsv, idbase, cossav, icbase, cosnpi
      common /cfiltr/ jbase, ibase
      common /cfilt1/  circle
      common /cfilt2/  ftarr

