c====================== include file "fdift.h" =========================
c
c     finite difference numerics for tracers
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#ifdef partial_cell
c
c     statement functions
c
      dhte(i,k,j)   = min(dht(i+1,k,j),dht(i,k,j))
      dhtn(i,k,j)   = min(dht(i,k,j+1),dht(i,k,j))
      dhwe(i,k,j)   = min(dhwt(i+1,k,j),dhwt(i,k,j))
      dhwn(i,k,j)   = min(dhwt(i,k,j+1),dhwt(i,k,j))
      T_i(i,k,j,n,ip) = t(i+ip,max(1,k-1),j,n,taum1) - dhwe(i,k-1,j)
     &            *(t(i+ip,max(1,k-1),j,n,taum1) - t(i+ip,k,j,n,taum1))
     &               /dhwt(i+ip,k-1,j)
      T_j(i,k,j,n,jp) = t(i,max(1,k-1),j+jp,n,taum1) - dhwn(i,k-1,j)
     &            *(t(i,max(1,k-1),j+jp,n,taum1) - t(i,k,j+jp,n,taum1))
     &               /dhwt(i,k-1,j+jp)
      dz_t2r(i,k,j) = p5/dht(i,k,j)
      dz_tr(i,k,j)  = c1/dht(i,k,j)
      dz_wtr(i,k,j) = c1/dhwt(i,k,j)
      dx_t2r(i,k,j) = cstdxt2r(i,j)/dht(i,k,j)
      dx_tr(i,k,j)  = cstdxtr(i,j)/dht(i,k,j)
      dy_t2r(i,k,j) = cstdyt2r(jrow)/dht(i,k,j)
      dy_tr(i,k,j)  = cstdytr(jrow)/dht(i,k,j)
#else
      T_i(i,k,j,n,ip) = t(i+ip,k,j,n,taum1)
      T_j(i,k,j,n,jp) = t(i,k,j+jp,n,taum1)
      dz_t2r(i,k,j) = dzt2r(k)
      dz_tr(i,k,j)  = dztr(k)
      dz_wtr(i,k,j) = dzwr(k)
      dx_t2r(i,k,j) = cstdxt2r(i,j)
      dx_tr(i,k,j)  = cstdxtr(i,j)
      dy_t2r(i,k,j) = cstdyt2r(jrow)
      dy_tr(i,k,j)  = cstdytr(jrow)
#endif
c
c-----------------------------------------------------------------------
c     advective terms
c-----------------------------------------------------------------------
c
#ifdef linearized_advection
      ADV_Tx(i,k,j) = 0.0
      ADV_Ty(i,k,j) = 0.0
      ADV_Tz(i,k,j) = (adv_fb(i,k-1,j) - adv_fb(i,k,j))*dz_t2r(i,k,j)
#endif
c
#if !defined linearized_advection
      ADV_Tx(i,k,j) = (adv_fe(i,k,j) - adv_fe(i-1,k,j))*dx_t2r(i,k,j)
# if defined fourth_order_tracer_advection || defined quicker
      ADV_Ty(i,k,j)=(adv_f4n(i,k,j,n)-adv_f4n(i,k,j-1,n))*dy_t2r(i,k,j)
# else
      ADV_Ty(i,k,j) = (adv_fn(i,k,j) - adv_fn(i,k,j-1))*dy_t2r(i,k,j)
# endif
      ADV_Tz(i,k,j) = (adv_fb(i,k-1,j) - adv_fb(i,k,j))*dz_t2r(i,k,j)
#endif
c
c-----------------------------------------------------------------------
c     diffusive terms
c-----------------------------------------------------------------------
c
      DIFF_Tx(i,k,j) = (diff_fe(i,  k,j)*tmask(i+1,k,j)
     &                - diff_fe(i-1,k,j)*tmask(i-1,k,j))*dx_tr(i,k,j)
      DIFF_Ty(i,k,j) = (diff_fn(i,k,j  )*tmask(i,k,j+1)
     &                - diff_fn(i,k,j-1)*tmask(i,k,j-1))*dy_tr(i,k,j)
      DIFF_Tz(i,k,j) = ((diff_fb(i,k-1,j) - diff_fb(i,k,j))
#if defined implicitvmix || defined redi_diffusion
     &                 *(c1-aidif)
#endif
#if defined isoneutralmix
     &               + (diff_fbiso(i,k-1,j) - diff_fbiso(i,k,j))
#endif
     &                   )*dz_tr(i,k,j)
