c====================== include file "fdifm.h" =========================
c
c     finite difference numerics for momentum
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#ifdef partial_cell
c
c     statement functions
c
      dhue(i,k,j)   = min(dhu(i+1,k,j),dhu(i,k,j))
      dhun(i,k,j)   = min(dhu(i,k,j+1),dhu(i,k,j))
c
      dx_ur(i,k,j)  = csudxur(i,j)/dhu(i,k,j)
      dy_ur(i,k,j)  = csudyur(jrow)/dhu(i,k,j)
      dz_ur(i,k,j)  = c1/dhu(i,k,j)
      dx_u2r(i,k,j) = csudxu2r(i,j)/dhu(i,k,j)
      dy_u2r(i,k,j) = csudyu2r(jrow)/dhu(i,k,j)
      dz_u2r(i,k,j) = p5/dhu(i,k,j)
#else
      dx_ur(i,k,j)  = csudxur(i,j)
      dy_ur(i,k,j)  = csudyur(jrow)
      dz_ur(i,k,j)  = dztr(k)
      dx_u2r(i,k,j) = csudxu2r(i,j)
      dy_u2r(i,k,j) = csudyu2r(jrow)
      dz_u2r(i,k,j) = dzt2r(k)
#endif
c
c-----------------------------------------------------------------------
c     advective terms
c-----------------------------------------------------------------------
c
#ifdef linearized_advection
      ADV_Ux(i,k,j) = 0.0
      ADV_Uy(i,k,j) = 0.0
      ADV_Uz(i,k,j) = 0.0
      ADV_metric(i,k,j,n)=0.0
#else
      ADV_Ux(i,k,j) = (adv_fe(i,k,j) - adv_fe(i-1,k,j))*dx_u2r(i,k,j)
      ADV_Uy(i,k,j) = (adv_fn(i,k,j) - adv_fn(i,k,j-1))*dy_u2r(i,k,j)
      ADV_Uz(i,k,j) = (adv_fb(i,k-1,j) - adv_fb(i,k,j))*dz_u2r(i,k,j)
      ADV_metric(i,k,j,n)=advmet(jrow,n)*u(i,k,j,1,tau)*u(i,k,j,3-n,tau)
#endif
c
c-----------------------------------------------------------------------
c     viscous terms
c-----------------------------------------------------------------------
c
      DIFF_Ux(i,k,j) = (diff_fe(i,k,j) - diff_fe(i-1,k,j))*dx_ur(i,k,j)
      DIFF_Uy(i,k,j) = (diff_fn(i,k,j) - diff_fn(i,k,j-1))*dy_ur(i,k,j)
      DIFF_Uz(i,k,j) = (diff_fb(i,k-1,j) - diff_fb(i,k,j))*dz_ur(i,k,j)
#ifdef implicitvmix
     &                *(c1-aidif)
#endif
c
#if defined velocity_horz_laplacian && !defined neptune
       DIFF_metric(i,k,j,n) = (am3(jrow)*u(i,k,j,n,taum1)
     &  		      +am4(jrow,n)*dxmetr(i)
     &                *(u(i+1,k,j,3-n,taum1)-u(i-1,k,j,3-n,taum1))
     &                        )
# if defined velocity_horz_mix_const
     &                        *am
# endif
# if defined velocity_horz_mix_var
     &                        *varh_cu(i,k,j)
     &   		      +varh_met(i,k,j)
# endif
# if defined partial_cell
     &                        +pc_sink(i,k,j)*u(i,k,j,n,taum1)
# endif
#endif
c
#if defined velocity_horz_biharmonic && !defined neptune
       DIFF_metric(i,k,j,n) = (am3(jrow)*del2_vel(i,k,j,n)
     &  		      +am4(jrow,n)*dxmetr(i)
     &                *(del2_vel(i+1,k,j,3-n)-del2_vel(i-1,k,j,3-n))
     &                        )
# if defined velocity_horz_mix_const
     &                        *abs(ambi)
# endif
# if defined velocity_horz_mix_var
     &                        *varh_cu(i,k,j)
     &   		      +varh_met(i,k,j)
# endif
# if defined partial_cell
     &                        +pc_sink(i,k,j)*del2_vel(i,k,j,n)
# endif
#endif
c
#if defined velocity_horz_laplacian && defined neptune
      DIFF_metric(i,k,j,n) = (am3(jrow)*(u(i,k,j,n,taum1) -
     &                        unep(i,jrow,n)*umask(i,k,j))
     &                        + am4(jrow,n)*dxmetr(i)*
     &        (u(i+1,k,j,3-n,taum1)-unep(i+1,jrow,3-n)*umask(i+1,k,j)
     &        -u(i-1,k,j,3-n,taum1)+unep(i-1,jrow,3-n)*umask(i-1,k,j))
     &                       )*am
#endif
c
c-----------------------------------------------------------------------
c     coriolis term
c-----------------------------------------------------------------------
c
#ifdef rot_grid
# if defined stream_function
#  if defined damp_inertial_oscillation
      CORIOLIS(i,k,j,n) = cori(i,jrow,n)*u(i,k,j,3-n,taum1)
#  else
      CORIOLIS(i,k,j,n) = cori(i,jrow,n)*u(i,k,j,3-n,tau) 
#  endif
# else
      CORIOLIS(i,k,j,n) = cori(i,jrow,n)*( 
     &             gcor*u(i,k,j,3-n,tau) + (c1-gcor)*u(i,k,j,3-n,taum1))
# endif
#else
# if defined stream_function
#  if defined damp_inertial_oscillation
      CORIOLIS(i,k,j,n) = cori(jrow,n)*u(i,k,j,3-n,taum1)
#  else
      CORIOLIS(i,k,j,n) = cori(jrow,n)*u(i,k,j,3-n,tau) 
#  endif
# else
      CORIOLIS(i,k,j,n) = cori(jrow,n)*( 
     &             gcor*u(i,k,j,3-n,tau) + (c1-gcor)*u(i,k,j,3-n,taum1))
# endif
#endif
c
