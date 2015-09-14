!======================= include file "hmixc.h" ========================
!
!                    horizontal mixing coefficients
!
!     visc_cnu = viscosity coeff for northern face of U cell
!     visc_ceu = viscosity coeff for eastern face of U cell
!     diff_cnt = diffusion coeff for northern face of T cell
!     diff_cet = diffusion coeff for eastern face of T cell
!     
!     am     = constant lateral viscosity coeff for momentum
!     ah     = constant lateral diffusion coeff for tracers
!     am3    = viscosity coeff for metric term on U cell
!     am4    = another viscosity coeff for metric term on U cell
!     ambi   = constant lateral biharmonic viscosity coeff for momentum
!     ahbi   = constant lateral biharmonic diffusion coeff for tracers
!
!     author:   r.c.pacanowski     e-mail  rcp@gfdl.gov
!=======================================================================
!
!-----------------------------------------------------------------------
!     Terms for velocity 
!-----------------------------------------------------------------------
!
      common /diffus0/ am3(jmt), am4(jmt,2)
      common /diffus0/ am, ambi, ah, ahbi
#if defined velocity_horz_mix_const
      common /diffus0/ visc_ceu, visc_cnu
#endif
!
#if defined velocity_horz_mix_var  
!
!     variable horizontal viscosity terms 
!
!     am_varh     = variable viscosity on NE corner of U cell
!     varh_cu     = variable viscosity averaged to U cell
!     varh_met    = extra metric term for U cell
!     
      common /diffus/ am_varh(imt,km,1:jmw-1)
      common /diffus/ varh_met(imt,km,2:jmw-1)
      common /diffus/ varh_cu(imt,km,2:jmw-1)
      common /diffus0/ varh1(jmt,2), varh2(jmt,2)
      common /diffus0/ varh3(jmt,2), varh4(jmt,2)
      common /diffus/ visc_ceu(imt,km,2:jmw-1)
      common /diffus/ visc_cnu(imt,km,1:jmw-1)
!
# if defined velocity_horz_mix_smag
!
!     Smagorinsky coefficients
!
!     k_smag      = overall scaling coeff
!     diff_c_back = background diffusivity for T cell (cm**2/sec)
!     visc_c_back = background viscosity   for U cell (cm**2/sec)
!     prandtl     = ratio (viscosity/diffusivity)
!
      real k_smag
      common /diffus0/  k_smag, diff_c_back, visc_c_back, prandtl
# endif
#endif
!
!-----------------------------------------------------------------------
!     Terms for tracers
!-----------------------------------------------------------------------
!
#if defined tracer_horz_mix_const
      common /diffus0/ diff_cnt, diff_cet
#endif
#if defined tracer_horz_mix_var && !defined isoneutralmix && !defined bryan_lewis_horizontal
      common /diffus/ diff_cet(imt,km,2:jmw-1)
      common /diffus/ diff_cnt(imt,km,1:jmw-1)
#endif
#if defined bryan_lewis_horizontal
      common /diffus0/ Ahh(km)
#endif
