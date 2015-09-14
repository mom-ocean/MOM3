c======================= include file "hmixc.h" ========================
c
c                    horizontal mixing coefficients
c
c     visc_cnu = viscosity coeff for northern face of U cell
c     visc_ceu = viscosity coeff for eastern face of U cell
c     diff_cnt = diffusion coeff for northern face of T cell
c     diff_cet = diffusion coeff for eastern face of T cell
c     
c     am     = constant lateral viscosity coeff for momentum
c     ah     = constant lateral diffusion coeff for tracers
c     am3    = viscosity coeff for metric term on U cell
c     am4    = another viscosity coeff for metric term on U cell
c     ambi   = constant lateral biharmonic viscosity coeff for momentum
c     ahbi   = constant lateral biharmonic diffusion coeff for tracers
c
c     author:   r.c.pacanowski     e-mail  rcp@gfdl.gov
c=======================================================================
c
c-----------------------------------------------------------------------
c     Terms for velocity 
c-----------------------------------------------------------------------
c
      common /diffus0/ am3(jmt), am4(jmt,2)
      common /diffus0/ am, ambi, ah, ahbi
#if defined velocity_horz_mix_const
      common /diffus0/ visc_ceu, visc_cnu
#endif
c
#if defined velocity_horz_mix_var  
c
c     variable horizontal viscosity terms 
c
c     am_varh     = variable viscosity on NE corner of U cell
c     varh_cu     = variable viscosity averaged to U cell
c     varh_met    = extra metric term for U cell
c     
      common /diffus/ am_varh(imt,km,1:jmw-1)
      common /diffus/ varh_met(imt,km,2:jmw-1)
      common /diffus/ varh_cu(imt,km,2:jmw-1)
      common /diffus0/ varh1(jmt,2), varh2(jmt,2)
      common /diffus0/ varh3(jmt,2), varh4(jmt,2)
      common /diffus/ visc_ceu(imt,km,2:jmw-1)
      common /diffus/ visc_cnu(imt,km,1:jmw-1)
c
# if defined velocity_horz_mix_smag
c
c     Smagorinsky coefficients
c
c     k_smag      = overall scaling coeff
c     diff_c_back = background diffusivity for T cell (cm**2/sec)
c     visc_c_back = background viscosity   for U cell (cm**2/sec)
c     prandtl     = ratio (viscosity/diffusivity)
c
      real k_smag
      common /diffus0/  k_smag, diff_c_back, visc_c_back, prandtl
# endif
#endif
c
c-----------------------------------------------------------------------
c     Terms for tracers
c-----------------------------------------------------------------------
c
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
