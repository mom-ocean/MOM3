c======================== include file "isopyc.h" ======================
c
c-----------------------------------------------------------------------
c Set default ifdefs
c-----------------------------------------------------------------------
c
#if defined redi_diffusion
# if !defined full_tensor
#  define small_tensor
# endif
#endif
c
#if defined gent_mcwilliams
# if !defined gm_advect
#  define gm_skew
# endif
#endif
c
#if !defined dm_taper
#  define gkw_taper
#endif
c
c
c-----------------------------------------------------------------------
c Arrays and scalars used for all isoneutral mixing schemes.
c
c Arrays
c     drho     = d(rho)/d(tracer), where tracer is temp and salt.
c                This is a quadratic polynomial defined in dens.h.
c     ddxt     = tracer gradients at eastern face of T cells
c     ddyt     = tracer gradients at northern face of T cells
c     ddyt     = tracer gradients at bottom face of T cells
c     fzisop   = vertical structure function for the mixing coefficients.
c     redi_gm  = array to help set the elements of the mixing tensor. 
c     agm      = array of gm diffusivity (cm**2/sec).  Must be 
c                defined as (imt,jmt) unless wish to compute
c                nonconstant diffusivity on every time step. 
c                This field is saved in restart when running with 
c                a nonconstant diffusivity scheme. 
c
c Scalars
c     slmxr  = reciprocal of maximum allowable slope of isoneutrals for
c              small angle Redi, for gent_mcwilliams, and for
c              biharmonic_rm.
c     del_dm  = parameter for the dm_taper steep slope scheme
c     s_dm    = parameter for the dm_taper steep slope scheme
c
c     ahisop  = isoneutral diffusivity (cm**2/sec) 
c     athkdf  = GM diffusivity (cm**2/sec)
c     abihrm  = Roberts and Marshall biharmonic diffusivity (cm**4/sec) 
c     ahsteep = horizontal diffusivity used in the regions of
c               steep slopes (cm**2/sec) 
c-----------------------------------------------------------------------
c
#ifdef isoneutralmix
      common /cisop/ drho(imt,km,jmw,2)
# ifdef full_tensor
      common /cisop/ ddxt(imt,km,jmw,2)
# else
      common /cisop/ ddxt(imt,km,2:jmw-1,2)
# endif
      common /cisop/ ddyt(imt,km,1:jmw-1,2)
      common /cisop/ ddzt(imt,0:km,jmw,2)
c
      common /cisopx/ fzisop(km), slmxr
      common /cisopx/ redi_gm(6) 
      common /cisopx/ del_dm, s_dmr 
      common /cisopx/ ahisop
      common /cisopx/ athkdf 
      common /cisopx/ agm(imt,jmt) 
      common /cisopx/ ahsteep
      common /cisopx/ abihrm
c
# if defined redi_diffusion || defined gent_mcwilliams
c
c-----------------------------------------------------------------------
c Arrays for either Redi diffusion or Gent McWilliams stirring
c     Ai_ez  = slope*diffusion coefficient on eastern face of T cell
c     Ai_nz  = slope*diffusion coefficient on northern face of T cell
c     Ai_bx  = slope*diffusion coefficient on bottom face of T cell
c     Ai_by  = slope*diffusion coefficient on bottom face of T cell
c-----------------------------------------------------------------------
c
      common /cisop/ Ai_ez(imt,km,2:jmw-1,0:1,0:1)
      common /cisop/ Ai_nz(imt,km,1:jmw-1,0:1,0:1)
      common /cisop/ Ai_bx(imt,km,2:jmw-1,0:1,0:1)
      common /cisop/ Ai_by(imt,km,2:jmw-1,0:1,0:1)
# endif
c
c
# ifdef redi_diffusion
c
c-----------------------------------------------------------------------
c Arrays and scalars for just Redi diffusion
c
c Arrays
c     K11,K22,K33 = diagonal components to the Redi diffusion tensor.
c
c Scalars
c     delta_iso1  = grid factor which determines when to scale the 
c                   full tensor's diffusivities. 
c     s_minus,s_plus = slopes which determine the range over which 
c                      need to scale the full tensor's diffusivities.
c-----------------------------------------------------------------------
c
      real K11, K22, K33
      common /cisop/ K11(imt,km,2:jmw-1)
      common /cisop/ K22(imt,km,1:jmw-1)
      common /cisop/ K33(imt,km,2:jmw-1)
      common /cisopx/ delta_iso1, s_minus, s_plus 
# endif
c
# ifdef gent_mcwilliams
c
c-----------------------------------------------------------------------
c Arrays for GM.  Note, the advection velocities are still needed 
c with the default gm_skew approach, since generally wish to map 
c these velocities in snapshots.
c
c     adv_vetiso = zonal isopycnal mixing velocity computed at the 
c                  center of the eastern face of the "t" cells
c     adv_vntiso = meridional isopycnal mixing velocity computed at
c                  the center of the northern face of the "t" cells
c                  (Note: this includes the cosine as in "adv_vnt")
c     adv_vbtiso = vertical isopycnal mixing velocity computed at the
c                  center of the top face of the "t" cells
c     adv_fbiso  = "adv_vbtiso" * (tracer) evaluated at the center of
c                  the bottom face of the "t" cells
c-----------------------------------------------------------------------
c
      common /cisop/ adv_vetiso(imt,km,2:jmw-1)
      common /cisop/ adv_vntiso(imt,km,1:jmw-1)
      common /cisop/ adv_vbtiso(imt,0:km,2:jmw-1)   
      common /cisop/ adv_fbiso(imt,0:km,2:jmw-1)   
# endif
c
# ifdef biharmonic_rm
c
c-----------------------------------------------------------------------
c Arrays for the Roberts and Marshall biharmonic operator.
c
c     Bi_ez  = slope*diffusivity for biharmonic scheme
c              on east face of T cell
c     Bi_nz  = slope*diffusivity for biharmonic scheme
c              on north  face of T cell
c     Bi_bx  = slope*diffusivity for biharmonic scheme
c              on bottom face of T cell
c     Bi_by  = slope*diffusivity for biharmonic scheme 
c              on bottom face of T cell
c-----------------------------------------------------------------------
c
      common /cisop/ Bi_ez(imt,km,2:jmw-1,0:1,0:1)
      common /cisop/ Bi_nz(imt,km,1:jmw-1,0:1,0:1)
      common /cisop/ Bi_bx(imt,km,2:jmw-1,0:1,0:1)
      common /cisop/ Bi_by(imt,km,2:jmw-1,0:1,0:1)
# endif
c
# ifdef save_density_terms
c
c     advrho   = contribution to local potential density tendency
c                due to advection.
c
      common /cisop/ advrho(imt,km,2:jmw-1)
# endif
c
# ifdef partial_cell
c
c     twice the thickness of quarter cells within T-cells
c     0 is for upper quarter cell and 1 is for lower quarter cell
c
      common /cisop/ delqc(imt,km,jmw,0:1)
# endif
c
# ifdef nonconst_diffusivity
c
c     parameters for nonconstant diffusivity schemes
c
c     depthbot         = Bottom depth for performing vertical average.
c     depthtop         = Top depth for performing vertical average.  
c     kncbot           = Depth k-level for bottom of vertical average.
c     knctop           = Depth k-level for top of vertical average.  
c     diffmin          = Minimum diffusivity (cm2/sec).
c     diffmax          = Maximum diffusivity (cm2/sec).
c     vmhs_alpha       = Dimensionless parameter for setting the 
c                        overall scale of the VMHS diffusivity.
c     vmhs_rate2_limit = minimum squared growth rate (units sec^-2)
c                        for defining a baroclinic zone in 
c                        vmhs_diffusivity.
c     ijvmhs_maxlen    = (1/2) the maximum number of grid points 
c                        to be accepted for the width of the
c                        baroclinic zone in vmhs_diffusivity
c     numfltrgrth      = number of FIR filter passes applied to the 
c                        growth rate.  default=0. 
c     growth2max       = maximum squared growth rate allowed for use in 
c                        computing diffusivities. 
c                        (1/day)**2 is reasonable.
c
      common /diffus1/ kncbot,knctop
      common /diffus1/ depthtop,depthbot,diffmin,diffmax
      common /diffus1/ vmhs_alpha,vmhs_rate2_limit,ijvmhs_maxlen
      common /diffus1/ numfltrgrth
      common /diffus1/ growth2max
c
c     variables for nonconstant diffusivity schemes
c
c     fznoncst   = vertical structure function for determining
c                  what vertical levels to be used in making 
c                  vertical average to compute time scale.
c                  fznoncst(k) ne 0 if include k-level.
c                  fznoncst(k) =  0 if do not include k-level.
c     betaeff    = effective beta = H |grad(f/H)|, based on FIR
c                  filtered bottom topog, defined in topog.F.
c     growth2a   = Old squared growth rate (sec^-2)
c                  This field is saved in the restart file.
c     growth2b   = New squared growth rate (sec^-2)
c     hldiff     = Held and Larichev diffusivity (cm**2/sec)
c     hllength   = Held and Larichev length scale (cm)
c     vmhsdiff   = Visbeck etal diffusivity (cm**2/sec)
c     vmhslength = VMHS length scale (cm)
c
      common /diffus1/ fznoncst(km),betaeff(imt,jmt)
      common /diffus1/ growth2a(imt,jmt),growth2b(imt,jmt)
      common /diffus1/ hllength(imt,jmw),hldiff(imt,jmw)
      common /diffus1/ vmhslength(imt,jmw),vmhsdiff(imt,jmw)
# endif
c
c
#endif
