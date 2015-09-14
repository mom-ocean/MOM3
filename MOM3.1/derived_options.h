#if defined explicit_free_surface
# if !defined explicit_eb && !defined explicit_efb
#  define time_averaged_surface
# endif
#endif

#if defined xlandmix || defined kppvmix || defined sponges 
# define source_term
#endif

#if !defined velocity_horz_laplacian && !defined velocity_horz_biharmonic
# define velocity_horz_laplacian
# define no_velocity_mixing_type
#endif
#if !defined tracer_horz_laplacian && !defined tracer_horz_biharmonic
# define tracer_horz_laplacian
# define no_tracer_mixing_type
#endif

#if defined parallel_1d
# define fourth_order_window
#endif

#if defined bbl_ag
# define fourth_order_window
#endif

#if defined fourth_order_tracer_advection
# define fourth_order_window
#endif

#if defined tracer_advection_fct
# define fourth_order_window
#endif

#if defined quicker
# define fourth_order_window
#endif

#if defined tracer_horz_biharmonic || defined velocity_horz_biharmonic
# define fourth_order_window
#endif

#if defined velocity_horz_mix_smag
# define velocity_horz_mix_var
#endif

#if defined tracer_horz_mix_smag
# define tracer_horz_mix_var
#endif

#if defined biharmonic_rm
# define fourth_order_window
# define isoneutralmix
#endif

#if defined redi_diffusion || defined gent_mcwilliams 
# define isoneutralmix
#endif

#if defined save_density_terms
# define isoneutralmix
#endif

#if defined isotropic_mixed
# define fourth_order_window
#endif

#if !defined tracer_advection_fct && !defined fourth_order_tracer_advection
# define no_fct_and_fourth
#endif

#if defined no_fct_and_fourth && !defined quicker
# define second_order_tracer_advection
#endif

#if defined implicit_free_surface || defined explicit_free_surface
# define free_surf_on
#endif

#if defined explicit_free_surface
# if defined explicit_eta_dpdc || defined explicit_eta_dpdc_all ||defined explicit_eta_dpdc_once 
#  if !defined explicit_eta_dpdc
#   define explicit_eta_dpdc
#  endif
#  if !defined explicit_eta_dpdc_all && !defined explicit_eta_dpdc_once
#   define explicit_eta_dpdc_once
#  endif
# endif
# if !defined simple_sbc && defined explicit_fresh_water_flux
#  if defined river_tracer_explicit || atmos_tracer_explicit || ice_tracer_explicit 
#   define need_sfft_expl
#  endif 
# endif
#endif

#if defined obc_south || defined obc_north
# define obc
#endif

#if defined obc_west || defined obc_east
# define obc
#endif

#if defined obc && !defined source_term
# define source_term
#endif

#if defined kppvmix
# define implicitvmix
#endif

#if defined hl_diffusivity || defined vmhs_diffusivity
# define nonconst_diffusivity
#endif


! window size must be appropriate for options with highest order numerics

#define second_order_window

#if defined fourth_order_window
# undef second_order_window
#endif

#if defined sixth_order_window
# undef second_order_window
# undef fourth_order_window
#endif

! specify lomem or himem on 3D diagnostics. default is himem
# ifdef diagnostics_lomem
#  undef diagnostics_himem
# else
#  define diagnostics_himem
# endif
# ifdef diagnostics_himem
#  undef diagnostics_lomem
# endif

!these are used to determine hardware/OS/compiler (V. Balaji vb@gfdl.gov)

#ifdef __sgi
#ifdef _COMPILER_VERSION
!the MIPSPro compiler defines _COMPILER_VERSION
#define sgi_mipspro
#else
#define sgi_generic
#endif
#endif

#if defined(_CRAY) || defined(sgi_mipspro)
#define SGICRAY
#endif

!compilers that support Cray pointers
#if defined(SGICRAY) || defined(__alpha)|| defined rs6000
#define use_CRI_pointers
#endif

!values of kind: double and long are 8-byte, float and int are 4-byte
#if defined(SGICRAY)
#define DOUBLE_KIND 8
#define FLOAT_KIND 4
#define LONG_KIND 8
#define INT_KIND 4
#else
!these might be different on non-SGICRAY, I believe
#define DOUBLE_KIND 8
#define FLOAT_KIND 4
#define LONG_KIND 8
#define INT_KIND 4
#endif

#ifdef sgi_generic
!this is for the Edinburgh n32/o32 compiler, which won't accept 8-byte ints
!at any price
#define LONG_KIND 4
#endif

!parallel machine types
#if defined(_CRAY) && !defined(_CRAYT3E) && !defined(_CRAYT3D)
#define CRAYPVP
#endif

#if defined(_CRAYT3E) || defined(_CRAYT3D) || defined(sgi_mipspro)
#define SGICRAY_MPP
#endif
