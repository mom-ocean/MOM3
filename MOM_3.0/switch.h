!====================== include file "switch.h" ========================
!
!     all time dependent decisions are made by time manager "tmngr.F"
!     and communicated elsewhere to the model via logical switches.
!
!     inputs: (defaulted in "blkdta.F", optionally reset via namelist)
!
!
!     runlen  = integration period (see rununits). note "runlen" should
!               be an integral number of density time steps. if not,
!               then "runlen" is automatically adjusted to insure this.
!               fractional days are supported but not fractional months
!               or years.
!     rununits= units of "runlen". may be "days", "months", or "years".
!               tmngr will convert "runlen" which is in "rununits"
!               to "rundays" in units of days.
!
!     segtim  = the integration time "runlen" is broken into a number of
!               segments each of length "segtim" days. updated surface
!               boundary conditions are applied to MOM every "segtim" 
!               days. this is useful when coupling to atmospheric models
!               in which case both models exchange surface boundary
!               conditions every "segtim" days where "segtim"
!               is 1/(coupling frequency). without an atmospheric model,
!               when getting surface boundary conditions from data, 
!               "segtim" is set to the time step (in days) by mom.F. in
!               either case, "runlen" (in days) should be an integral
!               number of "segtim".     
!
!     nmix    = number of time steps between mixing timesteps. used
!               to damp timestep splitting due to centered leapfrog.
!
!     init    = (true,false)  indicates that this run is a
!               (start from initial conditions, restart)
!
!     restrt  = (true,false) = (do,don`t) write a restart at the end
!               of the run
!
!     eb      = (true,false) configures for the use of a 
!               (euler backward,forward) type mixing timestep
!
!-----------------------------------------------------------------------
!     inputs to tmngr.F: diagnostic intervals
!-----------------------------------------------------------------------
!
!     note: switches are used to control the interval between doing
!           diagnostics. units for all switches are in days. 
!           setting a switch < 0.0 disables its diagnostic 
!           setting a switch = 0.0 causes the diagnostic to be 
!           done every time step, and setting a switch > 0.0 causes the
!           diagnostic to be done repeatedly at the specified interval
!           from a specified starting date. 
!
!     cmixint = number of days between writing estimated mixing coeffs
!               on faces of T cells and U cells
!
!     crossint = number of days between writing diapycnal and isopycnal
!               components of flow
!
!     fctint = number of days between writing difference between 
!              FCT and leapfrog advection
!
!     pressint = number of days between writing pressure field 
!
!     exconvint = number of days between writing temperature rate of
!                 change due to explicit convection
!
!     glenint =  number of days between global energetics integrals.
!
!     trmbint =  number of days between momentum and tracer term
!                balances (global and regional).
!
!     itrmb   = (true,false) = (do,don`t) write regional mask info for
!               the term balance diagnostic. Typically set true
!               at the beginning of a run; otherwise false since it is
!               not necessary to keep writing a time independent field
!               particularly when it may be a significant part of the
!               time dependent part of the diagnostic.
!
!     gyreint =  number of days between calculation of tracer northward
!                transport.
!     igyre   = (true,false) = (do,don`t) write regional mask info for
!               the gyre diagnostic. Typically set true
!               at the beginning of a run; otherwise false since it is
!               not necessary to keep writing a time independent field
!               particularly when it may be a significant part of the
!               time dependent part of the diagnostic.
!
!     vmsfint =  number of days between calculation of vertical and
!                meridional stream function.
!
!     tyzint  =  number of days between calculation of zonally averaged
!                tracer components.
!
!     rhoint  = number of days between writing instantaneous terms 
!               in the locally referenced potential density equation.
!               Relevant only when 
!               "save_density_terms"
!               is enabled. see "iounit.h" for more details.
!               all data is at time level "tau" so it is one time step
!               before the date shown in "stamp".
!
!     diffint   = number of days between computing nonconstant 
!               diffusivity.
!
!     diffsnapint = number of days between saving nonconstant 
!               diffusivity information to netcdf file.
!
!     prxzint =  number of days between printouts of x-z data.
!
!     extint  =  number of days between printouts of external mode.
!
!     dspint  =  number of days between surface pressure calculation.
!                Note: only when "diagnostic_surface_height" is enabled.
!     dspper  = averaging period for "diagnostic_surface_height"
!
!     tavgint = number of days between regional tracer averages (under
!               horizontal regions).
!
!     itavg   = (true,false) = (do,don`t) write regional mask info for
!               the tracer average diagnostic. Typically set true
!               at the beginning of a run; otherwise false since it is
!               not necessary to keep writing a time independent field
!               particularly when it may be a significant part of the
!               time dependent part of the diagnostic.
!
!     tmbint  = number of days over which tracer equation in averaged
!               in depth and longitude to determine the meridional 
!               balance among storage, divergence, dissipation and 
!               forcing.
!     tmbper  = averaging period for "meridional_tracer_balance"
!
!     itmb    = (true,false) = (do,don`t) write "msktmb" for tracer
!               the meridional balance diagnostic. Typically set true
!               at the beginning of a run; otherwise false since it is
!               not necessary to keep writing a time independent field
!               particularly when it may be a significant part of the
!               time dependent part of the diagnostic.
!
!
!     tsiint  = number of days between printing of time step integrals.
!
!     stabint = number of days between sampling for various stability
!               criteria.
!
!     snapint = number of days between writing instantaneous samples
!               of data from the MOM grid. note: only when "snapshots"
!               is enabled. see "iounit.h" for more details.
!               all data is at time level "tau" so it is one time step
!               before the date shown in "stamp".
!
!     timavgint= interval (days) for writing time mean data
!               (only when "time_averages" is enabled). 
!               if "timavgint" is not an integral number of
!               density time steps,"timavgint" is automatically adjusted
!               to insure this. if the number of days to integrate is
!               not an integral number of "timavgint" then the last
!               averaging period will be less than "timavgint" days.this
!               may lead to one more averaging period than expected. 
!               see "iounit.h" for more details.
!     timavgper= averaging period for "time_averages"
!
!     xbtint  = averaging period (days) for writing XBT data (only when
!               "xbts" is enabled). if "xbtint" is not an integral
!               number of density time steps, "xbtint" is automatically
!               adjusted to insure this. if the number of days to 
!               integrate is not an integral number of "xbtint" then the
!               last averaging period will be less than "xbtint" days.
!               this may lead to one more averaging period than  
!               expected. see "iounit.h" for more details.
!     xbtper  = averaging period for "xbts"
!
!     zmbcint = number of days between calculation of zonal mean
!               surface boundary conditions (and related  quantities)
!
!     trajint = number of days between writing particle trajectories.
!               particle trajectories are always written on the last
!               time step of the run to allow for restarting.
!
!-----------------------------------------------------------------------
!     outputs from tmngr.F: logical switches
!-----------------------------------------------------------------------
!
!     rundays = integration time in days (from "runlen")
!
!     the following are logical counterparts to the above switches and
!     are set every time step. logical switches control all decisions
!     about when to do things in MOM.
!
!     cmixts  = (false,true) = (don`t, do) do write estimated mixing
!               coefficients on this time step.
!               based on "cmixint".
!
!     crossts  = (false,true) = (don`t, do) write diapycnal and 
!               isopycnal components of flow on this time step.
!               based on "crossint".
!
!     fctts    = (false,true) = (don`t, do) write difference between 
!               FCT and leapfrog advection on this time step.
!               based on "fctint".
!
!     pressts  = (false,true) = (don`t, do) write pressure on this time
!               step. based on "pressint".
!               
!
!     exconvts  = (false,true) = (don`t, do) do write temperature change
!               due to explicit convection on this time step.
!               based on "exconvint".
!
!     glents  = (false,true) = (don`t, do) do calculation of global
!               energy integrals on this time step. based on "glenint".
!
!     trmbts  = (false,true) = (don`t, do) do calculation of momentum &
!               tracer term balance on this timestep. based on "trmbint"
!
!     gyrets  = (false,true) = (don`t, do) do calculation of tracer
!               northward transport on this timestep. based on "gyreint"
!
!     vmsfts  = (false,true) = (don`t, do) do calculation of vertical
!               and meridional stream function on this time step.
!               based on "vmsfint"
!
!     tyzts   = (false,true) = (don`t, do) do calculation of zonally
!               averaged tracer components on this time step.
!               based on "tyzint"
!
!     prxzts  = (false,true) = (don`t, do) do printouts of x-z data
!               on this time step. based on "prxzint"
!
!     extts  = (false,true) = (don`t, do) do printout of external mode
!               on this time step. based on "extint"
!
!     dspts  = (false,true) = (don`t, do) do calculation of diagnosti! 
!              surface pressure on this time step. based on "dspint"
!               
!     rhots  = (false,true) =(don`t, do) do calculation of diagnostic
!              for locally referenced potential density terms.
!              Based on "rhoint".
!
!     diffts    = (false,true) =(don`t, do) do calculation of new 
!              nonconstant diffusivity.
!              Based on "diffint".
!
!     diffsnapts = (false,true) =(don`t, do) save fields related 
!              to nonconstant diffusivity calculation.
!              Based on "diffsnapint".
!
!     stabts  = (false,true) = (don`t, do) test for stability on this
!               time step. based on "stabint"
!
!     tavgts  = (false,true) = (don`t do) do tracer averages on this
!               time step. based on "tavgint"
!
!     tmbts   = (false,true) = (don`t, do) write out tracer meridional .
!               balance on this time step. based on "tmbint"
!
!     tsits   = (false,true) = (don`t, do) print time step integrals
!               on this time step. based on "tsiint"
!
!     zmbcts  = (false,true) = (don`t, do) print zonal mean boundary
!               conditions on this time step.  based on "zmbcint"
!
!     trajts  = (false,true) = (don`t, do) print particle trajectories
!               on this time step. based on "trajint"
!
!     snapts  = (false,true) = (don`t, do) save a snapshot of the data
!               on this time step. based on "snapint"
!
!     timats  = (false,true) = (don`t, do) write time mean data
!               on this time step. based on "timavgint"
!                                          
!     xbtts   = (false,true) = (don`t, do) write averaged XBT data on
!               this time step based on "xbtint"
!
!
!     leapfrog= (false,true) on a (mixing, normal leapfrog) time step
!                based on "nmix"
!                                          
!     euler1  = true on the 1st pass of an euler backward time step
!               otherwise false. (applies when "eb" = true)
!     euler2  = true on the 2nd pass of an euler backward time step
!               otherwise false. (applies when "eb" = true)
!     forward = true on a forward time step. otherwise false
!		(applies when "eb" = false)		  
!                                          
!
!     the following logical switches are based on the model time step
!     and not on a specified interval in days. 
!    
!     first   = (true,false) =  when it`s (the first, not the first)
!                               time step of a run
!     eots    = end of a time step. always true except for first
!               pass of an euler backward time step
!     eorun   = last time step of a run. always false except during the
!               last time step of the run.
!
!     eoday   = true when within 1/2 time step of the end of a day
!               else ... false
!     midmon  = true when within 1/2 time step of the middle of a month
!               else ... false
!     eomon   = true when within 1/2 time step of the end of a month
!               else ... false
!     eoyear  = true when within 1/2 time step of the end of a year
!               else ... false
!     osegs   = true on the 1st time step of an ocean segment in mom.F
!               otherwise false.
!     osege  =  true on the last time step of an ocean segment in mom.F
!               otherwise false.
!
!     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
!
!
      logical eb, leapfrog, euler1, euler2, forward, eots
      logical init, first, restrt
      logical itavg, itmb, itrmb, igyre
      character *8 rununits
!
      common /switcc/ rununits
      common /switcr/ runlen, rundays
      common /switci/ nmix
      common /switcl/ eb, leapfrog, euler1, euler2, forward, eots
      common /switcl/ init, first, restrt
      common /switcl/ itavg, itmb, itrmb, igyre
!
      logical          eoday, eomon, midmon, eoyear, eorun
      common /switcl/  eoday, eomon, midmon, eoyear, eorun
!
!-----------------------------------------------------------------------
!
!     S W I T C H E S    B A S E D    O N    A N    I N T E R V A L
!
!     each interval switch needs two variables in common. The 
!     following naming convention is used.
!
!         1) an interval (real) for diagnosti! output (e.g,.  glenint)
!         2) a switch (logical) for the interval (e.g.,  glents )
!
!-----------------------------------------------------------------------
!
      logical           tavgts,    glents,    trmbts,   prxzts
      common /switcr/   tavgint,   glenint,   trmbint,  prxzint
      common /switcl/   tavgts,    glents,    trmbts,   prxzts
!
      logical          extts,   trajts,   exconvts,   cmixts
      common /switcr/  extint,  trajint,  exconvint,  cmixint
      common /switcl/  extts,   trajts,   exconvts,   cmixts
!
      logical          vmsfts,   gyrets,    tsits,   snapts,   fctts
      common /switcr/  vmsfint,  gyreint,   tsiint,  snapint,  fctint
      common /switcl/  vmsfts,   gyrets,    tsits,   snapts,   fctts
!
      logical          tyzts
      common /switcr/  tyzint
      common /switcl/  tyzts
!
      logical          stabts,   zmbcts,   crossts,   pressts
      common /switcr/  stabint,  zmbcint,  crossint,  pressint
      common /switcl/  stabts,   zmbcts,   crossts,   pressts
!
      logical          osegs,  osege
      common /switcr/          segtim
      common /switcl/  osegs,  osege
!
      logical          rhots
      common /switcr/  rhoint
      common /switcl/  rhots
!
      logical          diffts,  diffsnapts
      common /switcr/  diffint, diffsnapint
      common /switcl/  diffts,  diffsnapts
!   
!-----------------------------------------------------------------------
!
!     S W I T C H E S    B A S E D    O N    A N    I N T E R V A L
!
!              A N D   A V E R A G I N G   P E R I O D
!
!     each averaging period switch needs four variables in common. The 
!     following naming convention is used.
!
!         1) an interval (real) for diagnostic output    (e.g. xbtint  )
!         2) a switch (logical) for the interval         (e.g. xbtts   )
!         3) an averaging period (real)                  (e.g. xbtper  )
!         4) a switch (logical) for accumulating         (e.g. xbtperts)
!
!
!     The user must specify the interval [e.g., xbtint] for diagnostic
!     output in units of days and the averaging period [e.g., xbtper] 
!     in units of days. The averaging period may be less than or equal
!     to the interval. For example, if the interval is 30.0 days and the
!     averaging period is 5.0 days, results will be averaged over all
!     time steps within days 26, 27, 28, 29, and 30. An averaging period
!     of 0.0 days averages over the last time step of the interval (as
!     does xbtper = dt), and an averaging period less than zero turns
!     the switches off for all time steps.  
!
!     The logical switch for writing output at the specified interval 
!     [e.g., xbtts] is set to true on the last time step of the 
!     averaging period. The logical switch for accumulating results
!     [e.g., xbtperts] is true for all time steps within the averaging
!     period, otherwise it is false. 
!
!-----------------------------------------------------------------------
!
      logical         xbtts,     dspts,    tmbts,    timavgts
      logical         xbtperts,  dspperts, tmbperts, timavgperts
      common /switcr/ xbtint,    dspint,   tmbint,   timavgint
      common /switcr/ xbtper,    dspper,   tmbper,   timavgper
      common /switcl/ xbtts,     dspts,    tmbts,    timavgts
      common /switcl/ xbtperts,  dspperts, tmbperts, timavgperts
