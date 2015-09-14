c====================== add-in "mom_obc.inc" ===========================
c
c-----------------------------------------------------------------------
c         set vertical boundary conditions for momentum and tracers
c         compute internal mode velocities and tracers
c         (first window: start calculation from first row)
c         (last window: continue calculation to last row)   
c-----------------------------------------------------------------------
c
              iobcnmw = num_loads
#ifdef obc_south
          if (mw .eq. 1         ) then
             call cobc       (joff, js_calc,           is, ie)
             call tracer     (joff, js_calc+1, je_calc, is, ie)
             call baroclinic (joff, js_calc,   je_calc, is, ie)
#endif
#ifdef obc_north
# ifdef obc_south
          else if (mw .eq. iobcnmw) then
# else
          if (mw .eq. iobcnmw) then
# endif
             call tracer     (joff, js_calc,  je_calc-1, is, ie)
             call cobc       (joff, je_calc,            is, ie)
             call baroclinic (joff, js_calc,  je_calc-1, is, ie)
#endif
          else
             call tracer     (joff, js_calc, je_calc, is, ie)
             call baroclinic (joff, js_calc, je_calc, is, ie)
          endif
