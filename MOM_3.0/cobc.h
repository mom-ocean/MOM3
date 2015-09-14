c====================== include file "cobc.h" ==========================
#if defined obc_south || defined obc_north
c
c     vad(i,k)   = advective vel. for tracers at open walls
c     var1(i,k,m)= auxiliary field to enhance vectorization
      common /ccobc/ vad(imt,km), var1(imt,km,nt)
# ifdef obc_south
c
c     c1s(i,k,m) = phase vel. at southern boundary for the m_th tracer
      common /ccobc/ c1s(imt,km,nt)
#  if defined orlanski || defined explicit_free_surface
c
c     c1ps(i)    = phase vel. at southern boundary for psi
      common /ccobc/ c1ps(imt)
#  endif
# endif
c
# ifdef obc_north
c
c     c1n(i,k,m) = phase vel. at northern boundary for the m-th tracer
      common /ccobc/ c1n(imt,km,nt)
#  if defined orlanski || defined explicit_free_surface
c
c     c1pn(i)    = phase vel. at northern boundary for psi
      common /ccobc/ c1pn(imt)
#  endif
# endif
#endif
#if defined obc_west || defined obc_east
c
c     uad(k,j)   = advective vel. for tracers at open walls
c     var2(k,j,m)= auxiliary field to enhance vectorization
      common /ccobc/ uad(km,jmw), var2(km,jmt,nt)
c
# ifdef obc_west
c
c     c1w(k,j,m) = phase vel. at western boundary for the m_th tracer
      common /ccobc/ c1w(km,jmw,nt)
#  if defined orlanski || defined explicit_free_surface
c
c     c1pw(j)    = phase vel. at western boundary for psi
      common /ccobc/ c1pw(jmt)
#  endif
# endif
c
# ifdef obc_east
c
c     c1e(k,j,m) = phase vel. at eastern boundary for the m_th tracer
      common /ccobc/ c1e(km,jmw,nt)
#  if defined orlanski || defined explicit_free_surface
c
c     c1pe(j)    = phase vel. at eastern boundary for psi
      common /ccobc/ c1pe(jmt)
#  endif
# endif
#endif
      common /ccobc/ var
#ifdef explicit_free_surface && defined orlanski
     &,              obc_etanull
#endif
