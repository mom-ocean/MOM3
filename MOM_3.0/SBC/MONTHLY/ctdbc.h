c======================= include file "ctdbc.h" ======================== 
c
# ifdef minimize_sbc_memory
      common/ctdbc/ obc(imt,jmw,ntdbc,2)
# else
      common/ctdbc/ obc(imt,jmt,ntdbc,2)
# endif
