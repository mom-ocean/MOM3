c====================== include file "bbl.h" ===========================
#ifdef bbl_ag
c                                                                       
c     variables for use within the bottom boundary layer
c
      real Kbbl,Kvbbltop
      common /bbl_com/ raydrag, Kvbbltop, cdbotbbl
      common /bbl_com/ entrain_kbbl, entrain_kbbl_constant
!
      common /bbl_com/ zt_bbl(imt,jmw), dzw_bbl(imt,jmw)
      common /bbl_com/ uz2(imt,jmw), rhoz(imt,jmw)
!
      common /bbl_com/ pbot(imt,jmw)
!
      common /bbl_com/ entrain(imt,2:jmw-1)
      common /bbl_com/ px_bbl(imt,1:jmw), oppsignx(imt,1:jmw)
!
      common /bbl_com/ Kbbl(imt,2:jmw)

      common /bbl_com/ adv_fdets(imt,2:jmw,nt)
      common /bbl_com/ adv_fdete(imt,2:jmw,nt)
      common /bbl_com/ adv_fdetw(imt,2:jmw,nt)
!
      common /bbl_com/ py_bbl(imt,1:jmw-1), oppsigny(imt,1:jmw-1)
      common /bbl_com/ grad_p_bbl(imt,1:jmt,2)
      common /bbl_com/ adv_fdetn(imt,1:jmw-1,nt)
# ifdef bbl_mask
      common /bbl_com/ mask_t_bbl(1:imt,1:jmt), mask_u_bbl(1:imt,1:jmt)
# endif
#endif
