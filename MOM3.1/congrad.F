#if !defined explicit_free_surface
      subroutine congr  (npt, variable, bc_symm
     &,                  guess, dpsi, forc, res
     &,                  cf, js, je
     &,                  max_iterations, iterations, epsilon
     &,                  imask, iperm, jperm, iofs, nisle, nippts
     &,                  converged
     &,                  estimated_error
     &                  )
c
c=======================================================================
c
c                            C O N G R A D
c
c      solve:
c
c             A * dpsi = forc
c
c      for "dpsi" with dirichlet boundary conditions (dpsi=const on
c      each component of the boundary) by a preconditioned conjugate
c      gradient algorithm.
c
c      inputs:
c              npt   = 5 or 9 (active coefficients)
c              variable = character string identifying solution variable
c              bc_symm = equatorial symmetry type (used only when the
c                        symmetry option is on. otherwise ignore it)
c              guess = initial approximation to solution
c              A     = linear operator (assumed symmetric)
c                      typically A is  grad{(1/h)*grad(dpsi)} -
c                      2dt*acor*{grad(f/h) x grad(dpsi)}
c                      using 5 or 9 pt discretizations
c              cf    = imt x jmt x 3 x 3 array of coefficients of A
c              forc  = the sum of all terms evaluated at times tau
c                      or tau-1
c              epsilon = convergence criterion
c              max_iterations = maximum number of iterations
c              imask = shows which land masses have perimeter equations
c              iperm = i coordinate of island perimeter points
c              jperm = j coordinate of island perimeter points
c              iofs  = offset in iperm, jperm for start of perimeter
c                      of land_mass(isle)
c              nisle = actual number of land_masses
c              nippts = number of perimeter ocean points for a land_mass
c      output:
c              dpsi   = answer
c              iterations = actual number of iterations performed
c              converged = logical value 
c              estimated_error = estimated maximum error in solution
c                          based on step sizes and convergence rate
c
c
c      author:      Charles Goldberg        e-mail: chg@gfdl.gov
c     (mpp items)   v. balaji               e-mail=> vb@gfdl.gov
c
c      based on the preconditioned conjugate gradient algorithm given
c      in:
c
c     A Reformulation and Implementation of the Bryan-Cox-Semtner
c     Ocean Model on the Connection Machine
c     J.K. Dukowicz, R.D. Smith, and R.C. Malone
c     Journal of Atmospheric and Oceanic Technology
c     Vol 10. No. 2 April 1993
c
c=======================================================================
c
c      more specifically, the equations to be solved are
c
c             sum (A(ij,i'j') * dpsi(i'j')) = forc(ij)
c
c      where the subscripts ij and i'j' range over all "free ocean"
c      T cells ij=(i,j) that are not adjacent to land T cells,
c      and one ij=isle for each boundary component of the ocean.
c
c      with this choice of variables, in the absence of coriolis terms
c      (acor=0), the operator A is symmetric, i.e.,
c
c             A(ij,i'j') = A(i'j',ij)
c
c
c      the algorithm (essentially executable in Fortran 90) is...
c
c      subroutine congrad (A, guess, forc, dpsi, iterations)
c
c      use matrix_module
c
c      intent (in)     :: A, guess, forc
c      intent (out)    :: dpsi, iterations
c
c      type(dpsi_type) :: guess, dpsi, Zres, s
c      type(res_type)  :: res, As, forc
c      type(operator)  :: A
c      type(inv_op)    :: Z
c      dimension (0:max_iterations) :: dpsi, res, s, As, beta, alpha
c
c      dpsi(0) = guess
c      res(0)  = forc - A * dpsi(0)
c      beta(0) = 1
c      s(0)    = zerovector()
c      do k = 1 to max_iterations
c        Zres(k-1) = Z * res(k-1)
c        beta(k)   = res(k-1) * Zres(k-1)
c        s(k)      = Zres(k-1) + (beta(k)/beta(k-1)) * s(k-1)
c        As(k)     = A * s(k)
c        alpha(k)  = beta(k) / (s(k) * As(k))
c        dpsi(k)   = dpsi(k-1) + alpha(k) * s(k)
c        res(k)    = res(k-1) - alpha(k) * As(k)
c        estimated_error = err_est(k, alpha(k), s(k))
c        if (estimated_error) < epsilon) exit
c      enddo
c      if (k > max_iterations) then
c        write (stdout,*) 'did not converge in ',k,' iterations'
c        stop '=>congrad'
c      endif
c
c      iterations = k
c      dpsi = dpsi(k)
c
c      end
c
c      where...
c
c      the "vector" and "operator" types used in conjugate gradient
c      are mapped to ordinary 2-dimensional fortran arrays as follows:
c
c      type(dpsi_type) :: guess, dpsi, Zres, s
c          if ij=(i,j) is a mid-ocean point, map dpsi(ij)-->dpsi(i,j)
c          if ij=isle is an ocean boundary subscript, replicate the
c          value dpsi(isle) in dpsi(i,j) for each (i,j) in the ocean
c          perimeter of land_mass(isle).  the arrays iperm(isle) and
c          jperm(isle), along with iofs(isle) locate these ocean
c          perimeter T cells.
c      type(res_type)  :: res, As, forc
c          if ij=(i,j) is a mid-ocean point, res(ij)-->res(i,j)
c          if ij=isle is an ocean boundary subscript, the value of
c          res(isle) = sum (res(i,j))
c          where the sum is taken over all (i,j) in the ocean perimeter
c          of land_mass(isle).  sometimes, the computed values
c          res(i,j) represent contributions of T cell (i,j) to the
c          component res(isle), and sometimes the values are balanced
c          so that res(i,j)=res(isle)/nippts(isle).  note that, even
c          when balanced, the relation between type(res_type) variables
c          res(isle) and res(i,j) differs from that of type(dpsi_type)
c          variables dpsi(isle) and dpsi(i,j) on T cells in the ocean
c          perimeter.
c      type(operator)  :: A
c          the nearly diagonal quality of the operators used
c          permits a representation as a small collection of
c          2-dimensional arrays.
c          the diagonal, A(ij,ij), is stored in an array cfdiag(i,j)
c          as follows:
c          if ij=(i,j) is a mid-ocean point, A(ij,ij) = cfdiag(i,j)
c          if ij=isle is an ocean boundary subscript,
c          A(isle, isle) = sum (cfdiag(i,j))
c          where the sum is taken over all (i,j) in the ocean perimeter
c          of land_mass(isle).  each cfdiag(i,j) is the contribution
c          of T cell (i,j) to the island variable diagonal coefficient.
c          the off-diagonal terms A(ij,i`j`) are stored in 4 arrays
c          cfn, cfs, cfe, and cfw if A is a 5-point operator, and in
c          these and 4 additional arrays, cfne, cfnw, cfse, cfsw, if
c          A is a 9-point operator.  For example, if i`=i and j`=j+1,
c          then A(ij,i`j`) is stored in cfn(i,j).
c          if ij=(i,j) is a mid-ocean point and i`j`=isle` is and ocean
c          perimeter subscript, with i`=i and j`=j+1, then
c          cfn(i,j)=A(ij,isle`) is the coefficient of the island
c          variable dpsi(isle`) in the equation for mid-ocean point
c          dpsi(ij)=dpsi(i,j).
c          if ij=isle is an ocean perimeter point and i`j`=(i`,j`) is
c          a mid-ocean point, with i`=i and j`=j-1, then
c          cfs(i,j)=A(isle,i`j`) is the coefficient of the mid-ocean
c          variable dpsi(i`j`)=dpsi(i,j) in the equation for the island
c          variable dpsi(isle).  note that equations for island
c          variables dpsi(isle) are "non-local" in the sense that
c          they usually contain more than 5 or 9 terms, some of which
c          involve values dpsi(i`j`) outside of a compact 5-point
c          or 9-point neighborhood.
c      type(inv_op)    :: Z
c          the approximate inverse operator Z used at present is a
c          diagonal operator Z(ij,ij) = 1/A(ij,ij).
c          if ij=(i,j) is a mid-ocean point,
c            then Z(i,j)=Z(ij)=1/A(ij)=1/cfdiag(i,j)
c          if ij=isle is an ocean perimeter point, then
c          Z(isle) is replicated at each ocean perimeter T cell
c          bordering land_mass(isle).
c            Z(i,j)=Z(isle)=1/A(isle)=1/sum(A(i,j))
c
c=======================================================================
c
      use mom3_mpp_mod
#include "size.h"
#include "stdunits.h"
#include "taskrows.h"
c
      logical imask, converged, diverging, istop
      character * 16 variable
      character *(*) bc_symm
      dimension guess(imt,js:je), dpsi(imt,js:je), Zres(imt,js:je)
      dimension s(imt,js:je)
      dimension res(imt,js:je), As(imt,js:je), forc(imt,js:je)
      dimension cf(imt,js:je,-1:1,-1:1)
      dimension Z(imt,js:je)
      dimension imask(-mnisle:mnisle)
      dimension iperm(maxipp), jperm(maxipp), iofs(mnisle)
      dimension nippts(mnisle)
      logical do_iterations, continue_iterations
c
#ifdef debug_congrad
c
c-----------------------------------------------------------------------
c     verify that domain(forc) is a subset of domain(A)
c-----------------------------------------------------------------------
c
      call subset (forc, cf(1,js,0,0), js, je, nerror)
      if (nerror .ne. 0) then
        write (stdout,'(a,a,a,i5,a,a)')
     &    'fatal error: ',trim(variable)
     &,   ' forcing is nonzero at ',nerror,' points '
     &,   ' with no equation, i.e., cf(i,j,0,0)=0.'
c       call abort()
      endif
#endif
c
#ifdef symmetry
c
c-----------------------------------------------------------------------
c     check for allowable symmetry conditions at equator
c-----------------------------------------------------------------------
c
      if (bc_symm .ne. 't even' .and. bc_symm .ne. 't odd') then
        write (stdout,*) 
     &'=>Error: type of bc_symm is ',bc_symm,' in congrad.F'
     &, '... Only allowed type is "t even" or "t odd"'
        call abort()
      endif
#endif
c
c-----------------------------------------------------------------------
c     impose boundary conditions on guess
c     dpsi(0) = guess
c-----------------------------------------------------------------------
c
      call border (guess, bc_symm, js, je)
c
      do j=js,je
        do i=1,imt
          dpsi(i,j) = guess(i,j)
        enddo
      enddo
#if defined debug_congrad_parallel
      write (stdout,*) 'ck2d dpsi=',ck2d(dpsi(1,js), js, je)
#endif
c
c-----------------------------------------------------------------------
c     make approximate inverse operator Z (always even symmetry)
c-----------------------------------------------------------------------
c
#if defined debug_congrad_parallel
      do j=js,je
        cf(1,j,0,0) = cf(imt-1,j,0,0)
        cf(imt,j,0,0) = cf(2,j,0,0)
      enddo
      write (stdout,*) 'ck2d cf00=',ck2d(cf(1,js,0,0), js, je)
#endif
      call make_inv (cf, Z, js, je
     &,              imask, iperm, jperm, iofs, nisle, nippts)
      call border (Z, 't even', js, je)
#if defined debug_congrad_parallel
      write (stdout,*) 'ck2d Z=',ck2d(Z(1,js), js, je)
#endif
c
c-----------------------------------------------------------------------
c     res(0)  = forc - A * dpsi(0)
c     impose cyclic and/or symmetry conditions on res(i,j)
c-----------------------------------------------------------------------
c
      if (npt .eq. 5) then
        call op5_vec(cf, dpsi, res, js, je)
      else
        call op9_vec(cf, dpsi, res, js, je)
      endif
      do j=jscomp,jecomp
        do i=2,imt-1
          res(i,j) = forc(i,j) - res(i,j)
        enddo
      enddo
c
      call border (res, bc_symm, js, je)
c
#ifdef parallel_1d
      call mpp_update_domains (res, mom_domain)
#endif
#if defined debug_congrad_parallel
      write (stdout,*) 'In congrad: res ...pn=',pn,' js=',js,' je=',je
      write (stdout,*) 'ck2d res=',ck2d(res(1,js), js, je)
#endif
c
c-----------------------------------------------------------------------
c     Zres(k-1) = Z * res(k-1)
c     see if guess is a solution, bail out to avoid division by zero
c-----------------------------------------------------------------------
c
        k = 0
        call inv_op (Z, res, Zres, js, je
     &,              imask, iperm, jperm, iofs, nisle, nippts)
c
c       set borders of Zres using cyclic/symmetry, if defined.
c
        call border (Zres, bc_symm, js, je)
#ifdef parallel_1d
        call mpp_update_domains (Zres, mom_domain)
#endif
#if defined debug_congrad_parallel
        write (stdout,*) 'congrad: Zres ...pn=',pn,' js=',js,' je=',je
        write (stdout,*) 'ck2d Zres=',ck2d(Zres(1,js), js, je)
#endif
#ifdef debug_congrad
        Zresmax = absvecmax (Zres, imax, jmax, js, je)
        write (stdout, '(a,a,i5,a,e14.7,a,i3,a,i3,a)')
     &  trim(variable), ' iteration',k,'  max Zres = ',Zresmax
     &,          ' at (',imax,',',jmax,')'
#else
        Zresmax = absmax(Zres, js, je)
#endif
c
c       Assume convergence rate of 0.99 to extrapolate error
c
        do_iterations = .true.
        if (100.0 * Zresmax .lt. epsilon) then
	  estimated_error = 100.0 * Zresmax 
          do_iterations = .false.
	endif
c
      if (do_iterations) then  
c-----------------------------------------------------------------------
c     beta(0) = 1
c     s(0)    = zerovector()
c-----------------------------------------------------------------------
c
      betakm1 = 1.0
      do j=js,je
        do i=1,imt
          s(i,j) = 0.0
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     begin iteration loop
c-----------------------------------------------------------------------
c
      diverging = .false.
      istop     = .false.
      continue_iterations = .true.
      k=0
      do while (k .le. max_iterations .and. continue_iterations)
         k=k+1
c
c-----------------------------------------------------------------------
c       Zres(k-1) = Z * res(k-1)
c-----------------------------------------------------------------------
c
        call inv_op (Z, res, Zres, js, je
     &,              imask, iperm, jperm, iofs, nisle, nippts)
c
c       set borders of Zres using cyclic/symmetry, if defined.
c
        call border (Zres, bc_symm, js, je)
#ifdef parallel_1d
        call mpp_update_domains (Zres, mom_domain)
#endif
#if defined debug_congrad_parallel
        write (stdout,*) 'congrad: Zres ...pn=',pn,' js=',js,' je=',je
        write (stdout,*) 'ck2d k=',k,' Zres=',ck2d(Zres(1,js), js, je)
#endif
c
c-----------------------------------------------------------------------
c       beta(k)   = res(k-1) * Zres(k-1)
c-----------------------------------------------------------------------
c
        betak = dot2(Zres, res, js, je)
        if (k .eq. 1) then
          betak_min = abs(betak)
	elseif (k .gt. 2) then
          betak_min = min(betak_min, abs(betak))
          if (abs(betak) .gt. 100.0*betak_min) then
	    write (stdout,'(/2(a/))')
     &      'WARNING: conjugate gradient solver terminated because'
     &,     '         correction steps are diverging.'
	    write (stdout,'(/7(a/))')
     &      'PROBABLE CAUSES:'
     &,     '         1. convergence criterion is too tight...'
     &,     '            roundoff error prevents convergence'
     &,     '     or  2. the solution is beginning to blow up...'
     &,     '            if so, it is extremely unlikely that usable'
     &,     '            results can be obtained in subsequent time'
     &,     '            steps.'
	    write (stdout,'(/3(a/))')
     &      'ERROR:   It is assumed that the solution is blowing up.'
     &,     '         It is extremely unlikely that usable results can'
     &,     '         be obtained in subsequent time steps.'
            if (variable .ne. 'surfpres') then
#ifdef timing
	      call showtimers
#endif
              call abort()
            endif
            diverging = .true.
            smax = absmax(s, js, je)
            step = abs(alpha) * smax
	    estimated_error=step*convergence_rate/(1.0-convergence_rate)
	    istop = .true.
	  endif
	endif
c
c-----------------------------------------------------------------------
c       s(k)      = Zres(k-1) + (beta(k)/beta(k-1)) * s(k-1)
c-----------------------------------------------------------------------
c
        betaquot = betak/betakm1
        do j=js,je
          do i=1,imt
            s(i,j) = Zres(i,j) + betaquot * s(i,j)
          enddo
        enddo
c
c-----------------------------------------------------------------------
c       As(k)     = A * s(k)
c-----------------------------------------------------------------------
c
        if (npt .eq. 5) then
          call op5_vec(cf, s, As, js, je)
        else
          call op9_vec(cf, s, As, js, je)
        endif
c
        call border (As, bc_symm, js, je)
c
#ifdef parallel_1d
        call mpp_update_domains (As, mom_domain) 
#endif
#if defined debug_congrad_parallel
        write (stdout,*) 'congrad: As ...pn=',pn,' js=',js,' je=',je
        write (stdout,*) 'ck2d As=',ck2d(As(1,js), js, je)
#endif
c
c-----------------------------------------------------------------------
c       If s=0 then the division for alpha(k) gives a float exception.
c       Assume convergence rate of 0.99 to extrapolate error.
c       Also assume alpha(k) ~ 1.
c-----------------------------------------------------------------------
c
        s_dot_As = dot2(s, As, js, je)
        if (abs(s_dot_As) .lt. abs(betak)*1.e-10) then
          smax = absmax(s, js, je)
#ifdef debug_congrad
          write (stdout,*) 'Congrad exited on smax=',smax
#endif
	  estimated_error = 100.0 * smax 
	  istop = .true.
	endif
#if defined debug_congrad_parallel
        write (stdout,*) 'congrad: s_dot_As=',s_dot_As
#endif
c
c-----------------------------------------------------------------------
c       alpha(k)  = beta(k) / (s(k) * As(k))
c-----------------------------------------------------------------------
c
        alpha = betak / s_dot_As
#ifdef debug_congrad
        write (stdout,*) 'alpha=',alpha,' beta(k)=',betak
#endif
c
c-----------------------------------------------------------------------
c       update values:
c       dpsi(k)   = dpsi(k-1) + alpha(k) * s(k)
c       res(k)    = res(k-1) - alpha(k) * As(k)
c-----------------------------------------------------------------------
c
        do j=js,je
          do i=1,imt
            dpsi (i,j) = dpsi(i,j) + alpha * s(i,j)
            res  (i,j) = res (i,j) - alpha * As(i,j)
          enddo
        enddo
#ifdef stream_function
        call avg_dist (res, imask, iperm, jperm, iofs, nisle, nippts)
#endif
        call border (res, bc_symm, js, je)
c
#ifdef debug_congrad
        smax = absvecmax (s, imax, jmax, js, je)
        write (stdout,'(a,i5,a,e14.7,a,i3,a,i3,a)')
     &           'iteration',k,'   max correct = ',alpha*smax,
     &           ' at (',imax,',',jmax,')'
#else
        smax = absmax(s, js, je)
#endif
c
c-----------------------------------------------------------------------
c       test for convergence
c       if (estimated_error) < epsilon) exit
c-----------------------------------------------------------------------
c
        step = abs(alpha) * smax
        if (k .eq. 1) then
          step1 = step
          estimated_error = step
          if (step .lt. epsilon) istop = .true.
        else if (step .lt. epsilon) then
          cfactor = log(step/step1)
          convergence_rate = exp(cfactor/(k-1))
	  estimated_error = step*convergence_rate/(1.0-convergence_rate)
#ifdef debug_congrad
          write (stdout,'(a,i6,4(a,e13.4))')
     &       'iteration=', k, ', correction=', step
     &,      ', convergence rate=', convergence_rate
     &,      ', estimated error=', estimated_error
#endif
          if (estimated_error .lt. epsilon) istop = .true.
        endif
c
        betakm1 = betak
c
#ifdef parallel_1d
        call mpp_update_domains (dpsi, mom_domain)
        call mpp_update_domains (res,  mom_domain)
        call mpp_update_domains (s,    mom_domain)
#endif
#if defined debug_congrad_parallel
        write (stdout,*) 'betakm1=',betakm1
        write (stdout,*)' afer updt k=',k,' pn=',pn,' js=',js,' je=',je
        write (stdout,*) 'ck2d k=',k,' dpsi=',ck2d(dpsi(1,js), js, je)
        write (stdout,*) 'ck2d k=',k,'  res=',ck2d(res(1,js), js, je)
        write (stdout,*) 'ck2d k=',k,'  s=',ck2d(s(1,js), js, je)
        if (k .eq. 2) stop
#endif
        if (istop) continue_iterations = .false.
c
      enddo
c
      endif
c-----------------------------------------------------------------------
c     end of iteration loop
c-----------------------------------------------------------------------
c
      if (k .gt. max_iterations) then
        cfactor = log(step/step1)
        convergence_rate = exp(cfactor/(k-1))
	estimated_error = step*convergence_rate/(1.0-convergence_rate)
#ifdef debug_congrad
        write (stdout,'(a,i4,a,a,e10.3,a,e10.3)')
     & '****> congrad did not converge in ',k,' iterations.'
     &,' maximum dpsi is ',absmax(dpsi, js, je)
     &,' estimated error is ',estimated_error
#endif
        converged = .false.
      else
        if (diverging) then
          converged = .false.
        else
          converged = .true.
        endif
      endif

      iterations = k
c
c-----------------------------------------------------------------------
c     return the last increment of dpsi in the argument res
c-----------------------------------------------------------------------
c
      if (iterations .eq. 0) then
        do j=js,je
          do i=1,imt
            res(i,j) = Zres(i,j)
          enddo
        enddo
      else
        do j=js,je
          do i=1,imt
            res(i,j) = alpha * s(i,j)
          enddo
        enddo
      endif
#ifdef parallel_1d
      call mpp_update_domains (res, mom_domain)
#endif
#ifdef debug_congrad
c
c     write out convergence information each call
c
        write (stdout,'(a,l1,a,i4,a,a,a,a,e10.3,a,e10.3)')
     & '      congrad converged = ',converged,' in ',k,' iterations.'
     &,' maximum ',trim(variable),' is '
     &, absmax(dpsi, js, je),' estimated error is ',estimated_error
#endif
c
      return
      end




      function dot2 (dp_vec, res_vec, js, je)
c
c     this dot product produces the correct answers because for
c     ocean perimeter subscripts, ij=isle, the value on a
c     type(dpsi_type) vector, dp_vec(isle)=dp_vec(i,j), i.e., the true
c     value is replicated, and for a type(res_type) vector,
c     res_vec(isle) = sum (res_vec(i,j)), i.e., the true value is the
c     accumulation of the distributed values.
c
      use mom3_mpp_mod
#include "size.h"
#include "taskrows.h"
#include "stdunits.h"
      dimension dp_vec(imt,js:je), res_vec(imt,js:je)
      common /rowsum_dot2/ rowsum(jmt), dot2t
c
      do j=jscomp,jecomp
        rowsum(j) = 0.0
        do i=2,imt-1
          rowsum(j) = rowsum(j) + dp_vec(i,j) * res_vec(i,j)
        enddo
      enddo
#ifdef parallel_1d
c
      call mpp_sync()
c
# if defined debug_congrad_parallel
      write (stdout,*) 'entering dot2: ...pn=',pn 
# endif
      dot2t = 0.0
      do j = 2,jmt-1
         if( pe_of_lat(j).ne.0 )then
             if( pe.EQ.pe_of_lat(j) )call mpp_transmit(
     &            rowsum(j), 1, 0, rowsum(j), 1, NULL_PE )
             if( pe.EQ.0 )call mpp_transmit(
     &            rowsum(j), 1, NULL_PE, rowsum(j), 1, pe_of_lat(j) )
         endif
         if( pe.EQ.0 )dot2t = dot2t + rowsum(j)
      enddo
      call mpp_transmit( dot2t, 1, ALL_PES, dot2t, 1, 0 )
      call mpp_sync()
c
      dot2 = dot2t
# if defined debug_congrad_parallel
      write (stdout,*) 'leaving ...pn=',pn,' dot2=',dot2 
# endif
#else
      dot2 = 0.0
      do j=jscomp,jecomp
        dot2 = dot2 + rowsum(j)
      enddo
#endif
      return
      end


      subroutine op5_vec(cf, dpsi, res, js, je)
c
c                       res = A * dpsi
c
c     this subroutine does not collect the terms of the true value
c     of res(isle) = sum (res(i,j)).  the contributions to the sum
c     remain distributed among the T cells (i,j) that form the
c     ocean perimeter of land_mass(isle).
c
c     at present, borders are not computed [i=1 or imt] [j=1 or jmt]
c
#include "size.h"
#include "taskrows.h"
      dimension cf(imt,js:je,-1:1,-1:1)
      dimension dpsi(imt,js:je)
      dimension res(imt,js:je)
c
      do j=jscomp,jecomp
        do i=2,imt-1
          res(i,j) = cf(i,j, 0, 0) * dpsi(i,j)   +
     &               cf(i,j, 0, 1) * dpsi(i,j+1) +
     &               cf(i,j, 0,-1) * dpsi(i,j-1) +
     &               cf(i,j, 1, 0) * dpsi(i+1,j) +
     &               cf(i,j,-1, 0) * dpsi(i-1,j)
        enddo
      enddo
c
      return
      end


      subroutine op9_vec(cf, dpsi, res, js, je)
c
c                       res = A * dpsi
c
c     this subroutine does not collect the terms of the true value
c     of res(isle) = sum (res(i,j)).  the contributions to the sum
c     remain distributed among the T cells (i,j) that form the
c     ocean perimeter of land_mass(isle).
c
c     at present, borders are not computed [i=1 or imt] [j=1 or jmt]
c
#include "size.h"
#include "taskrows.h"
      dimension cf(imt,js:je,-1:1,-1:1)
      dimension dpsi(imt,js:je)
      dimension res(imt,js:je)
c
      do j=jscomp,jecomp
        do i=2,imt-1
          res(i,j) = cf(i,j, 0, 0) * dpsi(i  ,j  ) +
     &               cf(i,j, 0, 1) * dpsi(i  ,j+1) +
     &               cf(i,j, 0,-1) * dpsi(i  ,j-1) +
     &               cf(i,j, 1, 0) * dpsi(i+1,j  ) +
     &               cf(i,j,-1, 0) * dpsi(i-1,j  ) +
     &               cf(i,j, 1, 1) * dpsi(i+1,j+1) +
     &               cf(i,j,-1, 1) * dpsi(i-1,j+1) +
     &               cf(i,j, 1,-1) * dpsi(i+1,j-1) +
     &               cf(i,j,-1,-1) * dpsi(i-1,j-1)
        enddo
      enddo
c
      return
      end


      subroutine subset (a, b, js, je, nerror)
c
c     verifies that the set of subscripts for which a(i,j) .ne. 0.0
c     is a subset of the set of subscripts for which b(i,j) .ne. 0.0
c
#include "size.h"
#include "stdunits.h"
#include "taskrows.h"
      dimension a(imt,js:je), b(imt,js:je)
c
      nerror = 0
      do j=jscomp,jecomp
        do i=2,imt-1
          if (a(i,j) .ne. 0.0 .and. b(i,j) .eq. 0.0) then
            nerror = nerror + 1
            write (stdout,'(a)') ' forcing is reset to zero'
c           set forcing (i.e., a(i,j)) to zero
            a(i,j) = 0.0
          endif
        enddo
      enddo
c
      return
      end


      subroutine inv_op (Z, res, Zres, js, je
     &,                  imask, iperm, jperm, iofs, nisle, nippts)
c
c
c     apply and approximate inverse Z or the operator A
c
c     res is type(res_type), i.e., perimeter values res(isle)
c         are the sum of the distributed contributions res(i,j)
c     Zres is type(dpsi_type), i.e., perimeter values Zres(isle)
c         must be replicated at each perimeter point Zres(i,j)
c
c     borders  of Zres [i=1 or imt] [j=1 or jmt] must be defined
c     and must satisfy cyclic and/or symmetry, if defined.
c
c     currently, Z is diagonal:  Z(ij) = 1/A(ij)
c     and is stored in type(dpsi_type) format, i.e., Z(isle) is
c     replicated and stored in each Z(i,j) in the perimeter of
c     land_mass(isle).
c
#include "size.h"
c
      logical imask
      dimension imask (-mnisle:mnisle)
      dimension iperm(maxipp), jperm(maxipp)
      dimension iofs (mnisle)
      dimension nippts(mnisle)
c
      dimension Z(imt,js:je)
      dimension res(imt,js:je)
      dimension Zres(imt,js:je)
c
      do j=js,je
        do i=1,imt
          Zres(i,j) = Z(i,j) * res(i,j)
        enddo
      enddo
#ifdef stream_function
c
c     sum contributions to Zres(isle)
c     distribute Zres(isle) to all perimeter points
c
      call sum_dist (Zres, imask, iperm, jperm, iofs, nisle, nippts)
c
#endif
      return
      end


      function absvecmax (res, imax, jmax, js, je)
c
#include "size.h"
#include "taskrows.h"
      common /absvecmaxcom/ big, c, rimax, rjmax
      dimension res(imt,js:je)
c
      big = 0.0
c
      do j=jscomp,jecomp
        do i=2,imt-1
          if (abs(res(i,j)) .gt. big) then
            big = abs(res(i,j))
            rimax = i
            rjmax = j
          endif
        enddo
      enddo
#ifdef parallel_1d
      c   = big
!      call shmem_barrier
# if defined debug_congrad_parallel
      write (stdout,*) 'entering absvecmax: ...pn=',pn 
# endif
      if (pn .eq. 1) then
        len = 1
        do np=1,num_processors-1
          call shmem_get8(big, big, len, np)
          if (c .gt. big) then
	    c = big
            call shmem_get8(rimax, rimax, len, np)
            call shmem_get8(rjmax, rjmax, len, np)
          endif
        enddo
        do np=1,num_processors-1
          call shmem_put8(c, c, len, np)
          call shmem_put8(rimax, rimax, len, np)
          call shmem_put8(rjmax, rjmax, len, np)
        enddo
      endif
!      call shmem_barrier
      big = c
#endif
      absvecmax = big
      imax = nint(rimax)
      jmax = nint(rjmax)
#if defined debug_congrad_parallel
      write (stdout,*) 
     &'leaving: ...pn=',pn, 'absvecmax=',absvecmax, imax, jmax 
#endif
      return
      end


      function absmax (f, js, je)
c
#include "size.h"
#include "taskrows.h"
#include "stdunits.h"
      common /absmaxcom/ big, c
      dimension f(imt,js:je)
c
      big = 0.0
c
      do j=js,je
        do i=1,imt
	  big = max(big, abs(f(i,j)))
        enddo
      enddo
#ifdef parallel_1d
      c   = big
!      call shmem_barrier
# if defined debug_congrad_parallel
      write (stdout,*) 'entering absmax: ...pn=',pn 
# endif
      if (pn .eq. 1) then
        len = 1
        do np=1,num_processors-1
          call shmem_get8(big, big, len, np)
	  c = max(c,big)
        enddo
        do np=1,num_processors-1
          call shmem_put8(c, c, len, np)
        enddo
      endif
!      call shmem_barrier
      big = c
#endif
      absmax = big
#if defined debug_congrad_parallel
      write (stdout,*) 'leaving  ...pn=',pn,' absmax=',absmax 
#endif
      return
      end


      subroutine make_inv (cf, Z, js, je
     &,              imask, iperm, jperm, iofs, nisle, nippts)

c     (mpp items)   v. balaji             e-mail=> vb@gfdl.gov
      use mom3_mpp_mod
c
c     construct an approximate inverse Z to A
c
c     Z will be diagonal:  Z(ij) = 1/A(ij)
c     and values for ocean perimeter entries Z(isle) will be replicated
c     at all T cells Z(i,j) in the ocean perimeter of land_mass(isle).
c
c     T cells (i,j) for which there is no diagonal coefficient
c     i.e., A(ij)=A(i,j)=0, are masked off by assigning Z(i,j)=0.
c     there are effectively no equations and no variables dpsi(i,j)
c     at these points.
c
#include "size.h"
#include "stdunits.h"
#include "taskrows.h"
      logical imask
      dimension imask (-mnisle:mnisle)
      dimension iperm(maxipp), jperm(maxipp)
      dimension iofs (mnisle)
      dimension nippts(mnisle)
c
      dimension cf(imt,js:je,-1:1,-1:1)
      dimension Z(imt,js:je)
c
c     copy diagonal coefficients of A to Z
c
      do j=jscomp,jecomp
        do i=2,imt-1
          Z(i,j) = cf(i,j,0,0)
        enddo
      enddo
#ifdef stream_function
c
c     for each land_mass(isle),
c     sum the contributions to cfdiag(isle)=A(isle,isle)
c     now stored in Z(i,j) at ocean perimeter T cells and replicate
c     the sum in all Z(i,j) for which (i,j) is in ocean perimeter
c     of land_mass(isle).
c
      call sum_dist (Z, imask, iperm, jperm, iofs, nisle, nippts)
#endif
c
c     now invert Z
c
      do j=jscomp,jecomp
        do i=2,imt-1
          if (Z(i,j) .ne. 0.0) then
            Z(i,j) = 1/Z(i,j)
          else
            Z(i,j) = 0.0
          endif
        enddo
      enddo
#ifdef stream_function
c
c     make inverse zero on island perimeters that are not integrated
c
      do isle=1,nisle
        if (.not. imask(isle)) then
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            Z(i,j) = 0.0
          enddo
        endif
      enddo
#endif
c
#if defined parallel_1d
      do j=jscomp,jecomp
        Z(1,j) = Z(imt-1,j)
	Z(imt,j) = Z(2,j)
      enddo

      call mpp_update_domains( Z, mom_domain )
#endif
      return
      end



#ifdef stream_function
      subroutine sum_dist (Zres
     &,                    imask, iperm, jperm, iofs, nisle, nippts)
c
c
c     sum contributions to Zres(isle)
c     distribute Zres(isle) to all perimeter points
c
c     this subroutine converts a type(res_type) vector with
c     distributed contributions to perimeter values
c        Zres(isle) = sum (Zres(i,j))
c     into a type (dpsi_type) vector with replicated values
c     for land_mass perimeters
c        Zres(isle) = Zres(i,j)
c     for all (i,j) in the ocean perimeter of land_mass(isle).
c
# include "size.h"
      logical imask
      dimension imask (-mnisle:mnisle)
      dimension iperm(maxipp), jperm(maxipp)
      dimension iofs (mnisle)
      dimension nippts(mnisle)
c
      dimension Zres(imt,jmt)
      dimension Zresisle(mnisle)
c
c     sum contributions to Zres(isle)
c
# ifdef timing
      call tic ('isle','loop1')
# endif
      do isle=1,nisle
        if (imask(isle)) then
          Zresisle(isle) = 0.0
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            Zresisle(isle) = Zresisle(isle) + Zres(i,j)
          enddo
        endif
      enddo
# ifdef timing
      call toc ('isle','loop1')
# endif
c
c     distribute Zres(isle) to all perimeter points
c
# ifdef timing
      call tic ('isle','loop2')
# endif
      do isle=1,nisle
        if (imask(isle)) then
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            Zres(i,j) = Zresisle(isle)
          enddo
        endif
      enddo
# ifdef timing
      call toc ('isle','loop2')
# endif
c
      return
      end

      subroutine avg_dist (Zres,
     &              imask, iperm, jperm, iofs, nisle, nippts)
c
c     avg contributions to Zres(isle)
c     distribute Zres(isle) to all perimeter points
c
c     this subroutine converts a type(res_type) vector with
c     distributed contributions to perimeter values
c        Zres(isle) = avg (Zres(i,j))
c     into a type (dpsi_type) vector with replicated values
c     for land_mass perimeters
c        Zres(isle) = Zres(i,j)
c     for all (i,j) in the ocean perimeter of land_mass(isle).
c
# include "size.h"
      logical imask
      dimension imask (-mnisle:mnisle)
      dimension iperm(maxipp), jperm(maxipp)
      dimension iofs (mnisle)
      dimension nippts(mnisle)
c
      dimension Zres(imt,jmt)
      dimension Zresisle(mnisle)
c
c     avg contributions to Zres(isle)
c
# ifdef timing
      call tic ('isle','loop3')
# endif
      do isle=1,nisle
        if (imask(isle)) then
          Zresisle(isle) = 0.0
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            Zresisle(isle) = Zresisle(isle) + Zres(i,j)
          enddo
        endif
      enddo
# ifdef timing
      call toc ('isle','loop3')
# endif
c
c     distribute Zres(isle) to all perimeter points
c
# ifdef timing
      call tic ('isle','loop4')
# endif
      do isle=1,nisle
        if (imask(isle)) then
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            Zres(i,j) = Zresisle(isle)/nippts(isle)
          enddo
        endif
      enddo
# ifdef timing
      call toc ('isle','loop4')
# endif
c
      return
      end
#endif
#else

      subroutine congr()
      return
      end
#endif      
