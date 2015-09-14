!-----------------------------------------------------------------------
!             Checksums for bitwise comparison of real data
!
! AUTHOR: V. Balaji (vb@gfdl.gov)
!         SGI/GFDL Princeton University
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.  
!-----------------------------------------------------------------------

!these are used to determine hardware/OS/compiler

#if defined(_CRAY) || defined(__sgi)
#define SGICRAY
#endif

#if defined(_CRAY) && !defined(_CRAYT3E) && !defined(_CRAYT3D)
#define CRAYPVP
#endif

#if defined(_CRAYT3E) || defined(_CRAYT3D) || defined(__sgi)
#define SGICRAY_MPP
#endif

!machines that support Cray pointers
#if defined(SGICRAY) || defined(__alpha)
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

module chksum_mod
!this module contains long integer function chksum.
!this does an exact sum of its argument as an integer
!(floating point checksums do not guarantee that all the bits match)
!chksum_mod uses Cray pointers to equivalence real/cmplx to int,
!then sums the resulting array as int.
!this works with IEEE words, where any bit sequence is a valid integer
  implicit none
  private
  character(len=256), private :: version='$Id$'

  interface chksum
     module procedure chksum_int8
     module procedure chksum_int8_scalar
     module procedure chksum_real8
     module procedure chksum_real8_scalar
     module procedure chksum_cmplx8
     module procedure chksum_cmplx8_scalar
  end interface

  public :: chksum
  
  contains
    
    integer(LONG_KIND) function chksum_int8( array, length )
      integer, intent(in) :: length
      integer(LONG_KIND), intent(in) :: array(length)
!sum is f90 intrinsic
      chksum_int8 = sum(array)
      return
    end function chksum_int8

    integer(LONG_KIND) function chksum_int8_scalar( scalar, length )
      integer, intent(in) :: length
      integer(LONG_KIND), intent(in) :: scalar
      integer(LONG_KIND) :: array(length)
#ifdef use_CRI_pointers
      pointer( ptr, array )
      ptr = LOC(scalar)

      chksum_int8_scalar = chksum_int8(array,length)
#else
      call error( FATAL, 'CHKSUM_INT8_SCALAR currently requires CRI pointers.' )
#endif
      return
    end function chksum_int8_scalar

    integer(LONG_KIND) function chksum_real8( array, length )
      integer, intent(in) :: length
      real(DOUBLE_KIND), intent(in) :: array(length)
      integer(LONG_KIND) :: int_array(length)
#ifdef use_CRI_pointers
      pointer( ptr, int_array )
      ptr = LOC(array)

      chksum_real8 = chksum_int8(int_array,length)
#else
      call error( FATAL, 'CHKSUM_REAL8 currently requires CRI pointers.' )
#endif
      return
    end function chksum_real8

    integer(LONG_KIND) function chksum_real8_scalar( scalar, length )
      integer, intent(in) :: length
      real(DOUBLE_KIND), intent(in) :: scalar
      integer(LONG_KIND) :: int_array(length)
#ifdef use_CRI_pointers
      pointer( ptr, int_array )
      ptr = LOC(scalar)

      chksum_real8_scalar = chksum_int8(int_array,length)
#else
      call error( FATAL, 'CHKSUM_REAL8_SCALAR currently requires CRI pointers.' )
#endif
      return
    end function chksum_real8_scalar

    integer(LONG_KIND) function chksum_cmplx8( array, length )
      integer, intent(in) :: length
      complex(DOUBLE_KIND), intent(in) :: array(length)
      integer(LONG_KIND) :: int_array(length*2)
#ifdef use_CRI_pointers
      pointer( ptr, int_array )
      ptr = LOC(array)

      chksum_cmplx8 = chksum_int8(int_array,length*2)
#else
      call error( FATAL, 'CHKSUM_CMPLX8 currently requires CRI pointers.' )
#endif
      return
    end function chksum_cmplx8

    integer(LONG_KIND) function chksum_cmplx8_scalar( scalar, length )
      integer, intent(in) :: length
      complex(DOUBLE_KIND), intent(in) :: scalar
      integer(LONG_KIND) :: int_array(length*2)
#ifdef use_CRI_pointers
      pointer( ptr, int_array )
      ptr = LOC(scalar)

      chksum_cmplx8_scalar = chksum_int8(int_array,length*2)
#else
      call error( FATAL, 'CHKSUM_CMPLX8_SCALAR currently requires CRI pointers.' )
#endif
      return
    end function chksum_cmplx8_scalar

  end module chksum_mod

#if defined test_chksum
!
! Reals smaller than epsilon\approx 10^{-16} will not
! figure in a real checksum. This is illustrated in the
! following test program.  What is shown is that a floating point add
! 1 + 0.5*epsilon gives the same bits as 1, whereas an integer add gives
! different answers, as it should, if you want your checksum to tell you
! the bitwise exact state of your answers. 
!
    program test
    use chksum_mod

    integer, parameter :: n=1024
    real :: a, b(n)

    call random_number(b)
!    print '(x,z18)', chksum(b, size(b))

    a = 1.
    print '(x,a,z18)', 'chksum(1)=', chksum(a,1)
    a = 1. + epsilon(1.)
    print '(x,a,z18)', 'chksum(1+epsilon)=', chksum(a,1)
    a = 1. + 0.5*epsilon(1.)
    print '(x,a,z18)', 'chksum(1+0.5*epsilon)=', chksum(a,1)
  end program test
#endif
