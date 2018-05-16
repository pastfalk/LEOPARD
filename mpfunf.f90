!*****************************************************************************

!  MPFUN-Fort: A thread-safe arbitrary precision computation package
!  Precision level declaration module (module MPFUNF)

!  Revision date:  19 Jul 2015

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired) and University of California, Davis
!     Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2015 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to arbitrarily high numeric precision, by making only relatively
!    minor changes to existing Fortran-90 programs.  All basic arithmetic
!    operations and transcendental functions are supported, together with several
!    special functions.

!    In addition to fast execution times, one key feature of this package is a
!    100% THREAD-SAFE design, which means that user-level applications can be
!    easily converted for parallel execution, say using a threaded parallel
!    environment such as OpenMP.  There are NO global shared variables (except
!    static compile-time data), and NO initialization is necessary unless
!    extremely high precision (> 19,500 digits) is required.

!  DOCUMENTATION:
!    A detailed description of this package, and instructions for compiling
!    and testing this program on various specific systems are included in the
!    README file accompanying this package, and, in more detail, in the
!    following technical paper:
   
!    David H. Bailey, "MPFUN2015: A thread-safe arbitrary precision package," 
!    available at http://www.davidhbailey.com/dhbpapers/mpfun2015.pdf. 

!  DESCRIPTION OF THIS MODULE (MPFUNF):
!    This module defines the default precision level in digits (mpipl) and
!    the equivalent precision level in words (mpwds), which is calculated
!    below according to the formula:
!       mpwds = int (mpipl / mpdpw + 2)
!    where mpdpw = 14.449439791871... is an approximation to log10(2^48).
!    This precision level is the working precision level for all operations,
!    unless one dynamically varies the working precision level, in which case
!    this is the maximum precision level for the entire calculation.

module mpfunf
use mpfuna
use mpfunb
use mpfunc
use mpfund
use mpfune
implicit none
integer, public:: mpipl

!  *** Set the default precision level (in digits) here.

parameter (mpipl = 1200)

!----------------------------------------------------------------------------

!  *** Do not change the following code (in normal usage).

integer, public:: mpwds, mpwds6
parameter (mpwds = int (mpipl / mpdpw + 2.d0), mpwds6 = mpwds + 6)

end module mpfunf
