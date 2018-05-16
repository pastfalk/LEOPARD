!  MPFUN2015: A thread-safe Fortran arbitrary precision computation package
!  Main module (module MPMODULE) -- references all other modules for user.
!  Version date:  24 May 2015

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired) and University of California, Davis
!     Email: dhbailey@lbl.gov

!  All software in this package (c) 2015 David H. Bailey

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to arbitrarily high numeric precision, by making only relatively
!    minor changes to existing Fortran-90 programs (mostly changes to type
!    statements).  All basic arithmetic operations and transcendental functions
!    are supported.  Advanced techniques, including FFT-based multiplication and
!    quadratically convergent transcendental algorithms, are employed.

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

!  NOTE ON THE PARAMETERS MPIPL AND MPWDS:
!    mpipl, the default (maximum) precision level in digits, mpipl, is set
!    by the user in module MPMODF.  mpwds is the equivalent in
!    words.  It is calculated in MPMODF according to the formula:
!       mpwds = int (mpipl / mpdpw + 2)
!    where mpdpw = 14.449439791871... is an approximation to log10(2^48).
!    This precision level is the working precision level for all operations,
!    unless one dynamically varies the working precision level, in which case
!    this is the maximum precision level for the entire calculation.
    
module mpmodule
use mpfuna
use mpfunb
use mpfunc
use mpfund
use mpfune
use mpfunf
use mpfung

end module mpmodule

