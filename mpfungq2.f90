!*****************************************************************************

!  MPFUN-Fort: A thread-safe arbitrary precision computation package
!  Language interface module (module MPFUNG)
!  Variant 2: Precision level specifications are *required*; real*16 support.
!  Search for !> for version differences.

!  Revision date:  8 Feb 2016

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

!  DESCRIPTION OF THIS MODULE (MPFUNG):
!    This module contains all high-level Fortran-90 language interfaces.
!    There are two variants of this module, which take different approaches
!    to dynamically changing working precision within an application.  Variant
!    1 allows mixed-mode assignments, and does not require precision level
!    specifications in certain functions, whereas variant 2 does not permit
!    mixed-mode assignments, and requires precision level specifications in
!    certain functions.  For each of these two variants, there are "Q"
!    versions that provide interfaces for real*16 on systems that support
!    real*16.  See documentation for details.
 
module mpfung
use mpfuna
use mpfunb
use mpfunc
use mpfund
use mpfune
use mpfunf
implicit none
  
!   The mp_real and mp_complex datatypes are defined here:

type mp_real
  sequence
  double precision mpr(0:mpwds6-1)
end type
type mp_complex
  sequence
  double precision mpc(0:2*mpwds6-1)
end type

!  Assignments and the five arithmetic operators:

private &
  mp_eqrr, mp_eqdr, mp_eqrd, mp_eqir, mp_eqri, mp_eqra, mp_eqrz, &
  mp_eqzr, mp_eqzz, mp_eqdcz, mp_eqzdc, &
  mp_addrr, mp_adddr, mp_addrd, mp_addir, mp_addri, mp_addzz, &
  mp_adddcz, mp_addzdc, mp_addrz, mp_addzr, &
  mp_subrr, mp_subdr, mp_subrd, mp_subir, mp_subri, mp_subzz, &
  mp_subdcz, mp_subzdc, mp_subrz, mp_subzr, mp_negr, mp_negz, &
  mp_mulrr, mp_muldr, mp_mulrd, mp_mulir, mp_mulri, mp_mulzz, &
  mp_muldcz, mp_mulzdc, mp_mulrz, mp_mulzr, mp_muldz, mp_mulzd, &
  mp_divrr, mp_divdr, mp_divrd, mp_divir, mp_divri, mp_divzz, &
  mp_divdcz, mp_divzdc, mp_divrz, mp_divzr, &
  mp_expri, mp_exprr, mp_expzi, mp_expzz, mp_exprz, mp_expzr

!  The six comparison tests:

private &
  mp_eqtrr, mp_eqtdr, mp_eqtrd, mp_eqtir, mp_eqtri, mp_eqtzz, &
  mp_eqtdcz, mp_eqtzdc, mp_eqtrz, mp_eqtzr, &
  mp_netrr, mp_netdr, mp_netrd, mp_netir, mp_netri, mp_netzz, &
  mp_netdcz, mp_netzdc, mp_netrz, mp_netzr, &
  mp_letrr, mp_letdr, mp_letrd, mp_letir, mp_letri, &
  mp_getrr, mp_getdr, mp_getrd, mp_getir, mp_getri, &
  mp_lttrr, mp_lttdr, mp_lttrd, mp_lttir, mp_lttri, &
  mp_gttrr, mp_gttdr, mp_gttrd, mp_gttir, mp_gttri

!  Algebraic, transcendental and type conversion functions:
  
private &
  mp_absr, mp_absz, mp_acos, mp_acosh, mp_agm, mp_aimag, mp_aint, &
  mp_anint, mp_asin, mp_asinh, mp_atan, mp_atan2, mp_atanh, mp_ator1, &
  mp_atorn, mp_berne, mp_besselj, mp_bessel_j0, mp_bessel_j1, &
  mp_bessel_jn, mp_ccos, mp_cexp, mp_clog, mp_conjg, mp_cos, mp_cosh, &
  mp_csin, mp_csqrt, mp_cssh, mp_cssn, mp_dctoz, mp_dctoz2, mp_mdi, &
  mp_dtor, mp_dtor2, mp_eform, mp_exp, mp_fform, mp_gamma, mp_hypot, &
  mp_incgamma, mp_init, mp_log, mp_log2, mp_max, mp_min, mp_nrt, &
  mp_pi, mp_prodd, mp_quotd, mp_readr1, mp_readr2, mp_readr3, &
  mp_readr4, mp_readr5, mp_readz1, mp_readz2, mp_readz3, mp_readz4, &
  mp_readz5, mp_rtod, mp_rtor, mp_rtoz, mp_setwp, mp_sign, mp_sin, &
  mp_sinh, mp_sqrt, mp_tan, mp_tanh, mp_wprec, mp_wprecz, mp_writer, &
  mp_writez, mp_zeta, mp_zetaem, mp_ztodc, mp_ztor, mp_ztoz

!>  This line is for real*16 support.  No need to change either way.

private mp_rtoq, mp_qtor

!  Operator extension interface blocks:

interface assignment (=)
  module procedure mp_eqrr
  module procedure mp_eqdr
  module procedure mp_eqir
  module procedure mp_eqrz
  module procedure mp_eqzr
  module procedure mp_eqzz
  module procedure mp_eqdcz

!>  In version #1, the next four module precedure lines are uncommented;
!>  In version #2 they are commented out.

!  module procedure mp_eqrd
!  module procedure mp_eqri
!  module procedure mp_eqra
!  module procedure mp_eqzdc
end interface

interface operator (+)
  module procedure mp_addrr
  module procedure mp_adddr
  module procedure mp_addrd
  module procedure mp_addir
  module procedure mp_addri
  module procedure mp_addzz
  module procedure mp_adddcz
  module procedure mp_addzdc
  module procedure mp_addrz
  module procedure mp_addzr
end interface

interface operator (-)
  module procedure mp_subrr
  module procedure mp_subdr
  module procedure mp_subrd
  module procedure mp_subir
  module procedure mp_subri
  module procedure mp_subzz
  module procedure mp_subdcz
  module procedure mp_subzdc
  module procedure mp_subrz
  module procedure mp_subzr
  module procedure mp_negr
  module procedure mp_negz
end interface

interface operator (*)
  module procedure mp_mulrr
  module procedure mp_muldr
  module procedure mp_mulrd
  module procedure mp_mulir
  module procedure mp_mulri
  module procedure mp_mulzz
  module procedure mp_muldcz
  module procedure mp_mulzdc
  module procedure mp_mulrz
  module procedure mp_mulzr
  module procedure mp_muldz
  module procedure mp_mulzd
end interface

interface operator (/)
  module procedure mp_divrr
  module procedure mp_divdr
  module procedure mp_divrd
  module procedure mp_divir
  module procedure mp_divri
  module procedure mp_divzz
  module procedure mp_divdcz
  module procedure mp_divzdc
  module procedure mp_divrz
  module procedure mp_divzr
end interface

interface operator (**)
   module procedure mp_expri
   module procedure mp_exprr
   module procedure mp_expzi
   module procedure mp_expzz
   module procedure mp_exprz
   module procedure mp_expzr
end interface

interface operator (==)
  module procedure mp_eqtrr
  module procedure mp_eqtdr
  module procedure mp_eqtrd
  module procedure mp_eqtir
  module procedure mp_eqtri
  module procedure mp_eqtzz
  module procedure mp_eqtdcz
  module procedure mp_eqtzdc
  module procedure mp_eqtrz
  module procedure mp_eqtzr
end interface

interface operator (/=)
  module procedure mp_netrr
  module procedure mp_netdr
  module procedure mp_netrd
  module procedure mp_netir
  module procedure mp_netri
  module procedure mp_netzz
  module procedure mp_netdcz
  module procedure mp_netzdc
  module procedure mp_netrz
  module procedure mp_netzr
end interface

interface operator (<=)
  module procedure mp_letrr
  module procedure mp_letdr
  module procedure mp_letrd
  module procedure mp_letir
  module procedure mp_letri
end interface

interface operator (.ge.)
  module procedure mp_getrr
  module procedure mp_getdr
  module procedure mp_getrd
  module procedure mp_getir
  module procedure mp_getri
end interface

interface operator (<)
  module procedure mp_lttrr
  module procedure mp_lttdr
  module procedure mp_lttrd
  module procedure mp_lttir
  module procedure mp_lttri
end interface

interface operator (>)
  module procedure mp_gttrr
  module procedure mp_gttdr
  module procedure mp_gttrd
  module procedure mp_gttir
  module procedure mp_gttri
end interface

!  MP generic function interface blogs, listed alphabetically by interface name:

interface abs
  module procedure mp_absr
  module procedure mp_absz
end interface

interface acos
  module procedure mp_acos
end interface

interface acosh
  module procedure mp_acosh
end interface

interface agm
  module procedure mp_agm
end interface

interface aimag
  module procedure mp_aimag
end interface

interface aint
  module procedure mp_aint
end interface

interface anint
  module procedure mp_anint
end interface

interface asin
  module procedure mp_asin
end interface

interface asinh
  module procedure mp_asinh
end interface

interface atan
  module procedure mp_atan
end interface

interface atan2
  module procedure mp_atan2
end interface

interface atanh
  module procedure mp_atanh
end interface

interface berne
  module procedure mp_berne
end interface

interface besselj
  module procedure mp_besselj
end interface

interface bessel_j0
  module procedure mp_bessel_j0
end interface

interface bessel_j1
  module procedure mp_bessel_j1
end interface

interface bessel_jn
  module procedure mp_bessel_jn
end interface

interface conjg
  module procedure mp_conjg
end interface

interface cos
  module procedure mp_cos
  module procedure mp_ccos
end interface

interface cosh
  module procedure mp_cosh
end interface

interface dble
  module procedure mp_rtod
end interface

interface dcmplx
  module procedure mp_ztodc
end interface

interface exp
  module procedure mp_exp 
  module procedure mp_cexp
end interface

interface gamma
  module procedure mp_gamma
end interface

interface hypot
  module procedure mp_hypot
end interface

interface log
  module procedure mp_log
  module procedure mp_clog
end interface

interface max
  module procedure mp_max
end interface

interface min
  module procedure mp_min
end interface

interface mpcmplx
  module procedure mp_dctoz
  module procedure mp_rtoz
  module procedure mp_ztoz
end interface

interface mpcmplxdc
  module procedure mp_dctoz2
end interface

interface mpcssh
  module procedure mp_cssh
end interface

interface mpcssn
  module procedure mp_cssn
end interface

interface mpeform
  module procedure mp_eform
end interface

interface mpfform
  module procedure mp_fform
end interface

interface incgamma
  module procedure mp_incgamma
end interface

interface mpinit
  module procedure mp_init
end interface

interface mplog2
  module procedure mp_log2
end interface

interface mpmdi
  module procedure mp_mdi
end interface

interface mpnrt
  module procedure mp_nrt
end interface

interface mppi
  module procedure mp_pi
end interface

interface mpprodd
  module procedure mp_prodd
end interface

interface mpquotd
  module procedure mp_quotd
end interface

interface mpread
  module procedure mp_readr1
  module procedure mp_readr2
  module procedure mp_readr3
  module procedure mp_readr4
  module procedure mp_readr5
  module procedure mp_readz1
  module procedure mp_readz2
  module procedure mp_readz3
  module procedure mp_readz4
  module procedure mp_readz5
end interface

interface mpreal
  module procedure mp_ator1
  module procedure mp_atorn
  module procedure mp_dtor
  module procedure mp_rtor
  module procedure mp_ztor

!>  If real*16 is supported, uncomment this line; otherwise commented.

  module procedure mp_qtor
end interface

interface mpreald
  module procedure mp_dtor2
end interface

interface mpwprec
  module procedure mp_wprec
  module procedure mp_wprecz
end interface

interface mpwrite
  module procedure mp_writer
  module procedure mp_writez
end interface

!>  If real*16 is supported, uncomment these three lines; otherwise commented.

interface qreal
  module procedure mp_rtoq
end interface

interface sign
  module procedure mp_sign
end interface

interface sin
  module procedure mp_sin
  module procedure mp_csin
end interface

interface sinh
  module procedure mp_sinh
end interface

interface sqrt
  module procedure mp_sqrt
  module procedure mp_csqrt
end interface

interface tan
  module procedure mp_tan
end interface

interface tanh
  module procedure mp_tanh
end interface

interface zeta
  module procedure mp_zeta
end interface

interface zetaem
  module procedure mp_zetaem
end interface

contains

!  This routine outputs an error message if iprec exceeds mpwds.

  function mp_setwp (iprec)
    integer mp_setwp          
    integer, intent (in):: iprec
    if (iprec > mpwds) then
      write (mpldb, 1)
1       format ( &
        '*** MP_SETWP: requested precision level exceeds default precision.'/ &
        'Increase default precision in module MPMODF.')
      call mpabrt (98)
    endif
    mp_setwp = iprec
  end function

!  Assignment routines:

  subroutine mp_eqrr (ra, rb)
    implicit none
    type (mp_real), intent (out):: ra
    type (mp_real), intent (in):: rb
    integer mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    ra%mpr(0) = mpwds6
    call mpeq (rb%mpr, ra%mpr, mpnw) 
    return
  end subroutine

  subroutine mp_eqdr (da, rb)
    implicit none
    double precision, intent (out):: da
    type (mp_real), intent (in):: rb
    integer ib, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    call mpmdc (rb%mpr, da, ib, mpnw)
    da = da * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqrd (ra, db)
    implicit none
    type (mp_real), intent (out):: ra
    double precision, intent (in):: db
    integer i1, mpnw
    mpnw = mpwds
    ra%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (db, i1, ra%mpr, mpnw)
    return
  end subroutine

  subroutine mp_eqir (ia, rb)
    implicit none
    integer, intent (out):: ia
    type (mp_real), intent (in):: rb
    double precision da
    integer ib, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    call mpmdc (rb%mpr, da, ib, mpnw)
    ia = da * 2.d0 ** ib
    return
  end subroutine

  subroutine mp_eqri (ra, ib)
    implicit none
    type (mp_real), intent (out):: ra
    integer, intent (in):: ib
    double precision db
    integer i1, mpnw
    mpnw = mpwds
    ra%mpr(0) = mpwds6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, ra%mpr, mpnw)
    return
  end subroutine

  subroutine mp_eqra (ra, ab)
    implicit none
    type (mp_real), intent (out):: ra
    character*(*), intent (in):: ab
    character*1 :: chr1(len(ab))
    integer i, l1, mpnw
    mpnw = mpwds
    l1 = len (ab)
    do i = 1, l1
      chr1(i) = ab(i:i)
    enddo
    ra%mpr(0) = mpwds6
    call mpctomp (chr1, l1, ra%mpr, mpnw) 
    return
  end subroutine

  subroutine mp_eqzz (za, zb)
    implicit none
    type (mp_complex), intent (out):: za
    type (mp_complex), intent (in):: zb
    integer l1, l2, mpnw
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    l2 = mpwds6
    za%mpc(0) = mpwds6
    za%mpc(l2) = mpwds6
    call mpceq (zb%mpc, za%mpc, mpnw) 
    return
  end subroutine

  subroutine mp_eqdcz (dca, zb)
    implicit none
    complex (kind (0.d0)), intent (out):: dca
    type (mp_complex), intent (in):: zb
    integer l1, mpnw, n1, n2
    double precision d1, d2
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    call mpmdc (zb%mpc, d1, n1, mpnw)
    d1 = d1 * 2.d0 ** n1
    call mpmdc (zb%mpc(l1), d2, n2, mpnw)
    d2 = d2 * 2.d0 ** n2
    dca = dcmplx (d1, d2)
    return
  end subroutine

  subroutine mp_eqzdc (za, dcb)
    implicit none
    type (mp_complex), intent (out):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer l1, mpnw
    mpnw = mpwds
    l1 = mpwds6
    za%mpc(0) = mpwds6
    za%mpc(l1) = mpwds6
    call mpdmc40 (dble (dcb), 0, za%mpc, mpnw)
    call mpdmc40 (aimag (dcb), 0, za%mpc(l1), mpnw)
    return
  end subroutine

  subroutine mp_eqrz (ra, zb)
    implicit none
    type (mp_real), intent (out):: ra
    type (mp_complex), intent (in):: zb
    integer l1, mpnw
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    ra%mpr(0) = mpwds6
    call mpeq (zb%mpc, ra%mpr, mpnw) 
    return
  end subroutine

  subroutine mp_eqzr (za, rb)
    implicit none
    type (mp_complex), intent (out):: za
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    integer l1, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    call mpdmc (0.d0, 0, r1%mpr, mpnw)
    l1 = mpwds6
    za%mpc(0) = mpwds6
    za%mpc(l1) = mpwds6
    call mpeq (rb%mpr, za%mpc, mpnw)
    call mpeq (r1%mpr, za%mpc(l1), mpnw)
    return
  end subroutine

!  Addition routines:

  function mp_addrr (ra, rb)
    implicit none
    type (mp_real):: mp_addrr
    type (mp_real), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    mp_addrr%mpr(0) = mpwds6
    call mpadd (ra%mpr, rb%mpr, mp_addrr%mpr, mpnw) 
    return
  end function

  function mp_adddr (da, rb)
    implicit none
    type (mp_real):: mp_adddr
    double precision, intent (in):: da
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    integer i1, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    mp_adddr%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpadd (r1%mpr, rb%mpr, mp_adddr%mpr, mpnw) 
    return
  end function

  function mp_addrd (ra, db)
    implicit none
    type (mp_real):: mp_addrd
    type (mp_real), intent (in):: ra
    double precision, intent (in):: db
    type (mp_real) r1
    integer i1, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    mp_addrd%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpadd (ra%mpr, r1%mpr, mp_addrd%mpr, mpnw) 
    return
  end function

  function mp_addir (ia, rb)
    implicit none
    type (mp_real):: mp_addir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    double precision da
    integer i1, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    mp_addir%mpr(0) = mpwds6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpadd (r1%mpr, rb%mpr, mp_addir%mpr, mpnw) 
    return
  end function

  function mp_addri (ra, ib)
    implicit none
    type (mp_real):: mp_addri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    double precision db
    type (mp_real) r1
    integer i1, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    mp_addri%mpr(0) = mpwds6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpadd (ra%mpr, r1%mpr, mp_addri%mpr, mpnw) 
    return
  end function

  function mp_addzz (za, zb)
    implicit none
    type (mp_complex):: mp_addzz
    type (mp_complex), intent (in):: za, zb
    integer l1, l2, l3, mpnw
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (zb%mpc(1)), &
      int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwds)      
    l3 = mpwds6
    mp_addzz%mpc(0) = mpwds6
    mp_addzz%mpc(l3) = mpwds6
    call mpcadd (za%mpc, zb%mpc, mp_addzz%mpc, mpnw) 
    return
  end function

  function mp_adddcz (dca, zb)
    implicit none
    type (mp_complex):: mp_adddcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complex), intent (in):: zb
    integer l1, l2, mpnw
    type (mp_complex) z1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)      
    l2 = mpwds6
    z1%mpc(0) = mpwds6
    z1%mpc(l2) = mpwds6
    mp_adddcz%mpc(0) = mpwds6
    mp_adddcz%mpc(l2) = mpwds6
    call mpdmc40 (dble (dca), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dca), 0, z1%mpc(l2), mpnw)
    call mpcadd (z1%mpc, zb%mpc, mp_adddcz%mpc, mpnw) 
    return
  end function

  function mp_addzdc (za, dcb)
    implicit none
    type (mp_complex):: mp_addzdc
    type (mp_complex), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer l1, l2, mpnw
    type (mp_complex) z1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)      
    l2 = mpwds6
    z1%mpc(0) = mpwds6
    z1%mpc(l2) = mpwds6
    mp_addzdc%mpc(0) = mpwds6
    mp_addzdc%mpc(l2) = mpwds6
    call mpdmc40 (dble (dcb), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dcb), 0, z1%mpc(l2), mpnw)
    call mpcadd (za%mpc, z1%mpc, mp_addzdc%mpc, mpnw) 
    return
  end function

  function mp_addrz (ra, zb)
    implicit none
    type (mp_complex):: mp_addrz
    type (mp_real), intent (in):: ra
    type (mp_complex), intent (in):: zb
    integer l2, l3, mpnw
    l2 = zb%mpc(0)
    mpnw = max (int (ra%mpr(1)), int (zb%mpc(1)), int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwds)      
    l3 = mpwds6
    mp_addrz%mpc(0) = mpwds6
    mp_addrz%mpc(l3) = mpwds6
    call mpadd (ra%mpr, zb%mpc, mp_addrz%mpc, mpnw)
    call mpeq (zb%mpc(l2), mp_addrz%mpc(l3), mpnw)
    return
  end function

  function mp_addzr (za, rb)
    implicit none
    type (mp_complex):: mp_addzr
    type (mp_complex), intent (in):: za
    type (mp_real), intent (in):: rb
    integer l1, l3, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)      
    l3 = mpwds6
    mp_addzr%mpc(0) = mpwds6
    mp_addzr%mpc(l3) = mpwds6
    call mpadd (za%mpc, rb%mpr, mp_addzr%mpc, mpnw)
    call mpeq (za%mpc(l1), mp_addzr%mpc(l3), mpnw)
    return
  end function

!  Subtraction routines:

  function mp_subrr (ra, rb)
    implicit none
    type (mp_real):: mp_subrr
    type (mp_real), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    mp_subrr%mpr(0) = mpwds6
    call mpsub (ra%mpr, rb%mpr, mp_subrr%mpr, mpnw) 
    return
  end function

  function mp_subdr (da, rb)
    implicit none
    type (mp_real):: mp_subdr
    double precision, intent (in):: da
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    integer i1, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    mp_subdr%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpsub (r1%mpr, rb%mpr, mp_subdr%mpr, mpnw) 
    return
  end function

  function mp_subrd (ra, db)
    implicit none
    type (mp_real):: mp_subrd
    type (mp_real), intent (in):: ra
    double precision, intent (in):: db
    type (mp_real) r1
    integer i1, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    mp_subrd%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpsub (ra%mpr, r1%mpr, mp_subrd%mpr, mpnw) 
    return
  end function

  function mp_subir (ia, rb)
    implicit none
    type (mp_real):: mp_subir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    double precision da
    integer i1, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    mp_subir%mpr(0) = mpwds6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpsub (r1%mpr, rb%mpr, mp_subir%mpr, mpnw) 
    return
  end function

  function mp_subri (ra, ib)
    implicit none
    type (mp_real):: mp_subri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    type (mp_real) r1
    double precision db
    integer i1, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    mp_subri%mpr(0) = mpwds6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpsub (ra%mpr, r1%mpr, mp_subri%mpr, mpnw) 
    return
  end function

  function mp_subzz (za, zb)
    implicit none
    type (mp_complex):: mp_subzz
    type (mp_complex), intent (in):: za, zb
    integer l1, l2, l3, mpnw
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (zb%mpc(1)), &
      int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwds)
    l3 = mpwds6
    mp_subzz%mpc(0) = mpwds6
    mp_subzz%mpc(l3) = mpwds6
    call mpcsub (za%mpc, zb%mpc, mp_subzz%mpc, mpnw) 
    return
  end function

  function mp_subdcz (dca, zb)
    implicit none
    type (mp_complex):: mp_subdcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complex), intent (in):: zb
    integer l1, l2, mpnw
    type (mp_complex) z1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)      
    l2 = mpwds6
    z1%mpc(0) = mpwds6
    z1%mpc(l2) = mpwds6
    mp_subdcz%mpc(0) = mpwds6
    mp_subdcz%mpc(l2) = mpwds6
    call mpdmc40 (dble (dca), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dca), 0, z1%mpc(l2), mpnw)
    call mpcsub (z1%mpc, zb%mpc, mp_subdcz%mpc, mpnw) 
    return
  end function

  function mp_subzdc (za, dcb)
    implicit none
    type (mp_complex):: mp_subzdc
    type (mp_complex), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer l1, l2, mpnw
    type (mp_complex) z1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)      
    l2 = mpwds6
    z1%mpc(0) = mpwds6
    z1%mpc(l2) = mpwds6
    mp_subzdc%mpc(0) = mpwds6
    mp_subzdc%mpc(l2) = mpwds6
    call mpdmc40 (dble (dcb), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dcb), 0, z1%mpc(l2), mpnw)
    call mpcsub (za%mpc, z1%mpc, mp_subzdc%mpc, mpnw) 
    return
  end function

  function mp_subrz (ra, zb)
    implicit none
    type (mp_complex):: mp_subrz
    type (mp_real), intent (in):: ra
    type (mp_complex), intent (in):: zb
    integer l2, l3, mpnw
    l2 = zb%mpc(0)
    mpnw = max (int (ra%mpr(1)), int (zb%mpc(1)), int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwds)      
    l3 = mpwds6
    mp_subrz%mpc(0) = mpwds6
    mp_subrz%mpc(l3) = mpwds6
    call mpsub (ra%mpr, zb%mpc, mp_subrz%mpc, mpnw)
    call mpeq (zb%mpc(l2), mp_subrz%mpc(l3), mpnw)
    mp_subrz%mpc(l3+2) = - mp_subrz%mpc(l3+2)
    return
  end function

  function mp_subzr (za, rb)
    implicit none
    type (mp_complex):: mp_subzr
    type (mp_complex), intent (in):: za
    type (mp_real), intent (in):: rb
    integer l1, l3, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)      
    l3 = mpwds6
    mp_subzr%mpc(0) = mpwds6
    mp_subzr%mpc(l3) = mpwds6
    call mpsub (za%mpc, rb%mpr, mp_subzr%mpc, mpnw)
    call mpeq (za%mpc(l1), mp_subzr%mpc(l3), mpnw)
    return
  end function

!  Negation routines:

  function mp_negr (ra)
    implicit none
    type (mp_real):: mp_negr
    type (mp_real), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_negr%mpr(0) = mpwds6
    call mpeq (ra%mpr, mp_negr%mpr, mpnw) 
    mp_negr%mpr(2) = - ra%mpr(2)
    return
  end function

  function mp_negz (za)
    implicit none
    type (mp_complex):: mp_negz
    type (mp_complex), intent (in):: za
    integer l1, l2, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    l2 = mpwds6
    mp_negz%mpc(0) = mpwds6
    mp_negz%mpc(l2) = mpwds6
    call mpceq (za%mpc, mp_negz%mpc, mpnw)
    mp_negz%mpc(2) = - mp_negz%mpc(2)
    mp_negz%mpc(l2+2) = - mp_negz%mpc(l2+2)
    return
  end function

!  Multiplication routines:

  function mp_mulrr (ra, rb)
    implicit none
    type (mp_real):: mp_mulrr
    type (mp_real), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    mp_mulrr%mpr(0) = mpwds6
    call mpmul (ra%mpr, rb%mpr, mp_mulrr%mpr, mpnw) 
    return
  end function

  function mp_muldr (da, rb)
    implicit none
    type (mp_real):: mp_muldr
    double precision, intent (in):: da
    type (mp_real), intent (in):: rb
    integer mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    mp_muldr%mpr(0) = mpwds6
    call mpmuld40 (rb%mpr, da, mp_muldr%mpr, mpnw) 
    return
  end function

  function mp_mulrd (ra, db)
    implicit none
    type (mp_real):: mp_mulrd
    type (mp_real), intent (in):: ra
    double precision, intent (in):: db
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_mulrd%mpr(0) = mpwds6
    call mpmuld40 (ra%mpr, db, mp_mulrd%mpr, mpnw) 
    return
  end function

  function mp_mulir (ia, rb)
    implicit none
    type (mp_real):: mp_mulir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    double precision da
    integer mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    mp_mulir%mpr(0) = mpwds6
    da = ia
    call mpmuld40 (rb%mpr, da, mp_mulir%mpr, mpnw) 
    return
  end function

  function mp_mulri (ra, ib)
    implicit none
    type (mp_real):: mp_mulri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    double precision db
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_mulri%mpr(0) = mpwds6
    db = ib
    call mpmuld40 (ra%mpr, db, mp_mulri%mpr, mpnw) 
    return
  end function

  function mp_mulzz (za, zb)
    implicit none
    type (mp_complex):: mp_mulzz
    type (mp_complex), intent (in):: za, zb
    integer l1, l2, l3, mpnw
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (zb%mpc(1)), &
      int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwds)
    l3 = mpwds6
    mp_mulzz%mpc(0) = mpwds6
    mp_mulzz%mpc(l3) = mpwds6
    call mpcmul (za%mpc, zb%mpc, mp_mulzz%mpc, mpnw) 
    return
  end function

  function mp_muldcz (dca, zb)
    implicit none
    type (mp_complex):: mp_muldcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complex), intent (in):: zb
    integer l1, l2, mpnw
    type (mp_complex) z1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)      
    l2 = mpwds6
    z1%mpc(0) = mpwds6
    z1%mpc(l2) = mpwds6
    mp_muldcz%mpc(0) = mpwds6
    mp_muldcz%mpc(l2) = mpwds6
    call mpdmc40 (dble (dca), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dca), 0, z1%mpc(l2), mpnw)
    call mpcmul (z1%mpc, zb%mpc, mp_muldcz%mpc, mpnw) 
    return
  end function

  function mp_mulzdc (za, dcb)
    implicit none
    type (mp_complex):: mp_mulzdc
    type (mp_complex), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer l1, l2, mpnw
    type (mp_complex) z1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)      
    l2 = mpwds6
    z1%mpc(0) = mpwds6
    z1%mpc(l2) = mpwds6
    mp_mulzdc%mpc(0) = mpwds6
    mp_mulzdc%mpc(l2) = mpwds6
    call mpdmc40 (dble (dcb), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dcb), 0, z1%mpc(l2), mpnw)
    call mpcmul (za%mpc, z1%mpc, mp_mulzdc%mpc, mpnw) 
    return
  end function

  function mp_mulrz (ra, zb)
    implicit none
    type (mp_complex):: mp_mulrz
    type (mp_real), intent (in):: ra
    type (mp_complex), intent (in):: zb
    integer l2, l3, mpnw
    l2 = zb%mpc(0)
    mpnw = max (int (ra%mpr(1)), int (zb%mpc(1)), int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwds)      
    l3 = mpwds6
    mp_mulrz%mpc(0) = mpwds6
    mp_mulrz%mpc(l3) = mpwds6
    call mpmul (ra%mpr, zb%mpc, mp_mulrz%mpc, mpnw)
    call mpmul (ra%mpr, zb%mpc(l2), mp_mulrz%mpc(l3), mpnw)
    return
  end function

  function mp_mulzr (za, rb)
    implicit none
    type (mp_complex):: mp_mulzr
    type (mp_complex), intent (in):: za
    type (mp_real), intent (in):: rb
    integer l1, l3, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)      
    l3 = mpwds6
    mp_mulzr%mpc(0) = mpwds6
    mp_mulzr%mpc(l3) = mpwds6
    call mpmul (za%mpc, rb%mpr, mp_mulzr%mpc, mpnw)
    call mpmul (za%mpc(l1), rb%mpr, mp_mulzr%mpc(l3), mpnw)
    return
  end function

  function mp_muldz (da, zb)
    implicit none
    type (mp_complex):: mp_muldz
    double precision, intent (in):: da
    type (mp_complex), intent (in):: zb
    integer l2, l3, mpnw
    l2 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwds)      
    l3 = mpwds6
    mp_muldz%mpc(0) = mpwds6
    mp_muldz%mpc(l3) = mpwds6
    call mpmuld (zb%mpc, da, mp_muldz%mpc, mpnw)
    call mpmuld (zb%mpc(l2), da, mp_muldz%mpc(l3), mpnw)
    return
  end function

  function mp_mulzd (za, db)
    implicit none
    type (mp_complex):: mp_mulzd
    type (mp_complex), intent (in):: za
    double precision, intent (in):: db
    integer l1, l3, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)      
    l3 = mpwds6
    mp_mulzd%mpc(0) = mpwds6
    mp_mulzd%mpc(l3) = mpwds6
    call mpmuld (za%mpc, db, mp_mulzd%mpc, mpnw)
    call mpmuld (za%mpc(l1), db, mp_mulzd%mpc(l3), mpnw)
    return
  end function

!  Division routines:

  function mp_divrr (ra, rb)
    implicit none
    type (mp_real):: mp_divrr
    type (mp_real), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    mp_divrr%mpr(0) = mpwds6
    call mpdiv (ra%mpr, rb%mpr, mp_divrr%mpr, mpnw) 
    return
  end function

  function mp_divdr (da, rb)
    implicit none
    type (mp_real):: mp_divdr
    double precision, intent (in):: da
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    integer i1, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    mp_divdr%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpdiv (r1%mpr, rb%mpr, mp_divdr%mpr, mpnw) 
    return
  end function

  function mp_divrd (ra, db)
    implicit none
    type (mp_real):: mp_divrd
    type (mp_real), intent (in):: ra
    double precision, intent (in):: db
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_divrd%mpr(0) = mpwds6
    call mpdivd40 (ra%mpr, db, mp_divrd%mpr, mpnw) 
    return
  end function

  function mp_divir (ia, rb)
    implicit none
    type (mp_real):: mp_divir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    double precision da
    integer i1, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    mp_divir%mpr(0) = mpwds6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpdiv (r1%mpr, rb%mpr, mp_divir%mpr, mpnw) 
    return
  end function

  function mp_divri (ra, ib)
    implicit none
    type (mp_real):: mp_divri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    double precision db
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_divri%mpr(0) = mpwds6
    db = ib
    call mpdivd40 (ra%mpr, db, mp_divri%mpr, mpnw) 
    return
  end function

  function mp_divzz (za, zb)
    implicit none
    type (mp_complex):: mp_divzz
    type (mp_complex), intent (in):: za, zb
    integer l1, l2, l3, mpnw
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (zb%mpc(1)), &
      int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwds)
    l3 = mpwds6
    mp_divzz%mpc(0) = mpwds6
    mp_divzz%mpc(l3) = mpwds6
    call mpcdiv (za%mpc, zb%mpc, mp_divzz%mpc, mpnw) 
    return
  end function

  function mp_divdcz (dca, zb)
    implicit none
    type (mp_complex):: mp_divdcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complex), intent (in):: zb
    integer l1, l2, mpnw
    type (mp_complex) z1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)      
    l2 = mpwds6
    z1%mpc(0) = mpwds6
    z1%mpc(l2) = mpwds6
    mp_divdcz%mpc(0) = mpwds6
    mp_divdcz%mpc(l2) = mpwds6
    call mpdmc40 (dble (dca), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dca), 0, z1%mpc(l2), mpnw)
    call mpcdiv (z1%mpc, zb%mpc, mp_divdcz%mpc, mpnw) 
    return
  end function

  function mp_divzdc (za, dcb)
    implicit none
    type (mp_complex):: mp_divzdc
    type (mp_complex), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer l1, l2, mpnw
    type (mp_complex) z1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)      
    l2 = mpwds6
    z1%mpc(0) = mpwds6
    z1%mpc(l2) = mpwds6
    mp_divzdc%mpc(0) = mpwds6
    mp_divzdc%mpc(l2) = mpwds6
    call mpdmc40 (dble (dcb), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dcb), 0, z1%mpc(l2), mpnw)
    call mpcdiv (za%mpc, z1%mpc, mp_divzdc%mpc, mpnw) 
    return
  end function

  function mp_divrz (ra, zb)
    implicit none
    type (mp_complex):: mp_divrz
    type (mp_real), intent (in):: ra
    type (mp_complex), intent (in):: zb
    type (mp_real):: r1, r2, r3
    integer l2, l3, mpnw
    l2 = zb%mpc(0)
    mpnw = max (int (ra%mpr(1)), int (zb%mpc(1)), int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwds)
    r1%mpr(0) = mpwds6
    r2%mpr(0) = mpwds6
    r3%mpr(0) = mpwds6     
    l3 = mpwds6
    mp_divrz%mpc(0) = mpwds6
    mp_divrz%mpc(l3) = mpwds6
    call mpmul (zb%mpc, zb%mpc, r1%mpr, mpnw)
    call mpmul (zb%mpc(l2), zb%mpc(l2), r2%mpr, mpnw)
    call mpadd (r1%mpr, r2%mpr, r3%mpr, mpnw)
    call mpmul (ra%mpr, zb%mpc, r1%mpr, mpnw)
    call mpdiv (r1%mpr, r3%mpr, mp_divrz%mpc, mpnw)
    call mpmul (ra%mpr, zb%mpc(l2), r1%mpr, mpnw)
    call mpdiv (r1%mpr, r3%mpr, mp_divrz%mpc(l3), mpnw)
    mp_divrz%mpc(l3+2) = - mp_divrz%mpc(l3+2)
    return
  end function

  function mp_divzr (za, rb)
    implicit none
    type (mp_complex):: mp_divzr
    type (mp_complex), intent (in):: za
    type (mp_real), intent (in):: rb
    integer l1, l3, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)      
    l3 = mpwds6
    mp_divzr%mpc(0) = mpwds6
    mp_divzr%mpc(l3) = mpwds6
    call mpdiv (za%mpc, rb%mpr, mp_divzr%mpc, mpnw)
    call mpdiv (za%mpc(l1), rb%mpr, mp_divzr%mpc(l3), mpnw)
    return
  end function

!  Exponentiation routines:

  function mp_expri (ra, ib)
    implicit none
    type (mp_real):: mp_expri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_expri%mpr(0) = mpwds6
    call mpnpwr (ra%mpr, ib, mp_expri%mpr, mpnw) 
    return
  end function

  function mp_exprr (ra, rb)
    implicit none
    type (mp_real):: mp_exprr
    type (mp_real), intent (in):: ra, rb
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_exprr%mpr(0) = mpwds6
    call mppower (ra%mpr, rb%mpr, mp_exprr%mpr, mpnw) 
    return
  end function

  function mp_expzi (za, ib)
    implicit none
    type (mp_complex):: mp_expzi
    type (mp_complex), intent (in):: za
    integer, intent (in):: ib
    integer l1, l2, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    l2 = mpwds6
    mp_expzi%mpc(0) = mpwds6
    mp_expzi%mpc(l2) = mpwds6
    call mpcnpwr (za%mpc, ib, mp_expzi%mpc, mpnw) 
    return
  end function

  function mp_expzz (za, zb)
    implicit none
    type (mp_complex):: mp_expzz
    type (mp_complex), intent (in):: za, zb
    integer l1, l2, l3, mpnw
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (zb%mpc(1)), &
      int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwds)  
    l3 = mpwds6
    mp_expzz%mpc(0) = mpwds6
    mp_expzz%mpc(l3) = mpwds6
    call mpcpowcc (za%mpc, zb%mpc, mp_expzz%mpc, mpnw)
    return
  end function

  function mp_exprz (ra, zb)
    implicit none
    type (mp_complex):: mp_exprz
    type (mp_real), intent (in):: ra
    type (mp_complex), intent (in):: zb
    integer l2, l3, mpnw
    l2 = zb%mpc(0)
    mpnw = max (int (ra%mpr(1)), int (zb%mpc(1)), int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwds)  
    l3 = mpwds6
    mp_exprz%mpc(0) = mpwds6
    mp_exprz%mpc(l3) = mpwds6
    call mpcpowrc (ra%mpr, zb%mpc, mp_exprz%mpc, mpnw)
    return
  end function

  function mp_expzr (za, rb)
    implicit none
    type (mp_complex):: mp_expzr
    type (mp_complex), intent (in):: za
    type (mp_real), intent (in):: rb
    type (mp_complex) z1, z2
    integer l1, l3, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)      
    l3 = mpwds6
    mp_expzr%mpc(0) = mpwds6
    mp_expzr%mpc(l3) = mpwds6
    call mpcpowcr (za%mpc, rb%mpr, mp_expzr%mpc, mpnw)
    return
  end function

!  Equality test routines:

  function mp_eqtrr (ra, rb)
    implicit none
    logical mp_eqtrr
    type (mp_real), intent (in):: ra, rb
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw) 
    if (ic == 0) then
      mp_eqtrr = .true.
    else
      mp_eqtrr = .false.
    endif
    return
  end function

  function mp_eqtdr (da, rb)
    implicit none
    logical mp_eqtdr
    double precision, intent (in):: da
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw) 
    if (ic == 0) then
      mp_eqtdr = .true.
    else
      mp_eqtdr = .false.
    endif
    return
  end function

  function mp_eqtrd (ra, db)
    implicit none
    logical mp_eqtrd
    type (mp_real), intent (in):: ra
    double precision, intent (in):: db
    type (mp_real) r1
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw) 
    if (ic == 0) then
      mp_eqtrd = .true.
    else
      mp_eqtrd = .false.
    endif
    return
  end function

  function mp_eqtir (ia, rb)
    implicit none
    logical mp_eqtir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    double precision da
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw) 
    if (ic == 0) then
      mp_eqtir = .true.
    else
      mp_eqtir = .false.
    endif
    return
  end function

  function mp_eqtri (ra, ib)
    implicit none
    logical mp_eqtri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    type (mp_real) r1
    double precision db
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw) 
    if (ic == 0) then
      mp_eqtri = .true.
    else
      mp_eqtri = .false.
    endif
    return
  end function

  function mp_eqtzz (za, zb)
    implicit none
    logical mp_eqtzz
    type (mp_complex), intent (in):: za, zb
    integer ic1, ic2, l1, l2, mpnw
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (zb%mpc(1)), &
      int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwds)
    call mpcpr (za%mpc, zb%mpc, ic1, mpnw)
    call mpcpr (za%mpc(l1), zb%mpc(l2), ic2, mpnw)
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtzz = .true.
    else
      mp_eqtzz = .false.
    endif 
    return
  end function

  function mp_eqtdcz (dca, zb)
    implicit none
    logical mp_eqtdcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complex), intent (in):: zb
    integer ic1, ic2, l1, l2, mpnw
    type (mp_complex) z1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    l2 = mpwds6
    z1%mpc(0) = mpwds6
    z1%mpc(l2) = mpwds6
    call mpdmc40 (dble (dca), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dca), 0, z1%mpc(l2), mpnw)
    call mpcpr (z1%mpc, zb%mpc, ic1, mpnw)
    call mpcpr (z1%mpc(l2), zb%mpc(l1), ic2, mpnw)
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtdcz = .true.
    else
      mp_eqtdcz = .false.
    endif 
    return
  end function

  function mp_eqtzdc (za, dcb)
    implicit none
    logical mp_eqtzdc
    type (mp_complex), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer ic1, ic2, l1, l2, mpnw
    type (mp_complex) z1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    l2 = mpwds6
    z1%mpc(0) = mpwds6
    z1%mpc(l2) = mpwds6
    call mpdmc40 (dble (dcb), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dcb), 0, z1%mpc(l2), mpnw)
    call mpcpr (za%mpc, z1%mpc, ic1, mpnw)
    call mpcpr (za%mpc(l1), z1%mpc(l2), ic2, mpnw)
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtzdc = .true.
    else
      mp_eqtzdc = .false.
    endif 
    return
  end function

  function mp_eqtrz (ra, zb)
    implicit none
    logical mp_eqtrz
    type (mp_real), intent (in):: ra
    type (mp_complex), intent (in):: zb
    integer ic1, ic2, l2, mpnw
    l2 = zb%mpc(0)
    mpnw = max (int (ra%mpr(1)), int (zb%mpc(1)), int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwds)
    call mpcpr (ra%mpr, zb%mpc, ic1, mpnw)
    ic2 = int (zb%mpc(l2+2))
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtrz = .true.
    else
      mp_eqtrz = .false.
    endif 
    return
  end function

  function mp_eqtzr (za, rb)
    implicit none
    logical mp_eqtzr
    type (mp_complex), intent (in):: za
    type (mp_real), intent (in):: rb
    integer ic1, ic2, l1, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    call mpcpr (za%mpc, rb%mpr, ic1, mpnw)
    ic2 = int (za%mpc(l1+2))
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtzr = .true.
    else
      mp_eqtzr = .false.
    endif 
    return
  end function

!  Non-equality test routines:

  function mp_netrr (ra, rb)
    implicit none
    logical mp_netrr
    type (mp_real), intent (in):: ra, rb
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw) 
    if (ic /= 0) then
      mp_netrr = .true.
    else
      mp_netrr = .false.
    endif
    return
  end function

  function mp_netdr (da, rb)
    implicit none
    logical mp_netdr
    double precision, intent (in):: da
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw) 
    if (ic /= 0) then
      mp_netdr = .true.
    else
      mp_netdr = .false.
    endif
    return
  end function

  function mp_netrd (ra, db)
    implicit none
    logical mp_netrd
    type (mp_real), intent (in):: ra
    double precision, intent (in):: db
    type (mp_real) r1
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw) 
    if (ic /= 0) then
      mp_netrd = .true.
    else
      mp_netrd = .false.
    endif
    return
  end function

  function mp_netir (ia, rb)
    implicit none
    logical mp_netir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    double precision da
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw) 
    if (ic /= 0) then
      mp_netir = .true.
    else
      mp_netir = .false.
    endif
    return
  end function

  function mp_netri (ra, ib)
    implicit none
    logical mp_netri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    type (mp_real) r1
    double precision db
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw) 
    if (ic /= 0) then
      mp_netri = .true.
    else
      mp_netri = .false.
    endif
    return
  end function

  function mp_netzz (za, zb)
    implicit none
    logical mp_netzz
    type (mp_complex), intent (in):: za, zb
    integer ic1, ic2, l1, l2, mpnw
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (zb%mpc(1)), &
      int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwds)
    call mpcpr (za%mpc, zb%mpc, ic1, mpnw)
    call mpcpr (za%mpc(l1), zb%mpc(l2), ic2, mpnw)
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netzz = .true.
    else
      mp_netzz = .false.
    endif 
    return
  end function

  function mp_netdcz (dca, zb)
    implicit none
    logical mp_netdcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complex), intent (in):: zb
    integer ic1, ic2, l1, l2, mpnw
    type (mp_complex) z1
    l1 = zb%mpc(0)
    mpnw = max (int (zb%mpc(1)), int (zb%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    l2 = mpwds6
    z1%mpc(0) = mpwds6
    z1%mpc(l2) = mpwds6
    call mpdmc40 (dble (dca), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dca), 0, z1%mpc(l2), mpnw)
    call mpcpr (z1%mpc, zb%mpc, ic1, mpnw)
    call mpcpr (z1%mpc(l2), zb%mpc(l1), ic2, mpnw)
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netdcz = .true.
    else
      mp_netdcz = .false.
    endif 
    return
  end function

  function mp_netzdc (za, dcb)
    implicit none
    logical mp_netzdc
    type (mp_complex), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer ic1, ic2, l1, l2, mpnw
    type (mp_complex) z1
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    l2 = mpwds6
    z1%mpc(0) = mpwds6
    z1%mpc(l2) = mpwds6
    call mpdmc40 (dble (dcb), 0, z1%mpc, mpnw)
    call mpdmc40 (aimag (dcb), 0, z1%mpc(l2), mpnw)
    call mpcpr (za%mpc, z1%mpc, ic1, mpnw)
    call mpcpr (za%mpc(l1), z1%mpc(l2), ic2, mpnw)
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netzdc = .true.
    else
      mp_netzdc = .false.
    endif 
    return
  end function

  function mp_netrz (ra, zb)
    implicit none
    logical mp_netrz
    type (mp_real), intent (in):: ra
    type (mp_complex), intent (in):: zb
    integer ic1, ic2, l2, mpnw
    l2 = zb%mpc(0)
    mpnw = max (int (ra%mpr(1)), int (zb%mpc(1)), int (zb%mpc(l2+1)))
    mpnw = min (mpnw, mpwds)
    call mpcpr (ra%mpr, zb%mpc, ic1, mpnw)
    ic2 = int (zb%mpc(l2+2))
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netrz = .true.
    else
      mp_netrz = .false.
    endif 
    return
  end function

  function mp_netzr (za, rb)
    implicit none
    logical mp_netzr
    type (mp_complex), intent (in):: za
    type (mp_real), intent (in):: rb
    integer ic1, ic2, l1, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    call mpcpr (za%mpc, rb%mpr, ic1, mpnw)
    ic2 = int (za%mpc(l1+2))
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netzr = .true.
    else
      mp_netzr = .false.
    endif 
    return
  end function

!  Less-than-or-equal test routines:

  function mp_letrr (ra, rb)
    implicit none
    logical mp_letrr
    type (mp_real), intent (in):: ra, rb
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw) 
    if (ic <= 0) then
      mp_letrr = .true.
    else
      mp_letrr = .false.
    endif
    return
  end function

  function mp_letdr (da, rb)
    implicit none
    logical mp_letdr
    double precision, intent (in):: da
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw) 
    if (ic <= 0) then
      mp_letdr = .true.
    else
      mp_letdr = .false.
    endif
    return
  end function

  function mp_letrd (ra, db)
    implicit none
    logical mp_letrd
    type (mp_real), intent (in):: ra
    double precision, intent (in):: db
    type (mp_real) r1
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw) 
    if (ic <= 0) then
      mp_letrd = .true.
    else
      mp_letrd = .false.
    endif
    return
  end function

  function mp_letir (ia, rb)
    implicit none
    logical mp_letir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    double precision da
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw) 
    if (ic <= 0) then
      mp_letir = .true.
    else
      mp_letir = .false.
    endif
    return
  end function

  function mp_letri (ra, ib)
    implicit none
    logical mp_letri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    type (mp_real) r1
    double precision db
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw) 
    if (ic <= 0) then
      mp_letri = .true.
    else
      mp_letri = .false.
    endif
    return
  end function

!  Greater-than-or-equal test routines:

  function mp_getrr (ra, rb)
    implicit none
    logical mp_getrr
    type (mp_real), intent (in):: ra, rb
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw) 
    if (ic .ge. 0) then
      mp_getrr = .true.
    else
      mp_getrr = .false.
    endif
    return
  end function

  function mp_getdr (da, rb)
    implicit none
    logical mp_getdr
    double precision, intent (in):: da
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw) 
    if (ic .ge. 0) then
      mp_getdr = .true.
    else
      mp_getdr = .false.
    endif
    return
  end function

  function mp_getrd (ra, db)
    implicit none
    logical mp_getrd
    type (mp_real), intent (in):: ra
    double precision, intent (in):: db
    type (mp_real) r1
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw) 
    if (ic .ge. 0) then
      mp_getrd = .true.
    else
      mp_getrd = .false.
    endif
    return
  end function

  function mp_getir (ia, rb)
    implicit none
    logical mp_getir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    double precision da
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw) 
    if (ic .ge. 0) then
      mp_getir = .true.
    else
      mp_getir = .false.
    endif
    return
  end function

  function mp_getri (ra, ib)
    implicit none
    logical mp_getri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    type (mp_real) r1
    double precision db
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw) 
    if (ic .ge. 0) then
      mp_getri = .true.
    else
      mp_getri = .false.
    endif
    return
  end function

!  Less-than test routines:

  function mp_lttrr (ra, rb)
    implicit none
    logical mp_lttrr
    type (mp_real), intent (in):: ra, rb
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw) 
    if (ic < 0) then
      mp_lttrr = .true.
    else
      mp_lttrr = .false.
    endif
    return
  end function

  function mp_lttdr (da, rb)
    implicit none
    logical mp_lttdr
    double precision, intent (in):: da
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw) 
    if (ic < 0) then
      mp_lttdr = .true.
    else
      mp_lttdr = .false.
    endif
    return
  end function

  function mp_lttrd (ra, db)
    implicit none
    logical mp_lttrd
    type (mp_real), intent (in):: ra
    double precision, intent (in):: db
    type (mp_real) r1
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw) 
    if (ic < 0) then
      mp_lttrd = .true.
    else
      mp_lttrd = .false.
    endif
    return
  end function

  function mp_lttir (ia, rb)
    implicit none
    logical mp_lttir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    double precision da
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw) 
    if (ic < 0) then
      mp_lttir = .true.
    else
      mp_lttir = .false.
    endif
    return
  end function

  function mp_lttri (ra, ib)
    implicit none
    logical mp_lttri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    type (mp_real) r1
    double precision db
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw) 
    if (ic < 0) then
      mp_lttri = .true.
    else
      mp_lttri = .false.
    endif
    return
  end function

!  Greater-than test routines:

  function mp_gttrr (ra, rb)
    implicit none
    logical mp_gttrr
    type (mp_real), intent (in):: ra, rb
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw) 
    if (ic > 0) then
      mp_gttrr = .true.
    else
      mp_gttrr = .false.
    endif
    return
  end function

  function mp_gttdr (da, rb)
    implicit none
    logical mp_gttdr
    double precision, intent (in):: da
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw) 
    if (ic > 0) then
      mp_gttdr = .true.
    else
      mp_gttdr = .false.
    endif
    return
  end function

  function mp_gttrd (ra, db)
    implicit none
    logical mp_gttrd
    type (mp_real), intent (in):: ra
    double precision, intent (in):: db
    type (mp_real) r1
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw) 
    if (ic > 0) then
      mp_gttrd = .true.
    else
      mp_gttrd = .false.
    endif
    return
  end function
  
  function mp_gttir (ia, rb)
    implicit none
    logical mp_gttir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    double precision da
    integer i1, ic, mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    da = ia
    i1 = 0
    call mpdmc40 (da, i1, r1%mpr, mpnw)
    call mpcpr (r1%mpr, rb%mpr, ic, mpnw) 
    if (ic > 0) then
      mp_gttir = .true.
    else
      mp_gttir = .false.
    endif
    return
  end function

  function mp_gttri (ra, ib)
    implicit none
    logical mp_gttri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    type (mp_real) r1
    double precision db
    integer i1, ic, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    db = ib
    i1 = 0
    call mpdmc40 (db, i1, r1%mpr, mpnw)
    call mpcpr (ra%mpr, r1%mpr, ic, mpnw) 
    if (ic > 0) then
      mp_gttri = .true.
    else
      mp_gttri = .false.
    endif
    return
  end function
  
!   Algebraic and transcendental function definitions, listed alphabetically:

  function mp_absr (ra)
    implicit none
    type (mp_real):: mp_absr
    type (mp_real), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_absr%mpr(0) = mpwds6
    call mpeq (ra%mpr, mp_absr%mpr, mpnw)
    mp_absr%mpr(2) = abs (ra%mpr(2))
    return
  end function
  
  function mp_absz (za)
    implicit none
    type (mp_real):: mp_absz
    type (mp_complex), intent (in):: za
    integer l1, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    mp_absz%mpr(0) = mpwds6
    call mpcabs (za%mpc, mp_absz%mpr, mpnw)
    return
  end function

  function mp_acos (ra)
    implicit none
    type (mp_real):: mp_acos
    type (mp_real), intent (in):: ra
    type (mp_real) r1, r2, r3
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    r2%mpr(0) = mpwds6
    r3%mpr(0) = mpwds6
    call mpmul (ra%mpr, ra%mpr, r1%mpr, mpnw)
    call mpdmc (1.d0, 0, r2%mpr, mpnw)
    call mpsub (r2%mpr, r1%mpr, r3%mpr, mpnw)
    if (r3%mpr(2) < 0.d0) then
      write (mpldb, 1)
1     format ('*** MP_ACOS: argument is not in (-1, 1).')
      call mpabrt (24)
    endif
    call mpsqrt (r3%mpr, r1%mpr, mpnw)
    mp_acos%mpr(0) = mpwds6
    call mpang (ra%mpr, r1%mpr, mp_acos%mpr, mpnw)
    return
  end function
  
  function mp_acosh (ra)
    implicit none
    type (mp_real):: mp_acosh
    type (mp_real), intent (in):: ra
    type (mp_real) r1, r2, r3
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    r2%mpr(0) = mpwds6
    r3%mpr(0) = mpwds6
    mp_acosh%mpr(0) = mpwds6
    call mpmul (ra%mpr, ra%mpr, r1%mpr, mpnw)
    call mpdmc (1.d0, 0, r2%mpr, mpnw)
    call mpsub (r1%mpr, r2%mpr, r3%mpr, mpnw)
    if (r3%mpr(2) < 0.d0) then
      write (mpldb, 1)
1     format ('*** MP_ACOSH: argument is not >= 1.')
      call mpabrt (24)
    endif
    call mpsqrt (r3%mpr, r1%mpr, mpnw)
    call mpadd (ra%mpr, r1%mpr, r2%mpr, mpnw)
    call mplog (r2%mpr, mp_acosh%mpr, mpnw)
    return
  end function
  
   function mp_agm (ra, rb)
    implicit none
    type (mp_real):: mp_agm
    type (mp_real), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    mp_agm%mpr(0) = mpwds6
    call mpagmr (rb%mpr, ra%mpr, mp_agm%mpr, mpnw)
    return
  end function

  function mp_aimag (za)
    implicit none
    type (mp_real):: mp_aimag
    type (mp_complex), intent (in):: za
    integer l1, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    mp_aimag%mpr(0) = mpwds6
    call mpeq (za%mpc(l1), mp_aimag%mpr, mpnw)
    return
  end function

  function mp_aint (ra)
    implicit none
    type (mp_real):: mp_aint
    type (mp_real), intent (in):: ra
    type (mp_real) r1
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    mp_aint%mpr(0) = mpwds6
    call mpinfr (ra%mpr, mp_aint%mpr, r1%mpr, mpnw)
    return
  end function

   function mp_anint (ra)
    implicit none
    type (mp_real):: mp_anint
    type (mp_real), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_anint%mpr(0) = mpwds6
    call mpnint (ra%mpr, mp_anint%mpr, mpnw)
    return
  end function

   function mp_asin (ra)
    implicit none
    type (mp_real):: mp_asin
    type (mp_real), intent (in):: ra
    type (mp_real) r1, r2, r3
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    r2%mpr(0) = mpwds6
    r3%mpr(0) = mpwds6
    call mpmul (ra%mpr, ra%mpr, r1%mpr, mpnw)
    call mpdmc (1.d0, 0, r2%mpr, mpnw)
    call mpsub (r2%mpr, r1%mpr, r3%mpr, mpnw)
    if (r3%mpr(2) < 0.d0) then
      write (mpldb, 1)
1     format ('*** MP_ASIN: argument is not in (-1, 1).')
      call mpabrt (25)
    endif
    call mpsqrt (r3%mpr, r1%mpr, mpnw)
    mp_asin%mpr(0) = mpwds6
    call mpang (r1%mpr, ra%mpr, mp_asin%mpr, mpnw)
    return
  end function

  function mp_asinh (ra)
    implicit none
    type (mp_real):: mp_asinh
    type (mp_real), intent (in):: ra
    type (mp_real) r1, r2, r3
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    r2%mpr(0) = mpwds6
    r3%mpr(0) = mpwds6
    mp_asinh%mpr(0) = mpwds6
    call mpmul (ra%mpr, ra%mpr, r1%mpr, mpnw)
    call mpdmc (1.d0, 0, r2%mpr, mpnw)
    call mpadd (r1%mpr, r2%mpr, r3%mpr, mpnw)
    call mpsqrt (r3%mpr, r1%mpr, mpnw)
    call mpadd (ra%mpr, r1%mpr, r2%mpr, mpnw)
    call mplog (r2%mpr, mp_asinh%mpr, mpnw)
    return
  end function
  
   function mp_atan (ra)
    implicit none
    type (mp_real):: mp_atan
    type (mp_real), intent (in):: ra
    type (mp_real) r1
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    call mpdmc (1.d0, 0, r1%mpr, mpnw)
    mp_atan%mpr(0) = mpwds6
    call mpang (r1%mpr, ra%mpr, mp_atan%mpr, mpnw)
    return
  end function

   function mp_atan2 (ra, rb)
    implicit none
    type (mp_real):: mp_atan2
    type (mp_real), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    mp_atan2%mpr(0) = mpwds6
    call mpang (rb%mpr, ra%mpr, mp_atan2%mpr, mpnw)
    return
  end function

  function mp_atanh (ra)
    implicit none
    type (mp_real):: mp_atanh
    type (mp_real), intent (in):: ra
    type (mp_real) r1, r2, r3
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    r2%mpr(0) = mpwds6
    r3%mpr(0) = mpwds6
    mp_atanh%mpr(0) = mpwds6
    call mpdmc (1.d0, 0, r1%mpr, mpnw)
    call mpadd (r1%mpr, ra%mpr, r2%mpr, mpnw)
    call mpsub (r1%mpr, ra%mpr, r3%mpr, mpnw)
    call mpdiv (r2%mpr, r3%mpr, r1%mpr, mpnw)
    if (r1%mpr(2) < 0.d0) then
      write (mpldb, 1)
1     format ('*** MP_ATANH: argument is not in (-1, 1).')
      call mpabrt (24)
    endif
    call mplog (r1%mpr, r2%mpr, mpnw)
    call mpmuld (r2%mpr, 0.5d0, mp_atanh%mpr, mpnw)
    return
  end function
  
  function mp_ator1 (a, ib, iprec)
    implicit none
    type (mp_real) mp_ator1
    integer, intent (in):: ib
    character*1, intent (in):: a(ib)
    integer mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    mp_ator1%mpr(0) = mpwds6
    call mpctomp (a, ib, mp_ator1%mpr, mpnw)
    return
  end function

  function mp_atorn (aa, iprec)
    implicit none
    character*(*), intent (in):: aa
    type (mp_real):: mp_atorn
    character*1 :: chr1(len(aa))
    integer i, l1, mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    l1 = len (aa)
    do i = 1, l1
      chr1(i) = aa(i:i)
    enddo
    mp_atorn%mpr(0) = mpwds6
    call mpctomp (chr1, l1, mp_atorn%mpr, mpnw) 
    return
  end function
  
  subroutine mp_berne (nb, rb, iprec)
    implicit none
    integer, intent (in):: nb
    type (mp_real), intent (out):: rb(nb)
    integer i, n1, mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    do i = 1, nb
      rb(i)%mpr(0) = mpwds6
    enddo
    n1 = mpwds
    call mpberner (n1, nb, rb(1)%mpr, mpnw)
    return
  end subroutine

  function mp_besselj (ra, rb)
    implicit none
    type (mp_real):: mp_besselj
    type (mp_real), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    mp_besselj%mpr(0) = mpwds6
    call mpbesseljr (ra%mpr, rb%mpr, mp_besselj%mpr, mpnw)
    return
  end function
  
  function mp_bessel_j0 (ra)
    implicit none
    type (mp_real):: mp_bessel_j0
    type (mp_real), intent (in):: ra
    type (mp_real) r1
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    call mpdmc (0.d0, 0, r1%mpr, mpnw)
    mp_bessel_j0%mpr(0) = mpwds6
    call mpbesseljr (r1%mpr, ra%mpr, mp_bessel_j0%mpr, mpnw)
    return
  end function
  
  function mp_bessel_j1 (ra)
    implicit none
    type (mp_real):: mp_bessel_j1
    type (mp_real), intent (in):: ra
    type (mp_real) r1
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    call mpdmc (1.d0, 0, r1%mpr, mpnw)
    mp_bessel_j1%mpr(0) = mpwds6
    call mpbesseljr (r1%mpr, ra%mpr, mp_bessel_j1%mpr, mpnw)
    return
  end function
    
  function mp_bessel_jn (ia, rb)
    implicit none
    type (mp_real):: mp_bessel_jn
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    integer mpnw
    mpnw = min (int (rb%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    call mpdmc (dble (ia), 0, r1%mpr, mpnw)
    mp_bessel_jn%mpr(0) = mpwds6
    call mpbesseljr (r1%mpr, rb%mpr, mp_bessel_jn%mpr, mpnw)
    return
  end function
  
  function mp_ccos (za)
    implicit none
    type (mp_complex):: mp_ccos
    type (mp_complex), intent (in):: za
    integer l1, l2, l3, mpnw
    type (mp_complex) z1, z2, z3
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    l2 = mpwds6
    z1%mpc(0) = mpwds6
    z1%mpc(l2) = mpwds6
    z2%mpc(0) = mpwds6
    z2%mpc(l2) = mpwds6
    z3%mpc(0) = mpwds6
    z3%mpc(l2) = mpwds6
    l3 = mpwds6
    mp_ccos%mpc(0) = mpwds6
    mp_ccos%mpc(l3) = mpwds6
    call mpdmc (1.d0, 0, z1%mpc, mpnw)
    call mpdmc (0.d0, 0, z1%mpc(l2), mpnw)
    call mpeq (za%mpc, z3%mpc(l2), mpnw)
    call mpeq (za%mpc(l1), z3%mpc, mpnw)
    z3%mpc(2) = - z3%mpc(2)
    call mpcexp (z3%mpc, z2%mpc, mpnw)
    call mpcdiv (z1%mpc, z2%mpc, z3%mpc, mpnw)
    call mpcadd (z2%mpc, z3%mpc, z1%mpc, mpnw)
    call mpmuld (z1%mpc, 0.5d0, mp_ccos%mpc, mpnw)
    call mpmuld (z1%mpc(l2), 0.5d0, mp_ccos%mpc(l3), mpnw)    
    return
  end function

  function mp_cexp (za)
    implicit none
    type (mp_complex):: mp_cexp
    type (mp_complex), intent (in):: za
    integer l1, l2, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    l2 = mpwds6
    mp_cexp%mpc(0) = mpwds6
    mp_cexp%mpc(l2) = mpwds6
    call mpcexp (za%mpc, mp_cexp%mpc, mpnw)
    return
  end function
  
  function mp_clog (za)
    implicit none
    type (mp_complex):: mp_clog
    type (mp_complex), intent (in):: za
    integer l1, l2, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    l2 = mpwds6
    mp_clog%mpc(0) = mpwds6
    mp_clog%mpc(l2) = mpwds6
    call mpclog (za%mpc, mp_clog%mpc, mpnw)
    return
  end function
  
  function mp_conjg (za)
    implicit none
    type (mp_complex):: mp_conjg
    type (mp_complex), intent (in):: za
    integer l1, l2, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    l2 = mpwds6
    mp_conjg%mpc(0) = mpwds6
    mp_conjg%mpc(l2) = mpwds6
    call mpconjg (za%mpc, mp_conjg%mpc, mpnw)
    return
  end function
  
  function mp_cos (ra)
    implicit none
    type (mp_real):: mp_cos
    type (mp_real), intent (in):: ra
    integer mpnw
    type (mp_real) r1
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_cos%mpr(0) = mpwds6
    r1%mpr(0) = mpwds6
    call mpcssnr (ra%mpr, mp_cos%mpr, r1%mpr, mpnw)
    return
  end function
  
  function mp_cosh (ra)
    implicit none
    type (mp_real):: mp_cosh
    type (mp_real), intent (in):: ra
    integer mpnw
    type (mp_real) r1
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_cosh%mpr(0) = mpwds6
    r1%mpr(0) = mpwds6
    call mpcsshr (ra%mpr, mp_cosh%mpr, r1%mpr, mpnw)
    return
  end function
   
  function mp_csin (za)
    implicit none
    type (mp_complex):: mp_csin
    type (mp_complex), intent (in):: za
    integer l1, l2, l3, mpnw
    type (mp_complex) z1, z2, z3
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    l2 = mpwds6
    z1%mpc(0) = mpwds6
    z1%mpc(l2) = mpwds6
    z2%mpc(0) = mpwds6
    z2%mpc(l2) = mpwds6
    z3%mpc(0) = mpwds6
    z3%mpc(l2) = mpwds6
    l3 = mpwds6
    mp_csin%mpc(0) = mpwds6
    mp_csin%mpc(l3) = mpwds6
    call mpdmc (1.d0, 0, z1%mpc, mpnw)
    call mpdmc (0.d0, 0, z1%mpc(l2), mpnw)
    call mpeq (za%mpc, z3%mpc(l2), mpnw)
    call mpeq (za%mpc(l1), z3%mpc, mpnw)
    z3%mpc(2) = - z3%mpc(2)
    call mpcexp (z3%mpc, z2%mpc, mpnw)
    call mpcdiv (z1%mpc, z2%mpc, z3%mpc, mpnw)
    call mpcsub (z2%mpc, z3%mpc, z1%mpc, mpnw)
    call mpmuld (z1%mpc, 0.5d0, mp_csin%mpc(l3), mpnw)
    call mpmuld (z1%mpc(l2), 0.5d0, mp_csin%mpc, mpnw)
    mp_csin%mpc(l3+2) = - mp_csin%mpc(l3+2)
    return
  end function

  function mp_csqrt (za)
    implicit none
    type (mp_complex):: mp_csqrt
    type (mp_complex), intent (in):: za
    integer l1, l2, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    l2 = mpwds6
    mp_csqrt%mpc(0) = mpwds6
    mp_csqrt%mpc(l2) = mpwds6
    call mpcsqrt (za%mpc, mp_csqrt%mpc, mpnw)
    return
  end function
  
  subroutine mp_cssh (ra, rb, rc)
    implicit none
    type (mp_real), intent (in):: ra
    type (mp_real), intent (out):: rb, rc
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    rb%mpr(0) = mpwds6
    rc%mpr(0) = mpwds6
    call mpcsshr (ra%mpr, rb%mpr, rc%mpr, mpnw)
    return
  end subroutine
  
  subroutine mp_cssn (ra, rb, rc)
    implicit none
    type (mp_real), intent (in):: ra
    type (mp_real), intent (out):: rb, rc
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    rb%mpr(0) = mpwds6
    rc%mpr(0) = mpwds6
    call mpcssnr (ra%mpr, rb%mpr, rc%mpr, mpnw)
    return
  end subroutine
  
  function mp_dctoz (dca, iprec)
    implicit none
    type (mp_complex):: mp_dctoz
    complex (kind(0.d0)), intent (in):: dca
    integer l1, mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    l1 = mpwds6
    mp_dctoz%mpc(0) = mpwds6
    mp_dctoz%mpc(l1) = mpwds6
    call mpdmc40 (dble (dca), 0, mp_dctoz%mpc, mpnw)
    call mpdmc40 (imag (dca), 0, mp_dctoz%mpc(l1), mpnw)
    return
  end function
  
  function mp_dctoz2 (dca, iprec)
    implicit none
    type (mp_complex):: mp_dctoz2
    complex (kind(0.d0)), intent (in):: dca
    integer l1, mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    l1 = mpwds6
    mp_dctoz2%mpc(0) = mpwds6
    mp_dctoz2%mpc(l1) = mpwds6
    call mpdmc (dble (dca), 0, mp_dctoz2%mpc, mpnw)
    call mpdmc (imag (dca), 0, mp_dctoz2%mpc(l1), mpnw)
    return
  end function
  
  function mp_dtor (da, iprec)
    implicit none
    type (mp_real):: mp_dtor
    double precision, intent (in):: da
    integer i1, mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    mp_dtor%mpr(0) = mpwds6
    i1 = 0
    call mpdmc40 (da, i1, mp_dtor%mpr, mpnw)
    return
  end function
 
  function mp_dtor2 (da, iprec)
    implicit none
    type (mp_real):: mp_dtor2
    double precision, intent (in):: da
    integer i1, mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    mp_dtor2%mpr(0) = mpwds6
    i1 = 0
    call mpdmc (da, i1, mp_dtor2%mpr, mpnw)
    return
  end function
 
  subroutine mp_eform (ra, nb, nd, b)
    implicit none
    type (mp_real), intent (in):: ra
    integer, intent (in):: nb, nd
    character*1, intent (out):: b(nb)
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    call mpeformat (ra%mpr, nb, nd, b, mpnw)
    return
  end subroutine
 
  function mp_exp (ra)
    implicit none
    type (mp_real):: mp_exp
    type (mp_real), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_exp%mpr(0) = mpwds6
    call mpexp (ra%mpr, mp_exp%mpr, mpnw)
    return
  end function

  subroutine mp_fform (ra, nb, nd, b)
    implicit none
    type (mp_real), intent (in):: ra
    integer, intent (in):: nb, nd
    character*1, intent (out):: b(nb)
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    call mpfformat (ra%mpr, nb, nd, b, mpnw)
    return
  end subroutine
 
  function mp_gamma (ra)
    implicit none
    type (mp_real):: mp_gamma
    type (mp_real), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_gamma%mpr(0) = mpwds6
    call mpgammar (ra%mpr, mp_gamma%mpr, mpnw)
    return
  end function
  
  function mp_hypot (ra, rb)
    implicit none
    type (mp_real):: mp_hypot
    type (mp_real), intent (in):: ra, rb
    type (mp_real) r1, r2, r3
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_hypot%mpr(0) = mpwds6
    r1%mpr(0) = mpwds6
    r2%mpr(0) = mpwds6
    r3%mpr(0) = mpwds6
    call mpmul (ra%mpr, ra%mpr, r1%mpr, mpnw)
    call mpmul (rb%mpr, rb%mpr, r2%mpr, mpnw)
    call mpadd (r1%mpr, r2%mpr, r3%mpr, mpnw)
    call mpsqrt (r3%mpr, mp_hypot%mpr, mpnw)
    return
  end function
  
  function mp_incgamma (ra, rb)
    implicit none
    type (mp_real):: mp_incgamma
    type (mp_real), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    mp_incgamma%mpr(0) = mpwds6
    call mpincgammar (ra%mpr, rb%mpr, mp_incgamma%mpr, mpnw)
    return
  end function
  
  subroutine mp_init (iprec)
    implicit none
    integer mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    call mpinifft (mpnw)
    call mpinitran (mpnw)
    return
  end subroutine
  
  function mp_log (ra)
    implicit none
    type (mp_real):: mp_log
    type (mp_real), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_log%mpr(0) = mpwds6
    call mplog (ra%mpr, mp_log%mpr, mpnw)
    return
  end function

  function mp_log2 (iprec)
    implicit none
    type (mp_real):: mp_log2
    type (mp_real) qpi
    integer mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    mp_log2%mpr(0) = mpwds6
    qpi%mpr(0) = mpwds6
    call mppiq (qpi%mpr, mpnw)
    call mplog2q (qpi%mpr, mp_log2%mpr, mpnw)
  end function

  function mp_max (ra, rb, rc)
    implicit none
    type (mp_real):: mp_max
    type (mp_real), intent (in):: ra, rb
    type (mp_real), intent (in), optional:: rc
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    if (present (rc)) mpnw = max (mpnw, int (rc%mpr(1)))
    mpnw = min (mpnw, mpwds)
    mp_max%mpr(0) = mpwds6
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw)
    if (ic >= 0) then
      call mpeq (ra%mpr, mp_max%mpr, mpnw)
    else
      call mpeq (rb%mpr, mp_max%mpr, mpnw)
    endif
    if (present (rc)) then
      call mpcpr (rc%mpr, mp_max%mpr, ic, mpnw)
      if (ic >= 0) call mpeq (rc%mpr, mp_max%mpr, mpnw)
    endif
    return
  end function
  
  subroutine mp_mdi (ra, db, ic)
    implicit none
    type (mp_real), intent (in):: ra
    double precision, intent (out):: db
    integer, intent (out):: ic
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    call mpmdc (ra%mpr, db, ic, mpnw)
    return
  end subroutine
    
  function mp_min (ra, rb, rc)
    implicit none
    type (mp_real):: mp_min
    type (mp_real), intent (in):: ra, rb
    type (mp_real), intent (in), optional:: rc
    integer ic, mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    if (present (rc)) mpnw = max (mpnw, int (rc%mpr(1)))
    mpnw = min (mpnw, mpwds)
    mp_min%mpr(0) = mpwds6
    call mpcpr (ra%mpr, rb%mpr, ic, mpnw)
    if (ic <= 0) then
      call mpeq (ra%mpr, mp_min%mpr, mpnw)
    else
      call mpeq (rb%mpr, mp_min%mpr, mpnw)
    endif
    if (present (rc)) then
      call mpcpr (rc%mpr, mp_min%mpr, ic, mpnw)
      if (ic <= 0) call mpeq (rc%mpr, mp_min%mpr, mpnw)
    endif
    return
  end function

  function mp_nrt (ra, ib)
    implicit none
    type (mp_real):: mp_nrt
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_nrt%mpr(0) = mpwds6
    call mpnrtr (ra%mpr, ib, mp_nrt%mpr, mpnw) 
    return
  end function
  
  function mp_pi (iprec)
    implicit none
    type (mp_real):: mp_pi
    integer mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    mp_pi%mpr(0) = mpwds6
    call mppiq (mp_pi%mpr, mpnw)
  end function

  function mp_prodd (ra, db)
    implicit none
    type (mp_real):: mp_prodd
    type (mp_real), intent (in):: ra
    double precision, intent (in):: db
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_prodd%mpr(0) = mpwds6
    call mpmuld (ra%mpr, db, mp_prodd%mpr, mpnw)
  end function
  
  function mp_qtor (qa, iprec)
    implicit none
    type (mp_real):: mp_qtor
    integer mpnw

!>  If real*16 is supported, uncomment this line:
    real (kind (0.q0)), intent (in):: qa
!>  Otherwise uncomment this line:
!    real (kind (0.d0)), intent (in):: qa

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    mp_qtor%mpr(0) = mpwds6
    call mpqmc (qa, 0, mp_qtor%mpr, mpnw)
    return
  end function
  
  function mp_quotd (ra, db)
    implicit none
    type (mp_real):: mp_quotd
    type (mp_real), intent (in):: ra
    double precision, intent (in):: db
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_quotd%mpr(0) = mpwds6
    call mpdivd (ra%mpr, db, mp_quotd%mpr, mpnw)
  end function
    
!   Five variations are necessary here due to Fortran rules about optional arguments.

  subroutine mp_readr1 (iu, r1, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_real), intent (out):: r1
    integer mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    call mpinp (iu, r1%mpr, mpnw)
    return
  end subroutine

  subroutine mp_readr2 (iu, r1, r2, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_real), intent (out):: r1, r2
    integer mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    call mpinp (iu, r1%mpr, mpnw)
    call mpinp (iu, r2%mpr, mpnw)
    return
  end subroutine

  subroutine mp_readr3 (iu, r1, r2, r3, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_real), intent (out):: r1, r2, r3
    integer mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    call mpinp (iu, r1%mpr, mpnw)
    call mpinp (iu, r2%mpr, mpnw)
    call mpinp (iu, r3%mpr, mpnw)
    return
  end subroutine

  subroutine mp_readr4 (iu, r1, r2, r3, r4, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_real), intent (out):: r1, r2, r3, r4
    integer mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    call mpinp (iu, r1%mpr, mpnw)
    call mpinp (iu, r2%mpr, mpnw)
    call mpinp (iu, r3%mpr, mpnw)
    call mpinp (iu, r4%mpr, mpnw)
    return
  end subroutine

  subroutine mp_readr5 (iu, r1, r2, r3, r4, r5, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_real), intent (out):: r1, r2, r3, r4, r5
    integer mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    call mpinp (iu, r1%mpr, mpnw)
    call mpinp (iu, r2%mpr, mpnw)
    call mpinp (iu, r3%mpr, mpnw)
    call mpinp (iu, r4%mpr, mpnw)
    call mpinp (iu, r5%mpr, mpnw)
    return
  end subroutine

  subroutine mp_readz1 (iu, z1, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_complex), intent (out):: z1
    integer l1, mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    l1 = z1%mpc(0)
    call mpinp (iu, z1%mpc, mpnw)
    call mpinp (iu, z1%mpc(l1), mpnw)
    return
  end subroutine

  subroutine mp_readz2 (iu, z1, z2, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_complex), intent (out):: z1, z2
    integer l1, l2, mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    l1 = z1%mpc(0)
    l2 = z2%mpc(0)
    call mpinp (iu, z1%mpc, mpnw)
    call mpinp (iu, z1%mpc(l1), mpnw)
    call mpinp (iu, z2%mpc, mpnw)
    call mpinp (iu, z2%mpc(l2), mpnw)
    return
  end subroutine

  subroutine mp_readz3 (iu, z1, z2, z3, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_complex), intent (out):: z1, z2, z3
    integer l1, l2, l3, mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    l1 = z1%mpc(0)
    l2 = z2%mpc(0)
    l3 = z3%mpc(0)
    call mpinp (iu, z1%mpc, mpnw)
    call mpinp (iu, z1%mpc(l1), mpnw)
    call mpinp (iu, z2%mpc, mpnw)
    call mpinp (iu, z2%mpc(l2), mpnw)
    call mpinp (iu, z3%mpc, mpnw)
    call mpinp (iu, z3%mpc(l3), mpnw)
    return
  end subroutine

  subroutine mp_readz4 (iu, z1, z2, z3, z4, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_complex), intent (out):: z1, z2, z3, z4
    integer l1, l2, l3, l4, mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    l1 = z1%mpc(0)
    l2 = z2%mpc(0)
    l3 = z3%mpc(0)
    l4 = z4%mpc(0)
    call mpinp (iu, z1%mpc, mpnw)
    call mpinp (iu, z1%mpc(l1), mpnw)
    call mpinp (iu, z2%mpc, mpnw)
    call mpinp (iu, z2%mpc(l2), mpnw)
    call mpinp (iu, z3%mpc, mpnw)
    call mpinp (iu, z3%mpc(l3), mpnw)
    call mpinp (iu, z4%mpc, mpnw)
    call mpinp (iu, z4%mpc(l4), mpnw)
    return
  end subroutine

  subroutine mp_readz5 (iu, z1, z2, z3, z4, z5, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_complex), intent (out):: z1, z2, z3, z4, z5
    integer l1, l2, l3, l4, l5, mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    l1 = z1%mpc(0)
    l2 = z2%mpc(0)
    l3 = z3%mpc(0)
    l4 = z4%mpc(0)
    l5 = z5%mpc(0)
    call mpinp (iu, z1%mpc, mpnw)
    call mpinp (iu, z1%mpc(l1), mpnw)
    call mpinp (iu, z2%mpc, mpnw)
    call mpinp (iu, z2%mpc(l2), mpnw)
    call mpinp (iu, z3%mpc, mpnw)
    call mpinp (iu, z3%mpc(l3), mpnw)
    call mpinp (iu, z4%mpc, mpnw)
    call mpinp (iu, z4%mpc(l4), mpnw)
    call mpinp (iu, z5%mpc, mpnw)
    call mpinp (iu, z5%mpc(l5), mpnw)
    return
  end subroutine

  function mp_rtod (ra)
    implicit none
    double precision:: mp_rtod
    type (mp_real), intent (in):: ra
    integer n1, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    call mpmdc (ra%mpr, mp_rtod, n1, mpnw)
    mp_rtod = mp_rtod * 2.d0**n1
    return
  end function
  
  function mp_rtoq (ra)
    implicit none

!>  If real*16 is supported, uncomment this line:
    real (kind (0.q0)):: mp_rtoq
!>  Otherwise uncomment this line:
!    real (kind (0.d0)):: mp_rtoq

    type (mp_real), intent (in):: ra
    integer n1, mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    call mpmqc (ra%mpr, mp_rtoq, n1, mpnw)
    mp_rtoq = mp_rtoq * 2.d0**n1
    return
  end function

  function mp_rtor (ra, iprec)
    implicit none
    type (mp_real):: mp_rtor
    type (mp_real), intent (in):: ra
    integer mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    mp_rtor%mpr(0) = mpwds6
    call mpeq (ra%mpr, mp_rtor%mpr, mpnw)
    return
  end function

  function mp_rtoz (ra, rb, iprec)
    implicit none
    type (mp_complex):: mp_rtoz
    type (mp_real), intent (in):: ra, rb
    integer l1, mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    l1 = mpwds6
    mp_rtoz%mpc(0) = mpwds6
    mp_rtoz%mpc(l1) = mpwds6
    call mpeq (ra%mpr, mp_rtoz%mpc, mpnw)
    call mpeq (rb%mpr, mp_rtoz%mpc(l1), mpnw)
    return
  end function

  function mp_sign (ra, rb)
    implicit none
    type (mp_real):: mp_sign
    type (mp_real), intent (in):: ra, rb
    integer mpnw
    mpnw = max (int (ra%mpr(1)), int (rb%mpr(1)))
    mpnw = min (mpnw, mpwds)
    mp_sign%mpr(0) = mpwds6
    call mpeq (ra%mpr, mp_sign%mpr, mpnw)
    mp_sign%mpr(2) = sign (mp_sign%mpr(2), rb%mpr(2))
    return
  end function
  
  function mp_sin (ra)
    implicit none
    type (mp_real):: mp_sin
    type (mp_real), intent (in):: ra
    integer mpnw
    type (mp_real) r1
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_sin%mpr(0) = mpwds6
    r1%mpr(0) = mpwds6
    call mpcssnr (ra%mpr, r1%mpr, mp_sin%mpr, mpnw)
    return
  end function
  
  function mp_sinh (ra)
    implicit none
    type (mp_real):: mp_sinh
    type (mp_real), intent (in):: ra
    integer mpnw
    type (mp_real) r1
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_sinh%mpr(0) = mpwds6
    r1%mpr(0) = mpwds6
    call mpcsshr (ra%mpr, r1%mpr, mp_sinh%mpr, mpnw)
    return
  end function
  
  function mp_sqrt (ra)
    implicit none
    type (mp_real):: mp_sqrt
    type (mp_real), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_sqrt%mpr(0) = mpwds6
    call mpsqrt (ra%mpr, mp_sqrt%mpr, mpnw)
    return
  end function
  
  function mp_tan (ra)
    implicit none
    type (mp_real):: mp_tan
    type (mp_real), intent (in):: ra
    type (mp_real) r1, r2
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    r2%mpr(0) = mpwds6
    mp_tan%mpr(0) = mpwds6
    call mpcssnr (ra%mpr, r1%mpr, r2%mpr, mpnw)
    if (r1%mpr(2) == 0.d0) then
      write (mpldb, 1)
1     format ('*** MP_TAN: Cos of argument is zero.')
      call mpabrt (26)
    endif
    call mpdiv (r2%mpr, r1%mpr, mp_tan%mpr, mpnw)
    return
  end function

  function mp_tanh (ra)
    implicit none
    type (mp_real):: mp_tanh
    type (mp_real), intent (in):: ra
    type (mp_real) r1, r2
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    r1%mpr(0) = mpwds6
    r2%mpr(0) = mpwds6
    mp_tanh%mpr(0) = mpwds6
    call mpcsshr (ra%mpr, r1%mpr, r2%mpr, mpnw)
    call mpdiv (r2%mpr, r1%mpr, mp_tanh%mpr, mpnw)
    return
  end function

  function mp_wprec (ra)
    implicit none
    integer mp_wprec
    type (mp_real), intent (in):: ra
    integer mpnw
    mp_wprec = min (int (ra%mpr(1)), mpwds)
    return
  end function

  function mp_wprecz (za)
    implicit none
    integer mp_wprecz
    type (mp_complex), intent (in):: za
    integer l1, mpnw
    l1 = za%mpc(0)
    mp_wprecz = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mp_wprecz = min (mp_wprecz, mpwds)
    return
  end function

  subroutine mp_writer (iu, ln, ld, r1, r2, r3, r4, r5)
    implicit none
    integer, intent (in):: iu, ln, ld
    type (mp_real), intent (in):: r1, r2, r3, r4, r5
    optional:: r2, r3, r4, r5
    integer mpnw

    mpnw = min (int (r1%mpr(1)), mpwds)
    call mpout (iu, ln, ld, r1%mpr, mpnw)
    if (present (r2)) then
      mpnw = min (int (r2%mpr(1)), mpwds)
      call mpout (iu, ln, ld, r2%mpr, mpnw)
    endif
    if (present (r3)) then
      mpnw = min (int (r3%mpr(1)), mpwds)
      call mpout (iu, ln, ld, r3%mpr, mpnw)
    endif
    if (present (r4)) then
      mpnw = min (int (r4%mpr(1)), mpwds)
      call mpout (iu, ln, ld, r4%mpr, mpnw)
    endif
    if (present (r5)) then
      mpnw = min (int (r5%mpr(1)), mpwds)
      call mpout (iu, ln, ld, r5%mpr, mpnw)
    endif
    
    return
  end subroutine
  
  subroutine mp_writez (iu, ln, ld, z1, z2, z3, z4, z5)
    implicit none
    integer, intent (in):: iu, ln, ld
    type (mp_complex), intent (in):: z1, z2, z3, z4, z5
    optional:: z2, z3, z4, z5
    integer l1, l2, l3, l4, l5, mpnw

    l1 = z1%mpc(0)
    mpnw = min (int (z1%mpc(1)), mpwds)
    call mpout (iu, ln, ld, z1%mpc, mpnw)
    call mpout (iu, ln, ld, z1%mpc(l1), mpnw)
    if (present (z2)) then
      l2 = z2%mpc(0)
      mpnw = min (int (z2%mpc(1)), mpwds)
      call mpout (iu, ln, ld, z2%mpc, mpnw)
      call mpout (iu, ln, ld, z2%mpc(l2), mpnw)
    endif
    if (present (z3)) then
      l3 = z3%mpc(0)
      mpnw = min (int (z3%mpc(1)), mpwds)
      call mpout (iu, ln, ld, z3%mpc, mpnw)
      call mpout (iu, ln, ld, z3%mpc(l3), mpnw)
    endif
    if (present (z4)) then
      l4 = z4%mpc(0)
      mpnw = min (int (z4%mpc(1)), mpwds)
      call mpout (iu, ln, ld, z4%mpc, mpnw)
      call mpout (iu, ln, ld, z4%mpc(l4), mpnw)
    endif
    if (present (z5)) then
      l5 = z5%mpc(0)
      mpnw = min (int (z5%mpc(1)), mpwds)
      call mpout (iu, ln, ld, z5%mpc, mpnw)
      call mpout (iu, ln, ld, z5%mpc(l5), mpnw)
    endif

    return
  end subroutine
  
  function mp_zeta (ra)
    implicit none
    type (mp_real):: mp_zeta
    type (mp_real), intent (in):: ra
    integer mpnw
    mpnw = min (int (ra%mpr(1)), mpwds)
    mp_zeta%mpr(0) = mpwds6
    call mpzetar (ra%mpr, mp_zeta%mpr, mpnw)
    return
  end function

  function mp_zetaem (nb, rb, rc)
    implicit none
    integer, intent (in):: nb
    type (mp_real):: mp_zetaem
    type (mp_real), intent (in):: rb(nb), rc
    integer n1, mpnw
    mpnw = max (int (rb(1)%mpr(1)), int (rc%mpr(1)))
    mpnw = min (mpnw, mpwds)
    mp_zetaem%mpr(0) = mpwds6
    n1 = mpwds
    call mpzetaemr (n1, nb, rb(1)%mpr, rc%mpr, mp_zetaem%mpr, mpnw)
    return
  end function
 
  function mp_ztodc (za)
    implicit none
    complex (kind(0.d0)):: mp_ztodc
    type (mp_complex), intent (in):: za
    integer l1, mpnw, n1, n2
    double precision d1, d2
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    call mpmdc (za%mpc, d1, n1, mpnw)
    d1 = d1 * 2.d0 ** n1
    call mpmdc (za%mpc(l1), d2, n2, mpnw)
    d2 = d2 * 2.d0 ** n2
    mp_ztodc = dcmplx (d1, d2)
    return
  end function

  function mp_ztor (za)
    implicit none
    type (mp_real):: mp_ztor
    type (mp_complex), intent (in):: za
    integer l1, mpnw
    l1 = za%mpc(0)
    mpnw = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mpnw = min (mpnw, mpwds)
    mp_ztor%mpr(0) = mpwds6
    call mpeq (za%mpc, mp_ztor%mpr, mpnw)
    return
  end function

  function mp_ztoz (za, iprec)
    implicit none
    type (mp_complex):: mp_ztoz
    type (mp_complex), intent (in):: za
    integer l1, l2, mpnw

!>  In version #1, uncomment these six lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnw = mp_setwp (iprec)
!    else
!      mpnw = mpwds
!    endif
!>  Otherwise in version #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnw = mp_setwp (iprec)

    l1 = za%mpc(0)
    l2 = mpwds6
    mp_ztoz%mpc(0) = mpwds6
    mp_ztoz%mpc(l2) = mpwds6
    call mpceq (za%mpc, mp_ztoz%mpc, mpnw)
    return
  end function

end module mpfung

