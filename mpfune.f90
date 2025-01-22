!*****************************************************************************

!  MPFUN-Fort: A thread-safe arbitrary precision computation package
!  Special functions module (module MPFUNE)

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
 
!  DESCRIPTION OF THIS MODULE (MPFUNE):
!    This module contains subroutines to perform special functions. Additional
!    functions will be added as they are completed.

module mpfune
use mpfuna
use mpfunb
use mpfunc
use mpfund
contains

subroutine mpberner (nb1, nb2, berne, mpnw)

!  This returns the even Bernouli numbers B(2*k), from B(2) = 1/6 up to
!  B(2*nb2).  The array berne must be dimensioned as shown below.

implicit none
integer k, nb1, nb2, mpnw, mpnw1
double precision berne(0:nb1+5,nb2), t1(0:mpnw+6), t2(0:mpnw+6), &
  t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6)

!  End of declaration

if (mpnw < 4 .or. berne(0,1) < mpnw + 4 .or. berne(0,nb2) < mpnw + 4) &
  then
  write (mpldb, 2) int (mpndpw * mpnw)
2 format ('*** MPBERN: unitialized or inadequately sized arrays')
  call mpabrt (62)
endif

mpnw1 = mpnw + 1
t1(0) = mpnw + 7
t2(0) = mpnw + 7
t3(0) = mpnw + 7
t4(0) = mpnw + 7
t5(0) = mpnw + 7
call mpmuld (mppicon, 2.d0, t1, mpnw1)
call mpmul (t1, t1, t2, mpnw1)
call mpdmc (-2.d0, 0, t1, mpnw1)

do k = 1, nb2
  call mpmuld (t1, dble (2*k - 1), t3, mpnw1)
  call mpmuld (t3, dble (2*k), t4, mpnw1)
  call mpdiv (t4, t2, t1, mpnw1)
  t1(2) = - t1(2)
  call mpdmc (2.d0 * dble (k), 0, t3, mpnw1)
  call mpzetar (t3, t4, mpnw1)
  call mpmul (t1, t4, t5, mpnw1)
  call mproun (t5, mpnw)
  call mpeq (t5, berne(0,k), mpnw)
enddo

return
end subroutine mpberner

subroutine mpbesseljr (anu, t, z, mpnw)

!   This evaluates the function BesselJ (ANU, T).  ANU must be nonnegative and
!   not greater than 10^6 (this limit can be adjusted below).  To compensate
!   for an unsually large amount of internal cancelation in these formulas, all
!   computations are performed to 3*mpnw/2 words precision.

!   In the parameter statement below:
!     itrmx = limit of number of iterations in series; default = 100000.
!     dasy = factor used to decide if asymptic series is used; default = 25.
!     anumx = upper limit of anu argument; default = 1000.

implicit none
integer i, itrmx, j, mpnw, mpnw1, ndp, nu, n1
double precision dasy, anu(0:mpnw+5), anumx, t(0:mpnw+5), z(0:mpnw+5), &
  t0(0:3*mpnw/2+5), t1(0:3*mpnw/2+5), t2(0:3*mpnw/2+5), t3(0:3*mpnw/2+5), &
  t4(0:3*mpnw/2+5), t5(0:3*mpnw/2+5), t6(0:3*mpnw/2+5)
parameter (itrmx = 100000, dasy = 25.d0, anumx = 1.d6)
  
! End of declaration

if (mpnw < 4 .or. anu(0) < mpnw + 4 .or. anu(0) < abs (anu(2)) + 4 .or. &
  t(0) < mpnw + 4 .or. t(0) < abs (t(2)) + 4 .or. z(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELJR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

if (anu(2) < 0.d0 .or. anu(3) > 0.d0 .or. &
  (anu(3) == 0.d0 .and. anu(4) > anumx)) then
  write (6, 2) anumx
2 format ('*** MPBESSELJR: First argument must be >= 0 and <=',f10.0)
  call mpabrt (65)
endif

mpnw1 = 3 * mpnw / 2
t0(0) = mpnw1 + 6
t1(0) = mpnw1 + 6
t2(0) = mpnw1 + 6
t3(0) = mpnw1 + 6
t4(0) = mpnw1 + 6
t5(0) = mpnw1 + 6
t6(0) = mpnw1 + 6

!   Select either the direct or the asymptotic series.

if (t(3) < 0.d0 .or. t(3) == 0.d0 .and. t(4) < dasy * (mpnw - 2)) then
  t2(1) = mpnw1
  t2(2) = 1.d0
  t2(3) = 0.d0
  t2(4) = 1.d0
  t2(5) = 0.d0
  t2(6) = 0.d0
  call mpadd (anu, t2, t0, mpnw1)
  call mpgammar (t0, t1, mpnw1)
  call mpdiv (t2, t1, t3, mpnw1)
  call mpeq (t3, t1, mpnw1)
  call mpeq (t1, t0, mpnw1)
  call mpmul (t, t, t3, mpnw1)
  call mpmuld (t3, 0.25d0, t2, mpnw1)

  do i = 1, itrmx
    call mpmul (t1, t2, t3, mpnw1)
    call mpdivd (t3, dble (i), t4, mpnw1)
    call mpdmc (dble (i), 0, t5, mpnw1)
    call mpadd (anu, t5, t6, mpnw1)
    call mpdiv (t4, t6, t1, mpnw1)
    t1(2) = - t1(2)
    call mpadd (t0, t1, t3, mpnw1)
    call mpeq (t3, t0, mpnw1)
    if (t1(2) == 0.d0 .or. t1(3) < t0(3) - mpnw1) goto 100
  enddo

  write (6, 3)
3 format ('*** MPBESSELJR: loop overflow 1')
  call mpabrt (66)

100 continue

  call mpmuld (t, 0.5d0, t1, mpnw1)
  call mppower (t1, anu, t2, mpnw1)  
  call mpmul (t0, t2, t3, mpnw1)
  call mpeq (t3, t0, mpnw1)
else
  t0(1) = mpnw1
  t0(2) = 1.d0
  t0(3) = 0.d0
  t0(4) = 1.d0
  t0(5) = 0.d0
  t0(6) = 0.d0
  t1(1) = mpnw1
  t1(2) = 0.d0
  t1(3) = 0.d0
  t1(4) = 0.d0
  t2(1) = mpnw1
  t2(2) = 1.d0
  t2(3) = 0.d0
  t2(4) = 1.d0
  t2(5) = 0.d0
  t2(6) = 0.d0
  call mpmul (anu, anu, t3, mpnw1)
  call mpmuld (t3, 4.d0, t5, mpnw1)
  
  do i = 1, itrmx
    call mpdmc (dble (2*i - 1), 0, t4, mpnw1)
    call mpmul (t4, t4, t6, mpnw1)
    call mpsub (t5, t6, t4, mpnw1)
    call mpmul (t2, t4, t3, mpnw1)
    call mpdivd (t3, 8.d0 * dble (i), t4, mpnw1)
    call mpdiv (t4, t, t2, mpnw1)
    if (mod (i, 2) == 0) then
      call mpeq (t2, t3, mpnw1)
      if (mod (i, 4) == 2)  t3(2) = - t3(2)
      call mpadd (t0, t3, t4, mpnw1)
      call mpeq (t4, t0, mpnw1)
    else
      call mpeq (t2, t3, mpnw1)
      if (mod (i, 4) == 3) t3(2) = - t3(2)
      call mpadd (t1, t3, t4, mpnw1)
      call mpeq (t4, t1, mpnw1)
    endif
    if (t2(2) == 0.d0 .or. (t2(3) < t0(3) - mpnw1 .and. t2(3) < t1(3) - mpnw1)) &
      goto 110
  enddo
  
  write (6, 4)
4 format ('*** MPBESSELJR: loop overflow 2')
  call mpabrt (66)
  
110 continue

  call mpeq (mppicon, t2, mpnw1)
  call mpmul (t2, anu, t4, mpnw1)
  call mpmuld (t4, 0.5d0, t3, mpnw1)
  call mpsub (t, t3, t4, mpnw1)
  call mpmuld (t2, 0.25d0, t3, mpnw1)
  call mpsub (t4, t3, t5, mpnw1)
  call mpcssnr (t5, t3, t4, mpnw1)
  call mpmul (t3, t0, t5, mpnw1)
  call mpmul (t4, t1, t6, mpnw1)
  call mpsub (t5, t6, t3, mpnw1)
  t4(1) = mpnw1
  t4(2) = 1.d0
  t4(3) = 0.d0
  t4(4) = 2.d0
  t4(5) = 0.d0
  t4(6) = 0.d0
  call mpmul (t2, t, t5, mpnw1)
  call mpdiv (t4, t5, t6, mpnw1)
  call mpsqrt (t6, t4, mpnw1)
  call mpmul (t4, t3, t0, mpnw1)
endif

call mproun (t0, mpnw)
call mpeq (t0, z, mpnw)

return
end subroutine mpbesseljr

subroutine mpgammar (t, z, mpnw)

!   This evaluates the gamma function, using an algorithm of R. W. Potter.
!   The argument t must not exceed 10^8 in size (this limit is set below),
!   must not be zero, and if negative must not be integer.

!   In the parameter statement below:
!     itrmx = limit of number of iterations in series; default = 100000.
!     con1 = 1/2 * log (10) to DP accuracy.
!     dmax = maximum size of input argument.

implicit none
integer i, itrmx, j, k, mpnw, mpnw1, ndp, neps, nt, n1, n2, n3
double precision alpha, al2, dmax, d1, d2, d3
parameter (al2 = 0.69314718055994530942d0, dmax = 1d8, itrmx = 100000)
double precision t(0:mpnw+5), z(0:mpnw+5), f1(0:8), sum1(0:mpnw+6), &
  sum2(0:mpnw+6), tn(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), &
  t4(0:mpnw+6), t5(0:mpnw+6), t6(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. t(0) < mpnw + 4 .or. t(0) < abs (t(2)) + 4 .or. &
  z(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPGAMMAR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

if (t(2) == 0 .or. t(3) > 0 .or. (t(3) == 0 .and. t(4) > dmax) .or. &
  (t(2) < 0.d0 .and. t(3) == 0.d0 .and. abs (t(2)) == 1.d0)) then
  write (6, 2) dmax
2 format ('*** MPGAMMAR: input argument must have absolute value <=',f10.0,','/ &
  'must not be zero, and if negative must not be an integer.')
  call mpabrt (65)
endif

mpnw1 = mpnw + 1
f1(0) = 9.d0
f1(1) = mpnw1
f1(2) = 1.d0
f1(3) = 0.d0
f1(4) = 1.d0
f1(5) = 0.d0
f1(6) = 0.d0
sum1(0) = mpnw + 7
sum2(0) = mpnw + 7
tn(0) = mpnw + 7
t1(0) = mpnw + 7
t2(0) = mpnw + 7
t3(0) = mpnw + 7
t4(0) = mpnw + 7
t5(0) = mpnw + 7
t6(0) = mpnw + 7

!   Find the integer and fractional parts of t.

call mpinfr (t, t2, t3, mpnw1)

if (t3(2) == 0.d0) then

!   If t is a positive integer, then apply the usual factorial recursion.

  call mpmdc (t2, d2, n2, mpnw1)
  nt = d2 * 2.d0 ** n2
  call mpeq (f1, t1, mpnw1)

  do i = 2, nt - 1
    call mpmuld (t1, dble (i), t2, mpnw1)
    call mpeq (t2, t1, mpnw1)
  enddo

  call mproun (t1, mpnw)
  call mpeq (t1, z, mpnw)
  goto 120
elseif (t(2) > 0.d0) then

!   Apply the identity Gamma[t+1] = t * Gamma[t] to reduce the input argument
!   to the unit interval.

  call mpmdc (t2, d2, n2, mpnw1)
  nt = d2 * 2.d0 ** n2
  call mpeq (f1, t1, mpnw1)
  call mpeq (t3, tn, mpnw1)
  
  do i = 1, nt
    call mpdmc (dble (i), 0, t4, mpnw1)
    call mpsub (t, t4, t5, mpnw1)
    call mpmul (t1, t5, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
  enddo
else

!   Apply the gamma identity to reduce a negative argument to the unit interval.

  call mpsub (f1, t, t4, mpnw1)
  call mpinfr (t4, t3, t5, mpnw1)
  call mpmdc (t3, d3, n3, mpnw1)
  nt = d3 * 2.d0 ** n3

  call mpeq (f1, t1, mpnw1)
  call mpsub (f1, t5, t2, mpnw1)
  call mpeq (t2, tn, mpnw1)
    
  do i = 0, nt - 1
!    t1 = t1 / (t + dble (i))
    call mpdmc (dble (i), 0, t4, mpnw1)
    call mpadd (t, t4, t5, mpnw1)
    call mpdiv (t1, t5, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
  enddo
endif

!   Calculate alpha = bits of precision * log(2) / 2, then take the nearest integer
!   value, so that d2 = 0.25 * alpha^2 can be calculated exactly in DP.

alpha = anint (0.5d0 * mpnbt * al2 * (mpnw1 + 1))
d2 = 0.25d0 * alpha**2

call mpeq (tn, t2, mpnw1)
call mpdiv (f1, t2, t3, mpnw1)
call mpeq (t3, sum1, mpnw1)

!   Evaluate the series with t.

do j = 1, itrmx
  call mpdmc (dble (j), 0, t6, mpnw1)
  call mpadd (t2, t6, t4, mpnw1)
  call mpmuld (t4, dble (j), t5, mpnw1)
  call mpdiv (t3, t5, t6, mpnw1)
  call mpmuld (t6, d2, t3, mpnw1)
  call mpadd (sum1, t3, t4, mpnw1)
  call mpeq (t4, sum1, mpnw1)
  if (t3(2) == 0.d0 .or. t3(3) < sum1(3) - mpnw1) goto 100
enddo

write (6, 3) 1, itrmx
3 format ('*** MPGAMMAR: iteration limit execeeded',2i10)
call mpabrt (67)

100 continue

call mpeq (tn, t2, mpnw1)
t2(2) = - t2(2)
call mpdiv (f1, t2, t3, mpnw1)
call mpeq (t3, sum2, mpnw1)

!   Evaluate the same series with -t.

do j = 1, itrmx
  call mpdmc (dble (j), 0, t6, mpnw1)
  call mpadd (t2, t6, t4, mpnw1)
  call mpmuld (t4, dble (j), t5, mpnw1)
  call mpdiv (t3, t5, t6, mpnw1)
  call mpmuld (t6, d2, t3, mpnw1)
  call mpadd (sum2, t3, t4, mpnw1)
  call mpeq (t4, sum2, mpnw1)
  if (t3(2) == 0.d0 .or. t3(3) < sum2(3) - mpnw1) goto 110
enddo

write (6, 3) 2, itrmx
call mpabrt (67)

110 continue

!   Compute sqrt (mppic * sum1 / (tn * sin (mppic * tn) * sum2)) 
!   and (alpha/2)^tn terms.

call mpeq (mppicon, t2, mpnw1)
call mpmul (t2, tn, t3, mpnw1)
call mpcssnr (t3, t4, t5, mpnw1)
call mpmul (t5, sum2, t6, mpnw1)
call mpmul (tn, t6, t5, mpnw1)
call mpmul (t2, sum1, t3, mpnw1)
call mpdiv (t3, t5, t6, mpnw1)
t6(2) = - t6(2)
call mpsqrt (t6, t2, mpnw1)

call mpdmc (0.5d0 * alpha, 0, t3, mpnw1)
call mplog (t3, t4, mpnw1)
call mpmul (tn, t4, t5, mpnw1)
call mpexp (t5, t6, mpnw1)
call mpmul (t2, t6, t3, mpnw1)

call mpmul (t1, t3, t4, mpnw1)

!   Round to mpnw words precision.

call mproun (t4, mpnw)
call mpeq (t4, z, mpnw)

120 continue

return
end subroutine mpgammar

subroutine mpincgammar (s, z, g, mpnw)

!  This returns the incomplete gamma function, using a combination of formula
!  8.7.3 of the DLMF (for modest-sized z) and formula 8.11.2 (for large z).

implicit none
integer i, itrmax, k, mpnw, mpnw1, n
double precision dmax
parameter (dmax = 40.d0, itrmax = 1000000)
double precision g(0:mpnw+5), s(0:mpnw+5), t0(0:mpnw+6), t1(0:mpnw+6), &
  t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), t6(0:mpnw+6), &
  z(0:mpnw+5), f1(0:8)

! End of declaration

if (mpnw < 4 .or. s(0) < mpnw + 4 .or. s(0) < abs (s(2)) + 4 .or. &
  z(0) < mpnw + 4 .or. z(0) < abs (z(2)) + 4 .or. g(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPINCGAMMAR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

mpnw1 = mpnw + 1
t0(0) = mpnw + 7
t1(0) = mpnw + 7
t2(0) = mpnw + 7
t3(0) = mpnw + 7
t4(0) = mpnw + 7
t5(0) = mpnw + 7
t6(0) = mpnw + 7
f1(0) = 9.d0
f1(1) = mpnw1
f1(2) = 1.d0
f1(3) = 0.d0
f1(4) = 1.d0
f1(5) = 0.d0
f1(6) = 0.d0

! if (abs (z) < dmax * mpnw) then

  if (z(3) < 0.d0 .or. (z(3) == 0.d0 .and. z(4) < dmax * mpnw)) then

!  t1 = gamma (s)

  call mpgammar (s, t1, mpnw1)

!  t2 = 1.d0 / (s * t1)

  call mpmul (s, t1, t3, mpnw1)
  call mpdiv (f1, t3, t2, mpnw1)

!   t0 = t2

  call mpeq (t2, t0, mpnw1)
  
  do k = 1, itrmax

!    t2 = t2 * z / (s + dble (k))

    call mpmul (t2, z, t5, mpnw1)
    call mpdmc (dble (k), 0, t3, mpnw1)
    call mpadd (s, t3, t4, mpnw1)
    call mpdiv (t5, t4, t2, mpnw1)
    
!    t0 = t0 + t2

    call mpadd (t0, t2, t3, mpnw1)
    call mpeq (t3, t0, mpnw1)

    if (t2(2) == 0.d0 .or. t2(3) < t0(3) - mpnw) goto 100
  enddo
  
  write (mpldb, 2) 1, itrmax
2 format ('*** MPINCGAMMAR: iteration limit exceeded:',2i10)
  call mpabrt (101)

100 continue

!   gammainc = t1 * (1.d0 - z ** s / exp (z) * t0)

  call mppower (z, s, t2, mpnw1)
  call mpexp (z, t3, mpnw1)
  call mpdiv (t2, t3, t4, mpnw1)
  call mpmul (t4, t0, t5, mpnw1)
  call mpsub (f1, t5, t2, mpnw1)
  call mpmul (t1, t2, t3, mpnw1)
  call mpeq (t3, t1, mpnw1)
else
!  t0 = mpreal (1.d0, mpnw)

  t0(2) = 1.d0
  t0(3) = 0.d0
  t0(4) = 1.d0
  t0(5) = 0.d0
  t0(6) = 0.d0
  
!  t1 = mpreal (1.d0, mpnw)
  
  t1(2) = 1.d0
  t1(3) = 0.d0
  t1(4) = 1.d0
  t1(5) = 0.d0
  t1(6) = 0.d0

  do k = 1, itrmax
!    t1 = t1 * (s - dble (k)) / z

    call mpdmc (dble (k), 0, t2, mpnw1)
    call mpsub (s, t2, t3, mpnw1)
    call mpmul (t1, t3, t4, mpnw1)
    call mpdiv (t4, z, t1, mpnw1)

!    t0 = t0 + t1

    call mpadd (t0, t1, t2, mpnw1)
    call mpeq (t2, t0, mpnw1)

    if (t1(2) == 0.d0 .or. t1(3) < t0(3) - mpnw) goto 110
  enddo

  write (mpldb, 2) 2, itrmax
  call mpabrt (101)

110 continue

!  gammainc = z ** (s - 1.d0) / exp (z) * t0

   call mpsub (s, f1, t2, mpnw1)
   call mppower (z, t2, t3, mpnw1)
   call mpexp (z, t4, mpnw1)
   call mpdiv (t3, t4, t2, mpnw1)
   call mpmul (t2, t0, t1, mpnw1)
endif

200 continue

call mproun (t1, mpnw)
call mpeq (t1, g, mpnw)

return
end subroutine mpincgammar

subroutine mpzetar (ss, zz, mpnw)

!   This returns the zeta function at positive real argument SS using an algorithm
!   due to Peter Borwein.

implicit none
integer i, is, itrmax, j, mpnw, mpnw1, n, nwds
double precision dfrac, dlogb, d1
parameter (itrmax = 1000000, dfrac = 16.d0, dlogb = 33.27106466687737d0)
double precision s(0:mpnw+6), ss(0:mpnw+5), zz(0:mpnw+5), &
  t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), &
  tn(0:mpnw+6), tt(0:mpnw+6), f1(0:8)
double precision sgn

!  End of declaration

if (mpnw < 4 .or. ss(0) < mpnw + 4 .or. ss(0) < abs (ss(2)) + 4 .or. &
  zz(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPZETAR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

!   Check if argument is 1 -- undefined.

if (ss(2) == 1.d0 .and. ss(3) == 0.d0 .and. ss(4) == 1.d0) then
  write (mpldb, 2)
2 format ('*** MPZETAR: argument is 1')
  call mpabrt (63)
endif 

mpnw1 = mpnw + 1
s(0) = mpnw + 7
t1(0) = mpnw + 7
t2(0) = mpnw + 7
t3(0) = mpnw + 7
t4(0) = mpnw + 7
t5(0) = mpnw + 7
tn(0) = mpnw + 7
tt(0) = mpnw + 7

!   Set f1 = 1.

f1(0) = 9.d0
f1(1) = mpnw1
f1(2) = 1.d0
f1(3) = 0.d0
f1(4) = 1.d0
f1(5) = 0.d0
f1(6) = 0.d0

!   Check if argument is zero.  If so, the result is -1/2.

if (ss(2) == 0.d0) then
  call mpdmc (-0.5d0, 0, t1, mpnw1)
  goto 200
endif

!   Check if argument is negative.

if (ss(2) < 0.d0) then

!   Check if argument is a negative even integer.  If so, the result is zero.

  call mpmuld (ss, 0.5d0, t1, mpnw1)
  call mpinfr (t1, t2, t3, mpnw1)
  if (t3(2) == 0.d0) then
    t1(1) = mpnw1
    t1(2) = 0.d0
    t1(3) = 0.d0
    t1(4) = 0.d0
    goto 200
  endif

!   Otherwise compute zeta(1-ss), and later apply the reflection formula.

  call mpsub (f1, ss, tt, mpnw1)
else
  call mpeq (ss, tt, mpnw1)
endif

!  Check if argument is large enough that computing with definition is faster.

! if (tt > mpreald (dlogb * mpnw / log (32.d0 * mpnw), mpnw)) then

d1 = dlogb * mpnw / log (32.d0 * mpnw)
if (tt(2) >= 1.d0 .and. (tt(3) > 1.d0 .or. tt(3) == 0.d0 .and. tt(4) > d1)) then

!  t1 = mpreal (1.d0, mpnw)

t1(1) = mpnw1
t1(2) = 1.d0
t1(3) = 0.d0
t1(4) = 1.d0
t1(5) = 0.d0
t1(6) = 0.d0

  do i = 2, itrmax

!    t2 = mpreal (dble (i), mpnw) ** tt

    call mpdmc (dble (i), 0, t4, mpnw1)
    call mppower (t4, tt, t2, mpnw1)
    
!    t3 = 1.d0 / t2

    call mpdiv (f1, t2, t3, mpnw1)

!    t1 = t1 + t3

    call mpadd (t1, t3, t2, mpnw1)
    call mpeq (t2, t1, mpnw1)

    if (t3(2) == 0.d0 .or. t3(3) < - mpnw) goto 200
  enddo
  
  write (mpldb, 3) 1, itrmax
3 format ('*** MPZETAR: iteration limit exceeded',2i10)
  call mpabrt (101)
endif  

n = dfrac * mpnw

! tn = mpreal (2.d0, mpnw) ** n

tn(0) = mpnw + 7
call mpdmc (1.d0, n, tn, mpnw1)

! t1 = - tn

call mpeq (tn, t1, mpnw1)
t1(2) = - t1(2)

! t2 = mpreal (0.d0, mpnw)

t2(2) = 0.d0
t2(3) = 0.d0
t2(4) = 0.d0

! s = mpreal (0.d0, mpnw)

s(1) = mpnw1
s(2) = 0.d0
s(3) = 0.d0
s(4) = 0.d0

sgn = 1.d0

do j = 0, 2 * n - 1
!  t3 = mpreal (dble (j + 1), mpnw) ** tt

  call mpdmc (dble (j + 1), 0, t4, mpnw1)
  call mppower (t4, tt, t3, mpnw1)
  
!  s = s + sgn * t1 / t3

  call mpdiv (t1, t3, t4, mpnw1)
  if (sgn < 0.d0) t4(2) = - t4(2)
  call mpadd (s, t4, t5, mpnw1)
  call mpeq (t5, s, mpnw1)

  sgn = - sgn

  if (j < n - 1) then
!    t2 = mpreal (0.d0, mpnw)

    t2(2) = 0.d0
    t2(3) = 0.d0
    t2(4) = 0.d0

  elseif (j == n - 1) then
!    t2 = mpreal (1.d0, mpnw)

    t2(2) = 1.d0
    t2(3) = 0.d0
    t2(4) = 1.d0
    t2(5) = 0.d0
    t2(6) = 0.d0
    
  else
!     t2 = t2 * dble (2 * n - j) / dble (j + 1 - n)

     call mpmuld (t2, dble (2 * n - j), t3, mpnw1)
     call mpdivd (t3, dble (j + 1 - n), t2, mpnw1)

  endif
!  t1 = t1 + t2

  call mpadd (t1, t2, t3, mpnw1)
  call mpeq (t3, t1, mpnw1)
enddo

! t1 = - s / (tn * (1.d0 - mpreal (2.d0, mpnw) ** (1.d0 - tt)))

call mpsub (f1, tt, t3, mpnw1)
t2(2) = 1.d0
t2(3) = 0.d0
t2(4) = 2.d0
t2(5) = 0.d0
t2(6) = 0.d0
call mppower (t2, t3, t4, mpnw1)
call mpsub (f1, t4, t2, mpnw1)
call mpmul (tn, t2, t3, mpnw1)
call mpdiv (s, t3, t1, mpnw1)
t1(2) = - t1(2)

!   If original argument was negative, apply the reflection formula.

if (ss(2) < 0.d0) then
  call mpgammar (tt, t3, mpnw1)
  call mpmul (t1, t3, t2, mpnw1)
  call mpmul (mppicon, tt, t1, mpnw1)
  call mpmuld (t1, 0.5d0, t3, mpnw1)
  call mpcssnr (t3, t4, t5, mpnw1)
  call mpmul (t2, t4, t1, mpnw1)
  call mpmuld (mppicon, 2.d0, t2, mpnw1)
  call mppower (t2, tt, t3, mpnw1)
  call mpdiv (t1, t3, t2, mpnw1)
  call mpmuld (t2, 2.d0, t1, mpnw1)
endif

200 continue

! zetapbr = t1

call mproun (t1, mpnw)
call mpeq (t1, zz, mpnw)
return
end subroutine mpzetar

subroutine mpzetaemr (nb1, nb2, berne, s, z, mpnw)

!  This evaluates the Riemann zeta function, using the combination of 
!  the definition formula (for large s), and an Euler-Maclaurin scheme
!  (see formula 25.2.9 of the DLMF.  The array berne contains precomputed
!  Bernoulli numbers.  Its dimensions must be as shown below. 

implicit none
integer i, is, itrmax, k, mpnw, mpnw1, nb1, nb2, nn, n1, n2
double precision dfrac, dlogb, d1, d2
parameter (itrmax = 1000000, dfrac = 8.5d0, dlogb = 33.27106466687737d0)
double precision s(0:mpnw+5), z(0:mpnw+5), t0(0:mpnw+6), t1(0:mpnw+6), &
  t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), t6(0:mpnw+6), &
  t7(0:mpnw+6), t8(0:mpnw+6), t9(0:mpnw+6), tt(0:mpnw+6), f1(0:8)
double precision berne(0:nb1+5,nb2)

! End of declaration

if (mpnw < 4 .or. s(0) < mpnw + 4 .or. s(0) < abs (s(2)) + 4 .or. &
  z(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPZETAEMR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

!   Check if argument is 1 -- undefined.

if (s(2) == 1.d0 .and. s(3) == 0.d0 .and. s(4) == 1.d0) then
  write (mpldb, 2)
2 format ('*** MPZETAEMR: argument is 1')
  call mpabrt (63)
endif

!   Check if berne array has been initialized.

if (berne(0,1) < mpnw + 4 .or. berne(0,1) < abs (berne(2,1)) + 4 .or. &
  berne(0,nb2) < mpnw + 4 .or. berne(0,nb2) < abs (berne(2,nb2)) + 4 .or. &
  nb2 < int (mpndpw * mpnw)) then
  write (mpldb, 3) int (mpndpw * mpnw)
3 format ('*** MPZETAEMR: Bernoulli coefficient array must be initialized'/ &
   'with at least',i8,' entries.')
  call mpabrt (62)
endif

i = 0
k = 0
mpnw1 = mpnw + 1
t0(0) = mpnw + 7
t1(0) = mpnw + 7
t2(0) = mpnw + 7
t3(0) = mpnw + 7
t4(0) = mpnw + 7
t5(0) = mpnw + 7
t6(0) = mpnw + 7
t7(0) = mpnw + 7
t8(0) = mpnw + 7
t9(0) = mpnw + 7
tt(0) = mpnw + 7

!   Set f1 = 1.

f1(0) = 9.d0
f1(1) = mpnw1
f1(2) = 1.d0
f1(3) = 0.d0
f1(4) = 1.d0
f1(5) = 0.d0
f1(6) = 0.d0

!   Check if argument is zero.  If so, result is - 1/2.

if (s(2) == 0.d0) then
  call mpdmc (-0.5d0, 0, t1, mpnw)
  goto 200
endif

!   Check if argument is negative.

if (s(2) < 0.d0) then

!   Check if argument is a negative even integer.  If so, the result is zero.

  call mpmuld (s, 0.5d0, t1, mpnw1)
  call mpinfr (t1, t2, t3, mpnw1)
  if (t3(2) == 0.d0) then
    t1(1) = mpnw1
    t1(2) = 0.d0
    t1(3) = 0.d0
    t1(4) = 0.d0
    goto 200
  endif

!   Otherwise compute zeta(1-s), and later apply the reflection formula.

  call mpsub (f1, s, tt, mpnw1)
else
  call mpeq (s, tt, mpnw1)
endif

!  Check if argument is large enough that computing with definition is faster.

! if (tt > mpreald (dlogb * mpnw / log (32.d0 * mpnw), mpnw)) then

d1 = dlogb * mpnw / log (32.d0 * mpnw)
if (tt(3) > 1.d0 .or. (tt(3) == 0.d0 .and. tt(4) > d1)) then

!  t1 = mpreal (1.d0, mpnw)

  t1(1) = mpnw1
  t1(2) = 1.d0
  t1(3) = 0.d0
  t1(4) = 1.d0
  t1(5) = 0.d0
  t1(6) = 0.d0

  do i = 2, itrmax

!    t2 = mpreal (dble (i), mpnw) ** tt

    call mpdmc (dble (i), 0, t4, mpnw1)
    call mppower (t4, tt, t2, mpnw1)
    
!    t3 = 1.d0 / t2

    call mpdiv (f1, t2, t3, mpnw1)

!    t1 = t1 + t3

    call mpadd (t1, t3, t2, mpnw1)
    call mpeq (t2, t1, mpnw1)

    if (t3(2) == 0.d0 .or. t3(3) < - mpnw) goto 200
  enddo
  
  write (mpldb, 4) 1, itrmax
4 format ('*** MPZETAEMR: iteration limit exceeded',2i10)
  call mpabrt (101)
endif  

! t0 = mpreal (1.d0, mpnw)

t0(1) = mpnw1
t0(2) = 1.d0
t0(3) = 0.d0
t0(4) = 1.d0
t0(5) = 0.d0
t0(6) = 0.d0

nn = dfrac * mpnw1

do k = 2, nn
!  t1 = mpreal (dble (k), mpnw) ** tt

  call mpdmc (dble (k), 0, t2, mpnw1)
  call mppower (t2, tt, t1, mpnw1)
  
!  t0 = t0 + 1.d0 / t1

  call mpdiv (f1, t1, t2, mpnw1)
  call mpadd (t0, t2, t3, mpnw1)
  call mpeq (t3, t0, mpnw1)
enddo

! t0 = t0 + dble (nn) / (t1 * (tt - 1.d0)) - 0.5d0 / t1

call mpdmc (dble (nn), 0, t2, mpnw1)
call mpsub (tt, f1, t3, mpnw1)
call mpmul (t1, t3, t4, mpnw1)
call mpdiv (t2, t4, t3, mpnw1)
call mpadd (t0, t3, t2, mpnw1)
call mpdmc (0.5d0, 0, t3, mpnw1)
call mpdiv (t3, t1, t4, mpnw1)
call mpsub (t2, t4, t0, mpnw1)

! t3 = tt

call mpeq (tt, t3, mpnw1)

! t2 = t3 / (12.d0 * dble (nn) * t1)

call mpmuld (t1, 12.d0 * dble (nn), t4, mpnw1)
call mpdiv (t3, t4, t2, mpnw1)

! t5 = dble (nn) * t1

call mpmuld (t1, dble (nn), t5, mpnw1)

! t9 = dble (nn) ** 2

call mpdmc (dble (nn), 0, t6, mpnw1)
call mpmul (t6, t6, t9, mpnw1)

do k = 2, min (nb2, itrmax)
!  t3 = t3 * (tt + dble (2*k - 2)) * (tt + dble (2*k - 3)) / &
!         (dble (2 * k - 1) * dble (2 * k - 2))

  call mpdmc (dble (2 * k - 2), 0, t4, mpnw1)
  call mpadd (tt, t4, t6, mpnw1)
  call mpdmc (dble (2 * k - 3), 0, t7, mpnw1)
  call mpadd (tt, t7, t8, mpnw1)
  call mpmul (t6, t8, t7, mpnw1)
  call mpmul (t3, t7, t4, mpnw1)
  call mpdmc (dble (2 * k - 1), 0, t6, mpnw1)
  call mpdmc (dble (2 * k - 2), 0, t7, mpnw1)
  call mpmul (t6, t7, t8, mpnw1)
  call mpdiv (t4, t8, t3, mpnw1)
  
!  t5 = t5 * t9

  call mpmul (t5, t9, t6, mpnw1)
  call mpeq (t6, t5, mpnw1)
  
!  t7 = t3 * berne(k) / (dble (2 * k) * t5)

  call mpmul (t3, berne(0,k), t4, mpnw1)
  call mpmuld (t5, dble (2 * k), t6, mpnw1)
  call mpdiv (t4, t6, t7, mpnw1)

!  t2 = t2 + t7

  call mpadd (t2, t7, t4, mpnw1)
  call mpeq (t4, t2, mpnw1)

  if (t7(2) == 0.d0 .or. t7(3) < t2(3) - mpnw) goto 110
enddo

write (mpldb, 3) 2, min (nb2, itrmax)
call mpabrt (101)

110 continue

! zetaem = t0 + t2

call mpadd (t0, t2, t1, mpnw1)

!   If original argument was negative, apply the reflection formula.

if (s(2) < 0.d0) then
  call mpgammar (tt, t3, mpnw1)
  call mpmul (t1, t3, t2, mpnw1)
  call mpmul (mppicon, tt, t1, mpnw1)
  call mpmuld (t1, 0.5d0, t3, mpnw1)
  call mpcssnr (t3, t4, t5, mpnw1)
  call mpmul (t2, t4, t1, mpnw1)
  call mpmuld (mppicon, 2.d0, t2, mpnw1)
  call mppower (t2, tt, t3, mpnw1)
  call mpdiv (t1, t3, t2, mpnw1)
  call mpmuld (t2, 2.d0, t1, mpnw1)
endif

200 continue

call mproun (t1, mpnw)
call mpeq (t1, z, mpnw)

return
end subroutine mpzetaemr

end module mpfune
