!*****************************************************************************

!  MPFUN-Fort: A thread-safe arbitrary precision computation package
!  Basic function module (module MPFUNB); includes real*16 support.
!  Search for !> for version differences.

!  Revision date:  5 Feb 2016

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

!  DESCRIPTION OF THIS MODULE (MPFUNB):
!    This module contains routines for: add, subtract, multiply, divide;
!    comparison; double/multi conversion; double/multi multiplication/division;
!    integer and fractional parts; nearest integer; nth power; nth root;
!    square root; and rounding and normalization.  Routines in this package
!    are not intended to be directly called by the user; instead, use the 
!    high-level language interface modules.

module mpfunb
use mpfuna

contains

subroutine mpabrt (ier)

!   This routine terminates execution.  Users may wish to replace the
!   default STOP with a call to a system routine that provides a traceback.

implicit none
integer ier

! End of declaration

write (mpldb, 1) ier
1 format ('*** MPABRT: Execution terminated, error code =',i4)
stop
end subroutine mpabrt

subroutine mpadd (a, b, c, mpnw)

!   This routine adds MPR numbers A and B to yield C.

implicit none
double precision a(0:mpnw+5), b(0:mpnw+5), c(0:mpnw+5), d(0:mpnw+6), db
integer i, ia, ib, ic, ish, ixa, ixb, ixd, mpnw, m1, m2, m3, m4, m5, na, nb, &
  nd, nsh

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < abs (b(2)) + 4 .or. &
  c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPADD: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

ia = sign (1.d0, a(2))
ib = sign (1.d0, b(2))
na = min (int (abs (a(2))), mpnw)
nb = min (int (abs (b(2))), mpnw)

!   Check for zero inputs.

if (na == 0) then

!   A is zero -- the result is B.

  c(1) = mpnw
  c(2) = sign (nb, ib)

  do i = 2, nb + 2
    c(i+1) = b(i+1)
  enddo

  goto 100
elseif (nb == 0) then

!   B is zero -- the result is A.

  c(1) = mpnw
  c(2) = sign (na, ia)

  do i = 2, na + 2
    c(i+1) = a(i+1)
  enddo

  goto 100
endif

if (ia == ib) then
  db = 1.d0
else
  db = -1.d0
endif
ixa = a(3)
ixb = b(3)
ish = ixa - ixb

if (ish >= 0) then

!   A has greater exponent than B, so B must be shifted to the right.

  m1 = min (na, ish)
  m2 = min (na, nb + ish)
  m3 = na
  m4 = min (max (na, ish), mpnw + 1)
  m5 = min (max (na, nb + ish), mpnw + 1)

  do i = 1, m1
    d(i+3) = a(i+3)
  enddo

  do i = m1 + 1, m2
    d(i+3) = a(i+3) + db * b(i+2-ish+1)
  enddo

  do i = m2 + 1, m3
    d(i+3) = a(i+3)
  enddo

  do i = m3 + 1, m4
    d(i+3) = 0.d0
  enddo

  do i = m4 + 1, m5
    d(i+3) = db * b(i+2-ish+1)
  enddo

  nd = m5
  ixd = ixa
  d(nd+4) = 0.d0
  d(nd+5) = 0.d0
else

!   B has greater exponent than A, so A must be shifted to the right.

  nsh = - ish
  m1 = min (nb, nsh)
  m2 = min (nb, na + nsh)
  m3 = nb
  m4 = min (max (nb, nsh), mpnw + 1)
  m5 = min (max (nb, na + nsh), mpnw + 1)

  do i = 1, m1
    d(i+3) = db * b(i+3)
  enddo

  do i = m1 + 1, m2
    d(i+3) = a(i+2-nsh+1) + db * b(i+3)
  enddo

  do i = m2 + 1, m3
    d(i+3) = db * b(i+3)
  enddo

  do i = m3 + 1, m4
    d(i+3) = 0.d0
  enddo

  do i = m4 + 1, m5
    d(i+3) = a(i+2-nsh+1)
  enddo

  nd = m5
  ixd = ixb
  d(nd+4) = 0.d0
  d(nd+5) = 0.d0
endif

!   Call mpnorm to fix up result and store in c.

d(0) = mpnw + 7
d(1) = mpnw
d(2) = sign (nd, ia)
d(3) = ixd

call mpnorm (d, c, mpnw)

100 continue

return
end subroutine mpadd

subroutine mpcabs (a, b, mpnw)

!   This routine returns the absolute value of the MPC argument A (the
!   result is of type MPR).

implicit none
integer la, lb, mpnw, mpnw1
double precision a(0:2*mpnw+11), b(0:2*mpnw+5), s0(0:mpnw+6), s1(0:mpnw+6), &
  s2(0:mpnw+6)

! End of declaration

la = a(0)
lb = b(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCABS: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

mpnw1 = mpnw + 1
s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
call mpmul (a, a, s0, mpnw1)
call mpmul (a(la), a(la), s1, mpnw1)
call mpadd (s0, s1, s2, mpnw1)
call mpsqrt (s2, s0, mpnw1)
call mproun (s0, mpnw)
call mpeq (s0, b, mpnw)

return
end subroutine mpcabs

subroutine mpcadd (a, b, c, mpnw)

!   This routine adds the MPC numbers A and B.

implicit none
integer la, lb, lc, mpnw
double precision a(0:2*mpnw+11), b(0:2*mpnw+11), c(0:2*mpnw+11)

! End of declaration

la = a(0)
lb = b(0)
lc = c(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < abs (b(2)) + 4 .or. b(lb) < abs (b(lb+2)) + 4 .or. &
  c(0) < mpnw + 6 .or. c(lc) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCADD: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

call mpadd (a, b, c, mpnw)
call mpadd (a(la), b(lb), c(lc), mpnw)
return
end subroutine mpcadd

subroutine mpcdiv (a, b, c, mpnw)

!   This routine divides the MPC numbers A and B.

implicit none
integer la, lb, lc, mpnw, mpnw1
double precision a(0:2*mpnw+11), b(0:2*mpnw+11), c(0:2*mpnw+11), &
  s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6), s4(0:mpnw+6)

! End of declaration

la = a(0)
lb = b(0)
lc = c(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < abs (b(2)) + 4 .or. b(lb) < abs (b(lb+2)) + 4 .or. &
  c(0) < mpnw + 6 .or. c(lc) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCDIV: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

mpnw1 = mpnw + 1
s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
s4(0) = mpnw + 7

call mpmul (a, b, s0, mpnw1)
call mpmul (a(la), b(lb), s1, mpnw1)
call mpadd (s0, s1, s2, mpnw1)
call mpmul (a, b(lb), s0, mpnw1)
s0(2) = - s0(2)
call mpmul (a(la), b, s1, mpnw1)
call mpadd (s0, s1, s3, mpnw1)

call mpmul (b, b, s0, mpnw1)
call mpmul (b(lb), b(lb), s1, mpnw1)
call mpadd (s0, s1, s4, mpnw1)
call mpdiv (s2, s4, s0, mpnw1)
call mpdiv (s3, s4, s1, mpnw1)


call mproun (s0, mpnw)
call mproun (s1, mpnw)
call mpeq (s0, c, mpnw)
call mpeq (s1, c(lc), mpnw)

return
end subroutine mpcdiv

subroutine mpceq (a, b, mpnw)

!   Sets the MPC number B equal to A.

implicit none
integer i, ia, la, lb, mpnw, na
double precision a(0:2*mpnw+11), b(0:2*mpnw+11)

! End of declaration

la = a(0)
lb = b(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < mpnw + 6 .or. b(lb) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCEQ: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

ia = sign (1.d0, a(2))
na = min (int (abs (a(2))), mpnw)
if (na == 0)  then
  b(1) = mpnw
  b(2) = 0.d0
  b(3) = 0.d0
  goto 110
endif
b(1) = mpnw
b(2) = sign (na, ia)

do i = 2, na + 2
  b(i+1) = a(i+1)
enddo

b(na+4) = 0.d0
b(na+5) = 0.d0

110 continue

ia = sign (1.d0, a(la+2))
na = min (int (abs (a(la+2))), mpnw)
if (na == 0)  then
  b(lb+1) = mpnw
  b(lb+2) = 0.d0
  b(lb+3) = 0.d0
  goto 120
endif
b(lb+1) = mpnw
b(lb+2) = sign (na, ia)

do i = 2, na + 2
  b(i+lb+1) = a(i+la+1)
enddo

b(na+lb+4) = 0.d0
b(na+lb+5) = 0.d0

120 continue

return
end subroutine mpceq

subroutine mpcmul (a, b, c, mpnw)

!   This routine multiplies the MPC numbers A and B.

implicit none
integer la, lb, lc, mpnw, mpnw1
double precision a(0:2*mpnw+11), b(0:2*mpnw+11), c(0:2*mpnw+11), &
  s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6)

! End of declaration

la = a(0)
lb = b(0)
lc = c(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < abs (b(2)) + 4 .or. b(lb) < abs (b(lb+2)) + 4 .or. &
  c(0) < mpnw + 6 .or. c(lc) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCMUL: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

mpnw1 = mpnw + 1
s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7

call mpmul (a, b, s0, mpnw1)
call mpmul (a(la), b(lb), s1, mpnw1)
call mpsub (s0, s1, s2, mpnw1)
call mpmul (a, b(lb), s0, mpnw1)
call mpmul (a(la), b, s1, mpnw1)
call mpadd (s0, s1, s3, mpnw1)

call mproun (s2, mpnw)
call mproun (s3, mpnw)
call mpeq (s2, c, mpnw)
call mpeq (s3, c(lc), mpnw)

return
end subroutine mpcmul

subroutine mpcnpwr (a, n, b, mpnw)

!   This computes the N-th power of the MPC number A and returns the MPC result
!   in B.  When N is zero, 1 is returned.  When N is negative, the reciprocal
!   of A ^ |N| is returned.

implicit none
integer i, j, k, kk, kn, k0, k1, k2, la, lb, lc, mn, mpnw, mpnw1, n, na, nn, nws
double precision cl2, t1, mprxx
parameter (cl2 = 1.4426950408889633d0, mprxx = 1d-14)
double precision a(0:mpnw+11), b(0:mpnw+11), s0(0:2*mpnw+13), s1(0:2*mpnw+13), &
  s2(0:2*mpnw+13)
  
! End of declaration

la = a(0)
lb = b(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 .or. &
  b(0) < mpnw + 6 .or. b(lb) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCNPWR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

na = min (int (abs (a(2))), mpnw)
if (na == 0) then
  if (n >= 0) then
  	b(1) = mpnw
    b(2) = 0.d0
    b(3) = 0.d0
    goto 120
  else
    write (mpldb, 2)
2   format ('*** MPCNPWR: Argument is zero and N is negative or zero.')
    call mpabrt (57)
  endif
endif

mpnw1 = mpnw + 1
lc = mpnw + 7
s0(0) = mpnw + 7
s0(lc) = mpnw + 7
s1(0) = mpnw + 7
s1(lc) = mpnw + 7
s2(0) = mpnw + 7
s2(lc) = mpnw + 7

nn = abs (n)
if (nn == 0) then
  call mpdmc (1.d0, 0, b, mpnw)
  call mpdmc (0.d0, 0, b(lb), mpnw)
  goto 120
elseif (nn == 1) then
  call mpceq (a, s2, mpnw1)
  goto 110
elseif (nn == 2) then
  call mpcmul (a, a, s2, mpnw1) 
  goto 110
endif

!   Determine the least integer MN such that 2 ^ MN .GT. NN.

t1 = nn
mn = cl2 * log (t1) + 1.d0 + mprxx
call mpdmc (1.d0, 0, s2, mpnw1)
call mpceq (a, s0, mpnw1)
kn = nn

!   Compute B ^ N using the binary rule for exponentiation.

do j = 1, mn
  kk = kn / 2
  if (kn /= 2 * kk) then
    call mpcmul (s2, s0, s1, mpnw1)
    call mpceq (s1, s2, mpnw1)
  endif
  kn = kk
  if (j < mn) then
    call mpcmul (s0, s0, s1, mpnw1)
    call mpceq (s1, s0, mpnw1)
  endif
enddo

!   Compute reciprocal if N is negative.

110 continue

if (n < 0) then
  call mpdmc (1.d0, 0, s1, mpnw1)
  call mpdmc (0.d0, 0, s1(lc), mpnw1)
  call mpcdiv (s1, s2, s0, mpnw1)
  call mpceq (s0, s2, mpnw1)
endif

!   Restore original precision level.

call mproun (s2, mpnw)
call mproun (s2(lc), mpnw)
call mpceq (s2, b, mpnw)

120 continue
return
end subroutine mpcnpwr

subroutine mpconjg (a, b, mpnw)

!   This routine returns the conjugate of the MPC argument A.

implicit none
integer la, lb, mpnw
double precision a(0:2*mpnw+11), b(0:2*mpnw+11)

! End of declaration

la = a(0)
lb = b(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < mpnw + 6 .or. b(lb) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCONJ: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

call mpceq (a, b, mpnw)
b(lb+2) = - b(lb+2)
return
end subroutine mpconjg

subroutine mpcsqrt (a, b, mpnw)

!   This routine returns the square root of the MPC argument A.
!   The formula is:

!   1/Sqrt[2] * (Sqrt[r + a1] + I * a2 / Sqrt[r + a1])  if a1 >= 0, or
!   1/Sqrt[2] * (|a2| / Sqrt[r - a1] + I * Sgn[a2] * Sqrt[r - a1]) if a1 < 0,

!   where r = Sqrt[a1^2 + a2^2], and a1 and a2 are the real and imaginary
!   parts of A.

implicit none
integer la, lb, mpnw, mpnw1
double precision a(0:2*mpnw+11), b(0:2*mpnw+11), s0(0:mpnw+6), &
  s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6), s4(0:mpnw+6)

! End of declaration

la = a(0)
lb = b(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < mpnw + 6 .or. b(lb) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCSQRT: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

mpnw1 = mpnw + 1 
s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
s4(0) = mpnw + 7

call mpmul (a, a, s0, mpnw1)
call mpmul (a(la), a(la), s1, mpnw1)
call mpadd (s0, s1, s2, mpnw1)
call mpsqrt (s2, s0, mpnw1)
  
if (a(2) >= 0.d0) then
  call mpadd (s0, a, s1, mpnw1)
  call mpsqrt (s1, s0, mpnw1) 
  call mpdiv (a(la), s0, s1, mpnw1)
else
  call mpsub (s0, a, s2, mpnw1)
  call mpsqrt (s2, s1, mpnw1)  
  call mpdiv (a(la), s1, s0, mpnw1)
  s0(2) = abs (s0(2))
  s1(2) = sign (s1(2), a(la+2))
endif

call mpeq (mpsqrt22con, s2, mpnw1)
call mpmul (s0, s2, s3, mpnw1)
call mpmul (s1, s2, s4, mpnw1)
  
call mproun (s3, mpnw)
call mproun (s4, mpnw)
call mpeq (s3, b, mpnw)
call mpeq (s4, b(lb), mpnw)

return
end subroutine mpcsqrt

subroutine mpcsub (a, b, c, mpnw)

!   This routine subtracts the MPC numbers A and B.

implicit none
integer la, lb, lc, mpnw
double precision a(0:2*mpnw+11), b(0:2*mpnw+11), c(0:2*mpnw+11)

! End of declaration

la = a(0)
lb = b(0)
lc = c(0)
if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(la) < abs (a(la+2)) + 4 &
  .or. b(0) < abs (b(2)) + 4 .or. b(lb) < abs (b(lb+2)) + 4 .or. &
  c(0) < mpnw + 6 .or. c(lc) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPCSUB: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

call mpsub (a, b, c, mpnw)
call mpsub (a(la), b(lb), c(lc), mpnw)
return
end subroutine mpcsub

subroutine mpcpr (a, b, ic, mpnw)

!   This routine compares the MPR numbers A and B and returns in IC the value
!   -1, 0, or 1 depending on whether A < B, A = B, or A > B.
!   Note that the first and second words do NOT need to be the same for the
!   result to be "equal".

implicit none
double precision a(0:mpnw+5), b(0:mpnw+5), s0(0:mpnw+5)
integer i, ia, ib, ic, ma, mb, mpnw, na, nb

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < abs (b(2)) + 4) then
  write (mpldb, 1)
1 format ('*** MPCPR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

s0(0) = mpnw + 6
call mpsub (a, b, s0, mpnw)
if (s0(2) < 0.d0) then
  ic = -1
elseif (s0(2) == 0.d0) then
  ic = 0
else
  ic = 1
endif

return
end subroutine mpcpr

subroutine mpdiv (a, b, c, mpnw)

!   This divides the MPR number A by the MP number B to yield C.

implicit none
integer i, i1, i2, i3, ia, ib, ij, is, j, j3, md, mpnw, mpnwx, na, nb, nc
parameter (mpnwx = 200)
double precision a1, a2, b1, b2, c1, c2, dc, rb, ss, t0, t1, t2
double precision a(0:mpnw+5), b(0:mpnw+5), c(0:mpnw+5), d(0:mpnw+5)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < abs (b(2)) + 4 .or. &
  c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPDIV: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

ia = sign (1.d0, a(2))
ib = sign (1.d0, b(2))
na = min (int (abs (a(2))), mpnw)
nb = min (int (abs (b(2))), mpnw)

!   Check if dividend is zero.                                                            

if (na .eq. 0) then
  c(1) = mpnw
  c(2) = 0.d0
  c(3) = 0.d0
  c(4) = 0.d0
  goto 190
endif

!   Check if divisor is zero.

if (nb == 0) then
  write (mpldb, 2)
2 format ('*** MPDIV: Divisor is zero.')
  call mpabrt (31)
endif

if (na > mpnwx .and. nb > mpnwx) then

!   Precision is very high, so call mpdivx.

  call mpdivx (a, b, c, mpnw)
  goto 200
endif

!   Compute double precision approximation to trial divisor and its reciprocal.

t0 = mpbdx * b(4)
if (nb >= 2) t0 = t0 + b(5)
if (nb >= 3) t0 = t0 + mprdx * b(6)
rb = 1.d0 / t0

md = min (na + nb, mpnw)
d(0) = mpnw + 6
d(1) = mpnw
d(2) = 0.d0

do i = 2, na + 1
  d(i+1) = a(i+2)
enddo

do i = na + 2, md + 4
  d(i+1) = 0.d0
enddo

!   Perform ordinary long division algorithm.

do j = 2, mpnw + 3

!   Compute trial dividend and trial quotient term.

  t1 = mpbdx**2 * d(j) + mpbdx * d(j+1) + d(j+2)
  if (j <= mpnw + 2) t1 = t1 + mprdx * d(j+3)
  t0 = anint (rb * t1)

!   Split trial quotient term into high and low order halves, 24-bits each.

  a1 = mpb24x * aint (mpr24x * t0)
  a2 = t0 - a1
  j3 = j - 3
  i2 = min (nb, mpnw + 2 - j3) + 2
  ij = i2 + j3

!   Multiply trial quotient term by each term of divisor, saving high and low-
!   order parts into appropriate terms of the running dividend.

  do i = 3, i2
    i3 = i + j3
    b1 = mpb24x * aint (mpr24x * b(i+1))
    b2 = b(i+1) - b1
    dc = a1 * b2 + a2 * b1
    c1 = mpbdx * aint (mprdx * dc)
    c2 = dc - c1
    d(i3) = d(i3) - mprdx * (a1 * b1 + c1)
    d(i3+1) = d(i3+1) - a2 * b2 - c2
  enddo

!   Release carries periodically to avoid overflowing the exact integer
!   capacity of double precision floating point words in D.

  if (mod (j - 1, mpnpr) .eq. 0) then
    do i = j + 1, ij
      t1 = d(i)
      t2 = int (mprdx * t1)
      d(i) = t1 - mpbdx * t2
      d(i-1) = d(i-1) + t2
    enddo
  endif

  d(j+1) = d(j+1) + mpbdx * d(j)
  d(j) = t0

!  Continue computing quotient terms past the end of the input dividend, until
!  trial running dividend is zero.

  if (j >= na + 2) then
    if (ij <= mpnw + 1) d(ij+4) = 0.d0
  endif
enddo

!   Final bookkeeping.

j = mpnw + 3

170 continue

d(j+1) = 0.d0
if (d(2) == 0.d0) then
  is = 1
else
  is = 2
endif
nc = min (j - 1, mpnw)
d(nc+4) = 0.d0
d(nc+5) = 0.d0

do i = j + 1, 3, -1
  d(i+1) = d(i-is+1)
enddo

d(0) = mpnw + 6
d(1) = mpnw
d(2) = sign (nc, ia * ib)
d(3) = a(3) - b(3) + is - 2
c(1) = mpnw

!   Call mpnorm to fix up any remaining bugs and perform rounding.

call mpnorm (d, c, mpnw)

190 continue
200 continue
return
end subroutine mpdiv

subroutine mpdivd (a, b, c, mpnw)

!   This routine divides the MPR number A by the DP number B to yield C.

!   NOTE however that if A = 0.1D0, for example, then C will NOT be the true 
!   multiprecision equivalent of the quotient, since 0.1d0 is not an exact
!   binary value.

!   Examples of exact binary values (good): 123456789.d0, 0.25d0, -5.3125d0.
!   Examples of inexact binary values (bad): 0.1d0, 1234567.8d0, -3333.3d0.

implicit none
integer i, ij, is, i1, i2, i3, ia, ib, j, j3, k, md, mpnw, n, na, nb, nc, n1
double precision a1, a2, b, bb, b1, b2, br, c1, c2, dc, dd, d1, d2, rb, &
  t0, t1
double precision a(0:mpnw+5), c(0:mpnw+5), d(0:mpnw+5)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPDIVD: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

ia = sign (1.d0, a(2))
na = min (int (abs (a(2))), mpnw)
ib = sign (1.d0, b)

!   Check if dividend is zero.

if (na == 0) then
  c(1) = mpnw
  c(2) = 0.d0
  c(3) = 0.d0
  goto 190
endif

!   Check if divisor is zero.

if (b == 0.d0) then
  write (mpldb, 3)
3 format ('*** MPDIVD: Divisor is zero.')
  call mpabrt (32)
endif

n1 = 0
bb = abs (b)

!   Reduce BB to within 1 and MPBDX.

if (bb >= mpbdx) then

  do k = 1, 100
    bb = mprdx * bb
    if (bb < mpbdx) then
      n1 = n1 + k
      goto 120
    endif
 enddo

elseif (bb < 1.d0) then

  do k = 1, 100
    bb = mpbdx * bb
    if (bb >= 1.d0) then
      n1 = n1 - k
      goto 120
    endif
 enddo

endif

120 continue

!   If B cannot be represented exactly in a single mantissa word, use MPDIV.

if (bb /= aint (bb)) then
  bb = sign (bb, b)
  d(0) = mpnw + 6
  d(1) = mpnw
  call mpdmc (bb, n1 * mpnbt, d, mpnw)
  call mpdiv (a, d, c, mpnw)
  goto 190
endif

!   Compute double precision approximation to trial divisor and its reciprocal.

t0 = mpbdx * bb
rb = 1.d0 / t0
b1 = mpb24x * aint (mpr24x * bb)
b2 = bb - b1

md = min (na + 1, mpnw)
d(0) = mpnw + 6
d(1) = mpnw
d(2) = 0.d0

do i = 2, na + 1
  d(i+1) = a(i+2)
enddo

do i = na + 2, md + 4
  d(i+1) = 0.d0
enddo

!   Perform ordinary short division algorithm.

do j = 2, mpnw + 3

!   Compute trial dividend and trial quotient term.

  t1 = mpbdx**2 * d(j) + mpbdx * d(j+1) + d(j+2)
  if (j <= mpnw + 2) t1 = t1 + mprdx * d(j+3)
  t0 = anint (rb * t1)

!   Split trial quotient term into high and low order halves, 24-bits each.

  a1 = mpb24x * aint (mpr24x * t0)
  a2 = t0 - a1

!   Multiply trial quotient term by each term of divisor, saving high and low-
!   order parts into appropriate terms of the running dividend.

  dc = a1 * b2 + a2 * b1
  c1 = mpbdx * aint (mprdx * dc)
  c2 = dc - c1
  d(j) = d(j) - mprdx * (a1 * b1 + c1)
  d(j+1) = d(j+1) - a2 * b2 - c2
  d(j+1) = d(j+1) + mpbdx * d(j)
  d(j) = t0

!  Continue computing quotient terms past the end of the input dividend, until
!  trial running dividend is zero.

  if (j >= na + 2) then
    if (j <= mpnw + 1) d(j+4) = 0.d0
  endif
enddo

!   Final bookkeeping.

j = mpnw + 3

170 continue

d(j+1) = 0.d0
if (d(2) == 0.d0) then
  is = 1
else
  is = 2
endif
nc = min (j - 1, mpnw)
d(nc+4) = 0.d0
d(nc+5) = 0.d0

do i = j + 1, 3, -1
  d(i+1) = d(i-is+1)
enddo

d(0) = mpnw + 6
d(1) = mpnw
d(2) = sign (nc, ia * ib)
d(3) = a(3) - n1 + is - 2
c(1) = mpnw

!   Call mpnorm to fix up any remaining bugs and perform rounding.

call mpnorm (d, c, mpnw)

190 continue
return
end subroutine mpdivd

subroutine mpdivd40 (a, b, c, mpnw)

!   This routine divides the MPR number A by the DP number B to yield C.
!   In contrast to mpdivd, this routine only allows 40 significant bits
!   (approximately 12 significant decimal digits) in B.  If more nonzero bits
!   are present, an error is flagged.

!   Examples of exact binary values (good): 123456789.d0, 0.25d0, -5.3125d0.
!   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.

implicit none
integer mpnw
double precision a(0:mpnw+5), b, c(0:mpnw+5), t1, t2

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. c(0) < mpnw + 6) then
 write (mpldb, 1)
1 format ('*** MPDIVD40: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

!   This convoluted-looking code tests whether B has more than 40
!   significant bits (actually whether the trailing 13 bits are zero).

t1 = mpb13x * abs (b)
t2 = abs (abs (b) + t1) - abs (t1)
if (t2 == abs (b)) then
  call mpdivd (a, b, c, mpnw)
else
  write (mpldb, 2) b
2 format ('*** MPDIVD40: DP value has more than 40 significant bits:', &
  1p,d25.15/'and thus very likely represents an unintended loss of accuracy.'/ &
  'Fix the issue, or else use functions mpprodd, mpquotd, mpreald or mpcmplxdc.'/ &
  'See documentation for details.')
  call mpabrt (81)
endif

return
end subroutine mpdivd40

subroutine mpdmc (a, n, b, mpnw)

!   This routine converts the DP number A * 2^N to MPR form in B.

!   NOTE however that if A = 0.1D0, for example, then B will NOT be the true 
!   multiprecision equivalent of 1/10, since 0.1d0 is not an exact binary value.

!   Examples of exact binary values (good): 123456789.d0, 0.25d0, -5.3125d0.
!   Examples of inexact binary values (bad): 0.1d0, 1234567.8d0, -3333.3d0.

implicit none
integer i, k, mpnw, n, n1, n2
double precision a, aa, b(0:*)

! End of declaration

if (mpnw < 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPDMC: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

!   Check for zero.

if (a == 0.d0) then
  b(1) = mpnw
  b(2) = 0.d0
  b(3) = 0.d0
  goto 150
endif
n1 = n / mpnbt
n2 = n - mpnbt * n1
aa = abs (a) * 2.d0 ** n2

!   Reduce AA to within 1 and MPBDX.

if (aa >= mpbdx) then

  do k = 1, 100
    aa = mprdx * aa
    if (aa < mpbdx) then
      n1 = n1 + k
      goto 120
    endif
 enddo

elseif (aa < 1.d0) then

  do k = 1, 100
    aa = mpbdx * aa
    if (aa >= 1.d0) then
      n1 = n1 - k
      goto 120
    endif
  enddo

endif

!   Store successive sections of AA into B.

120  continue

b(3) = n1
b(4) = aint (aa)
aa = mpbdx * (aa - b(3+1))
b(5) = aint (aa)
aa = mpbdx * (aa - b(4+1))
b(6) = aint (aa)
b(7) = 0.d0
b(8) = 0.d0

do i = 6, 3, -1
  if (b(i+1) /= 0.d0) goto 140
enddo

140  continue

b(1) = mpnw
aa = i - 2
b(2) = sign (aa, a)

150 continue
return
end subroutine mpdmc

subroutine mpdmc40 (a, n, b, mpnw)

!   This routine converts the DP number A * 2^N to MPR form in B.  In contrast
!   to mpdmc, this routine only allows 40 significant bits (approximately
!   12 significant decimal digits) in A.  If more nonzero bits are present,
!   an error is flagged.

!   Examples of exact binary values (good): 123456789.d0, 0.25d0, -5.3125d0.
!   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.

implicit none
integer mpnw, n
double precision a, b(0:*), t1, t2

! End of declaration

if (mpnw < 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPDMC40: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

!   This convoluted-looking code tests whether A has more than 40
!   significant bits (actually whether the trailing 13 bits are zero).

t1 = mpb13x * abs (a)
t2 = abs (abs (a) + t1) - abs (t1)
if (t2 == abs (a)) then
  call mpdmc (a, n, b, mpnw)
else
  write (mpldb, 2) a
2 format ('*** MPDMC40: DP value has more than 40 significant bits:', &
  1p,d25.15/'and thus very likely represents an unintended loss of accuracy.'/ &
  'Fix the issue, or else use functions mpprodd, mpquotd, mpreald or mpcmplxdc.'/ &
  'See documentation for details.')
  call mpabrt (82)
endif

return
end subroutine mpdmc40

subroutine mpeq (a, b, mpnw)

!   Sets the MPR number B equal to the MPR number A.

implicit none
integer i, ia, mpnw, na
double precision a(0:mpnw+5), b(0:mpnw+5)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPEQ: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

ia = sign (1.d0, a(2))
na = min (int (abs (a(2))), mpnw)
if (na == 0)  then
  b(1) = mpnw
  b(2) = 0.d0
  b(3) = 0.d0
  goto 110
endif
b(1) = mpnw
b(2) = sign (na, ia)

do i = 2, na + 2
  b(i+1) = a(i+1)
enddo

b(na+4) = 0.d0
b(na+5) = 0.d0

110 continue

return
end subroutine mpeq

subroutine mpinfr (a, b, c, mpnw)

!   Sets B to the integer part of the MPR number A and sets C equal to the
!   fractional part of A.  Note this is NOT the quite same as the greatest
!   integer function as often defined in some mathematical books and papers.
!   Examples:  If A = 1.95, then B = 1., C = 0.95.
!     If A = -3.25, then B = -3., C = -0.25.

implicit none
integer i, ia, ma, mpnw, na, nb, nc
double precision a(0:mpnw+5), b(0:mpnw+5), c(0:mpnw+5)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < mpnw + 6 .or. &
  c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPINFR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

!   Check if  A  is zero.

ia = sign (1.d0, a(2))
na = min (int (abs (a(2))), mpnw)
ma = a(3)
if (na == 0)  then
  b(1) = mpnw
  b(2) = 0.d0
  b(3) = 0.d0
  c(1) = mpnw
  c(2) = 0.d0
  c(3) = 0.d0
  goto 120
endif

if (ma >= mpnw - 1) then
  write (mpldb, 2)
2 format ('*** MPINFR: Argument is too large.')
  call mpabrt (40)
endif

!   Place integer part in  B.

nb = min (max (ma + 1, 0), na)
if (nb == 0) then
  b(1) = mpnw
  b(2) = 0.d0
  b(3) = 0.d0
else
  b(1) = mpnw
  b(2) = sign (nb, ia)
  b(3) = ma
  b(nb+4) = 0.d0
  b(nb+5) = 0.d0

  do i = 3, nb + 2
    b(i+1) = a(i+1)
  enddo
endif

!   Place fractional part in C.

nc = na - nb
if (nc <= 0) then
  c(1) = mpnw
  c(2) = 0.d0
  c(3) = 0.d0
else
  c(1) = mpnw
  c(2) = sign (nc, ia)
  c(3) = ma - nb
  c(nc+4) = 0.d0
  c(nc+5) = 0.d0

  do i = 3, nc + 2
    c(i+1) = a(i+nb+1)
  enddo
endif

!   Fix up results.  B may have trailing zeros and C may have leading zeros.

call mproun (b, mpnw)
call mproun (c, mpnw)

120  continue
return
end subroutine mpinfr

subroutine mpmdc (a, b, n, mpnw)

!   This returns a DP approximation the MPR number A in the form B * 2^n.

implicit none
integer i, mpnw, n, na
double precision aa, b
double precision a(0:mpnw+5)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4) then
  write (mpldb, 1)
1 format ('*** MPMDC: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

if (a(2) == 0.d0)  then
  b = 0.d0
  n = 0
  goto 100
endif

na = abs (a(2))
aa = a(4)
if (na >= 2) aa = aa + mprdx * a(5)
if (na >= 3) aa = aa + mprx2 * a(6)
if (na >= 4) aa = aa + mprdx * mprx2 * a(7)

n = mpnbt * a(3)
b = sign (aa, dble (a(2)))

!   Reduce b to within 1 and 2.

na = log (abs (b)) / log (2.d0) + mprdx
b = b / 2.d0**na
n = n + na
if (abs (b) < 1.d0) then
  b = 2.d0 * b
  n = n - 1
elseif (abs (b) > 2.d0) then
  b = 0.5d0 * b
  n = n + 1
endif

100  continue
return
end subroutine mpmdc

subroutine mpmul (a, b, c, mpnw)

!   This routine multiplies MPR numbers A and B to yield C.

!   This routine returns up to MPNW mantissa words of the product.  If the
!   complete double-long product of A and B is desired (for example in large
!   integer applications), then MPNW must be at least as large as the sum of
!   the mantissa lengths of A and B.  In other words, if the precision levels
!   of A and B are both 64 words, then MPNW must be at least 128 words to
!   produce the complete double-long product in C.

implicit none
integer i, i1, i2, j, j3, ia, ib, mpnw, mpnwx, na, nb, nc, n2
parameter (mpnwx = 200)
double precision a1, a2, c1, c2, dc, d1, d2, t1, t2, t3
double precision a(0:mpnw+5), b(0:mpnw+5), c(0:mpnw+5), d(0:mpnw+5), &
  b1(0:mpnw+5), b2(0:mpnw+5)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < abs (b(2)) + 4 .or. &
  c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPMUL: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

ia = sign (1.d0, a(2))
ib = sign (1.d0, b(2))
na = min (int (abs (a(2))), mpnw)
nb = min (int (abs (b(2))), mpnw)
nc = min (na + nb, mpnw)

if (na == 0 .or. nb == 0) then

!   One of the inputs is zero -- result is zero.

  c(1) = mpnw
  c(2) = 0.d0
  c(3) = 0.d0
  goto 200
endif

if (na == 1 .and. a(4) == 1.d0) then

!   A is 1 or -1 -- result is B or -B.

  c(1) = mpnw
  c(2) = sign (nb, ia * ib)
  c(3) = a(3) + b(3)

  do i = 3, nb + 2
    c(i+1) = b(i+1)
  enddo

  goto 200
elseif (nb == 1 .and. b(4) == 1.d0) then

!   B is 1 or -1 -- result is A or -A.

  c(1) = mpnw
  c(2) = sign (na, ia * ib)
  c(3) = a(3) + b(3)

  do i = 3, na + 2
    c(i+1) = a(i+1)
  enddo

  goto 200
endif

if (na > mpnwx .and. nb > mpnwx) then

!   Precision levels of both arguments are higher than mpnwx, so call mpmulx.

  call mpmulx (a, b, c, mpnw)
  goto 200
endif

d2 = a(3) + b(3)
d(0) = mpnw + 6
d(1) = mpnw

do i = 1, nc + 4
  d(i+1) = 0.d0
enddo

do i = 0, nb + 4
  b1(i) = mpb24x * aint (mpr24x * b(i+1))
  b2(i) = b(i+1) - b1(i)
enddo

!   Perform ordinary long multiplication algorithm.  Accumulate at most MPNW+4
!   mantissa words of the product.

do j = 3, na + 2
  a1 = mpb24x * aint (mpr24x * a(j+1))
  a2 = a(j+1) - a1
  j3 = j - 3
  n2 = min (nb + 2, mpnw + 4 - j3)

  do i = 3, n2
    dc = a1 * b2(i) + a2 * b1(i)
    c1 = mpbdx * aint (mprdx * dc)
    c2 = dc - c1
    d(i+j3) = d(i+j3) + mprdx * (a1 * b1(i) + c1)    
    d(i+j3+1) = d(i+j3+1) + a2 * b2(i) + c2
  enddo

!   Release carries periodically to avoid overflowing the exact integer
!   capacity of double precision floating point words in D.

  if (mod (j - 2, mpnpr) == 0) then
    i1 = max (3, j - mpnpr)
    i2 = n2 + j3

    do i = i1, i2
      t1 = d(i+1)
      t2 = int (mprdx * t1)
      d(i+1) = t1 - mpbdx * t2
      d(i) = d(i) + t2
    enddo
  endif
enddo

!   If D(3) is nonzero, shift the result one cell right.

if (d(3) /= 0.d0) then
  d2 = d2 + 1.d0

  do i = nc + 4, 3, -1
    d(i+1) = d(i)
  enddo
endif
d(2) = sign (nc, ia * ib)
d(3) = d2

!   Fix up result, since some words may be negative or exceed MPBDX.

call mpnorm (d, c, mpnw)

200 continue

return
end subroutine mpmul

subroutine mpmuld (a, b, c, mpnw)

!   This routine multiplies the MPR number A by the DP number B to yield C.

!   Note, however, that if B = 0.1D0 for example (or any other value that 
!   is not either a whole number or exact binary fraction), then C will NOT
!   be the true multiprecision product A * B.  This is because the double
!   precision value 0.1d0 is only a 15-digit approximation of 1/10, and thus
!   the product C will only be good to 15 digits or so.

!   Examples of exact binary values (good): 35.0, 395.5, 0.125, -5.3125.
!   Examples of inexact binary values (bad):  0.1, 0.8, -2.95, 33.3.

implicit none
integer i, ia, ib, k, mpnw, na, n1
double precision a1, a2, b, bb, b1, b2, c1, c2, dc
double precision a(0:mpnw+5), c(0:mpnw+5), d(0:mpnw+5)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPMULD: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

!   Check for zero inputs.

ia = sign (1.d0, a(2))
na = min (int (abs (a(2))), mpnw)
ib = sign (1.d0, b)
if (na == 0 .or. b == 0.d0) then
  c(1) = mpnw
  c(2) = 0.d0
  c(3) = 0.d0
  goto 140
endif
bb = abs (b)
n1 = 0

!   Reduce BB to within 1 and MPBDX.

if (bb >= mpbdx) then
  do k = 1, 100
    bb = mprdx * bb
    if (bb < mpbdx) then
      n1 = n1 + k
      goto 120
    endif
  enddo
elseif (bb < 1.d0) then
  do k = 1, 100
    bb = mpbdx * bb
    if (bb >= 1.d0) then
      n1 = n1 - k
      goto 120
    endif
  enddo
endif

!   If B cannot be represented exactly in a single mantissa word, use MPMUL.

120  continue

if (bb /= aint (bb)) then
  bb = sign (bb, b)
  d(0) = mpnw + 6
  d(1) = mpnw
  call mpdmc (bb, n1 * mpnbt, d, mpnw)
  call mpmul (a, d, c, mpnw)
  goto 140
endif

!   Perform short multiply operation.

b1 = mpb24x * aint (mpr24x * bb)
b2 = bb - b1
d(0) = mpnw + 6
d(1) = mpnw
d(2) = sign (na + 1, ia * ib)

do i = 2, min (na + 5, mpnw + 4)
  d(i+1) = 0.d0
enddo

do i = 3, na + 2
    a1 = mpb24x * aint (mpr24x * a(i+1))
    a2 = a(i+1) - a1
    dc = a1 * b2 + a2 * b1
    c1 = mpbdx * aint (mprdx * dc)
    c2 = dc - c1
    d(i) = d(i) + mprdx * (a1 * b1 + c1)
    d(i+1) = a2 * b2 + c2
enddo

!   If d(3) is nonzero, shift all words right by one.

if (d(3) > 0.d0) then
  do i = na + 3, 3, -1
    d(i+1) = d(i)
  enddo
  
  d(3) = a(3) + n1 + 1
else
  d(3) = a(3) + n1
endif

!   Fix up the result.

call mpnorm (d, c, mpnw)

d(3) = a(3) + n1

140 continue

return
end subroutine mpmuld

subroutine mpmuld40 (a, b, c, mpnw)

!   This routine multiples the MP number A by the DP number B to yield C.
!   In contrast to mpmuld, this routine only allows 40 significant bits
!   (approximately 12 significant decimal digits) in B.  If more nonzero bits
!   are present, an error is flagged.

!   Examples of exact binary values (good): 123456789.d0, 0.25d0, -5.3125d0.
!   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.

implicit none
integer mpnw
double precision a(0:mpnw+5), b, c(0:mpnw+5), t1, t2

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. c(0) < mpnw + 6) then
 write (mpldb, 1)
1 format ('*** MPMULD40: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

!   This convoluted-looking code tests whether B has more than 40
!   significant bits (actually whether the trailing 13 bits are zero).

t1 = mpb13x * abs (b)
t2 = abs (abs (b) + t1) - abs (t1)
if (t2 == abs (b)) then
  call mpmuld (a, b, c, mpnw)
else
  write (mpldb, 2) b
2 format ('*** MPMULD40: DP value has more than 40 significant bits:', &
  1p,d25.15/'and thus very likely represents an unintended loss of accuracy.'/ &
  'Fix the issue, or else use functions mpprodd, mpquotd, mpreald or mpcmplxdc.'/ &
  'See documentation for details.')
  call mpabrt (83)
endif

return
end subroutine mpmuld40

subroutine mpnint (a, b, mpnw)

!   This sets B to the nearest integer to the MPR number A.
!   Examples:  If A = 1.49, B = 1.; if A = 3.5, B = 4; if A = -2.5, B = -3.

implicit none
integer i, ia, ic, k0, k1, ma, na, mpnw
double precision a(0:mpnw+5), b(0:mpnw+5), s0(0:mpnw+5), s1(0:mpnw+5)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPNINT: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

ia = sign (1.d0, a(2))
na = min (int (abs (a(2))), mpnw)
ma = a(3)
if (na == 0)  then

!   A is zero -- result is zero.

  b(1) = mpnw
  b(2) = 0.d0
  b(3) = 0.d0
  goto 110
endif

if (ma >= mpnw) then

!   A cannot be represented exactly as an integer.

  write (mpldb, 2)
2 format ('*** MPNINT: Argument is too large.')
  call mpabrt (56)
endif

!   Add or subtract 1/2 from the input, depending on its sign, then
!   return the greatest integer.

s0(0) = mpnw + 6
s1(0) = mpnw + 6

call mpdmc (0.5d0, 0, s0, mpnw)
if (ia == 1) then
  call mpadd (a, s0, s1, mpnw)
else
  call mpsub (a, s0, s1, mpnw)
endif
call mpinfr (s1, b, s0, mpnw)

110 continue
return
end subroutine mpnint

subroutine mpnorm (d, a, mpnw)

!   This converts the MP number in array D to the standard normalized form
!   in A.

!   MPNORM assumes that two extra mantissa words are input at the end of D.
!   This reduces precision loss when it is necessary to shift the result to
!   the left. All words up to index A(2) + 5 in A *must* have data, even if 0.

implicit none
integer i, ia, k, mpnw, na, n4
double precision a2, s1, t1, t2, t3
double precision d(0:mpnw+5), a(0:mpnw+5)

! End of declaration

if (mpnw < 4 .or. d(0) < abs (d(2)) + 4 .or. a(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPNORM: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

ia = sign (1.d0, d(2))
na = min (int (abs (d(2))), mpnw)
if (na == 0)  then
  a(1) = mpnw
  a(2) = 0.d0
  a(3) = 0.d0
  goto 170
endif
n4 = na + 4
a2 = d(3)
d(3) = 0.d0

110 continue

t1 = 0.d0

do i = n4, 3, -1
  t3 = t1 + d(i+1)
  t2 = mprdx * (t3)
  t1 = int (t2)
  if (t2 < 0.d0 .and. t1 /= t2) t1 = t1 - 1.d0
  d(i+1) = t3 - t1 * mpbdx
enddo

d(3) = d(3) + t1

140  continue

if (d(3) < 0.d0) then

!   D(2) is negative -- negate all words and re-normalize.

  ia = - ia
  d(4) = d(4) + mpbdx * d(3)
  d(3) = 0.d0

  do i = 2, n4
    d(i+1) = - d(i+1)
  enddo

  goto 110
elseif (d(3) > 0.d0) then

!   The fixup loops above "spilled" a nonzero number into D(2).  Shift the
!   entire number right one cell.  The exponent and length of the result
!   are increased by one.

  do i = n4, 3, -1
    a(i+1) = d(i)
  enddo

  na = min (na + 1, mpnw)
  a2 = a2 + 1.
else
  do i = 3, n4
    a(i+1) = d(i+1)
  enddo
endif

!   Perform rounding and truncation.

a(1) = mpnw
a(2) = sign (na, ia)
a(3) = a2

call mproun (a, mpnw)

170 continue
return
end subroutine mpnorm

subroutine mpnpwr (a, n, b, mpnw)

!   This computes the N-th power of the MPR number A and returns the result
!   in B.  When N is zero, 1 is returned.  When N is negative, the reciprocal
!   of A ^ |N| is returned.

implicit none
integer i, j, k, kk, kn, k0, k1, k2, mn, mpnw, mpnw1, n, na, nn, nws
double precision cl2, t1, mprxx
parameter (cl2 = 1.4426950408889633d0, mprxx = 1d-14)
double precision a(0:mpnw+5), b(0:mpnw+5), s0(0:mpnw+6), s1(0:mpnw+6), &
  s2(0:mpnw+6)
  
! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPNPWR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

na = min (int (abs (a(2))), mpnw)
if (na == 0) then
  if (n >= 0) then
  	b(1) = mpnw
    b(2) = 0.d0
    b(3) = 0.d0
    goto 120
  else
    write (mpldb, 2)
2   format ('*** MPNPWR: Argument is zero and N is negative or zero.')
    call mpabrt (57)
  endif
endif

mpnw1 = mpnw + 1
s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7

nn = abs (n)
if (nn == 0) then
  call mpdmc (1.d0, 0, b, mpnw)
  goto 120
elseif (nn == 1) then
  call mpeq (a, s2, mpnw1)
  goto 110
elseif (nn == 2) then
  call mpmul (a, a, s2, mpnw1) 
  goto 110
endif

!   Determine the least integer MN such that 2 ^ MN .GT. NN.

t1 = nn
mn = cl2 * log (t1) + 1.d0 + mprxx
call mpdmc (1.d0, 0, s2, mpnw1)
call mpeq (a, s0, mpnw1)
kn = nn

!   Compute B ^ N using the binary rule for exponentiation.

do j = 1, mn
  kk = kn / 2
  if (kn /= 2 * kk) then
    call mpmul (s2, s0, s1, mpnw1)
    call mpeq (s1, s2, mpnw1)
  endif
  kn = kk
  if (j < mn) then
    call mpmul (s0, s0, s1, mpnw1)
    call mpeq (s1, s0, mpnw1)
  endif
enddo

!   Compute reciprocal if N is negative.

110 continue

if (n < 0) then
  call mpdmc (1.d0, 0, s1, mpnw1)
  call mpdiv (s1, s2, s0, mpnw1)
  call mpeq (s0, s2, mpnw1)
endif

!   Restore original precision level.

call mproun (s2, mpnw)
call mpeq (s2, b, mpnw)

120 continue
return
end subroutine mpnpwr

subroutine mpnrtr (a, n, b, mpnw)

!   This computes the N-th root of the MPR number A and returns result inB.
!   N must be at least one and must not exceed 2 ^ 30.

!   This subroutine employs the following Newton-Raphson iteration, which
!   converges to A ^ (-1/N):

!    X_{k+1} = X_k + (X_k / N) * (1 - A * X_k^N)

!   The reciprocal of the final approximation to A ^ (-1/N) is the N-th root.
!   These iterations are performed with a maximum precision level MPNW that
!   is dynamically changed, approximately doubling with each iteration.

!   When N is large and A is very near one, the following binomial series is
!   employed instead of the Newton scheme:

!   (1 + x)^(1/N)  =  1  +  x / N  +  x^2 * (1 - N) / (2! N^2)  +  ...

!   See the comment about the parameter NIT in MPDIVX.

implicit none
integer i, ia, iq, k, k0, k1, k2, k3, mpnw, mpnw1, mq, n, na, nit, &
  n1, n2, n3, n5, n30
double precision alt, cl2, t1, t2, tn, mprxx
parameter (alt = 0.693147180559945309d0, cl2 = 1.4426950408889633d0, &
  nit = 3, n30 = 2 ** 30, mprxx = 1d-14)
double precision a(0:mpnw+5), b(0:mpnw+5), f1(0:8), s0(0:mpnw+6), &
  s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPNRTR: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

ia = sign (1.d0, a(2))
na = min (int (abs (a(2))), mpnw)

if (na == 0) then
  b(2) = 0.
  b(3) = 0.
  goto 140
endif
if (ia < 0) then
  write (mpldb, 2)
2 format ('*** MPNRTR: Argument is negative.')
  call mpabrt (59)
endif

if (n <= 0 .or. n > n30) then
  write (mpldb, 3) n
3 format ('*** MPNRTR: Improper value of N',i10)
  call mpabrt (60)
endif

!   If N = 1 or 2, call MPEQ or MPSQRT instead.

if (n == 1) then
  call mpeq (a, b, mpnw)
  goto 140
elseif (n == 2) then
  call mpsqrt (a, b, mpnw)
  goto 140
endif

s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7
mpnw1 = mpnw + 1

!   Set f1 = 1.

f1(0) = 9.d0
f1(1) = mpnw1
f1(2) = 1.d0
f1(3) = 0.d0
f1(4) = 1.d0
f1(5) = 0.d0
f1(6) = 0.d0

!   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.

t1 = mpnw
mq = cl2 * log (t1) + 1.d0 - mprxx

!   Check how close A is to 1.

call mpsub (a, f1, s0, mpnw1)
if (s0(2) == 0.d0) then
  call mpeq (f1, b, mpnw)
  goto 140
endif
call mpmdc (s0, t1, n1, mpnw1)
n2 = cl2 * log (abs (t1))
t1 = t1 * 0.5d0 ** n2
n1 = n1 + n2

if (n1 <= -30) then
  t2 = n
  n2 = cl2 * log (t2) + 1.d0 + mprxx
  n3 = - mpnbt * mpnw1 / n1
  if (n3 < 1.25d0 * n2) then

!   A is so close to 1 that it is cheaper to use the binomial series.

    call mpdivd (s0, t2, s1, mpnw1)
    call mpadd (f1, s1, s2, mpnw1)
    k = 0

100 continue

	k = k + 1
    t1 = 1 - k * n
    t2 = (k + 1) * n
    call mpmuld (s1, t1, s3, mpnw1)
    call mpdivd (s3, t2, s1, mpnw1)
    call mpmul (s0, s1, s3, mpnw1)
    call mpeq (s3, s1, mpnw1)
    call mpadd (s1, s2, s3, mpnw1)
    call mpeq (s3, s2, mpnw1)
    if (s1(2) /= 0.d0 .and. s1(3) >= - mpnw1) then
    	  goto 100
    else
      call mpeq (s2, s1, mpnw1)
      goto 130
    endif
  endif
endif

!   Compute the initial approximation of A ^ (-1/N).

tn = n
call mpmdc (a, t1, n1, mpnw1)
n2 = - n1 / tn
t2 = exp (-1.d0 / tn * (log (t1) + (n1 + tn * n2) * alt))
call mpdmc (t2, n2, s2, mpnw1)
mpnw1 = 5
iq = 0

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (one greater than powers of two).

do k = 1, mq
  if (k > 2) mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1

110  continue

  call mpnpwr (s2, n, s0, mpnw1)
  call mpmul (a, s0, s1, mpnw1)
  call mpsub (f1, s1, s0, mpnw1)
  call mpmul (s2, s0, s1, mpnw1)
  call mpdivd (s1, tn, s0, mpnw1)
  call mpadd (s2, s0, s1, mpnw1)
  call mpeq (s1, s2, mpnw1)
  if (k == mq - nit .and. iq == 0) then
    iq = 1
    goto 110
  endif
enddo

!   Take the reciprocal to give final result.

call mpdiv (f1, s2, s1, mpnw)

!   Restore original precision level.

130 continue

call mproun (s1, mpnw)
call mpeq (s1, b, mpnw)

140 continue
return
end subroutine mpnrtr

subroutine mpoutw (iu, anam, a, mpnw)

!   This outputs the words of A up to the end of the active mantissa.
!   This is for internal debugging only; it should not be called by user.

implicit none
integer i, iu, mpnw, na
character*(*) anam 
double precision a(0:mpnw+5)

! End of declaration

na = min (int (abs (a(2))), mpnw)
write (iu, '(a)') anam
write (iu, '(4f18.0)') (a(i), i = 0, na + 5)
return
end subroutine mpoutw

subroutine mproun (a, mpnw)

!   This performs rounding and truncation of the MPR number A.  It is called
!   by MPNORM, and also by other subroutines when the precision level is
!   modified.  It is not intended to be directly called by the user.
!   The parameter MPEXPMX is the absolute value of the largest exponent word
!   allowed for MP numbers (see system parameters at start of this module).

implicit none
integer i, ia, k, mpnw, na, n4
double precision a(0:mpnw+5), a2

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. a(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPROUN: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

!   Check for initial zeroes.

a2 = a(3)
a(3) = 0.d0
ia = sign (1.d0, a(2))
na = min (int (abs (a(2))), mpnw)
n4 = na + 4

if (a(4) == 0.d0) then

!   Find the first nonzero word and shift the entire number left.  The length
!   of the result is reduced by the length of the shift.

  do i = 4, n4
    if (a(i+1) /= 0.d0) goto 110
  enddo

  a(2) = 0.d0
  a(3) = 0.d0
  goto 170

110 continue

  k = i - 3

  do i = 3, n4 - k
    a(i+1) = a(i+k+1)
  enddo

  a2 = a2 - k
  na = na - max (k - 2, 0)
  if (k == 2) a(na+4) = 0.d0
endif

!   Perform rounding.

if (na == mpnw) then
  if (a(na+4) >= 0.5d0 * mpbdx) a(na+3) = a(na+3) + 1.d0

!   Release carries as far as necessary due to rounding.

  do i = na + 2, 3, -1
    if (a(i+1) < mpbdx) goto 140
    a(i+1) = a(i+1) - mpbdx
    a(i) = a(i) + 1.d0
  enddo

!   Release of carries due to rounding continued all the way to the start --
!   i.e. number was entirely 9's.

  a(4) = a(3)
  na = 1
  a2 = a2 + 1.d0
endif

140 continue

  if (a(na+3) == 0.d0) then

!   At least the last mantissa word is zero.  Find the last nonzero word
!   and adjust the length of the result accordingly.

  do i = na + 2, 3, -1
    if (a(i+1) /= 0.d0) goto 160
  enddo

  a(2) = 0.d0
  a(3) = 0.d0
  goto 170

160  continue

  na = i - 2
endif

!   Check for overflow and underflow.

if (a2 < - mpexpmx) then
  write (mpldb, 2)
2 format ('*** MPROUN: Exponent underflow.')
  call mpabrt (68)
elseif (a2 > mpexpmx) then
  write (mpldb, 3)
3 format ('*** MPROUN: Exponent overflow.')
  call mpabrt (69)
endif

!   Check for zero.

if (a(4) == 0.d0) then
  a(1) = mpnw
  a(2) = 0.d0
  a(3) = 0.d0
else
  a(1) = mpnw
  a(2) = sign (na, ia)
  a(3) = a2
  a(na+4) = 0.d0
  a(na+5) = 0.d0
endif

170  continue

return
end subroutine mproun

subroutine mpsqrt (a, b, mpnw)

!   This computes the square root of the MPR number A and returns the result in B.

!   This subroutine employs the following Newton-Raphson iteration, which
!   converges to 1 / Sqrt(A):

!    X_{k+1} = X_k + 0.5 * (1 - X_k^2 * A) * X_k

!   where the multiplication () * X_k is performed with only half of the
!   normal level of precision.  These iterations are performed with a
!   working precision level MPNW that is dynamically changed, approximately 
!   doubling with each iteration (except that at iteration NIT before the final 
!   iteration, the iteration is repeated without doubling the precision, in order to
!   enhance accuracy) .  The final iteration is performed as follows
!   (this is due to A. Karp):

!    Sqrt(A) = (A * X_n) + 0.5 * [A - (A * X_n)^2] * X_n  (approx.)

!   where the multiplications A * X_n and [] * X_n are performed with only
!   half of the final level of precision.

implicit none
integer i, ia, iq, k, k0, k1, k2, k3, mpnw, mpnw1, mq, n, na, nit, nws, &
  nw1, nw2, n2, n7
double precision cl2, t1, t2, mprxx
parameter (cl2 = 1.4426950408889633d0, mprxx = 1d-14,  nit = 3)
double precision a(0:mpnw+5), b(0:mpnw+5), s0(0:mpnw+6), s1(0:mpnw+6), &
  s2(0:mpnw+6), s3(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPSQRT: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

ia = sign (1.d0, a(2))
na = min (int (abs (a(2))), mpnw)

if (na == 0) then
  b(1) = mpnw
  b(2) = 0.d0
  b(3) = 0.d0
  goto 120
endif
if (ia < 0.d0) then
  write (mpldb, 2)
2 format ('*** MPSQRT: Argument is negative.')
  call mpabrt (70)
  return
endif

s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7

!   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.

t1 = mpnw
mq = cl2 * log (t1) + 1.d0 - mprxx

!   Compute the initial approximation of 1 / Sqrt(A).

call mpmdc (a, t1, n, mpnw)
n2 = - n / 2
t2 = sqrt (t1 * 2.d0 ** (n + 2 * n2))
t1 = 1.d0 / t2
call mpdmc (t1, n2, s2, mpnw)
call mpdmc (1.d0, 0, s3, mpnw)

mpnw1 = 5
iq = 0
nw1 = mpnw1
nw2 = mpnw1

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (one greater than powers of two).

do k = 1, mq - 1
  if (k > 2) then
    nw1 = mpnw1
    mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
    nw2 = mpnw1
  endif
  
100  continue

  call mpmul (s2, s2, s0, nw2)
  call mpmul (a, s0, s1, nw2)
  call mpsub (s3, s1, s0, nw2)
  call mpmul (s2, s0, s1, nw1)
  call mpmuld (s1, 0.5d0, s0, nw1)  
  call mpadd (s2, s0, s1, nw2)
  call mpeq (s1, s2, nw2)

  if (k == mq - nit .and. iq == 0) then
    iq = 1
    goto 100
  endif
enddo

!   Perform last iteration using Karp's trick.

nw1 = mpnw1
mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
nw2 = mpnw1

call mpmul (a, s2, s0, nw1)
call mpmul (s0, s0, s1, nw2)
call mpsub (a, s1, s3, nw2)
call mpmul (s3, s2, s1, nw1)
call mpmuld (s1, 0.5d0, s3, nw1)
call mpadd (s0, s3, s2, nw2)

!   Restore original precision level.

call mproun (s2, mpnw)
call mpeq (s2, b, mpnw)

120 continue
return
end subroutine mpsqrt

subroutine mpsub (a, b, c, mpnw)

!   This routine subtracts MPR numbers A and B to yield C.

implicit none
integer i, nb, mpnw
double precision a(0:mpnw+5), b(0:mpnw+5), c(0:mpnw+5), s(0:mpnw+5)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < abs (b(2)) + 4 .or. &
  c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPSUB: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

nb = min (abs (int (b(2))), mpnw)
s(0) = mpnw + 6
s(1) = mpnw
if (b(2) == 0.d0) then
  s(2) = 0.d0
elseif (b(2) > 0.d0) then
  s(2) = - nb
else
  s(2) = nb
endif

do i = 3, nb + 5
  s(i) = b(i)
enddo

call mpadd (a, s, c, mpnw)

return
end subroutine mpsub

! ***  The following are the extra-high precision routines:

subroutine mpdivx (a, b, c, mpnw)

!   This divides A by B and returns the result in C.

!   This subroutine employs the following Newton-Raphson iteration, which
!   converges to 1 / B:

!    X_{k+1} = X_k + (1 - X_k * B) * X_k

!   where the multiplication () * X_k is performed with only half of the
!   normal level of precision.  These iterations are performed with a
!   working precision level MPNW that is dynamically changed, approximately 
!   doubling with each iteration (except that at iteration NIT before the
!   final iteration, the iteration is repeated without doubling the
!   precision, in order to enhance accuracy).  The final iteration is
!   performed as follows (this is due to A. Karp):

!    A / B = (A * X_n) + [A - (A * X_n) * B] * X_n  (approx.)

!   where the multiplications A * X_n and [] * X_n are performed with only
!   half of the final level of precision.

implicit none
integer i, ia, ib, iq, k, k0, k1, k2, k3, mpnw, mpnw1, mq, n, na, nb, &
  nit, nws, nw1, nw2, n2, n7
double precision cl2, t1, t2, mprxx
parameter (cl2 = 1.4426950408889633d0, mprxx = 1d-14,  nit = 3)
double precision a(0:mpnw+5), b(0:mpnw+5), c(0:mpnw+5), s0(0:mpnw+6), &
  s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < abs (b(2)) + 4 .or. &
  c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPDIVX: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

ia = sign (1.d0, a(2))
ib = sign (1.d0, b(2))
na = min (int (abs (a(2))), mpnw)
nb = min (int (abs (b(2))), mpnw)

if (na == 0) then
  b(1) = mpnw
  b(2) = 0.d0
  b(3) = 0.d0
  goto 120
endif
if (nb == 0.d0) then
  write (mpldb, 2)
2 format ('*** MPDIVX: Divisor is zero.')
  call mpabrt (33)
  return
endif

s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
s3(0) = mpnw + 7

!   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.

t1 = mpnw
mq = cl2 * log (t1) + 1.d0 - mprxx

!   Compute the initial approximation of 1 / B.

call mpmdc (b, t1, n, mpnw)
t2 = 1.d0 / t1
call mpdmc (t2, -n, s2, mpnw)
call mpdmc (1.d0, 0, s3, mpnw)

mpnw1 = mpnw + 1
mpnw1 = 5
iq = 0
nw1 = mpnw1
nw2 = mpnw1

!   Perform the Newton-Raphson iteration described above with a dynamically
!   changing precision level MPNW (one greater than powers of two).

do k = 1, mq - 1
  if (k > 2) then
    nw1 = mpnw1
    mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
    nw2 = mpnw1
  endif
  
100  continue

  call mpmul (b, s2, s1, nw2)
  call mpsub (s3, s1, s0, nw2)
  call mpmul (s2, s0, s1, nw1)
  call mpadd (s2, s1, s0, nw2)
  call mpeq (s0, s2, nw2)
  if (k == mq - nit .and. iq == 0) then
    iq = 1
    goto 100
  endif
enddo

!   Perform last iteration using Karp's trick.

nw1 = mpnw1
mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
nw2 = mpnw1

call mpmul (a, s2, s0, nw1)
call mpmul (s0, b, s1, nw2)
call mpsub (a, s1, s3, nw2)
call mpmul (s3, s2, s1, nw1)
call mpadd (s0, s1, s2, nw2)

!   Restore original precision level.

call mproun (s2, mpnw)
call mpeq (s2, c, mpnw)

120 continue
return
end subroutine mpdivx

subroutine mpfftcr (is, m, n, nsq, x, y)

!   This performs an N-point complex-to-real FFT, where N = 2^M.  X is the
!   double complex input array, and Y is the double precision output array.
!   The array X is used as a scratch array in MPFFT1, and so is overwritten.
!   X and Y must be dimensioned as shown below.  IS is the sign of the FFT.
!   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
!   This routine is not intended to be called directly by the user.

implicit none
integer i, is, k, ku, m, mx, n, nsq, n1, n2, n21, n4
real (mprknd) y(n)
complex (mprknd) dc1(n/2), x(n/2+nsq*mpnsp1+1), ai, a1, a2, x1, x2

mx = mpuu1(1)

!   Check if input parameters are invalid.

if ((is .ne. 1 .and. is .ne. -1) .or. m .lt. 3 .or. m .gt. mx) then
  write (mpldb, 1)  is, m, mx
1 format ('*** MPFFTCR: Either the UU arrays have not been initialized'/ &
  'or else one of the input parameters is invalid', 3i5)
  call mpabrt (677)
endif

n1 = 2 ** (m / 2)
n2 = n / 2
n21 = n2 + 1
n4 = n / 4
ai = cmplx (0.d0, 1.d0, mprknd)

!   Construct the input to MPFFT1.

dc1(1) = 0.5d0 * cmplx (real (x(1) + x(n2+1), mprknd), &
  real (x(1) - x(n2+1), mprknd), mprknd)
if (is .eq. 1) then
  dc1(n4+1) = conjg (x(n4+1))
else
  dc1(n4+1) = x(n4+1)
endif
ku = n2

if (is .eq. 1) then
  do k = 2, n4
    x1 = x(k)
    x2 = conjg (x(n2+2-k))
    a1 = x1 + x2
    a2 = ai * mpuu1(k+ku) * (x1 - x2)
    dc1(k) = 0.5d0 * (a1 + a2)
    dc1(n2+2-k) = 0.5d0 * conjg (a1 - a2)
  enddo
else
  do k = 2, n4
    x1 = x(k)
    x2 = conjg (x(n2+2-k))
    a1 = x1 + x2
    a2 = ai * conjg (mpuu1(k+ku)) * (x1 - x2)
    dc1(k) = 0.5d0 * (a1 + a2)
    dc1(n2+2-k) = 0.5d0 * conjg (a1 - a2)
  enddo
endif

!   Perform a normal N/2-point FFT on DC1.

call mpfft1 (is, m - 1, n1, n2 / n1, dc1, x)

!   Copy DC1 to Y such that DC1(k) = Y(2k-1) + i Y(2k).

do k = 1, n / 2
  y(2*k-1) = real (dc1(k), mprknd)
  y(2*k) = aimag (dc1(k))
enddo

return
end subroutine mpfftcr

subroutine mpfftrc (is, m, n, nsq, x, y)

!   This performs an N-point real-to-complex FFT, where N = 2^M.  X is the
!   dobule precision input array, and Y is the double complex output array.
!   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
!   This routine is not intended to be called directly by the user.

implicit none
integer i, is, k, ku, m, mx, n, nsq, n1, n2, n21, n4
real (mprknd) x(n)
complex (mprknd) dc1(n/2), y(n/2+nsq*mpnsp1+1), ai, a1, a2, z1, z2

mx = mpuu1(1)

!   Check if input parameters are invalid.

if ((is .ne. 1 .and. is .ne. -1) .or. m .lt. 3 .or. m .gt. mx) then
  write (mpldb, 1)  is, m, mx
1 format ('*** MPFFTRC: either the UU arrays have not been initialized'/ &
  'or else one of the input parameters is invalid',3i5)
  call mpabrt (677)
endif

n1 = 2 ** (m / 2)
n2 = n / 2
n21 = n2 + 1
n4 = n / 4
ai = cmplx (0.d0, -1.d0, mprknd)

!   Copy X to DC1 such that DC1(k) = X(2k-1) + i X(2k).

do k = 1, n2
  dc1(k) = cmplx (x(2*k-1), x(2*k), mprknd)
enddo

!   Perform a normal N/2-point FFT on DC1.

call mpfft1 (is, m - 1, n1, n2 / n1, dc1, y)

!   Reconstruct the FFT of X.

y(1) = cmplx (2.d0 * (real (dc1(1), mprknd) + aimag (dc1(1))), &
  0.d0, mprknd)
if (is .eq. 1) then
  y(n4+1) = 2.d0 * dc1(n4+1)
else
  y(n4+1) = 2.d0 * conjg (dc1(n4+1))
endif
y(n2+1) = cmplx (2.d0 * (real (dc1(1), mprknd) - aimag (dc1(1))), &
  0.d0, mprknd)
ku = n2

if (is .eq. 1) then
  do k = 2, n4
    z1 = dc1(k)
    z2 = conjg (dc1(n2+2-k))
    a1 = z1 + z2
    a2 = ai * mpuu1(k+ku) * (z1 - z2)
    y(k) = a1 + a2
    y(n2+2-k) = conjg (a1 - a2)
  enddo
else
  do k = 2, n4
    z1 = dc1(k)
    z2 = conjg (dc1(n2+2-k))
    a1 = z1 + z2
    a2 = ai * conjg (mpuu1(k+ku)) * (z1 - z2)
    y(k) = a1 + a2
    y(n2+2-k) = conjg (a1 - a2)
  enddo
endif

return
end subroutine mpfftrc

subroutine mpfft1 (is, m, n1, n2, x, y)

!   This routine performs a complex-to-complex FFT.  IS is the sign of the
!   transform, N = 2^M is the size of the transform.  N1 = 2^M1 and N2 = 2^M2,
!   where M1 and M2 are defined as below.  X is the input and output array,
!   and Y is a scratch array.  X must have at N, and Y at least N + N1*MPNSP1,
!   double complex cells.  The arrays MPUU1 and MPUU2 must have been
!   initialized by calling MPINIFFT.  This routine is not intended to be called
!   directly by the user.

!   This employs the two-pass variant of the "four-step" FFT.  See the
!   article by David H. Bailey in J. of Supercomputing, March 1990, p. 23-35.

implicit none
integer i, is, iu, j, j2, k, ku, m, m1, m2, n, n1, n2, nr1, nr2
complex (mprknd) x(n1,n2), y(n2+mpnsp1,n1), z1(mpnrow+mpnsp1,n1), &
  z2(mpnrow+mpnsp1,n1)

n = 2 ** m
m1 = (m + 1) / 2
m2 = m - m1
nr1 = min (n1, mpnrow)
nr2 = min (n2, mpnrow)
ku = mpuu2(m)

do i = 0, n1 - 1, nr1

!   Copy NR1 rows of X (treated as a N1 x N2 complex array) into Z1.

  do j = 1, n2
    do k = 1, nr1
      z1(k,j) = x(i+k,j)
    enddo
  enddo

!   Perform NR1 FFTs, each of length N2.

  call mpfft2 (is, nr1, m2, n2, z1, z2)

!   Multiply the resulting NR1 x N2 complex block by roots of unity and
!   store transposed into the appropriate section of Y.

  iu = i + ku - n1 - 1
  if (is .eq. 1) then
    do j = 1, n2
      do k = 1, nr1
        y(j,i+k) = mpuu2(iu+k+j*n1) * z1(k,j)
      enddo
    enddo
  else
    do j = 1, n2
      do k = 1, nr1
        y(j,i+k) = conjg (mpuu2(iu+k+j*n1)) * z1(k,j)
      enddo
    enddo
  endif
enddo

do i = 0, n2 - 1, nr2

!   Copy NR2 rows of the Y array into Z2.

  do j = 1, n1
    do k = 1, nr2
      z2(k,j) = y(i+k,j)
    enddo
  enddo

!   Perform NR2 FFTs, each of length N1.

  call mpfft2 (is, nr2, m1, n1, z2, z1)

!   Copy NR2 x N1 complex block back into X array.  It's a little more
!   complicated if M is odd.

  if (mod (m, 2) .eq. 0) then
    do j = 1, n1
      do k = 1, nr2
        x(i+k,j) = z2(k,j)
      enddo
    enddo
  else
    do j = 1, n1 / 2
      j2 = 2 * j - 1

      do k = 1, nr2
        x(i+k,j) = z2(k,j2)
        x(i+k+n2,j) = z2(k,j2+1)
      enddo
    enddo
  endif
enddo

return
end subroutine mpfft1

subroutine mpfft2 (is, ns, m, n, x, y)

!   This performs NS simultaneous N-point complex-to-complex FFTs, where
!   N = 2^M.  X is the input and output array, and Y is a scratch array.
!   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
!   This routine is not intended to be called directly by the user.

implicit none
integer i, is, j, l, m, n, ns
complex (mprknd) x(mpnrow+mpnsp1,n), y(mpnrow+mpnsp1,n)

!   Perform the second variant of the Stockham FFT.

do l = 1, m, 2
  call mpfft3 (is, l, ns, m, n, x, y)
  if (l .eq. m) goto 100
  call mpfft3 (is, l + 1, ns, m, n, y, x)
enddo

goto 110

!   Copy Y to X.

100 continue

do j = 1, n
  do i = 1, ns
    x(i,j) = y(i,j)
  enddo
enddo

110 continue

return
end subroutine mpfft2

subroutine mpfft3 (is, l, ns, m, n, x, y)

!   This performs the L-th iteration of the second variant of the Stockham FFT
!   on the NS vectors in X.  X is input/output, and Y is a scratch array.
!   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
!   This routine is not intended to be called directly by the user.

implicit none
integer i, is, i11, i12, i21, i22, j, k, l, li, lj, lk, ku, m, n, n1, ns
complex (mprknd) x(mpnrow+mpnsp1,n), y(mpnrow+mpnsp1,n), u1, x1, x2

!   Set initial parameters.

n1 = n / 2
lk = 2 ** (l - 1)
li = 2 ** (m - l)
lj = 2 * lk
ku = li + 1

do i = 0, li - 1
  i11 = i * lk + 1
  i12 = i11 + n1
  i21 = i * lj + 1
  i22 = i21 + lk
  if (is .eq. 1) then
    u1 = mpuu1(i+ku)
  else
    u1 = conjg (mpuu1(i+ku))
  endif

  do k = 0, lk - 1
    do j = 1, ns
      x1 = x(j,i11+k)
      x2 = x(j,i12+k)
      y(j,i21+k) = x1 + x2
      y(j,i22+k) = u1 * (x1 - x2)
    enddo
  enddo
enddo

return
end subroutine mpfft3

subroutine mpinifft (mpnw)

!   This computes the root of unity arrays UU1 and UU2, which are required by
!   the FFT routines, and places this data in the proper arrays defined in 
!   module MPFUNA.  MPNW is the largest precision level (in words) that will be 
!   subsequently used for this run.

implicit none
integer i, iu, j, k, ku, ln, m, mm, mm1, mm2, mq, nn, nn1, nn2, nq, mpnw, nwds
double precision cl2, d1, mprxx
parameter (cl2 = 1.4426950408889633d0, mprxx = 1d-14)
real (mprknd) pi, t1, ti, tpn

!  Determine sizes for FFT arrays.  Three words are added to mpnw, since many
!  routines in MPFUND in particular increase the working precision upon entry.

nwds = mpnw + 3
d1 = 1.5d0 * (nwds + 1)
m = cl2 * log (d1) + 1.d0 - mprxx
mq = m + 2
nq = 2 ** mq

if (mq + nq > mplfftx) then
  write (6, 1) mq + nq
1 format ('*** MPINIFFT: Insufficient space for arrays mpuu1 and mpuu2.'/ &
  'At least',i12,' double complex cells must be allocated for each of'/ &
  'these arrays in module mpfuna. See documentation for details.')
  call mpabrt (91)
endif

mpuu1(1) = mq
ku = 2
ln = 1
pi = acos (-1.d0)

do j = 1, mq
  t1 = pi / ln

  do i = 0, ln - 1
    ti = i * t1
    mpuu1(i+ku) = cmplx (cos (ti), sin (ti), mprknd)
  enddo

  ku = ku + ln
  ln = 2 * ln
enddo

! write (6, 2) ku - 1
! 2 format ('MPINIFFT: Size of table mpuu1 =',i10)

ku = mq + 1
mpuu2(1) = mq

do k = 2, mq
  mpuu2(k) = cmplx (0.d0, 0.d0, mprknd)
enddo

do k = 2, mq - 1
  mpuu2(k) = ku
  mm = k
  nn = 2 ** mm
  mm1 = (mm + 1) / 2
  mm2 = mm - mm1
  nn1 = 2 ** mm1
  nn2 = 2 ** mm2
  tpn = 2.d0 * pi / nn

  do j = 0, nn2 - 1
    do i = 0, nn1 - 1
      iu = ku + i + j * nn1
      t1 = tpn * i * j
      mpuu2(iu) = cmplx (cos (t1), sin (t1), mprknd)
    enddo
  enddo

  ku = ku + nn
enddo

! write (6, 3) ku - 1
! 3 format ('MPINIFFT: Size of table mpuu2 =',i10)

return
end subroutine mpinifft

subroutine mplconv (iq, n, nsq, a, b, c)

!   This computes the linear convolution of A and B, returning the result
!   in C.  If IQ is 1, then it is presumed B = A; if IQ = 2, then A /= B.
!   NSQ is a spacing parameter, which should be set to more than sqrt (3*n).

implicit none
integer i, iq, m1, m2, n, n1, n2, n4, nm, nsq
double precision cl2, c0, mprxx, mpffterrmx
parameter (cl2 = 1.4426950408889633d0, mprxx = 1d-14, mpffterrmx = 0.375d0)
real (mprknd) a(n), an, b(n), c(2*n), d1(6*n+2), d2(6*n+2), d3(6*n+2), t1, t2
complex (mprknd) dc1(3*n+nsq*mpnsp1+3), dc2(3*n+nsq*mpnsp1+3)

t1 = n
m1 = cl2 * log (t1) + 1.d0 - mprxx
n1 = 2 ** m1
m2 = m1 + 1
n2 = 2 * n1
n4 = 2 * n2
nm = min (2 * n, n2)

if (abs (iq) .eq. 1) then

!   Compute the square of a -- only one forward FFT is needed.

  do i = 1, n
    d1(i) = a(i)
  enddo

  do i = n + 1, n2
    d1(i) = 0.d0
  enddo

!   Perform a forward real-to-complex FFT on the vector in a.

  call mpfftrc (1, m2, n2, nsq, d1, dc1)

!   Square the resulting complex vector.

  do i = 1, n1 + 1
    dc1(i) = dc1(i) ** 2
  enddo
else

!   Compute the product of a and b -- two forward FFTs are needed.
  do i = 1, n
    d1(i) = a(i)
    d2(i) = b(i)
  enddo

  do i = n + 1, n2
    d1(i) = 0.d0
    d2(i) = 0.d0
  enddo

!   Perform forward real-to-complex FFTs on the vectors in a and b.

  call mpfftrc (1, m2, n2, nsq, d1, dc1)
  call mpfftrc (1, m2, n2, nsq, d2, dc2)

!   Multiply the resulting complex vectors.

  do i = 1, n1 + 1
    dc1(i) = dc1(i) * dc2(i)
  enddo
endif

!   Perform an inverse complex-to-real FFT on the resulting data.

call mpfftcr (-1, m2, n2, nsq, dc1, d3)

!   Divide by N4.

an = 1.d0 / n4
c0 = 0.d0

do i = 1, nm
  t1 = an * d3(i)
!  t2 = anint (t1)
  if (d3(i) >= 0.d0) then
    t2 = aint (t1 + 0.5d0)
  else
    t2 = aint (t1 - 0.5d0)
  endif
  c(i) = t2
  c0 = max (c0, abs (dble ((t2 - t1))))
enddo

mpffterr = max (c0, mpffterr)
if (c0 > mpffterrmx) then
  write (6, 1) c0
1 format ('*** MPLCONV: excessive rounding error =',f12.6)
  call mpabrt (55)
endif

return
end subroutine mplconv

subroutine mpmulx (a, b, c, mpnw)

!   This routine multiplies MP numbers A and B to yield the MP product C,
!   using a FFT-convolution technique.  Before calling MPMULX, the arrays
!   UU1 and UU2 must be initialized by calling MPINIFFT.  For modest levels
!   of precision, use MPMUL.

implicit none
integer i, i1, ia, ib, j, k, mpnw, na, nb, nc, nn, nx
double precision mprxx
parameter (mprxx = 1d-14)
double precision a(0:mpnw+5), b(0:mpnw+5), c(0:mpnw+5), d(0:mpnw+7), &
  t0, t1, t2, t3, t4, t5, r16, t16, r32, t32, r48, t48
parameter (r16 = 0.5d0**16, t16 = 2.d0**16, r32 = 0.5d0**32, t32 = 2.d0**32, &
  r48 = 0.5d0**48, t48 = 2.d0**48)
real (mprknd) d1(0:3*mpnw+16), d2(0:3*mpnw+16), d3(0:6*mpnw+32)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. b(0) < abs (b(2)) + 4 .or. &
  c(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPMULX: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

ia = sign (1.d0, a(2))
ib = sign (1.d0, b(2))
na = min (int (abs (a(2))), mpnw)
nb = min (int (abs (b(2))), mpnw)
nc = min (na + nb, mpnw)
nn = 3 * max (na, nb)
nx = sqrt (3.d0 * nn) + mprxx

!   Divide each word into three 16-bit chunks.

do i = 0, na - 1
  t1 = a(i+4)
  d1(3*i) = aint (r32 * t1)
  t1 = t1 - t32 * d1(3*i)
  d1(3*i+1) = aint (r16 * t1)
  d1(3*i+2) = t1 - t16 * d1(3*i+1)
enddo

do i = 3 * na, nn - 1
  d1(i) = 0.d0
enddo

!   If A is the same array as B, then there is no need to deal with B.

if (loc (a) == loc (b)) then
  call mplconv (1, nn, nx, d1, d2, d3)
else
  do i = 0, nb - 1
    t1 = b(i+4)
    d2(3*i) = aint (r32 * t1)
    t1 = t1 - t32 * d2(3*i)
    d2(3*i+1) = aint (r16 * t1)
    d2(3*i+2) = t1 - t16 * d2(3*i+1)
  enddo

  do i = 3 * nb, nn - 1
    d2(i) = 0.d0
  enddo

  call mplconv (2, nn, nx, d1, d2, d3)
endif

!   Release carries.

do i = 0, 3 * nc + 13
  d1(i) = 0.d0
enddo

do i = 0, min (3 * nc + 9, 2 * nn - 1)
  t0 = d3(i)
  t1 = aint (r48 * t0)
  t2 = t0 - t48 * t1
  t3 = aint (r32 * t2)
  t4 = t2 - t32 * t3
  t5 = aint (r16 * t4)
  d1(i+3) = t4 - t16 * t5
  d1(i+2) = d1(i+2) + t5
  d1(i+1) = d1(i+1) + t3
  d1(i) = d1(i) + t1
enddo

!  Recombine words, with proper offset.

d(0) = 0.d0
d(1) = 0.d0
d(2) = 0.d0
d(3) = 0.d0

do i = 0, nc + 3
  d(i+4) = t32 * d1(3*i+2) + t16 * d1(3*i+3) + d1(3*i+4)
enddo

d(0) = mpnw + 6
d(1) = mpnw
d(2) = sign (nc, ia * ib)
d(3) = a(3) + b(3) + 1

!   Fix up the result.

d1(0) = mpnw + 6
call mpnorm (d, c, mpnw)

190 continue

return
end subroutine mpmulx

!>   These two subroutines are for real*16 support.  See "Uncomment" below
!>   for differences.

subroutine mpmqc (a, b, n, mpnw)

!   This returns a quad (real*16) approximation the MPR number A in the form B * 2^n.

implicit none
integer i, knd, mpnw, n, na

!>  Uncomment this line if real*16 is supported.
 parameter (knd = kind (0.q0))
!>  Otherwise uncomment this line.
! parameter (knd = kind (0.d0))

real (knd) aa, b
double precision a(0:mpnw+5)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4) then
  write (mpldb, 1)
1 format ('*** MPMQC: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

if (a(2) == 0.d0)  then
  b = 0.d0
  n = 0
  goto 100
endif

na = abs (a(2))
aa = a(4)
if (na >= 2) aa = aa + mprdx * a(5)
if (na >= 3) aa = aa + mprx2 * a(6)
if (na >= 4) aa = aa + mprdx * mprx2 * a(7)

n = mpnbt * a(3)
b = sign (aa, real (a(2), knd))

!   Reduce b to within 1 and 2.

na = log (abs (b)) / log (2.d0) + mprdx
b = b / 2.d0**na
n = n + na
if (abs (b) < 1.d0) then
  b = 2.d0 * b
  n = n - 1
elseif (abs (b) > 2.d0) then
  b = 0.5d0 * b
  n = n + 1
endif

100  continue
return
end subroutine mpmqc

subroutine mpqmc (a, n, b, mpnw)

!   This routine converts the quad (real*16) number A * 2^N to MPR form in B.

!   NOTE however that if A = 0.1q0, for example, then B will NOT be the true 
!   multiprecision equivalent of 1/10, since 0.1q0 is not an exact binary value.

!   Examples of exact binary values (good): 123456789.d0, 0.25d0, -5.3125d0.
!   Examples of inexact binary values (bad): 0.1d0, 1234567.8d0, -3333.3d0.

implicit none
integer i, k, knd, mpnw, n, n1, n2

!>  Uncomment this line if real*16 is supported.
 parameter (knd = kind (0.q0))
!>  Otherwise uncomment this line.
! parameter (knd = kind (0.d0))

real (knd) a, aa
double precision b(0:*)

! End of declaration

if (mpnw < 4 .or. b(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPQMC: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

!   Check for zero.

if (a == 0.d0) then
  b(1) = mpnw
  b(2) = 0.d0
  b(3) = 0.d0
  goto 150
endif
n1 = n / mpnbt
n2 = n - mpnbt * n1
aa = abs (a) * 2.d0 ** n2

!   Reduce AA to within 1 and MPBDX.

if (aa >= mpbdx) then

  do k = 1, 100
    aa = mprdx * aa
    if (aa < mpbdx) then
      n1 = n1 + k
      goto 120
    endif
 enddo

elseif (aa < 1.d0) then

  do k = 1, 100
    aa = mpbdx * aa
    if (aa >= 1.d0) then
      n1 = n1 - k
      goto 120
    endif
  enddo

endif

!   Store successive sections of AA into B.

120  continue

b(3) = n1
b(4) = aint (aa)
aa = mpbdx * (aa - b(3+1))
b(5) = aint (aa)
aa = mpbdx * (aa - b(4+1))
b(6) = aint (aa)
aa = mpbdx * (aa - b(5+1))
b(7) = aint (aa)
b(8) = 0.d0
b(9) = 0.d0

do i = 7, 3, -1
  if (b(i+1) /= 0.d0) goto 140
enddo

140  continue

b(1) = mpnw
aa = i - 2
b(2) = sign (aa, a)

150 continue
return
end subroutine mpqmc

end module mpfunb

