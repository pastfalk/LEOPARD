!*****************************************************************************

!  MPFUN-Fort: A thread-safe arbitrary precision computation package
!  Binary-decimal, decimal-binary and I/O functions (module MPFUNC)

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

!  DESCRIPTION OF THIS MODULE (MPFUNC):
!    This module contains subroutines for binary-decimal and decimal-binary
!    conversion, together with low-level E-format and F-format conversion, and
!    basic input and output.
 
module mpfunc
use mpfuna
use mpfunb

contains

subroutine mpctomp (a, n, b, mpnw)

!  Converts the CHARACTER*1 array A of length N into the MPR number B.
!  Restrictions: (a) no embedded blanks; (b) a leading digit (possibly
!  zero) must be present; and (c) a period must be present.  An exponent
!  (with "d" or "e") may optionally follow the numeric value.
 
implicit none
integer i, i1, i2, iexp, isgn, iexpsgn, ix, j, kde, kend, kexpend, kexpst, &
  kexpsgn, knumend1, knumend2, knumst1, knumst2, kper, ksgn, kstart, lexp, &
  lexpmx, lnum, lnum1, lnum2, mpnw, mpnw1, n, n1, n2
double precision d10w, t1, t2
character*1 a(n), ai
character*10 digits
character*32 ca
parameter (lexpmx = 9, digits = '0123456789', d10w = 10.d0**mpndpw)
double precision b(0:mpnw+5), f(0:8), s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6)

! write (6, *) 'mpctomp: a, n, mpnw =', n, mpnw
! write (6, '(100a1)') 'X',(a(i),i=1,n),'X'

! End of declaration

if (mpnw < 4 .or. b(0) < mpnw + 6) then
 write (mpldb, 1)
1 format ('*** MPCTOMP: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
f(0) = 9.d0
f(1) = mpnw1

do i = 2, 8
  f(i) = 0.d0
enddo

mpnw1 = mpnw + 1
kde = 0
kend = 0
kexpend = 0
kexpst = 0
kexpsgn = 0
knumend1 = 0
knumend2 = 0
knumst1 = 0
knumst2 = 0
kper = 0
ksgn = 0
kstart = 0

!   Locate:
!     kstart = index of first nonblank character.
!     kend = index of last nonblank character.

do i = 1, n
  if (a(i) /= ' ') goto 100
enddo

!   Input is completely blank.

write (6, 2)
2 format ('*** MPCTOMP: Syntax error in input string.'/ &
  'Restrictions: (a) no embedded blanks; (b) a leading digit (possibly'/ &
  'zero) must be present; and (c) a period must be present.  An exponent'/ &
  '(with "d" or "e") may optionally follow the numeric value.')
call mpabrt (41)

100 continue

kstart = i

do i = n, kstart, -1
  if (a(i) /= ' ') goto 110
enddo

i = kstart

110 continue

kend = i

!   Scan input for:
!     kde = index of 'd' or 'e'.
!     kexpend = index of end of exponent.
!     kexpst = index of start of exponent.
!     kespsgn = index of sign of exponent.
!     knumend1 = index of end of numeric part prior to period.
!     knumend2 = index of end of numeric part after period.
!     knumst1 = index of start of numeric part prior to period.
!     knumst2 = index of start of numeric part after period.
!     kper = index of period.
!     ksgn = index of sign of number.

do i = kstart, kend
  if (a(i) == ' ') then
    write (6, 2)
    call mpabrt (41)
  elseif (a(i) == '+' .or. a(i) == '-') then
    if (i == kstart) then
      ksgn = i
    elseif (kde > 0 .and. kexpsgn == 0 .and. kexpst == 0 .and. i < kend) then
      kexpsgn = i
    else
      write (6, 2)
      call mpabrt (41)
    endif
  elseif (a(i) == 'e' .or. a(i) == 'E' .or. a(i) == 'd' .or. a(i) == 'D') then
    if (kde == 0 .and. kper > 0 .and. i < kend) then
      kde = i
      knumend2 = i - 1
    else
      write (6, 2)
      call mpabrt (41)
    endif
  elseif (a(i) == '.') then
    if (kper == 0 .and. kde == 0 .and. knumst1 > 0 .and. knumst2 == 0) then
      kper = i
      knumend1 = i - 1
    else
      write (6, 2)
      call mpabrt (41)
    endif
  elseif (index (digits, a(i)) > 0) then
    if (knumst1 == 0) then
      knumst1 = i
    elseif (kper > 0 .and. knumst2 == 0 .and. kde ==  0) then
      knumst2 = i
    elseif (kde > 0 .and. kexpst == 0) then
      kexpst = i
    endif
    if (i == kend) then
      if (knumst2 > 0 .and. kde == 0) then
        knumend2 = i
      elseif (kexpst > 0) then
        kexpend = i
      else
        write (6, 2)
        call mpabrt (41)
      endif
    endif
  else
    write (6, 2)
    call mpabrt (41)
  endif
enddo

! write (6, *) 'kde, kend, kexpend, kexpst =', kde, kend, kexpend, kexpst
! write (6, *) 'kexpsgn, numend1, knumend2, knumst1 =', kexpsgn, knumend1, knumend2, knumst1
! write (6, *) 'knumst2, kper, ksgn, kstart =', knumst2, kper, ksgn, kstart

!   Decode exponent.

if (kexpst > 0) then
  lexp = kexpend - kexpst + 1
  if (lexp > lexpmx) then
    write (6, 3)
3   format ('*** MPCTOMP: exponent string is too long.')
    call mpabrt (85)
  endif

  do i = 1, lexp
    ca(i:i) = a(i+kexpst-1)
  enddo

  iexp = mpdigin (ca, lexp)
  if (a(kexpsgn) == '-') iexp = -iexp
else
  iexp = 0
endif

!   Determine sign of number.

if (ksgn == 0 .or. a(ksgn) == '+') then
  isgn = 1
elseif (a(ksgn) == '-') then
  isgn = -1
endif

!   Determine lengths of two sections of number.

lnum1 = knumend1 - knumst1 + 1
if (knumst2 > 0) then
  lnum2 = knumend2 - knumst2 + 1
else
  lnum2 = 0
endif
lnum = lnum1 + lnum2

! write (6, *) 'iexp, lnum1, lnum2 =', iexp, lnum1, lnum2

!   Determine the number of chunks of digits and the left-over.

n1 = lnum / mpndpw
n2 = mod (lnum, mpndpw)

!   Construct first (left-over) portion, right-justified in CA.

ca(1:mpndpw - n2) = ' '
ix = knumst1 - 1

do i = 1, n2
  ix = ix + 1
  if (ix == kper) ix = ix + 1
  ca(i+mpndpw-n2:i+mpndpw-n2) = a(ix)
enddo

t1 = mpdigin (ca, mpndpw)
if (t1 > 0) then
  f(2) = 1.d0
  f(3) = 0.d0
  f(4) = t1
else
  f(2) = 0.d0
  f(3) = 0.d0
  f(4) = 0.d0
endif
call mpeq (f, s0, mpnw1)

!   Process remaining chunks of digits.

do j = 1, n1
  do i = 1, mpndpw
    ix = ix + 1
    if (ix == kper) ix = ix + 1
    ca(i:i) = a(ix)
  enddo
  
  t1 = mpdigin (ca, mpndpw)
  if (t1 > 0) then
    f(2) = 1.d0
    f(3) = 0.d0
    f(4) = t1
  else
    f(2) = 0.d0
    f(3) = 0.d0
    f(4) = 0.d0
  endif
  
  call mpmuld (s0, d10w, s1, mpnw1)
  call mpadd (s1, f, s0, mpnw1)
enddo

!  Correct exponent.

iexp = iexp - lnum2
f(2) = 1.d0
f(3) = 0.d0
f(4) = 10.d0
call mpnpwr (f, iexp, s1, mpnw1)
call mpmul (s0, s1, s2, mpnw1)
if (a(ksgn) == '-') s2(2) = -s2(2)

!   Restore original precision and exit.

call mproun (s2, mpnw)
call mpeq (s2, b, mpnw)

! write (6, *) 'mpctomp: output ='
! call mpout (6, 420, 400, b, mpnw)

return
end subroutine mpctomp

double precision function mpdigin (ca, n)

!   This converts the string CA of nonblank length N to double precision.
!   CA may only be modest length and may only contain digits.  Blanks are ignored.
!   This is intended for internal use only.  

  implicit none
  double precision d1
  character*(*), ca
  character*10 digits
  integer i, k, n
  parameter (digits = '0123456789')

! End of declaration

  d1 = 0.d0

  do i = 1, n
    if (ca(i:i) /= ' ') then
      k = index (digits, ca(i:i)) - 1
      if (k < 0) then
        write (mpldb, 1) ca(i:i)
1       format ('*** MPDIGIN: non-digit in character string = ',a)
        call mpabrt (86)
      elseif (k <= 9) then
        d1 = 10.d0 * d1 + k
      endif
    endif
  enddo

  mpdigin = d1
end function mpdigin

character*32 function mpdigout (a, n)

!   This converts the double precision input A to a character*32 string of 
!   nonblank length N.  A must be a whole number, and N must be sufficient
!   to hold it.  This is intended for internal use only.

  implicit none
  double precision a, d1, d2
  character*32 ca
  character*10 digits
  parameter (digits = '0123456789')
  integer i, k, n

! End of declaration

  ca = ' '
  d1 = abs (a)

  do i = n, 1, -1
    d2 = aint (d1 / 10.d0)
    k = 1.d0 + (d1 - 10.d0 * d2)
    d1 = d2
    ca(i:i) = digits(k:k)
  enddo

  mpdigout = ca
  return
end function mpdigout

subroutine mpeformat (a, nb, nd, b, mpnw)

!   Converts the MPR number A into character form in the CHARACTER*1 array B.
!   NB (input) is the length of the output string, and ND (input) is the
!   number of digits after the decimal point.  The format is analogous to 
!   Fortran E format.  The result is left-justified among the NB cells of B.
!   The condition NB >= ND + 10 must hold or an error message will result.
!   NB cells must be available in array B.

implicit none
integer i, ia, ix, ixp, i1, i2, j, k, mpnw, mpnw1, na, nb, nd, nexp, nl
character*1 b(nb), b2(nb+50)
character*10 digits
parameter (digits = '0123456789')
character*32 ca
double precision aa, an, t1, t2
double precision a(0:mpnw+5), f(0:8), s0(0:mpnw+6), s1(0:mpnw+6), d10w

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. nb < nd + 10) then
  write (mpldb, 1)
1 format ('*** MPEFORMAT: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

d10w = 10.d0**mpndpw
ia = sign (1.d0, a(2))
na = min (int (abs (a(2))), mpnw)
s0(0) = mpnw + 7
s1(0) = mpnw + 7
mpnw1 = mpnw + 1

!   Set f = 10.

f(0) = 9.d0
f(1) = mpnw1
f(2) = 1.d0
f(3) = 0.d0
f(4) = 10.d0
f(5) = 0.d0
f(6) = 0.d0

!   Determine power of ten for exponent, and scale input to within 1 and 10.

if (na > 0) then
  aa = a(4)
  if (na .ge. 2) aa = aa + mprdx * a(5)
  if (na .ge. 3) aa = aa + mprx2 * a(6)
  t1 = log10 (2.d0) * mpnbt * a(3) + log10 (aa)
  
  if (t1 .ge. 0.d0) then
    nexp = t1
  else
    nexp = t1 - 1.d0
  endif

  if (nexp == 0) then
    call mpeq (a, s1, mpnw1)
  elseif (nexp > 0) then
    call mpnpwr (f, nexp, s0, mpnw1)
    call mpdiv (a, s0, s1, mpnw1)
  elseif (nexp < 0) then
    call mpnpwr (f, - nexp, s0, mpnw1)
    call mpmul (a, s0, s1, mpnw1)
  endif

!   If we didn't quite get it exactly right, multiply or divide by 10 to fix.

100 continue
   
  if (s1(3) < 0.d0) then
    nexp = nexp - 1
    call mpmuld (s1, 10.d0, s0, mpnw1)
    call mpeq (s0, s1, mpnw1)
    goto 100
  elseif (s1(4) >= 10.d0) then
    nexp = nexp + 1
    call mpdivd (s1, 10.d0, s0, mpnw1)
    call mpeq (s0, s1, mpnw1)
    goto 100
  endif
  
  s1(2) = abs (s1(2))
else
  nexp = 0
  call mpeq (a, s1, mpnw1)
endif

!   Insert sign and first digit.

ix = 0
if (ia == -1) then
  ix = ix + 1
  b2(ix) = '-'
endif
if (na > 0) then
  an = s1(4)
else
  an = 0.d0
endif
ca = mpdigout (an, 1)
ix = ix + 1
b2(ix) = ca(1:1)
ix = ix + 1
b2(ix) = '.'
ixp = ix

!   Set f = an.

f(0) = 9.d0
f(1) = mpnw1
f(2) = 1.d0
f(3) = 0.d0
f(4) = an
f(5) = 0.d0
f(6) = 0.d0
call mpsub (s1, f, s0, mpnw1)
call mpmuld (s0, d10w, s1, mpnw1)

!   Calculate the number of remaining chunks.

nl = nd / mpndpw + 1

!   Insert the digits of the remaining words.

do j = 1, nl
  if (s1(2) /= 0.d0 .and. s1(3) == 0.d0) then
    an = s1(4)
    f(2) = 1.d0
    f(3) = 0.d0
    f(4) = an
  else
    f(2) = 0.d0
    f(3) = 0.d0
    f(4) = 0.d0
    an = 0.d0
  endif

  ca = mpdigout (an, mpndpw)
  
  do i = 1, mpndpw
    ix = ix + 1
    if (ix > nb + 50) then
      write (6, 2)
2     format ('MPEFORMAT: Insufficient space in B2 array.')
      call mpabrt (84)
    endif
    b2(ix) = ca(i:i)
  enddo

  call mpsub (s1, f, s0, mpnw1)
  call mpmuld (s0, d10w, s1, mpnw1)
enddo

!   Round the result.

if (ix >= nd + 1) then
  i1 = index (digits, b2(nd+1)) - 1
  if (i1 >= 5) then

!   Perform rounding, beginning at the last digit (position ND).  If the rounded
!   digit is 9, set to 0, then repeat at position one digit to left.  Continue
!   rounding if necessary until the decimal point is reached.

    do i = nd, ixp + 1, -1
      i2 = index (digits, b2(i)) - 1
      if (i2 <= 8) then
        b2(i) = digits(i2+2:i2+2)
        goto 180
      else
        b2(i) = '0'
      endif
    enddo

!   We have rounded up all digits to the right of the decimal point.  If the 
!   digit to the left of the decimal point is a 9, then set that digit to 1
!   and increase the exponent by one; otherwise increase that digit by one.

    if (b2(ixp-1) == '9') then
      b2(ixp-1) = '1'
      nexp = nexp + 1
    else
      i1 = index (digits, b2(ixp-1)) - 1
      b2(ixp-1) = digits(i1+2:i1+2)
    endif
  endif
endif

180 continue

!   Done with mantissa.  Insert exponent.

ix = nd + 1
b2(ix) = 'e'
if (nexp < 0) then
  ix = ix + 1
  b2(ix) = '-'
endif
ca = mpdigout (dble (abs (nexp)), 10)

do k = 1, 10
  if (ca(k:k) /= '0') goto 190
enddo

k = 10

190 continue

do i = k, 10
  ix = ix + 1
  b2(ix) = ca(i:i)
enddo

do i = ix + 1, nb
  b2(i) = ' '
enddo

!   Copy entire b2 array to B.

do i = 1, nb
  b(i) = b2(i)
enddo

return
end subroutine mpeformat

subroutine mpfformat (a, nb, nd, b, mpnw)

!   Converts the MPR number A into character form in the CHARACTER*1 array B.
!   NB (input) is the length of the output string, and ND (input) is the
!   number of digits after the decimal point.  The format is analogous to 
!   Fortran F format; the result is right-justified among the NB cells of B.
!   The condition NB >= ND + 10 must hold or an error message will result.
!   However, if it is found during execution that there is not sufficient space,
!   to hold all digits, the entire output field will be filled with asterisks.
!   NB cells of type CHARACTER*1 must be available in B.

implicit none
integer i, ia, ix, ixp, i1, i2, j, k, mpnw, mpnw1, na, nb, nb2, nd, &
  nexp, nl, n2
character*1 b(nb), b2(nb+20)
character*16 ca
double precision aa, an, an1, t1
double precision a(0:mpnw+5), f(0:8), s0(0:mpnw+6), s1(0:mpnw+6)
! character*16 mpdigout

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. nb < nd + 10) then
  write (mpldb, 1)
1 format ('*** MPFFORMAT: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

ia = sign (1.d0, a(2))
na = min (int (abs (a(2))), mpnw)
s0(0) = mpnw + 7
s1(0) = mpnw + 7
mpnw1 = mpnw + 1

!   Set f = 10.

f(0) = 9.d0
f(1) = mpnw1
f(2) = 1.d0
f(3) = 0.d0
f(4) = 10.d0
f(5) = 0.d0
f(6) = 0.d0

nb2 = nb + 20
call mpeformat (a, nb2, nb, b2, mpnw+1)

!   Trim off trailing blanks.

do i = nb2, 1, -1
  if (b2(i) /= ' ') goto 90
enddo

90 continue

nb2 = i

!   Look for the 'e' in B2.

do k = 1, nb2
  if (b2(k) == 'e') goto 100
enddo

write (6, 2)
2 format ('*** MPFFORMAT: Syntax error in output of mpeformat')
call mpabrt (84)

100 continue

!   Check the sign of the exponent.

k = k + 1
if (b2(k) == '-') then
  ixp = -1
  k = k + 1
else
  ixp = 1
endif
j = 0
ca = ' '

!   Copy the exponent into CA.

do i = k, nb2
  j = j + 1
  if (j <= 16) ca(j:j) = b2(i)
enddo

t1 = mpdigin (ca, j)

!   Check if there is enough space in the ouput array for all digits.

if (t1 + nd + 3 > nb) then
  do i = 1, nb
    b(i) = '*'
  enddo

  goto 200
endif
nexp = ixp * t1

!   Insert the sign of the number, if any.

i1 = 0
i2 = 0
if (b2(1) == '-') then
  i1 = i1 + 1
  b(i1) = '-'
  i2 = i2 + 1
endif

if (nexp == 0) then

!   Exponent is zero.  Copy first digit, period and ND more digits. 
  
  do i = 1, nd + 2
    i1 = i1 + 1
    i2 = i2 + 1
    b(i1) = b2(i2)
  enddo

  goto 200
elseif (nexp > 0) then

!   Exponent is positive.  Copy first digit, skip the period, then copy 
!   nexp digits.

  i1 = i1 + 1
  i2 = i2 + 1
  b(i1) = b2(i2)
  i2 = i2 + 1

  do i = 1, nexp
    i1 = i1 + 1
    i2 = i2 + 1
    b(i1) = b2(i2)
  enddo
  
!   Insert the period.

  i1 = i1 + 1
  b(i1) = '.'

!   Copy nd more digits.

  do i = 1, nd
    i1 = i1 + 1
    i2 = i2 + 1
    b(i1) = b2(i2)
  enddo
  
  goto 200
else

!   Exponent is negative.  Insert a zero, then a period, then nexp - 1
!   zeroes, then the first digit, then the remaining digits up to ND total 
!   fractional digits.

  i1 = i1 + 1
  b(i1) = '0'
  i1 = i1 + 1
  b(i1) = '.'

  do i = 1, nexp - 1
    i1 = i1 + 1
    b(i1) = '0'
  enddo

  i1 = i1 + 1
  i2 = i2 + 1
  b(i1) = b2(i2)
  i2 = i2 + 1

  do i = nexp, nd
    i1 = i1 + 1
    i2 = i2 + 1
    b(i1) = b2(i2)
  enddo
endif

200 continue

!   Right-justify in field.

k = nb - i1

do i = 1, i1
  b(nb-i+1) = b(nb-i-k+1)
enddo

do i = 1, k
  b(i) = ' '
enddo

return
end subroutine mpfformat

subroutine mpinp (iu, a, mpnw)

!   This routine reads the MPR number A from logical unit IU.  The digits of A 
!   may span more than one line, provided that a "\" appears at the end of 
!   a line to be continued (any characters after the "\" on the same line
!   are ignored).  Individual input lines may not exceed 2048 characters in
!   length (although this can be changed here and in the system parameters
!   (parameter mpnstr) at the beginning of this module.  Embedded blanks are
!   allowed anywhere. An exponent with "e" or "d" may optionally follow the
!   numeric value.

!   A scratch array below (CHR1) holds character data for input to mpctomp.
!   It is dimensioned MPNW * (MPNDPW + 1) + 1000 (see below).  If more nonblank
!   input characters than this are input, they are ignored.

implicit none
integer i, i1, iu, lnc1, lncx, ln1, mpnw
character*2048 line1
character*18 validc
parameter (validc = ' 0123456789+-.dDeE')
character*1 chr1(mpnw*(mpndpw+1)+1000)
double precision a(0:mpnw+5)

! End of declaration

if (mpnw < 4 .or. a(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPINP: uninitialized or inadequately sized arrays')
  call mpabrt (99)
endif

lnc1 = 0
lncx = mpnw * (mpndpw + 1) + 1000

100 continue

read (iu, '(a)', end = 200) line1

!   Find the last nonblank character.

do i = mpnstr, 1, -1
  if (line1(i:i) /= ' ') goto 110
enddo

!   Input line is blank -- ignore.

goto 100

110 continue

ln1 = i

!   Scan input line, looking for valid characters.

do i = 1, ln1
  if (line1(i:i) == '\') goto 100
  i1 = index (validc, line1(i:i))
  if (i1 == 0 .and. line1(i:i) /= ' ') then
      write (6, 2) line1(i:i)
2     format ('*** MPINP: Invalid input character = ',a)
      call mpabrt (87)
  elseif (line1(i:i) /= ' ') then
    if (lnc1 < lncx) then
      lnc1 = lnc1 + 1
      chr1(i) = line1(i:i)
    endif
  endif
enddo

call mpctomp (chr1, lnc1, a, mpnw)
goto 300

200  continue

write (mpldb, 4)
4 format ('*** MPINP: End-of-file encountered.')
call mpabrt (72)

300 return
end subroutine mpinp

subroutine mpout (iu, ln, nd, a, mpnw)

!   This routine writes MPR number A to logical unit IU in E(LN,ND) format.
!   This is output on MPOUTL characters per line.  The value of MPOUTL is set 
!   in the system parameters at the start of module MPFUNA.

implicit none
integer i, iu, ln, ln1, mpnw, nd
character*1 chr1(ln)
character*32 cform1, cform2
double precision a(0:mpnw+5)

! End of declaration

call mpeformat (a, ln, nd, chr1, mpnw)

write (cform1, 1) mpoutl
1 format ('(',i8,'a1)')
write (cform2, 2) mpoutl
2 format ('(',i8,'a1,"\")')

if (ln <= mpoutl) then
  write (iu, fmt = cform1) (chr1(i), i = 1, ln)
elseif (mod (ln, mpoutl) == 0) then
  ln1 = mpoutl * (ln / mpoutl) - mpoutl
  write (iu, fmt = cform2) (chr1(i), i = 1, ln1)
  write (iu, fmt = cform1) (chr1(i), i = ln1 + 1, ln)
else
  ln1 = mpoutl * (ln / mpoutl)
  write (iu, fmt = cform2) (chr1(i), i = 1, ln1)
  write (iu, fmt = cform1) (chr1(i), i = ln1 + 1, ln)
endif

return
end subroutine mpout

end module mpfunc

