!> Computes the (exponentially-scaled) modified Bessel function of the first kind using its series expansion, allowing for big arguments 
!! \param n index of the modified Bessel function
!! \param z argument of the modified Bessel function
function exp_Bessel_In_mpfun(n,z_in,ndp)
  use mpmodule
  use param_mod
  implicit none
  real :: exp_Bessel_In_mpfun, error
  real(kind=16) :: z
  real :: z_in
  integer :: n, j
  integer :: m
  real(kind=16) :: gamma_func
  type(mp_real) :: expBes, old
  type(mp_real) :: fac1, fac2
  integer :: ndp, ndws

  external gamma_func

  z=z_in

  ndws = int (ndp / mpdpw + 2)

  if(n<170) then
     expBes=mpreal(z/2.0_16,ndws)**abs(n) /mpreal(gamma_func(abs(n)+1),ndws)*exp(-mpreal(z,ndws))
  else
     write(*,*) 'n too big in gamma_func'
     stop
  endif
  
  m=1

  do while(.true.)

     fac1=mpreal(1.0_16,ndws)

     do j=1,m

        fac1=fac1*mpreal(z/(2*j),ndws)

     enddo

     fac2=mpreal(1.0_16,ndws)

     do j=1,abs(n)

        fac2=fac2*mpreal(z/(2*(m+j)),ndws)

     enddo

     old=expBes

     expBes=expBes+fac1*fac1*fac2


     if(expBes==mpreal(0.0_16,ndws)) then
        exit
     endif

     error=abs(mpreal(1.0,ndws)-old/expBes)

     if (error<10.0_16**(-15)) then
        expBes=expBes*exp(-mpreal(z,ndws))
        exit
     endif

     m=m+1

  enddo

  exp_Bessel_In_mpfun=expBes


end function exp_Bessel_In_mpfun
