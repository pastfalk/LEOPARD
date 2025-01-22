!> Computes the (exponentially-scaled) first derivative of the modified Bessel function of the first kind using the relation dI_n/dz = 1/2*(I_(n-1) + I_(n+1))
!! \param n index of the modified Bessel function
!! \param z argument of the modified Bessel function
function exp_dBessel_In(n,z)
  implicit none
  integer :: n
  real :: z
  real :: exp_dBessel_In, exp_Bessel_In, exp_Bessel_In_mpfun

  external exp_Bessel_In
  external exp_Bessel_In_mpfun

  if(z<3400.0) then
     exp_dBessel_In=0.5*(exp_Bessel_In(n-1,z) + exp_Bessel_In(n+1,z))
  else
     exp_dBessel_In=0.5*(exp_Bessel_In_mpfun(n-1,z,40) + exp_Bessel_In_mpfun(n+1,z,40))
  endif


end function exp_dBessel_In
