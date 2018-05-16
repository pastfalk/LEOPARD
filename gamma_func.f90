!> Computes the gamma function value for a given integer argument x
!! \param x (integer) argument of the gamma function
function gamma_func(x)
  implicit none
  real(kind=16) :: gamma_func
  integer :: n,x

  gamma_func=1.0_16

  do n=1,x-1
     gamma_func=gamma_func*n
  enddo

end function gamma_func
