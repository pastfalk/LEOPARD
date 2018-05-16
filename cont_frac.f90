!> Computes the plasma dispersion function from a continued fraction formula including the first three terms only
!! \param zeta (complex) argument of the plasma dispersion function
!! \param Z_func value of the plasma dispersion function
subroutine cont_frac(zeta,Z_func)
  implicit none
  complex :: zeta
  complex :: Z_func
  real :: a

  a=0.5

  Z_func=zeta/(-zeta**2 + a/(1.0+(a/2.0)/(-zeta**2 +(a/4.0))))    


end subroutine cont_frac
