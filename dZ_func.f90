!> Computes the first derivative of the plasma dispersion function according to Fried & Conte (1961), Z'(zeta)=-2 (1+zeta*Z(zeta))
!! \param zeta (complex) argument of the plasma dispersion function
function dZ_func(zeta)
  implicit none

  complex :: zeta, dZ_func, Z_func
  external Z_func

  dZ_func=-2*(1.0+zeta*Z_func(zeta))

end function dZ_func
