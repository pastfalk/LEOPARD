!> Computes the plasma dispersion function using its relation to the complex error function or - in case of high zeta - a continued fraction formula
!! \param zeta (complex) argument of the plasma dispersion function
function Z_func(zeta)
  use param_mod
  implicit none
  complex :: Z_func, zeta
  complex :: error_func

  if(abs(real(zeta))<26.0) then

     call cerror(i*zeta,error_func)
     Z_func=i*sqrt(pi)*exp(-zeta*zeta)*(1.0+error_func)

  else

     call cont_frac(zeta,Z_func)

  endif


end function Z_func
