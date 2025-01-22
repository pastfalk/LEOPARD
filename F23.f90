!> Computes the generalized hypergeometric function 2F3 (a1,a2;b1,b2,b3;z) from a continued fraction formula
!! \param a array containing the arguments a1 and a2 of 2F3
!! \param b array containing the arguments b1, b2, and b3 of 2F3
!! \param z argument of 2F3
!! \param nfrac required number of terms to be included in continued fraction formula for an accurate evaluation of 2F3
!! \param F23_val estimated value of 2F3
subroutine F23(a,b,z,nfrac,F23_val)
  implicit none
  real(kind=16), dimension(2) :: a
  real(kind=16), dimension(3) :: b
  real(kind=16) :: z, sol
  integer :: nfrac,j
  real(kind=16) :: F23_val
  real(kind=16) :: fac1, fac2
  
  sol=-z/(1.0_16*(nfrac+2)) * (1.0_16*(nfrac+1)+a(1))*(1.0_16*(nfrac+1)+a(2))/&
       & ((1.0_16*(nfrac+1)+b(1))*(1.0_16*(nfrac+1)+b(2))*(1.0_16*(nfrac+1)+b(3)))

  do j=nfrac,0,-1

     if(j==0) then
        sol=z *a(1)*a(2)/(b(1)*b(2)*b(3))/&
             & (1.0_16+sol)
     else

        fac1=(1.0_16*j+a(1))*(1.0_16*j+a(2))
        fac2=(1.0_16*j+b(1))*(1.0_16*j+b(2))*(1.0_16*j+b(3))
        sol=   -z/(1.0_16*(j+1)) * fac1/fac2/&  
             & (1.0_16+z/(1.0_16*(j+1)) * fac1/fac2+sol)

     endif

  enddo
  
  sol=sol+1.0_16

  F23_val=sol
  

end subroutine F23
