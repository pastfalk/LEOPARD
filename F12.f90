!> Computes the generalized hypergeometric function 1F2 (a1;b1,b2;z) from a continued fraction formula
!! \param a argument a1 of 1F2
!! \param b array containing the arguments b1 and b2 of 1F2
!! \param z argument of 1F2
!! \param nfrac required number of terms to be included in continued fraction formula for an accurate evaluation of 1F2
!! \param F12_val estimated value of 1F2
subroutine F12(a,b,z,nfrac,F12_val)
  implicit none
  real(kind=16), dimension(1) :: a
  real(kind=16), dimension(2) :: b
  real(kind=16) :: z, sol 
  integer :: nfrac,j
  real(kind=16) :: F12_val
  real(kind=16) :: fac1,fac2

  sol=-z/(1.0_16*(nfrac+2)) * (1.0_16*(nfrac+1)+a(1))/&
       & ((1.0_16*(nfrac+1)+b(1))*(1.0_16*(nfrac+1)+b(2)))

  do j=nfrac,0,-1

     if(j==0) then
        sol=z *a(1)/(b(1)*b(2))/&
             & (1.0_16+sol)
     else

        fac1=1.0_16*j+a(1)
        fac2=(1.0_16*j+b(1))*(1.0_16*j+b(2))

        sol= -z/(1.0_16*(j+1)) * fac1/fac2/&  
             & (1.0_16  +z/(1.0_16*(j+1)) * fac1/fac2 +sol)

     endif

  enddo

  sol=sol+1.0_16
  F12_val=sol

  
end subroutine F12
