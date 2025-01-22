!> Computes the generalized hypergeometric function 1F2 (a1;b1,b2;z) from a continued fraction formula using arbitrary precision arithmetic for better accuracy
!! \param a argument a1 of 1F2
!! \param b array containing the arguments b1 and b2 of 1F2
!! \param z argument of 1F2
!! \param ndp measure for the number of significant digits to be included in the arbitrary precision arithmetic
!! \param nfrac required number of terms to be included in continued fraction formula for an accurate evaluation of 1F2
!! \param F12_val estimated value of 1F2
subroutine F12_mpfun(a,b,z,ndp,nfrac,F12_val)                       
  use mpmodule
  implicit none
  real(kind=16), dimension(1) :: a
  real(kind=16), dimension(2) :: b
  real(kind=16) :: z
  real(kind=16) :: F12_val
  real :: F12_dummy
  integer :: j
  integer :: ndp,ndws,nfrac
  type(mp_real) :: fac1, fac2
  type(mp_real) :: sol

  ndws = int (ndp / mpdpw + 2)

  sol=-mpreal(z,ndws)/( mpreal(1.0_16,ndws)*(nfrac+2)) *&
       & (mpreal(1.0_16,ndws)*(nfrac+1)+mpreal(a(1),ndws))/&
       &((mpreal(1.0_16,ndws)*(nfrac+1)+mpreal(b(1),ndws))*(mpreal(1.0_16,ndws)*(nfrac+1)+mpreal(b(2),ndws)))

  do j=nfrac,0,-1

     if(j==0) then
        sol=     mpreal(z,ndws)   *mpreal(a(1),ndws)/&
             &  (mpreal(b(1),ndws)*mpreal(b(2),ndws))/&
             &  (mpreal(1.0_16,ndws)+sol)
     else
        fac1=   mpreal(1.0_16*j,ndws)+mpreal(a(1),ndws)
        fac2=  (mpreal(1.0_16*j,ndws)+mpreal(b(1),ndws))*(mpreal(1.0_16*j,ndws)+mpreal(b(2),ndws))
        
        sol=   -mpreal(z,ndws)/(mpreal(1.0_16*(j+1),ndws)) * fac1/fac2/&  
             & (mpreal(1.0_16,ndws)+mpreal(z,ndws)/(mpreal(1.0_16*(j+1),ndws)) * fac1/fac2+sol)
     endif

  enddo

  sol=sol+mpreal(1.0_16,ndws)

  F12_dummy=sol
  F12_val=F12_dummy


end subroutine F12_mpfun
