!> Computes the generalized hypergeometric function 2F3 (a1,a2;b1,b2,b3;z) from a continued fraction formula using arbitrary precision arithmetic for better accuracy
!! \param a array containing the arguments a1 and a2 of 2F3
!! \param b array containing the arguments b1, b2, and b3 of 2F3
!! \param z argument of 2F3
!! \param ndp measure for the number of significant digits to be included in the arbitrary precision arithmetic
!! \param nfrac required number of terms to be included in continued fraction formula for an accurate evaluation of 2F3
!! \param F23_val estimated value of 2F3
subroutine F23_mpfun(a,b,z,ndp,nfrac,F23_val)
  use mpmodule
  implicit none
  real(kind=16), dimension(2) :: a
  real(kind=16), dimension(3) :: b
  real(kind=16) :: z
  real(kind=16) :: F23_val
  real :: F23_dummy
  integer :: j
  integer :: ndp, ndws,nfrac
  type(mp_real) :: fac1, fac2
  type(mp_real) :: sol
  
  ndws = int (ndp / mpdpw + 2)

  sol=-mpreal(z,ndws)/( mpreal(1.0_16,ndws)*(nfrac+2)) *&
       & (mpreal(1.0_16,ndws)*(nfrac+1)+mpreal(a(1),ndws))*(mpreal(1.0_16,ndws)*(nfrac+1)+mpreal(a(2),ndws))/&
       &((mpreal(1.0_16,ndws)*(nfrac+1)+mpreal(b(1),ndws))*(mpreal(1.0_16,ndws)*(nfrac+1)+mpreal(b(2),ndws))*&
       & (mpreal(1.0_16,ndws)*(nfrac+1)+mpreal(b(3),ndws)))

  do j=nfrac,0,-1
     if(j==0) then
        sol=     mpreal(z,ndws)   *mpreal(a(1),ndws)*mpreal(a(2),ndws)&
             & /(mpreal(b(1),ndws)*mpreal(b(2),ndws)*mpreal(b(3),ndws))/&
             &  (mpreal(1.0_16,ndws)+sol)
     else
        fac1=  (mpreal(1.0_16*j,ndws)+mpreal(a(1),ndws))*(mpreal(1.0_16*j,ndws)+mpreal(a(2),ndws))
        fac2=  (mpreal(1.0_16*j,ndws)+mpreal(b(1),ndws))*(mpreal(1.0_16*j,ndws)+mpreal(b(2),ndws))*&
             & (mpreal(1.0_16*j,ndws)+mpreal(b(3),ndws))
        
        sol=   -mpreal(z,ndws)/(mpreal(1.0_16*(j+1),ndws)) * fac1/fac2/&  
             & (mpreal(1.0_16,ndws)+mpreal(z,ndws)/(mpreal(1.0_16*(j+1),ndws)) * fac1/fac2+sol)
     endif
  enddo
  
  sol=sol+mpreal(1.0_16,ndws)

  F23_dummy=sol
  F23_val=F23_dummy


end subroutine F23_mpfun
