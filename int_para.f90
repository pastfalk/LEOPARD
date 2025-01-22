!> Computes the parallel velocity integrals over the interval [vara(ipara),vpara(ipara+1)] required for the dielectric tensor components
!! \param iarb index of species among the particle species with arbitrary velocity distribution
!! \param ipara index of parallel velocity gridpoint of the interval (lower limit) which is to be integrated over
!! \param k wavenumber
!! \param zeta (complex) parameter showing up in the integral
!! \param Kvpa array containing the values of the parallel velocity integrals for the given velocity interval
subroutine int_para(iarb,ipara,k,zeta,Kvpa)
  use param_mod
  implicit none
  real :: k
  integer :: iarb
  integer :: ipara
  complex(kind=16), dimension(6) :: Kvpa
  complex(kind=16) :: zeta
  integer :: sigma,l
  real(kind=16) :: vpa1, vpa2

  vpa2=vpara(ipara+1,iarb)
  vpa1=vpara(ipara,iarb)

  !compute the integrals
  Kvpa(1)=log((vpa2-zeta)/(vpa1-zeta))
  Kvpa(2)=(vpa2-vpa1)+zeta*Kvpa(1)
  Kvpa(3)=(vpa2**2 -vpa1**2)/2.0_16 + zeta*Kvpa(2)
  Kvpa(4)=(vpa2**3 -vpa1**3)/3.0_16 + zeta*Kvpa(3)
  Kvpa(5)=(vpa2**4 -vpa1**4)/4.0_16 + zeta*Kvpa(4)
  Kvpa(6)=(vpa2**5 -vpa1**5)/5.0_16 + zeta*Kvpa(5)

  !account for poles that lie within the integration interval
  !according to Landau's prescription
  if((vpa2.ge.real(zeta)).and.(vpa1<real(zeta))) then

     if(aimag(zeta)>0.0) then
        sigma=0
     else if(aimag(zeta)==0.0) then
        sigma=1
     else
        sigma=2
     endif

     do l=1,6 
        Kvpa(l)=Kvpa(l)+i*pi* k/abs(k) * sigma *zeta**(l-1)
     enddo

  endif


end subroutine int_para
