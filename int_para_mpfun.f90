!> Computes all parallel velocity integrals over the interval [vara(ipara),vpara(ipara+1)] required for the epsilon tensor components using arbitrary precision arithmetic for better accuracy
!! \param iarb index of species among the particle species with arbitrary velocity distribution
!! \param ipara index of parallel velocity gridpoint of the considered interval (lower limit) which is to be integrated over
!! \param ndp measure for the number of significant digits to be included in the arbitrary precision arithmetic
!! \param k wavenumber
!! \param zeta complex parameter showing up in the integral
!! \param Kvpa array containing the values of the parallel velocity integrals for the given velocity interval
subroutine int_para_mpfun(iarb,ipara,ndp,k,zeta,Kvpa)
  use param_mod
  use mpmodule
  implicit none
  integer :: iarb
  integer :: ipara
  real :: k
  complex(kind=16), dimension(6) :: Kvpa
  complex(kind=16) :: zeta
  integer :: l
  integer :: ndp,ndws
  type(mp_complex), dimension(6) :: sol
  type(mp_complex) :: zeta_mp
  type(mp_real) :: vpa1,vpa2

  type(mp_real) :: zeta_r_mp
  type(mp_real) :: zeta_i_mp

  ndws = int (ndp / mpdpw + 2)

  vpa1=mpreald(vpara(ipara,iarb),ndws)
  vpa2=mpreald(vpara(ipara+1,iarb),ndws)

  zeta_mp=mpcmplxdc(cmplx(real(zeta),aimag(zeta)),ndws)

  zeta_r_mp=mpreal(real(zeta),ndws)
  zeta_i_mp=mpreal(aimag(zeta),ndws)


  !compute the integrals

!  sol(1)=log((vpa2-zeta_mp)/(vpa1-zeta_mp))
  sol(1)=log(sqrt( ( ((vpa2-zeta_r_mp)*(vpa1-zeta_r_mp)+zeta_i_mp**2)/&
       &             ( (vpa1-zeta_r_mp)**2 +zeta_i_mp**2)     )**2 +&
       &           (zeta_i_mp*(vpa2-vpa1)/((vpa1-zeta_r_mp)**2+zeta_i_mp**2)   )**2 )  )+&
       & atan2( zeta_i_mp*(vpa2-vpa1)/((vpa1-zeta_r_mp)**2+zeta_i_mp**2),&
       &        ((vpa2-zeta_r_mp)*(vpa1-zeta_r_mp)+zeta_i_mp**2)/&
       &             ( (vpa1-zeta_r_mp)**2 +zeta_i_mp**2) )*mpcmplx(cmplx(0.0_16,aimag(i)),ndws)

  sol(2)=(vpa2-vpa1)+zeta_mp*sol(1)
  sol(3)=(vpa2**2-vpa1**2)/mpreal(2.0_16,ndws) + zeta_mp*sol(2)
  sol(4)=(vpa2**3 -vpa1**3)/mpreal(3.0_16,ndws) + zeta_mp*sol(3)
  sol(5)=(vpa2**4 -vpa1**4)/mpreal(4.0_16,ndws) + zeta_mp*sol(4)
  sol(6)=(vpa2**5 -vpa1**5)/mpreal(5.0_16,ndws) + zeta_mp*sol(5)


  !account for poles that lie within the integration interval
  !according to Landau's prescription
  if((vpara(ipara+1,iarb).ge.real(zeta)).and.(vpara(ipara,iarb)<real(zeta))) then

     if(aimag(zeta)==0.0) then

        if(k>0.0) then

           sol(1)=sol(1)+mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(pi,ndws)
           sol(2)=sol(2)+mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(pi,ndws)*zeta_mp
           sol(3)=sol(3)+mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(pi,ndws)*zeta_mp*zeta_mp
           sol(4)=sol(4)+mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(pi,ndws)*&
                & zeta_mp*zeta_mp*zeta_mp
           sol(5)=sol(5)+mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(pi,ndws)*&
                & zeta_mp*zeta_mp*zeta_mp*zeta_mp
           sol(6)=sol(6)+mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(pi,ndws)*&
                & zeta_mp*zeta_mp*zeta_mp*zeta_mp*zeta_mp

        else

           sol(1)=sol(1)-mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(pi,ndws)
           sol(2)=sol(2)-mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(pi,ndws)*zeta_mp
           sol(3)=sol(3)-mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(pi,ndws)*zeta_mp*zeta_mp
           sol(4)=sol(4)-mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(pi,ndws)*&
                & zeta_mp*zeta_mp*zeta_mp
           sol(5)=sol(5)-mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(pi,ndws)*&
                & zeta_mp*zeta_mp*zeta_mp*zeta_mp
           sol(6)=sol(6)-mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(pi,ndws)*&
                & zeta_mp*zeta_mp*zeta_mp*zeta_mp*zeta_mp

        endif
   
     else if(aimag(zeta)<0.0) then

        if(k>0.0) then

           sol(1)=sol(1)+mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(2*pi,ndws)
           sol(2)=sol(2)+mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(2*pi,ndws)*zeta_mp
           sol(3)=sol(3)+mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(2*pi,ndws)*zeta_mp*zeta_mp
           sol(4)=sol(4)+mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(2*pi,ndws)*&
                & zeta_mp*zeta_mp*zeta_mp
           sol(5)=sol(5)+mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(2*pi,ndws)*&
                & zeta_mp*zeta_mp*zeta_mp*zeta_mp
           sol(6)=sol(6)+mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(2*pi,ndws)*&
                & zeta_mp*zeta_mp*zeta_mp*zeta_mp*zeta_mp

        else

           sol(1)=sol(1)-mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(2*pi,ndws)
           sol(2)=sol(2)-mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(2*pi,ndws)*zeta_mp
           sol(3)=sol(3)-mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(2*pi,ndws)*zeta_mp*zeta_mp
           sol(4)=sol(4)-mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(2*pi,ndws)*&
                & zeta_mp*zeta_mp*zeta_mp
           sol(5)=sol(5)-mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(2*pi,ndws)*&
                & zeta_mp*zeta_mp*zeta_mp*zeta_mp
           sol(6)=sol(6)-mpcmplx(cmplx(0.0_16,aimag(i)),ndws)*mpreald(2*pi,ndws)*&
                & zeta_mp*zeta_mp*zeta_mp*zeta_mp*zeta_mp

        endif

     endif

  endif

  do l=1,6
     Kvpa(l)=qreal(mpreal(sol(l)))+i*qreal(aimag(sol(l)))
  enddo


end subroutine int_para_mpfun
