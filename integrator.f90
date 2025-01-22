!> Computes the total values of the 2d velocity integrals required for the dielectric tensor components by performing the parallel velocity integration and summing up the contributions of both the parallel and the perpendicular integrals
!! \param omega complex wave frequency
!! \param k wavenumber
!! \param m index of particle species
!! \param iarb index of species among the particle species with arbitrary velocity distribution
!! \param n Bessel index
!! \param Ivpe array containing the contributions of the velocity integrals from each interval in the perpendicular velocity grid
!! \param intgrl array containing the total values of the 2d velocity integrals
!! \param splcoeff1 array of coefficients obtained from cubic spline interpolation of the provided velocity distribution data
!! \param splcoeff2 array of coefficients obtained from cubic spline interpolation of the provided velocity distribution data
subroutine integrator(omega,k,m,iarb,n,Ivpe,intgrl,splcoeff1,splcoeff2)
  use param_mod
  implicit none
  complex :: omega
  real :: k
  integer :: m,n,iarb
  complex, dimension(3,3,2,2) :: intgrl
  integer :: pe, pa
  integer :: ipara,iperp
  complex(kind=16) :: zeta
  real :: zeta_r

  real, dimension (npara(iarb)-1,nperp(iarb)-1,4,3) :: splcoeff1
  real, dimension (npara(iarb)-1,nperp(iarb)-1,4,3) :: splcoeff2

  complex, dimension(3) :: sum0
  real, dimension(3) :: sum1, sum2, sum3

  real, dimension (3,5,nperp(iarb)) :: Ivpe

  complex, dimension (6,npara(iarb)-1) :: Kvpa
  complex(kind=16), dimension (6) :: Kvpa_dummy

  real :: dvpa

  integer :: ndp


  !perform velocity integration in parallel direction

  zeta=(omega-n*mu(m)*q(m))/(k*cos(theta))
  zeta_r=real(zeta)

  dvpa=abs(vpara(2,iarb)-vpara(1,iarb))

  do ipara=1,npara(iarb)-1


     !determine whether or not the parallel velocity integration requires arbitrary precision arithmetic
     call acc_Kvpa(dvpa,zeta_r,ndp)

     if(ndp==0) then
        call int_para(iarb,ipara,k,zeta,Kvpa_dummy)
     else
        call int_para_mpfun(iarb,ipara,ndp,k,zeta,Kvpa_dummy)
     endif


     Kvpa(1,ipara)=Kvpa_dummy(1)
     Kvpa(2,ipara)=Kvpa_dummy(2)
     Kvpa(3,ipara)=Kvpa_dummy(3)
     Kvpa(4,ipara)=Kvpa_dummy(4)
     Kvpa(5,ipara)=Kvpa_dummy(5)
     Kvpa(6,ipara)=Kvpa_dummy(6)

  enddo


  Kvpa=-Kvpa/(k*cos(theta))


  !sum up the contributions of each grid interval in both parallel and perpendicular direction
  !to compute the total values of all 2d velocity integrals required for the dispersion tensor elements

  intgrl=(0.0,0.0)

  do iperp=1,nperp(iarb)-1

     do pa=1,3
        do pe=1,2

           sum0=(0.0,0.0)

           do ipara=1,npara(iarb)-1

              sum0=sum0+splcoeff1(ipara,iperp,1,:)*Kvpa(3+pa,ipara)+splcoeff1(ipara,iperp,2,:)*Kvpa(2+pa,ipara)+&
                   &    splcoeff1(ipara,iperp,3,:)*Kvpa(1+pa,ipara)+splcoeff1(ipara,iperp,4,:)*Kvpa(  pa,ipara)

           enddo

           intgrl(1,pa,pe,1)=intgrl(1,pa,pe,1)+&
                & (Ivpe(1,2+pe,iperp+1)-Ivpe(1,2+pe,iperp))*sum0(1)+&
                & (Ivpe(1,1+pe,iperp+1)-Ivpe(1,1+pe,iperp))*sum0(2)+&
                & (Ivpe(1,  pe,iperp+1)-Ivpe(1,  pe,iperp))*sum0(3)

           intgrl(2,pa,pe,1)=intgrl(2,pa,pe,1)+&
                & (Ivpe(2,2+pe,iperp+1)-Ivpe(2,2+pe,iperp))*sum0(1)+&
                & (Ivpe(2,1+pe,iperp+1)-Ivpe(2,1+pe,iperp))*sum0(2)+&
                & (Ivpe(2,  pe,iperp+1)-Ivpe(2,  pe,iperp))*sum0(3)

           intgrl(3,pa,pe,1)=intgrl(3,pa,pe,1)+&
                & (Ivpe(3,2+pe,iperp+1)-Ivpe(3,2+pe,iperp))*sum0(1)+&
                & (Ivpe(3,1+pe,iperp+1)-Ivpe(3,1+pe,iperp))*sum0(2)+&
                & (Ivpe(3,  pe,iperp+1)-Ivpe(3,  pe,iperp))*sum0(3)


        enddo

     enddo

  enddo



  do ipara=1,npara(iarb)-1

     do pa=1,3
        do pe=1,2


           sum1=0.0
           sum2=0.0
           sum3=0.0


           do iperp=1,nperp(iarb)-1

              sum1=sum1+splcoeff2(ipara,iperp,1,:)*(Ivpe(1,3+pe,iperp+1)-Ivpe(1,3+pe,iperp))+&
                   &    splcoeff2(ipara,iperp,2,:)*(Ivpe(1,2+pe,iperp+1)-Ivpe(1,2+pe,iperp))+&
                   &    splcoeff2(ipara,iperp,3,:)*(Ivpe(1,1+pe,iperp+1)-Ivpe(1,1+pe,iperp))+&
                   &    splcoeff2(ipara,iperp,4,:)*(Ivpe(1,  pe,iperp+1)-Ivpe(1,  pe,iperp))

              sum2=sum2+splcoeff2(ipara,iperp,1,:)*(Ivpe(2,3+pe,iperp+1)-Ivpe(2,3+pe,iperp))+&
                   &    splcoeff2(ipara,iperp,2,:)*(Ivpe(2,2+pe,iperp+1)-Ivpe(2,2+pe,iperp))+&
                   &    splcoeff2(ipara,iperp,3,:)*(Ivpe(2,1+pe,iperp+1)-Ivpe(2,1+pe,iperp))+&
                   &    splcoeff2(ipara,iperp,4,:)*(Ivpe(2,  pe,iperp+1)-Ivpe(2,  pe,iperp))

              sum3=sum3+splcoeff2(ipara,iperp,1,:)*(Ivpe(3,3+pe,iperp+1)-Ivpe(3,3+pe,iperp))+&
                   &    splcoeff2(ipara,iperp,2,:)*(Ivpe(3,2+pe,iperp+1)-Ivpe(3,2+pe,iperp))+&
                   &    splcoeff2(ipara,iperp,3,:)*(Ivpe(3,1+pe,iperp+1)-Ivpe(3,1+pe,iperp))+&
                   &    splcoeff2(ipara,iperp,4,:)*(Ivpe(3,  pe,iperp+1)-Ivpe(3,  pe,iperp))


           enddo


           intgrl(1,pa,pe,2)=intgrl(1,pa,pe,2)+&
                & Kvpa(2+pa,ipara)*sum1(1)+&
                & Kvpa(1+pa,ipara)*sum1(2)+&
                & Kvpa(  pa,ipara)*sum1(3)

           intgrl(2,pa,pe,2)=intgrl(2,pa,pe,2)+&
                & Kvpa(2+pa,ipara)*sum2(1)+&
                & Kvpa(1+pa,ipara)*sum2(2)+&
                & Kvpa(  pa,ipara)*sum2(3)

           intgrl(3,pa,pe,2)=intgrl(3,pa,pe,2)+&
                & Kvpa(2+pa,ipara)*sum3(1)+&
                & Kvpa(1+pa,ipara)*sum3(2)+&
                & Kvpa(  pa,ipara)*sum3(3)

        enddo

     enddo

  enddo


end subroutine integrator
