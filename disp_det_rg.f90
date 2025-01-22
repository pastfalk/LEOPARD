!> Computes the determinant of the dispersion tensor D for a plasma with bi-Maxwellian particle species and/or particles with arbitrary gyrotropic velocity distribution
!! \param omega complex wave frequency
!! \param k wavenumber
!! \param splcoeff1 array of coefficients obtained from cubic spline interpolation of the provided velocity distribution data 
!! \param splcoeff2 array of coefficients obtained from cubic spline interpolation of the provided velocity distribution data 
function disp_det_rg(omega,k,splcoeff1,splcoeff2) result (disp_det)
  use param_mod
  implicit none
  complex :: epsilon_xx, epsilon_yy, epsilon_zz
  complex :: epsilon_xy, epsilon_xz, epsilon_yz
  complex :: del_xx, del_yy, del_zz
  complex :: del_xy, del_xz, del_yz

  logical, dimension(6) :: esc

  real, dimension(npara_max-1,nperp_max-1,4,3,narb) :: splcoeff1, splcoeff2

  complex :: omega
  real :: k
  complex :: disp_det
  integer :: m, n

  complex, dimension(3,3,2,2) :: intgrl
  complex, dimension(3,3,2,2) :: intgrl2

  real, allocatable, dimension(:,:,:) :: Ivpe
  real, dimension (3,5) :: Ivpe_dummy
  integer :: iperp
  integer :: iarb

  real :: lambda
  complex :: zeta1, zeta2

  complex :: Z_func, dZ_func
  real :: exp_Bessel_In, exp_dBessel_In, exp_Bessel_In_mpfun
  real :: expBes, expdBes
  
  real :: tdelr, tdeli

  external exp_Bessel_In
  external exp_Bessel_In_mpfun
  external exp_dBessel_In
  external Z_func
  external dZ_func


  epsilon_xx=0.0
  epsilon_yy=0.0
  epsilon_zz=0.0
  epsilon_xy=0.0
  epsilon_xz=0.0
  epsilon_yz=0.0

  !note that the dielectric tensor is renormalized according to epsilon_ij -> epsilon_ij * (v_A ^2/c^2 * omega^2)  

  epsilon_xx=epsilon_xx+(delta*omega)**2 
  epsilon_yy=epsilon_yy+(delta*omega)**2 
  epsilon_zz=epsilon_zz+(delta*omega)**2 


  iarb=0

  do m=1,Nspecies


     if(mode(m)==0) then

        !for bi-Maxwellian scenario

        esc=.true.

        epsilon_xx=epsilon_xx+ (beta_ratio(m)-1.0) *mu(m)*dens(m)*q(m)**2
        epsilon_yy=epsilon_yy+ (beta_ratio(m)-1.0) *mu(m)*dens(m)*q(m)**2

        lambda=(k*sin(theta))**2 * beta_perp(m) / (2.0 * q(m)**2 * mu(m) *dens(m))

        n=0

        do while(esc(1).or.esc(2).or.esc(3).or.esc(4).or.esc(5).or.esc(6))


           !compute exponentially-scaled modified Bessel functions required for the evaluation of the dielectric tensor

           if(lambda<3400.0) then
              expBes=exp_Bessel_In(n,lambda)
           else 
              expBes=exp_Bessel_In_mpfun(n,lambda,40)
           endif

           expdBes=exp_dBessel_In(n,lambda)

           
           !compute the dielectric tensor components

           if(n==0) then
              zeta1=(omega-k*cos(theta)*drift(m))/sqrt(beta_para(m))/(k*cos(theta))/sqrt(mu(m)) *sqrt(dens(m))

              epsilon_yy = epsilon_yy - sqrt(mu(m)) *dens(m)**1.5 * q(m)**2 / sqrt(beta_para(m)) / (k*cos(theta)) *&
                   & (2*lambda*(expdBes-expBes)) *(beta_ratio(m)*(omega-k*cos(theta)*drift(m)))*Z_func(zeta1)

              epsilon_yz = epsilon_yz + i/2 *dens(m)* q(m) * ( k*sin(theta))/(k*cos(theta)) *&
                   &  (beta_ratio(m)*omega)* (expdBes-expBes) *dZ_func(zeta1)


              epsilon_zz = epsilon_zz - dens(m)**2 * q(m)**2  / beta_perp(m) /(k*cos(theta))**2 * expBes *&
                   & omega*omega*beta_ratio(m)* dZ_func(zeta1)    

           else

              zeta1=(omega-k*cos(theta)*drift(m)-n*mu(m)*q(m))/sqrt(beta_para(m))/(k*cos(theta))/sqrt(mu(m)) *sqrt(dens(m))
              zeta2=(omega-k*cos(theta)*drift(m)+n*mu(m)*q(m))/sqrt(beta_para(m))/(k*cos(theta))/sqrt(mu(m)) *sqrt(dens(m))


              if(esc(1)) then

                 del_xx=sqrt(mu(m))*dens(m)**1.5 *q(m)**2 /sqrt(beta_para(m))/(k*cos(theta))*&
                      & n**2 *expBes / lambda *(beta_ratio(m) * (omega-k*cos(theta)*drift(m)) - &
                      & (beta_ratio(m)-1.0)*n*mu(m)*q(m))*Z_func(zeta1)+&
                      
                      & sqrt(mu(m))*dens(m)**1.5 *q(m)**2 /sqrt(beta_para(m))/(k*cos(theta))*&
                      & n**2 *expBes / lambda *(beta_ratio(m) * (omega-k*cos(theta)*drift(m)) + &
                      & (beta_ratio(m)-1.0)*n*mu(m)*q(m))*Z_func(zeta2)

                 if((n<=4).or.((n>4).and.((abs(real(del_xx)/ real(epsilon_xx))>eps_error).or. &
                      & (abs(aimag(del_xx)/aimag(epsilon_xx))>eps_error)))) then

                    epsilon_xx=epsilon_xx+ del_xx

                 else
                    esc(1)=.false.
                 endif
              endif

              if(esc(2)) then


                 del_yy = sqrt(mu(m)) *dens(m)**1.5 * q(m)**2 / sqrt(beta_para(m)) / (k*cos(theta)) *&
                      & (n**2 * expBes / lambda - 2*lambda*(expdBes-expBes)) *&
                      & (beta_ratio(m)*(omega-k*cos(theta)*drift(m)) - (beta_ratio(m)-1.0)*n*mu(m)*q(m))*Z_func(zeta1)+&
                      
                      & sqrt(mu(m)) *dens(m)**1.5 * q(m)**2 / sqrt(beta_para(m)) / (k*cos(theta)) *&
                      & (n**2 * expBes / lambda - 2*lambda*(expdBes-expBes)) *&
                      & (beta_ratio(m)*(omega-k*cos(theta)*drift(m)) + (beta_ratio(m)-1.0)*n*mu(m)*q(m))*Z_func(zeta2)

                 if((n<=4).or.((n>4).and.((abs(real(del_yy)/ real(epsilon_yy))>eps_error).or. &
                      & (abs(aimag(del_yy)/aimag(epsilon_yy))>eps_error)))) then

                    epsilon_yy=epsilon_yy+ del_yy

                 else
                    esc(2)=.false.
                 endif
              endif


              if(esc(3)) then

                 del_zz = - dens(m)**2 * q(m)**2  / beta_perp(m)  /(k*cos(theta))**2 * expBes*&
                      & (   (omega-n*mu(m)*q(m)) *(omega*beta_ratio(m)-&
                      &     (beta_ratio(m)-1.0)*n*mu(m)*q(m))* dZ_func(zeta1)+&
                      &      k*cos(theta)*drift(m)*n*mu(m)*q(m) * dZ_func(zeta1)  -&
                      &     2*n*q(m)*sqrt(mu(m)*dens(m)/beta_para(m))*k*cos(theta)*drift(m)**2 *Z_func(zeta1))-&
                      
                      &     dens(m)**2 * q(m)**2  / beta_perp(m)  /(k*cos(theta))**2 * expBes*&
                      & (   (omega+n*mu(m)*q(m)) *(omega*beta_ratio(m)+&
                      &     (beta_ratio(m)-1.0)*n*mu(m)*q(m))* dZ_func(zeta2)-&
                      &      k*cos(theta)*drift(m)*n*mu(m)*q(m) * dZ_func(zeta2)  +&
                      &     2*n*q(m)*sqrt(mu(m)*dens(m)/beta_para(m))*k*cos(theta)*drift(m)**2 *Z_func(zeta2))

                 if((n<=4).or.((n>4).and.((abs(real(del_zz)/ real(epsilon_zz))>eps_error).or. &
                      & (abs(aimag(del_zz)/aimag(epsilon_zz))>eps_error)))) then

                    epsilon_zz=epsilon_zz+ del_zz

                 else
                    esc(3)=.false.
                 endif
              endif



              if(esc(4)) then

                 del_xy = i*sqrt(mu(m))*dens(m)**1.5 *q(m)**2 /(k*cos(theta)) / sqrt(beta_para(m)) *n*&
                      & (expdBes-expBes) *(beta_ratio(m) *(omega-k*cos(theta)*drift(m))-&
                      & (beta_ratio(m) -1.0)*n*mu(m)*q(m))* Z_func(zeta1) -&
                      
                      &   i*sqrt(mu(m))*dens(m)**1.5 *q(m)**2 /(k*cos(theta)) / sqrt(beta_para(m)) *n*&
                      & (expdBes-expBes) *(beta_ratio(m) *(omega-k*cos(theta)*drift(m))+&
                      & (beta_ratio(m) -1.0)*n*mu(m)*q(m))* Z_func(zeta2)

                 if((n<=4).or.((n>4).and.((abs(real(del_xy)/ real(epsilon_xy))>eps_error).or. &
                      & (abs(aimag(del_xy)/aimag(epsilon_xy))>eps_error)))) then

                    epsilon_xy=epsilon_xy+ del_xy

                 else
                    esc(4)=.false.
                 endif
              endif



              if(esc(5)) then

                 del_xz = -mu(m)*dens(m)**2 *q(m)**3 / beta_perp(m)/ (k*sin(theta)) / (k*cos(theta))* n*expBes*&
                      &(    (beta_ratio(m) * omega-n*mu(m)*q(m)*(beta_ratio(m)-1.0))*dZ_func(zeta1)-&
                      &     2*n*sqrt(dens(m)*mu(m)/beta_para(m))*q(m)*drift(m)*Z_func(zeta1) )+&
                      
                      
                      & mu(m)*dens(m)**2 *q(m)**3 / beta_perp(m)/ (k*sin(theta)) / (k*cos(theta))* n*expBes*&
                      & (    (beta_ratio(m) * omega+n*mu(m)*q(m)*(beta_ratio(m)-1.0))*dZ_func(zeta2)+&           
                      &     2*n*sqrt(dens(m)*mu(m)/beta_para(m))*q(m)*drift(m)*Z_func(zeta2) )

                 if((n<=4).or.((n>4).and.((abs(real(del_xz)/ real(epsilon_xz))>eps_error).or. &
                      & (abs(aimag(del_xz)/aimag(epsilon_xz))>eps_error)))) then

                    epsilon_xz=epsilon_xz+ del_xz

                 else
                    esc(5)=.false.
                 endif
              endif




              if(esc(6)) then

                 del_yz = 0.5*i *dens(m)* q(m) * ( k*sin(theta))/(k*cos(theta))* (expdBes-expBes)  *&
                      &  (   (beta_ratio(m)*omega-n*mu(m)*q(m)*(beta_ratio(m)-1.0))*dZ_func(zeta1)-&
                      &      2*n*sqrt(dens(m)*mu(m)/beta_para(m))*q(m)*drift(m)*Z_func(zeta1) )+&    
                      
                      &   0.5*i *dens(m)* q(m) * ( k*sin(theta))/(k*cos(theta))* (expdBes-expBes) *&
                      &  (   (beta_ratio(m)*omega+n*mu(m)*q(m)*(beta_ratio(m)-1.0)) *dZ_func(zeta2)+&
                      &     2*n*sqrt(dens(m)*mu(m)/beta_para(m))*q(m)*drift(m)*Z_func(zeta2) )

                 if((n<=4).or.((n>4).and.((abs(real(del_yz)/ real(epsilon_yz))>eps_error).or. &
                      & (abs(aimag(del_yz)/aimag(epsilon_yz))>eps_error)))) then

                    epsilon_yz=epsilon_yz+ del_yz

                 else
                    esc(6)=.false.
                 endif
              endif

           endif

           n=n+1


        enddo


     else if(mode(m)==1) then
   
       !for particles with arbitrary gyrotropic velocity distributions

        iarb=iarb+1

        allocate(Ivpe(3,5,nperp(iarb)))

        esc=.true.

        n=0

        do while(esc(1).or.esc(2).or.esc(3).or.esc(4).or.esc(5).or.esc(6))


           !perform integration over perpendicular velocity

           do iperp=1,nperp(iarb)

              call Bessel_int(n,m,iarb,k,iperp,esc,Ivpe_dummy)

              Ivpe(1,1,iperp)=Ivpe_dummy(1,1)
              Ivpe(1,2,iperp)=Ivpe_dummy(1,2)
              Ivpe(1,3,iperp)=Ivpe_dummy(1,3)
              Ivpe(1,4,iperp)=Ivpe_dummy(1,4)
              Ivpe(1,5,iperp)=Ivpe_dummy(1,5)

              Ivpe(2,1,iperp)=Ivpe_dummy(2,1)
              Ivpe(2,2,iperp)=Ivpe_dummy(2,2)
              Ivpe(2,3,iperp)=Ivpe_dummy(2,3)
              Ivpe(2,4,iperp)=Ivpe_dummy(2,4)
              Ivpe(2,5,iperp)=Ivpe_dummy(2,5)

              Ivpe(3,1,iperp)=Ivpe_dummy(3,1)
              Ivpe(3,2,iperp)=Ivpe_dummy(3,2)
              Ivpe(3,3,iperp)=Ivpe_dummy(3,3)
              Ivpe(3,4,iperp)=Ivpe_dummy(3,4)
              Ivpe(3,5,iperp)=Ivpe_dummy(3,5)

           enddo



           !perform integration over parallel velocity and compute total value of the 2d velocity integrals required 
           !for the evaluation of the dielectric tensor

           call integrator(omega,k,m,iarb,n,Ivpe,intgrl,splcoeff1(:,:,:,:,iarb),splcoeff2(:,:,:,:,iarb))

           if(n/=0) then
              call integrator(omega,k,m,iarb,-n,Ivpe,intgrl2,splcoeff1(:,:,:,:,iarb),splcoeff2(:,:,:,:,iarb))
           endif

           

           !compute dielectric tensor components

           if (n==0) then

              epsilon_xx=epsilon_xx+ &
                   & 2*pi *mu(m)**3 *q(m)**4 *&
                   & dens(m) * omega/(k*sin(theta))**2 *n*n*&
                   & (  intgrl(1,1,1,1)-&
                   &    (k*cos(theta))/omega*intgrl(1,2,1,1)+&
                   &    (k*cos(theta))/omega*intgrl(1,1,2,2))

              epsilon_yy=epsilon_yy+ &
                   & 2*pi *mu(m)*q(m)**2 *dens(m) * omega *&
                   & (  intgrl(2,1,1,1)-&
                   & (k*cos(theta))/omega*intgrl(2,2,1,1)+&
                   & (k*cos(theta))/omega*intgrl(2,1,2,2))

              epsilon_zz=epsilon_zz+ &
                   & 2*pi *mu(m) *q(m)**2 * dens(m)* omega *&
                   & (  (1-n*mu(m)*q(m)/omega)*intgrl(1,2,2,2)+&
                   & n*mu(m)*q(m)/omega*intgrl(1,3,1,1))

              epsilon_xy=epsilon_xy+&
                   & 2 *pi*i * mu(m)**2 * q(m)**3 *&
                   & dens(m)*n *omega/(k*sin(theta))*&
                   & ( intgrl(3,1,1,1)-&
                   &    (k*cos(theta))/omega*intgrl(3,2,1,1)+&
                   &    (k*cos(theta))/omega*intgrl(3,1,2,2))

              epsilon_xz=epsilon_xz+&
                   & 2*pi * mu(m)**2 * q(m)**3 *&
                   & dens(m) *n* omega/(k*sin(theta))*&
                   & (  (1-n*mu(m)*q(m)/omega)*intgrl(1,1,2,2)+&
                   & n*mu(m)*q(m)/omega*intgrl(1,2,1,1))

              epsilon_yz=epsilon_yz+&
                   & (-2)*pi*i *mu(m)*q(m)**2 * dens(m)*omega*&
                   & (  (1-n*mu(m)*q(m)/omega)*intgrl(3,1,2,2)+&
                   & n*mu(m)*q(m)/omega*intgrl(3,2,1,1))

           else

              if(esc(1)) then
                 del_xx=2*pi *mu(m)**3 *q(m)**4 *&
                      & dens(m) * omega/(k*sin(theta))**2 *n*n*&
                      & (  intgrl(1,1,1,1)-&
                      &    (k*cos(theta))/omega*intgrl(1,2,1,1)+&
                      &    (k*cos(theta))/omega*intgrl(1,1,2,2))+&
                      
                      & 2*pi *mu(m)**3 *q(m)**4 *&
                      & dens(m) * omega/(k*sin(theta))**2 *n*n*&
                      & (  intgrl2(1,1,1,1)-&
                      &    (k*cos(theta))/omega*intgrl2(1,2,1,1)+&
                      &    (k*cos(theta))/omega*intgrl2(1,1,2,2))


                 if((n<=4).or.((n>4).and.((abs(real(del_xx)/ real(epsilon_xx))>eps_error).or. &
                      & (abs(aimag(del_xx)/aimag(epsilon_xx))>eps_error)))) then

                    epsilon_xx=epsilon_xx+ del_xx

                 else
                    esc(1)=.false.
                 endif
              endif

              if(esc(2)) then
                 del_yy=2*pi *mu(m)*q(m)**2 *dens(m) * omega *&
                      & (  intgrl(2,1,1,1)-&
                      & (k*cos(theta))/omega*intgrl(2,2,1,1)+&
                      & (k*cos(theta))/omega*intgrl(2,1,2,2))+&
                      
                      & 2*pi *mu(m)*q(m)**2 *dens(m) * omega *&
                      & (  intgrl2(2,1,1,1)-&
                      & (k*cos(theta))/omega*intgrl2(2,2,1,1)+&
                      & (k*cos(theta))/omega*intgrl2(2,1,2,2))


                 if((n<=4).or.((n>4).and.((abs(real(del_yy)/ real(epsilon_yy))>eps_error).or. &
                      & (abs(aimag(del_yy)/aimag(epsilon_yy))>eps_error)))) then
                    epsilon_yy=epsilon_yy+ del_yy
                 else
                    esc(2)=.false.
                 endif
              endif

              if(esc(3)) then
                 del_zz=2*pi *mu(m) *q(m)**2 * dens(m)* omega *&
                      & (  (1-n*mu(m)*q(m)/omega)*intgrl(1,2,2,2)+&
                      & n*mu(m)*q(m)/omega*intgrl(1,3,1,1))+&
                      
                      & 2*pi *mu(m) *q(m)**2 * dens(m)* omega *&
                      & (  (1+n*mu(m)*q(m)/omega)*intgrl2(1,2,2,2)-&
                      & n*mu(m)*q(m)/omega*intgrl2(1,3,1,1))


                 if((n<=4).or.((n>4).and.((abs(real(del_zz)/ real(epsilon_zz))>eps_error).or. &
                      & (abs(aimag(del_zz)/aimag(epsilon_zz))>eps_error)))) then
                    epsilon_zz=epsilon_zz+ del_zz
                 else
                    esc(3)=.false.
                 endif
              endif

              if(esc(4)) then
                 del_xy=2 *pi*i * mu(m)**2 * q(m)**3 *&
                      & dens(m)*n *omega/(k*sin(theta))*&
                      & ( intgrl(3,1,1,1)-&
                      &    (k*cos(theta))/omega*intgrl(3,2,1,1)+&
                      &    (k*cos(theta))/omega*intgrl(3,1,2,2))+& 
                      
                      & (-2) *pi*i * mu(m)**2 * q(m)**3 *&
                      & dens(m)*n *omega/(k*sin(theta))*&
                      & ( intgrl2(3,1,1,1)-&
                      &    (k*cos(theta))/omega*intgrl2(3,2,1,1)+&
                      &    (k*cos(theta))/omega*intgrl2(3,1,2,2)) 

                 if((n<=4).or.((n>4).and.((abs( real(del_xy)/ real(epsilon_xy))>eps_error).or. &
                      & (abs(aimag(del_xy)/aimag(epsilon_xy))>eps_error)))) then
                    epsilon_xy=epsilon_xy+ del_xy
                 else
                    esc(4)=.false.
                 endif
              endif

              if(esc(5)) then
                 del_xz=2*pi * mu(m)**2 * q(m)**3 *&
                      & dens(m) *n* omega/(k*sin(theta))*&
                      & (  (1-n*mu(m)*q(m)/omega)*intgrl(1,1,2,2)+&
                      & n*mu(m)*q(m)/omega*intgrl(1,2,1,1))+&
                      
                      & (-2)*pi * mu(m)**2 * q(m)**3 *&
                      & dens(m) *n* omega/(k*sin(theta))*&
                      & (  (1+n*mu(m)*q(m)/omega)*intgrl2(1,1,2,2)-&
                      & n*mu(m)*q(m)/omega*intgrl2(1,2,1,1))

                 if(abs(del_xz%re) >= huge(1.0)*min(1.0, abs(epsilon_xz%re)))then
                    tdelr= huge(1.0)
                 else
                    tdelr= abs(real(del_xz)/ real(epsilon_xz))
                 end if
                 if(abs(del_xz%im) >= huge(1.0)*min(1.0, abs(epsilon_xz%im)))then
                    tdeli= huge(1.0)
                 else
                    tdeli= abs(aimag(del_xz)/aimag(epsilon_xz))
                 end if
                 if((n<=4).or.((n>4) .and. ((tdelr>eps_error) .or. (tdeli>eps_error)))) then
                    epsilon_xz=epsilon_xz+ del_xz
                 else
                    esc(5)=.false.
                 endif
              endif

              if(esc(6)) then
                 del_yz=(-2)*pi*i *mu(m)*q(m)**2 * dens(m)*omega*&
                      & (  (1-n*mu(m)*q(m)/omega)*intgrl(3,1,2,2)+&
                      & n*mu(m)*q(m)/omega*intgrl(3,2,1,1))+&
                      
                      & (-2)*pi*i *mu(m)*q(m)**2 * dens(m)*omega*&
                      & (  (1+n*mu(m)*q(m)/omega)*intgrl2(3,1,2,2)-&
                      & n*mu(m)*q(m)/omega*intgrl2(3,2,1,1))

                 if((n<=4).or.((n>4).and.((abs(real(del_yz)/ real(epsilon_yz))>eps_error).or. &
                      & (abs(aimag(del_yz)/aimag(epsilon_yz))>eps_error)))) then
                    epsilon_yz=epsilon_yz+ del_yz
                 else
                    esc(6)=.false.
                 endif
              endif

           endif

           n=n+1

        enddo

        deallocate(Ivpe)

     else
        write(*,*) 'Check dist_in!'
        stop
     endif

  enddo


  !the dispersion relation for a collisionless plasma can be generally written as:
  !0=det(kk - 1k^2 + omega^2/c^2 * epsilon)

  !to rewrite this in the form given below, use  k_x = k_para = k*cos(theta), k_y = 0, k_z = k_perp = k*sin(theta) 
  !and the symmetry relations of the dielectric tensor


  disp_det=(epsilon_xx-(k*cos(theta))**2)  *  (epsilon_yy-k**2)  *  (epsilon_zz-(k*sin(theta))**2 ) +&
       &    2*epsilon_xy*epsilon_yz*(epsilon_xz+k**2 * sin(theta)*cos(theta) ) -&
       &   (epsilon_yy -k**2 )  *  (epsilon_xz + k**2 * sin(theta)*cos(theta) )**2 +&
       &   (epsilon_xx-(k*cos(theta))**2 )* epsilon_yz**2 +&
       &   (epsilon_zz - (k*sin(theta))**2  )*epsilon_xy**2


!  write(*,'(E20.10,E20.10,E20.10,E20.10)') real(disp_det), aimag(disp_det), real(omega), aimag(omega)


end function disp_det_rg
