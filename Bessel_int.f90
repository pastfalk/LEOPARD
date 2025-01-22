!> Computes perpendicular velocity integrals over the interval [vperp(iperp),vperp(iperp+1)] using generalized hypergeomtric function / Bessel function identities
!! \param n Bessel index
!! \param m index of particle species
!! \param iarb index of species among the particle species with arbitrary velocity distribution
!! \param k wavenumber
!! \param iperp index of perpendicular velocity gridpoint of the interval (lower limit) which is to be integrated over
!! \param esc boolean array indicating which dielectric tensor components are still not converged with respect to summation over the Bessel index 
!! \param Ivpe array containing the values of the perpendicular velocity integrals for the given velocity interval
subroutine Bessel_int(n,m,iarb,k,iperp,esc,Ivpe)
  use param_mod
  implicit none
  real :: k
  integer :: m,n,ndp,nfrac,iarb
  integer :: iperp
  real, dimension(3,5) :: Ivpe
  real :: lambda
  real(kind=16) :: arg
  real(kind=16) :: h1, x
  integer :: j
  real(kind=16), dimension(4) :: HGF
  integer :: il


  real(kind=16), allocatable, dimension(:) :: a
  real(kind=16), allocatable, dimension(:) :: b

  integer :: na,nb
  logical, dimension(6) :: esc

  real :: Bes_term1, Bes_term2


  if(vperp(iperp,iarb)==0.0) then
     Ivpe=0.0
     return
  endif

  lambda=(k*sin(theta))/(mu(m)*q(m))
  arg=lambda*lambda*vperp(iperp,iarb)*vperp(iperp,iarb)


  !determine how many terms have to be included in the continued fraction formula that is used for 
  !the computation of the generalized hypergeometric function needed for the evaluation of the perpendicular velocity integral
  !and determine whether arbitrary precision arithmetic has to be employed or not
  call acc_F(arg,ndp,nfrac)

  !compute the Bessel function terms required for the estimation of some of the perpendicular velocity integrals
  call fort_Bes(n,lambda*vperp(iperp,iarb),Bes_term1,Bes_term2)


  if ( n/=0) then


     !compute a prefactor required for the evaluation of the perpendicular velocity integrals
     x=lambda*vperp(iperp,iarb)

     h1=1.0_16

     do j=1,n
        h1=h1*(x/(2.0_16 *j))
     enddo
     h1=h1*h1


     !compute generalized hypergeometric functions 2F3 and 1F2, and use them for evaluating the perpendicular velocity integrals

     na=2
     nb=3
     allocate(a(na),b(nb))

     if(esc(1).or.esc(3).or.esc(5).or.esc(4).or.esc(6)) then

        do il=1,5

           if(il/=2) then

              a(1)=1.0_16 *n+0.5_16
              a(2)=1.0_16 *n+il*0.5_16
              b(1)=1.0_16 *n+1.0_16
              b(2)=1.0_16 *n+(il+2)*0.5_16
              b(3)=2.0_16 *n+1.0_16

              if(ndp==0) then
                 call F23(a,b,-arg,nfrac,HGF(1))
              else
                 call F23_mpfun(a,b,-arg,ndp,nfrac,HGF(1))
              endif

              Ivpe(1,il)=vperp(iperp,iarb)**il *h1/(2.0*n+1.0*il) *HGF(1)

           else

              if ((vperp(iperp,iarb)/=0.0).and.((Bes_term1-Bes_term2)/=0.0).and.&
                   & ((log10(abs(0.5*vperp(iperp,iarb)**2))+log10(abs(Bes_term1-Bes_term2)))>-299.0)) then

                 Ivpe(1,2)=0.5*vperp(iperp,iarb)*vperp(iperp,iarb) *&
                      & (Bes_term1-Bes_term2)
              else
                 Ivpe(1,2)=0.0
              endif

           endif

        enddo

     else

        Ivpe(1,1)=0.0
        Ivpe(1,2)=0.0
        Ivpe(1,3)=0.0
        Ivpe(1,4)=0.0
        Ivpe(1,5)=0.0

     endif


     if(esc(2)) then

        a(1)=1.0_16 *n+0.5_16
        a(2)=1.0_16 *n+1.0_16
        b(1)=1.0_16 *n
        b(2)=1.0_16 *n+2.0_16
        b(3)=2.0_16 *n+1.0_16

        if(ndp==0) then
           call F23(a,b,-arg,nfrac,HGF(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,HGF(1))
        endif

        do il=1,5

           a(1)=1.0_16 *n+0.5_16
           a(2)=1.0_16 *n+(il+2)*0.5_16
           b(1)=1.0_16 *n
           b(2)=1.0_16 *n+(il+4)*0.5_16
           b(3)=2.0_16 *n+1.0_16

           if(ndp==0) then
              call F23(a,b,-arg,nfrac,HGF(2))
           else
              call F23_mpfun(a,b,-arg,ndp,nfrac,HGF(2))
           endif

           a(1)=1.0_16 *n-0.5_16
           a(2)=1.0_16 *n+il*0.5_16
           b(1)=1.0_16 *n
           b(2)=1.0_16 *n+(il+2)*0.5_16
           b(3)=2.0_16 *n-1.0_16

           if(ndp==0) then
              call F23(a,b,-arg,nfrac,HGF(3))
           else
              call F23_mpfun(a,b,-arg,ndp,nfrac,HGF(3))
           endif

           a(1)=1.0_16 *n+1.5_16
           a(2)=1.0_16 *n+(il+4)*0.5_16
           b(1)=1.0_16 *n+2.0_16
           b(2)=1.0_16 *n+(il+6)*0.5_16
           b(3)=2.0_16 *n+3.0_16

           if(ndp==0) then
              call F23(a,b,-arg,nfrac,HGF(4))
           else
              call F23_mpfun(a,b,-arg,ndp,nfrac,HGF(4))
           endif

           Ivpe(2,il)=-vperp(iperp,iarb)**(il+2) /(2.0*il)*n*h1 /(1.0*(n+1)) * HGF(1)+&
                &      vperp(iperp,iarb)**(il+2)/(1.0*il)*n*h1 /(1.0*(2*n+il+2)) * HGF(2)+&
                &      vperp(iperp,iarb)**(il+2)*n*n*h1/arg /(1.0*(2*n+il)) * HGF(3)+&
                &      vperp(iperp,iarb)**(il+2)*h1*arg/(16.0*(2*n+il+4)*(n+1)*(n+1)) *HGF(4)

        enddo

     else

        Ivpe(2,1)=0.0
        Ivpe(2,2)=0.0
        Ivpe(2,3)=0.0
        Ivpe(2,4)=0.0
        Ivpe(2,5)=0.0

     endif


     if(esc(4).or.esc(6)) then

        do il=1,5

           if ((vperp(iperp,iarb)/=0.0).and.(Bes_term1/=0.0).and.&
                & ((log10(abs(vperp(iperp,iarb)**il/(2*lambda)))+log10(abs(Bes_term1)))>-299.0)) then

              Ivpe(3,il)=vperp(iperp,iarb)**il/(2 *lambda)* Bes_term1-  il*Ivpe(1,il)/(2*lambda)

           else
              Ivpe(3,il)=0.0
           endif

        enddo

     else

        Ivpe(3,1)=0.0
        Ivpe(3,2)=0.0
        Ivpe(3,3)=0.0
        Ivpe(3,4)=0.0
        Ivpe(3,5)=0.0

     endif

     deallocate(a,b)

  else

     if(esc(1).or.esc(3).or.esc(4).or.esc(5).or.esc(6)) then

        Ivpe(1,1)=0.0 !not needed

        if ((vperp(iperp,iarb)/=0.0).and.((Bes_term1-Bes_term2)/=0.0).and.&
             & ((log10(abs(0.5*vperp(iperp,iarb)*vperp(iperp,iarb)))+log10(abs(Bes_term1-Bes_term2)))>-299.0)) then
           Ivpe(1,2)= 0.5*vperp(iperp,iarb)*vperp(iperp,iarb) *&
                &(Bes_term1-Bes_term2)  
        else
           Ivpe(1,2)=0.0
        endif

        na=2
        nb=3

        allocate(a(na),b(nb))

        do il=3,5

           a(1)=0.5_16 
           a(2)=il*0.5_16 
           b(1)=1.0_16 
           b(2)=1.0_16 
           b(3)=(il+2)*0.5_16 

           if(ndp==0) then
              call F23(a,b,-arg,nfrac,HGF(1))
           else
              call F23_mpfun(a,b,-arg,ndp,nfrac,HGF(1))
           endif

           Ivpe(1,il)=vperp(iperp,iarb)**il  /(1.0*il) *HGF(1)

        enddo

        deallocate(a,b)

     else

        Ivpe(1,1)=0.0
        Ivpe(1,2)=0.0
        Ivpe(1,3)=0.0
        Ivpe(1,4)=0.0
        Ivpe(1,5)=0.0

     endif


     if(esc(2)) then

        na=2
        nb=3
        allocate(a(na),b(nb))

        do il=1,5

           if(il/=2) then

              a(1)=1.5_16
              a(2)=(il+4)*0.5_16
              b(1)=2.0_16
              b(2)=3.0_16
              b(3)=(il+6)*0.5_16

              if(ndp==0) then
                 call F23(a,b,-arg,nfrac,HGF(1))
              else
                 call F23_mpfun(a,b,-arg,ndp,nfrac,HGF(1))
              endif

              Ivpe(2,il)= arg*vperp(iperp,iarb)**(il+2) /(4.0*(il+4)) *HGF(1)

           else

              deallocate(a,b)
              na=1
              nb=2
              allocate(a(na),b(nb))

              a(1)=1.5_16
              b(1)=2.0_16
              b(2)=4.0_16

              if(ndp==0) then
                 call F12(a,b,-arg,nfrac,HGF(1))
              else
                 call F12_mpfun(a,b,-arg,ndp,nfrac,HGF(1))
              endif

              Ivpe(2,il)= arg*vperp(iperp,iarb)**(il+2) /(4.0*(il+4)) *HGF(1)

              deallocate(a,b)
              na=2
              nb=3
              allocate(a(na),b(nb))

           endif

        enddo

        deallocate(a,b)

     else

        Ivpe(2,1)=0.0
        Ivpe(2,2)=0.0
        Ivpe(2,3)=0.0
        Ivpe(2,4)=0.0
        Ivpe(2,5)=0.0

     endif


     if(esc(4).or.esc(6)) then

        Ivpe(3,1)=0.0 !not needed

        do il=2,5

           if ((vperp(iperp,iarb)/=0.0).and.(Bes_term1/=0.0).and.&
                & ((log10(abs(vperp(iperp,iarb)**il/(2*lambda)))+log10(abs(Bes_term1)))>-299.0)) then

              Ivpe(3,il)=vperp(iperp,iarb)**il /(2*lambda)*Bes_term1- il*1.0/(2*lambda) * Ivpe(1,il)
           else
              Ivpe(3,il)=0.0
           endif

        enddo

     else

        Ivpe(3,1)=0.0
        Ivpe(3,2)=0.0
        Ivpe(3,3)=0.0
        Ivpe(3,4)=0.0
        Ivpe(3,5)=0.0

     endif

  endif


end subroutine Bessel_int
