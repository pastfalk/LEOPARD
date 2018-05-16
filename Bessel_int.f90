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
  real(kind=16), dimension(3) :: F23_3

  real(kind=16), allocatable, dimension(:) :: a
  real(kind=16), allocatable, dimension(:) :: b

  integer :: na,nb
  logical, dimension(6) :: esc

  real :: Bes_term1, Bes_term2


  if(vperp(iperp,iarb).eq.0.0) then
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


  if ( n.ne.0) then


     !compute a prefactor required for the evaluation of the perpendicular velocity integrals
     x=lambda*vperp(iperp,iarb)
     
     h1=1.0_16

     do j=1,n
        h1=h1*(x/(1.0_16 *j))**2
     enddo



     !compute generalized hypergeometric functions 2F3 and 1F2, and use them for evaluating the perpendicular velocity integrals

     na=2
     nb=3
     allocate(a(na),b(nb))


     if(esc(1).or.esc(3).or.esc(5).or.esc(4).or.esc(6)) then

        a(1)=1.0_16 *n+0.5_16
        a(2)=1.0_16 *n+0.5_16
        b(1)=1.0_16 *n+1.0_16
        b(2)=1.0_16 *n+1.5_16
        b(3)=2.0_16 *n+1.0_16

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        a(1)=1.0_16 *n+0.5_16
        a(2)=1.0_16 *n+1.5_16
        b(1)=1.0_16 *n+1.0_16
        b(2)=1.0_16 *n+2.5_16
        b(3)=2.0_16 *n+1.0_16

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(2))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(2))
        endif

        a(1)=1.0_16 *n+0.5_16
        a(2)=1.0_16 *n+2.5_16
        b(1)=1.0_16 *n+1.0_16
        b(2)=1.0_16 *n+3.5_16
        b(3)=2.0_16 *n+1.0_16

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(3))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(3))
        endif


        Ivpe(1,1)=vperp(iperp,iarb) *h1 /&
             & (2.0**(2*n) * (2*n+1)) *&
             & F23_3(1)


        if ((vperp(iperp,iarb).ne.0.0).and.((Bes_term1-Bes_term2).ne.0.0).and.&
             & ((log10(abs(0.5*vperp(iperp,iarb)**2))+log10(abs(Bes_term1-Bes_term2))).gt.-299.0)) then

           Ivpe(1,2)=0.5*vperp(iperp,iarb)*vperp(iperp,iarb) *&
                & (Bes_term1-Bes_term2)
        else
           Ivpe(1,2)=0.0
        endif


        Ivpe(1,3)=vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)*h1 /&
             & (2.0**(2*n) * (2*n+3)) *&
             & F23_3(2)

        Ivpe(1,5)=vperp(iperp,iarb)**5 * h1 /&
             & (2.0**(2*n) * (2*n+5) ) *&
             & F23_3(3)


        a(1)=1.0_16 *n+0.5_16 
        a(2)=1.0_16 *n+2.0_16 
        b(1)=1.0_16 *n+1.0_16 
        b(2)=1.0_16 *n+3.0_16 
        b(3)=2.0_16 *n+1.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif


        Ivpe(1,4)=vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)*h1 /&
             & (2.0**(2*n) * (2*n+4)) *&
             & F23_3(1)


     else
        Ivpe(1,1)=0.0
        Ivpe(1,2)=0.0
        Ivpe(1,3)=0.0
        Ivpe(1,4)=0.0
        Ivpe(1,5)=0.0
     endif





     if(esc(2)) then

        a(1)=1.0_16 *n+0.5_16 
        a(2)=1.0_16 *n+1.5_16 
        b(1)=1.0_16 *n 
        b(2)=1.0_16 *n+2.5_16 
        b(3)=2.0_16 *n+1.0_16

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        a(1)=1.0_16 *n+0.5_16 
        a(2)=1.0_16 *n+2.5_16 
        b(1)=1.0_16 *n 
        b(2)=1.0_16 *n+3.5_16 
        b(3)=2.0_16 *n+1.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(2))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(2))
        endif

        a(1)=1.0_16 *n+0.5_16 
        a(2)=1.0_16 *n+3.5_16 
        b(1)=1.0_16 *n 
        b(2)=1.0_16 *n+4.5_16 
        b(3)=2.0_16 *n+1.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(3))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(3))
        endif

        Ivpe(2,1)=-vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)*n*h1 /&
             & (2.0**(2*n) * (2*n+3)) *&
             & F23_3(1)

        Ivpe(2,3)=vperp(iperp,iarb)**5 /3 *n *h1 /&
             & (2.0**(2*n) * (2*n+5)) *&
             & F23_3(2)

        Ivpe(2,5)=vperp(iperp,iarb)**7 /5 *n*h1 /&
             & (2.0**(2*n) *(2*n+7)) *&
             & F23_3(3)



        a(1)=1.0_16 *n-0.5_16 
        a(2)=1.0_16 *n+0.5_16 
        b(1)=1.0_16 *n
        b(2)=1.0_16 *n+1.5_16 
        b(3)=2.0_16 *n-1.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        a(1)=1.0_16 *n-0.5_16 
        a(2)=1.0_16 *n+1.5_16 
        b(1)=1.0_16 *n
        b(2)=1.0_16 *n+2.5_16 
        b(3)=2.0_16 *n-1.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(2))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(2))
        endif

        a(1)=1.0_16 *n-0.5_16 
        a(2)=1.0_16 *n+2.5_16 
        b(1)=1.0_16 *n
        b(2)=1.0_16 *n+3.5_16 
        b(3)=2.0_16 *n-1.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(3))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(3))
        endif


        Ivpe(2,1)=Ivpe(2,1)+vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)*n*n*h1/arg /&
             & (2.0**(2*n) * (2*n+1)) *&
             & F23_3(1)

        Ivpe(2,3)=Ivpe(2,3)+vperp(iperp,iarb)**5 *n*n*h1/arg /&
             & (2.0**(2*n) * (2*n+3)) *&
             & F23_3(2)

        Ivpe(2,5)=Ivpe(2,5)+vperp(iperp,iarb)**7 *n*n*h1/arg /&
             & (2.0**(2*n) * (2*n+5)) *&
             & F23_3(3)



        a(1)=1.0_16 *n+1.5_16 
        a(2)=1.0_16 *n+2.5_16 
        b(1)=1.0_16 *n+2.0_16 
        b(2)=1.0_16 *n+3.5_16 
        b(3)=2.0_16 *n+3.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        a(1)=1.0_16 *n+1.5_16 
        a(2)=1.0_16 *n+3.5_16 
        b(1)=1.0_16 *n+2.0_16 
        b(2)=1.0_16 *n+4.5_16 
        b(3)=2.0_16 *n+3.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(2))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(2))
        endif

        a(1)=1.0_16 *n+1.5_16 
        a(2)=1.0_16 *n+4.5_16 
        b(1)=1.0_16 *n+2.0_16 
        b(2)=1.0_16 *n+5.5_16 
        b(3)=2.0_16 *n+3.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(3))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(3))
        endif



        Ivpe(2,1)=Ivpe(2,1)+vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)*h1*arg /&
             & (2.0**(2*n+4) * (2*n+5)*(n+1)*&
             & (n+1)) *&
             & F23_3(1)

        Ivpe(2,3)=Ivpe(2,3)+vperp(iperp,iarb)**5 * h1*arg /&
             & (2.0**(2*n+4) * (2*n+7)*&
             & (n+1)*(n+1)) *&
             & F23_3(2)

        Ivpe(2,5)=Ivpe(2,5)+vperp(iperp,iarb)**7 * h1*arg /&
             & (2.0**(2*n+4) * (2*n+9)*&
             & (n+1)*(n+1)) *&
             & F23_3(3)




        a(1)=1.0_16 *n+1.5_16 
        a(2)=1.0_16 *n+3.0_16 
        b(1)=1.0_16 *n+2.0_16 
        b(2)=1.0_16 *n+4.0_16 
        b(3)=2.0_16 *n+3.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        a(1)=1.0_16 *n+1.5_16 
        a(2)=1.0_16 *n+4.0_16 
        b(1)=1.0_16 *n+2.0_16 
        b(2)=1.0_16 *n+5.0_16 
        b(3)=2.0_16 *n+3.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(2))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(2))
        endif

        Ivpe(2,2)=vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)*h1*arg /&
             & (2.0**(2*n+4) * (2*n+6)*(n+1)*&
             & (n+1)) *&
             & F23_3(1)

        Ivpe(2,4)=vperp(iperp,iarb)**6 * h1*arg /&
             & (2.0**(2*n+4) * (2*n+8)*&
             & (n+1)*(n+1)) *&
             & F23_3(2)



        a(1)=1.0_16 *n+0.5_16 
        a(2)=1.0_16 *n+1.0_16 
        b(1)=1.0_16 *n
        b(2)=1.0_16 *n+2.0_16 
        b(3)=2.0_16 *n+1.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        a(1)=1.0_16 *n+0.5_16 
        a(2)=1.0_16 *n+1.0_16 
        b(1)=1.0_16 *n
        b(2)=1.0_16 *n+3.0_16 
        b(3)=2.0_16 *n+1.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(2))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(2))
        endif

        Ivpe(2,3)=Ivpe(2,3)-vperp(iperp,iarb)**5 /3 *h1 *n /&
             & (2.0**(2*n+1) *(n+1)) *&
             & F23_3(1)

        Ivpe(2,4)=Ivpe(2,4)-vperp(iperp,iarb)**6 *n* h1 /&
             & (2.0**(2*n+2) *(2*n+2)) *&
             & F23_3(1)

        Ivpe(2,5)=Ivpe(2,5)-vperp(iperp,iarb)**7 /5 *n* h1 /&
             & (2.0**(2*n) *(2*n+2)) *&
             & F23_3(1)

        Ivpe(2,2)=Ivpe(2,2)-vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)*n*h1 /&
             & (2.0**(2*n+1) * (2*n+4)*(n+1)) *&
             & F23_3(2)



        a(1)=1.0_16 *n-0.5_16 
        a(2)=1.0_16 *n+1.0_16 
        b(1)=1.0_16 *n
        b(2)=1.0_16 *n+2.0_16 
        b(3)=2.0_16 *n-1.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        a(1)=1.0_16 *n-0.5_16 
        a(2)=1.0_16 *n+2.0_16 
        b(1)=1.0_16 *n
        b(2)=1.0_16 *n+3.0_16 
        b(3)=2.0_16 *n-1.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(2))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(2))
        endif

        Ivpe(2,2)=Ivpe(2,2)+vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)*n*n*h1/arg /&
             & (2.0**(2*n) * (2*n+2)) *&
             & F23_3(1)

        Ivpe(2,4)=Ivpe(2,4)+vperp(iperp,iarb)**6 *n*n*h1/arg /&
             & (2.0**(2*n) * (2*n+4)) *&
             & F23_3(2)


        a(1)=1.0_16 *n+1.0_16 
        a(2)=1.0_16 *n+1.5_16 
        b(1)=1.0_16 *n
        b(2)=1.0_16 *n+2.0_16 
        b(3)=2.0_16 *n+1.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        Ivpe(2,1)=Ivpe(2,1)+vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)*&
             & (2*n+1)*n*h1 /&
             & (2.0**(2*n+1) * (n+1) ) *&
             & F23_3(1)

        a(1)=1.0_16 *n+1.5_16 
        a(2)=1.0_16 *n+1.5_16 
        b(1)=1.0_16 *n
        b(2)=1.0_16 *n+2.5_16 
        b(3)=2.0_16 *n+1.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        Ivpe(2,1)=Ivpe(2,1)-vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)*&
             & (2*n+1)*n*h1 /&
             & (2.0**(2*n) * (2*n+3)) *&
             & F23_3(1)

        a(1)=1.0_16 *n+0.5_16 
        a(2)=1.0_16 *n+3.0_16 
        b(1)=1.0_16 *n
        b(2)=1.0_16 *n+4.0_16 
        b(3)=2.0_16 *n+1.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        Ivpe(2,4)=Ivpe(2,4)+vperp(iperp,iarb)**6 *n*h1 /&
             & (2.0**(2*n+2) *(2*n+6)) *&
             & F23_3(1)


     else
        Ivpe(2,1)=0.0
        Ivpe(2,2)=0.0
        Ivpe(2,3)=0.0
        Ivpe(2,4)=0.0
        Ivpe(2,5)=0.0

     endif





     if(esc(4).or.esc(6)) then


        if ((vperp(iperp,iarb).ne.0.0).and.(Bes_term1.ne.0.0).and.&
             & ((log10(abs(vperp(iperp,iarb)/(2*lambda)))+log10(abs(Bes_term1))).gt.-299.0)) then

           Ivpe(3,1)=vperp(iperp,iarb)/(2 *lambda)*&
                & Bes_term1-&
                & 1.0 /(2*lambda) * Ivpe(1,1)          

        else
           Ivpe(3,1)=0.0
        endif


        if ((vperp(iperp,iarb).ne.0.0).and.(Bes_term1.ne.0.0).and.&
             & ((log10(abs(vperp(iperp,iarb)**2/(2*lambda)))+log10(abs(Bes_term1))).gt.-299.0)) then

           Ivpe(3,2)=vperp(iperp,iarb)*vperp(iperp,iarb)/(2*lambda)*&
                & Bes_term1-&
                & 2.0 /(2*lambda) * Ivpe(1,2)          

        else
           Ivpe(3,2)=0.0
        endif




        if ((vperp(iperp,iarb).ne.0.0).and.(Bes_term1.ne.0.0).and.&
             & ((log10(abs(vperp(iperp,iarb)**3/(2*lambda)))+log10(abs(Bes_term1))).gt.-299.0)) then

           Ivpe(3,3)=vperp(iperp,iarb)**3 /(2 *lambda)*&
                & Bes_term1-&
                & 3.0 /(2*lambda) * Ivpe(1,3)          
        else
           Ivpe(3,3)=0.0
        endif



        if ((vperp(iperp,iarb).ne.0.0).and.(Bes_term1.ne.0.0).and.&
             & ((log10(abs(vperp(iperp,iarb)**4/(2*lambda)))+log10(abs(Bes_term1))).gt.-299.0)) then

           Ivpe(3,4)=vperp(iperp,iarb)**4 /(2 *lambda)*&
                & Bes_term1-&
                & 4.0 /(2*lambda) * Ivpe(1,4)          
        else
           Ivpe(3,4)=0.0
        endif


        if ((vperp(iperp,iarb).ne.0.0).and.(Bes_term1.ne.0.0).and.&
             & ((log10(abs(vperp(iperp,iarb)**5/(2*lambda)))+log10(abs(Bes_term1))).gt.-299.0)) then
           Ivpe(3,5)=vperp(iperp,iarb)**5 /(2 *lambda)*&
                & Bes_term1-&
                & 5.0 /(2*lambda) * Ivpe(1,5)   
        else
           Ivpe(3,5)=0.0
        endif



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

        Ivpe(1,1)=0.0


        if ((vperp(iperp,iarb).ne.0.0).and.((Bes_term1-Bes_term2).ne.0.0).and.&
             & ((log10(abs(0.5*vperp(iperp,iarb)*vperp(iperp,iarb)))+log10(abs(Bes_term1-Bes_term2))).gt.-299.0)) then
           Ivpe(1,2)= 0.5*vperp(iperp,iarb)*vperp(iperp,iarb) *&
                &(Bes_term1-Bes_term2)  
        else
           Ivpe(1,2)=0.0
        endif



        na=2
        nb=3
        allocate(a(na),b(nb))

        a(1)=0.5_16 
        a(2)=1.5_16 
        b(1)=1.0_16 
        b(2)=1.0_16 
        b(3)=2.5_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        a(1)=0.5_16 
        a(2)=2.5_16 
        b(1)=1.0_16 
        b(2)=1.0_16 
        b(3)=3.5_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(2))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(2))
        endif


        Ivpe(1,3)= vperp(iperp,iarb)*vperp(iperp,iarb)* vperp(iperp,iarb)/3.0 *F23_3(1)

        Ivpe(1,5)= vperp(iperp,iarb)**5 /5.0 * F23_3(2)

        a(1)=0.5_16 
        a(2)=2.0_16 
        b(1)=1.0_16 
        b(2)=1.0_16 
        b(3)=3.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        Ivpe(1,4)= vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)/4.0 *&
             & F23_3(1)

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



        a(1)=1.5_16 
        a(2)=2.5_16 
        b(1)=2.0_16 
        b(2)=3.0_16 
        b(3)=3.5_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        a(1)=1.5_16 
        a(2)=3.5_16 
        b(1)=2.0_16 
        b(2)=3.0_16 
        b(3)=4.5_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(2))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(2))
        endif

        a(1)=1.5_16 
        a(2)=4.5_16 
        b(1)=2.0_16 
        b(2)=3.0_16 
        b(3)=5.5_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(3))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(3))
        endif

        Ivpe(2,1)= arg*vperp(iperp,iarb)*vperp(iperp,iarb)*vperp(iperp,iarb)/20.0 *F23_3(1)

        Ivpe(2,3)= arg*vperp(iperp,iarb)**5 /28.0 *F23_3(2)

        Ivpe(2,5)= arg *vperp(iperp,iarb)**7 /36.0 *F23_3(3)



        a(1)=0.5_16 
        a(2)=3.0_16 
        b(1)=1.0_16 
        b(2)=2.0_16 
        b(3)=4.0_16 

        if(ndp.eq.0) then
           call F23(a,b,-arg,nfrac,F23_3(1))
        else
           call F23_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        Ivpe(2,4)= 0.5*vperp(iperp,iarb)**6 * F23_3(1)


        deallocate(a,b)

        na=1
        nb=2
        allocate(a(na),b(nb))

        a(1)=1.5_16 
        b(1)=2.0_16 
        b(2)=4.0_16 

        if(ndp.eq.0) then
           call F12(a,b,-arg,nfrac,F23_3(1))
        else
           call F12_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        Ivpe(2,2)= arg*vperp(iperp,iarb)*vperp(iperp,iarb)*&
             & vperp(iperp,iarb)*vperp(iperp,iarb)/24.0 *&
             & F23_3(1)

        a(1)=0.5_16
        b(1)=1.0_16
        b(2)=2.0_16

        if(ndp.eq.0) then
           call F12(a,b,-arg,nfrac,F23_3(1))
        else
           call F12_mpfun(a,b,-arg,ndp,nfrac,F23_3(1))
        endif

        Ivpe(2,4)=Ivpe(2,4)- 0.5*vperp(iperp,iarb)**6  *F23_3(1)

        deallocate(a,b)



     else
        Ivpe(2,1)=0.0
        Ivpe(2,2)=0.0
        Ivpe(2,3)=0.0
        Ivpe(2,4)=0.0
        Ivpe(2,5)=0.0

     endif

     if(esc(4).or.esc(6)) then

        Ivpe(3,1)=0.0


        if ((vperp(iperp,iarb).ne.0.0).and.(Bes_term1.ne.0.0).and.&
             & ((log10(abs(vperp(iperp,iarb)*vperp(iperp,iarb)/(2*lambda)))+log10(abs(Bes_term1))).gt.-299.0)) then

           Ivpe(3,2)=vperp(iperp,iarb)*vperp(iperp,iarb)/(2*lambda)*&
                & Bes_term1-&
                & 2.0/(2*lambda) * Ivpe(1,2)

        else
           Ivpe(3,2)=0.0
        endif


        if ((vperp(iperp,iarb).ne.0.0).and.(Bes_term1.ne.0.0).and.&
             & ((log10(abs(vperp(iperp,iarb)**3/(2*lambda)))+log10(abs(Bes_term1))).gt.-299.0)) then

           Ivpe(3,3)=vperp(iperp,iarb)**3 /(2*lambda)*&
                & Bes_term1-&
                & 3.0/(2*lambda) * Ivpe(1,3)
        else
           Ivpe(3,3)=0.0
        endif


        if ((vperp(iperp,iarb).ne.0.0).and.(Bes_term1.ne.0.0).and.&
             & ((log10(abs(vperp(iperp,iarb)**4/(2*lambda)))+log10(abs(Bes_term1))).gt.-299.0)) then

           Ivpe(3,4)=vperp(iperp,iarb)**4 /(2*lambda)*&
                & Bes_term1-&
                & 4.0/(2*lambda) * Ivpe(1,4)
        else
           Ivpe(3,4)=0.0
        endif


        if ((vperp(iperp,iarb).ne.0.0).and.(Bes_term1.ne.0.0).and.&
             & ((log10(abs(vperp(iperp,iarb)**5/(2*lambda)))+log10(abs(Bes_term1))).gt.-299.0)) then

           Ivpe(3,5)=vperp(iperp,iarb)**5 /(2 *lambda)*&
                & Bes_term1-&
                & 5.0/(2 *lambda) * Ivpe(1,5)

        else
           Ivpe(3,5)=0.0
        endif

     else

        Ivpe(3,1)=0.0
        Ivpe(3,2)=0.0
        Ivpe(3,3)=0.0
        Ivpe(3,4)=0.0
        Ivpe(3,5)=0.0

     endif

  endif

end subroutine Bessel_int
