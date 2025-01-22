!> Computes the (exponentially-scaled) modified Bessel function of the first kind using its series expansion
!! \param n index of the modified Bessel function
!! \param z argument of the modified Bessel function
function exp_Bessel_In(n,z_in)
  use param_mod
  implicit none
  real(kind=16) :: expBes, old, z
  real :: error_old, error, z_in
  real :: exp_Bessel_In
  integer :: n, j
  real(kind=16) :: fac1, fac2
  integer :: m
  integer :: zfrac

  z=z_in

  expBes=0.0_16
  error_old=1.0

  do j=0,10000

     old=expBes

     if(abs(n)==j) then
        expBes=expBes+2.0_16 **abs(n)

     else if(abs(n)>j) then

        fac1=1.0_16

        do m=j+1,abs(n)
           fac1=fac1*(m*1.0_16/z)
        enddo
        expBes=expBes+fac1*2.0_16 **abs(n)

     else if(j>abs(n)) then

        fac1=1.0_16

        do m=abs(n)+1,j
           fac1=fac1*(z/(1.0_16*m))
        enddo
        expBes=expBes+fac1*2.0_16 **abs(n)

     endif

     error=abs(1.0-old/expBes)

     if (error<10.0**(-15)) then
        expBes=1.0_16/expBes
        exit
     else if(expBes>10.0_16**300) then
        expBes=0.0_16
        exit
     endif

  enddo


  do m=1,10000 

     zfrac=2*m+abs(n)


     fac1=1.0_16

     do j=1,m

        fac1=fac1*(z*exp(-z/zfrac)/(2.0_16*j))

     enddo

     fac2=1.0_16

     do j=1,abs(n)

        fac2=fac2*(z*exp(-z/zfrac)/(2*(m+j)))
     enddo

     if(fac1*fac1*fac2<10.0_16**(-300)) then
        cycle
     else
        old=expBes
        expBes=expBes+fac1*fac1*fac2
     endif

     error=abs(1.0-old/expBes)

     if (error<10.0**(-15)) then
        exit
     endif

     error_old=error

  enddo


  exp_Bessel_In=expBes


end function exp_Bessel_In

