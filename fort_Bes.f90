!> Computes the Bessel function terms required for the estimation of some of the perpendicular velocity integrals
!! \param n Bessel index
!! \param Bes_arg argument of the Bessel function
!! \param Bes_term1 estimate of J_n(x)**2
!! \param Bes_term2 estimate of J_n-1(x)*J_n+1(x)
subroutine fort_Bes(n,Bes_arg,Bes_term1,Bes_term2)
  implicit none
  integer :: n
  real :: Bes_arg
  real :: Bes_term1, Bes_term2

  if(abs(Bes_arg).lt.10.0 **(-14)) then
     
     write(*,*) 'abs(Bes_arg) very low - check fort_Bes...'
     stop

  endif

  if(n.le.10) then
     Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
     Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))  

  else if((n.gt.10).and.(n.le.12)) then
     if(abs(Bes_arg).lt.10.0 **(-11)) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.12).and.(n.le.15)) then

     if(abs(Bes_arg).lt.10.0**(-8)) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.15).and.(n.le.20)) then

     if(abs(Bes_arg).lt.10.0**(-6)) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.20).and.(n.le.25)) then

     if(abs(Bes_arg).lt.10.0**(-4)) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.25).and.(n.le.30)) then

     if(abs(Bes_arg).lt.10.0**(-3)) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.30).and.(n.le.40)) then

     if(abs(Bes_arg).lt.10.0**(-2)) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.40).and.(n.le.50)) then

     if(abs(Bes_arg).lt.10.0**(-1)) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif


  else if((n.gt.50).and.(n.le.70)) then

     if(abs(Bes_arg).lt.0.5) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.70).and.(n.le.80)) then

     if(abs(Bes_arg).lt.1.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.80).and.(n.le.110)) then

     if(abs(Bes_arg).lt.4.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.110).and.(n.le.140)) then

     if(abs(Bes_arg).lt.10.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.140).and.(n.le.180)) then

     if(abs(Bes_arg).lt.20.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.180).and.(n.le.230)) then

     if(abs(Bes_arg).lt.40.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.230).and.(n.le.270)) then

     if(abs(Bes_arg).lt.60.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.270).and.(n.le.310)) then

     if(abs(Bes_arg).lt.80.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.310).and.(n.le.350)) then

     if(abs(Bes_arg).lt.100.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.350).and.(n.le.430)) then

     if(abs(Bes_arg).lt.150.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif


  else if((n.gt.430).and.(n.le.510)) then

     if(abs(Bes_arg).lt.200.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.510).and.(n.le.580)) then

     if(abs(Bes_arg).lt.250.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.580).and.(n.le.650)) then

     if(abs(Bes_arg).lt.300.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.650).and.(n.le.730)) then

     if(abs(Bes_arg).lt.360.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n.gt.730).and.(n.le.840)) then

     if(abs(Bes_arg).lt.450.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else
     write(*,*) 'Bessel n very large - check fort_Bes'
     stop
  endif

end subroutine fort_Bes
