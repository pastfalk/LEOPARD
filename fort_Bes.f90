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

  if(abs(Bes_arg)<10.0 **(-14)) then
     
     write(*,*) 'abs(Bes_arg) very low - check fort_Bes...'
     stop

  endif

  if(n<=10) then
     Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
     Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))  

  else if((n>10).and.(n<=12)) then
     if(abs(Bes_arg)<10.0 **(-11)) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>12).and.(n<=15)) then

     if(abs(Bes_arg)<10.0**(-8)) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>15).and.(n<=20)) then

     if(abs(Bes_arg)<10.0**(-6)) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>20).and.(n<=25)) then

     if(abs(Bes_arg)<10.0**(-4)) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>25).and.(n<=30)) then

     if(abs(Bes_arg)<10.0**(-3)) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>30).and.(n<=40)) then

     if(abs(Bes_arg)<10.0**(-2)) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>40).and.(n<=50)) then

     if(abs(Bes_arg)<10.0**(-1)) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif


  else if((n>50).and.(n<=70)) then

     if(abs(Bes_arg)<0.5) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>70).and.(n<=80)) then

     if(abs(Bes_arg)<1.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>80).and.(n<=110)) then

     if(abs(Bes_arg)<4.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>110).and.(n<=140)) then

     if(abs(Bes_arg)<10.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>140).and.(n<=180)) then

     if(abs(Bes_arg)<20.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>180).and.(n<=230)) then

     if(abs(Bes_arg)<40.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>230).and.(n<=270)) then

     if(abs(Bes_arg)<60.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>270).and.(n<=310)) then

     if(abs(Bes_arg)<80.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>310).and.(n<=350)) then

     if(abs(Bes_arg)<100.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>350).and.(n<=430)) then

     if(abs(Bes_arg)<150.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif


  else if((n>430).and.(n<=510)) then

     if(abs(Bes_arg)<200.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>510).and.(n<=580)) then

     if(abs(Bes_arg)<250.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>580).and.(n<=650)) then

     if(abs(Bes_arg)<300.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>650).and.(n<=730)) then

     if(abs(Bes_arg)<360.0) then
        Bes_term1=0.0
        Bes_term2=0.0
     else
        Bes_term1=BesJn(n,abs(Bes_arg))*BesJn(n,abs(Bes_arg))
        Bes_term2=BesJn(n-1,abs(Bes_arg))*BesJn(n+1,abs(Bes_arg))

     endif

  else if((n>730).and.(n<=840)) then

     if(abs(Bes_arg)<450.0) then
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
