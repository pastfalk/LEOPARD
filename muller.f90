!> Rootfinding algorithm based on the Muller method presented, e.g., in Gerald & Wheatley (2003)
!! \param omega_start initial frequency guess as starting value for the iteration
!! \param k wavenumber
!! \param sol approximate root of the dispersion relation obtained from Muller iteration
!! \param splcoeff1 array of coefficients obtained from cubic spline interpolation of the provided velocity distribution data 
!! \param splcoeff2 array of coefficients obtained from cubic spline interpolation of the provided velocity distribution data 
subroutine muller(omega_start,k,sol,splcoeff1,splcoeff2)
  use param_mod
  implicit none

  real :: k
  complex :: omega_start
  complex :: sol
  complex :: disp_det
  complex, dimension(1:4) :: fx
  complex, dimension(1:4) :: omega
  complex :: a, b, c, d1, d2
  integer :: n, j

  real :: om, ga, om_old, ga_old

  real, dimension(npara_max-1,nperp_max-1,4,3,narb) :: splcoeff1, splcoeff2

  external disp_det


  !Muller's method requires three starting points
  !one point is chosen by the user - the other two points are taken to be slightly left and right of it 
  omega(1)=0.999*omega_start
  omega(2)=omega_start
  omega(3)=1.001*omega_start

  fx(1)=disp_det(omega(1),k,splcoeff1, splcoeff2)
  fx(2)=disp_det(omega(2),k,splcoeff1, splcoeff2)
  fx(3)=disp_det(omega(3),k,splcoeff1, splcoeff2)

  !perform Muller iteration

  n=0

  do while(.true.)

     !determine coefficients from preceding three points

     a=((omega(2)-omega(3))*(fx(1)-fx(3))-(omega(1)-omega(3))*(fx(2)-fx(3)))/&
          & ((omega(1)-omega(3))*(omega(2)-omega(3))*(omega(1)-omega(2)))
     b=((omega(1)-omega(3))**2 * (fx(2)-fx(3)) - (omega(2) - omega(3))**2 *&
          & (fx(1)-fx(3)))/((omega(1)-omega(3))*(omega(2)-omega(3))*(omega(1)-omega(2)))
     c=fx(3)

     d1=b+sqrt(b**2 - 4.0*a*c)
     d2=b-sqrt(b**2 - 4.0*a*c)

     !compute new root from coefficients

     if  (abs(d1) .GE. abs(d2)) then

        omega(4)=omega(3)-2.0*c/d1

     else

        omega(4)=omega(3)-2.0*c/d2

     endif

     fx(4)=disp_det(omega(4),k,splcoeff1, splcoeff2)

     !measure the accuracy of iterated root and check exit-condition

     om=real(omega(4))
     ga=aimag(omega(4))

     om_old=real(omega(3))
     ga_old=aimag(omega(3))


     if(   ( ( ((om.ge.om_old).and.(abs(1.0-abs(om_old/om))<rf_error)).or.&
          &    ((om<om_old).and.(abs(1.0-abs(om/om_old))<rf_error))).and.&
          &  ( ((ga.ge.ga_old).and.(abs(1.0-abs(ga_old/ga))<rf_error)).or.&
          &    ((ga<ga_old).and.(abs(1.0-abs(ga/ga_old))<rf_error)))).or.&
          &( ( ((om.ge.om_old).and.(abs(1.0-abs(om_old/om))<rf_error)).or.&
          &    ((om<om_old).and.(abs(1.0-abs(om/om_old))<rf_error))).and.&
          &  ( (abs(ga)<10.0**(-10)).and.(abs(ga_old)<10.0**(-10)))).or.&
          &( ( ((ga.ge.ga_old).and.(abs(1.0-abs(ga_old/ga))<rf_error)).or.&
          &    ((ga<ga_old).and.(abs(1.0-abs(ga/ga_old))<rf_error))).and.&
          &  ( (abs(om)<10.0**(-10)).and.(abs(om_old)<10.0**(-10)))).or.&
          &( ( (abs(om)<10.0**(-10)).and.(abs(ga)<10.0**(-10)))    )  ) exit


     !stop iteration if last step was ineffective
     if( (abs((real(fx(4))-real(fx(3)))/real(fx(4))) < 10.0 **(-12)).and. &
          ( abs((aimag(fx(4))-aimag(fx(3)))/aimag(fx(4))) < 10.0 **(-12))) then
        write(*,*) 'Last step in Muller iteration was ineffective'
        exit
     endif

     do j=1,3
        omega(j)=omega(j+1)
        fx(j)=fx(j+1)
     enddo

     if(n>40) then
        write(*,*) 'Error: Muller method did not converge'
        exit
     endif
     
     n=n+1

  enddo


  !solution of root finding procedure
  sol=omega(4)

end subroutine muller
