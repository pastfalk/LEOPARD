!> Rootfinding algorithm based on the Muller method presented, e.g., in Gerald & Wheatley (2003)
!! \param omega_start initial frequency guess as starting value for the iteration
!! \param k wavenumber
!! \param solsk approximate roots of the dispersion relation obtained from Muller iteration
!! \param splcoeff1 array of coefficients obtained from cubic spline interpolation of the provided velocity distribution data 
!! \param splcoeff2 array of coefficients obtained from cubic spline interpolation of the provided velocity distribution data 

! Modified in July/2022 by Rudi Gaelzer

subroutine muller_rg(logu, nroots, omega_start, k, solsk, splcoeff1, splcoeff2)
use param_mod
implicit none
integer, intent(in) :: logu, nroots
real :: k
complex, dimension(nroots) :: omega_start
complex, dimension(nroots) :: solsk
complex :: omega_ini, disp_det_rg
real :: tfxr, tfxi
complex, dimension(1:4) :: fx
complex, dimension(1:4) :: omega
complex :: a, b, c, d1, d2, fzero
integer :: n, j, izero

real :: om, ga, om_old, ga_old

real, dimension(npara_max-1,nperp_max-1,4,3,narb) :: splcoeff1, splcoeff2

integer, dimension(:), allocatable :: vmit

external disp_det_rg

allocate(vmit(0))

do izero= 1, nroots

!Muller's method requires three starting points
!one point is chosen by the user - the other two points are taken to be slightly left and right of it 
   omega_ini= omega_start(izero)
   omega(1)=0.999*omega_ini
   omega(2)=omega_ini
   omega(3)=1.001*omega_ini

!    fx(1)=disp_det(omega(1),k,splcoeff1, splcoeff2)
!    fx(2)=disp_det(omega(2),k,splcoeff1, splcoeff2)
!    fx(3)=disp_det(omega(3),k,splcoeff1, splcoeff2)
   call dflate(izero, omega(1), solsk, fzero, fx(1))! ; PRINT*, '1:', FX(1)
   call dflate(izero, omega(2), solsk, fzero, fx(2))! ; PRINT*, '2:', FX(2)
   call dflate(izero, omega(3), solsk, fzero, fx(3))! ; PRINT*, '3:', FX(3)

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

!       fx(4)=disp_det(omega(4),k,splcoeff1, splcoeff2)
      call dflate(izero, omega(4), solsk, fzero, fx(4))! ; PRINT*, '4:', OMEGA(4), FX(4)

      !measure the accuracy of iterated root and check exit-condition

      om=real(omega(4))
      ga=aimag(omega(4))

      om_old=real(omega(3))
      ga_old=aimag(omega(3))


      if(   ( ( ((om.ge.om_old).and.(abs(1.0-abs(om_old/om)).lt.rf_error)).or.&
            &    ((om.lt.om_old).and.(abs(1.0-abs(om/om_old)).lt.rf_error))).and.&
            &  ( ((ga.ge.ga_old).and.(abs(1.0-abs(ga_old/ga)).lt.rf_error)).or.&
            &    ((ga.lt.ga_old).and.(abs(1.0-abs(ga/ga_old)).lt.rf_error)))).or.&
            &( ( ((om.ge.om_old).and.(abs(1.0-abs(om_old/om)).lt.rf_error)).or.&
            &    ((om.lt.om_old).and.(abs(1.0-abs(om/om_old)).lt.rf_error))).and.&
            &  ( (abs(ga).lt.10.0**(-10)).and.(abs(ga_old).lt.10.0**(-10)))).or.&
            &( ( ((ga.ge.ga_old).and.(abs(1.0-abs(ga_old/ga)).lt.rf_error)).or.&
            &    ((ga.lt.ga_old).and.(abs(1.0-abs(ga/ga_old)).lt.rf_error))).and.&
            &  ( (abs(om).lt.10.0**(-10)).and.(abs(om_old).lt.10.0**(-10)))).or.&
            &( ( (abs(om).lt.10.0**(-10)).and.(abs(ga).lt.10.0**(-10)))    )  ) exit


      !stop iteration if last step was ineffective
! PRINT*, 'TESTE:', abs((real(fx(4))-real(fx(3)))/real(fx(4))), abs((aimag(fx(4))-aimag(fx(3)))/aimag(fx(4)))
      
      if(abs(fx(3)%re) >= huge(1.0)*min(1.0, abs(fx(4)%re)))then
         tfxr= huge(1.0)
      else
         tfxr= abs((real(fx(4))-real(fx(3)))/real(fx(4)))
      end if
      if(abs(fx(3)%im) >= huge(1.0)*min(1.0, abs(fx(4)%im)))then
         tfxi= huge(1.0)
      else
         tfxi= abs((aimag(fx(4))-aimag(fx(3)))/aimag(fx(4)))
      end if
      
      if((tfxr .lt. 1.0e-12) .and. (tfxi .lt. 1.0e-12)) then
         write(logu,*) 'Last step in Muller iteration was ineffective'
         exit
      endif
!       if( (abs((real(fx(4))-real(fx(3)))/real(fx(4))) .lt. 1.0e-12).and. &
!             ( abs((aimag(fx(4))-aimag(fx(3)))/aimag(fx(4))) .lt. 1.0e-12)) then
!          write(logu,*) 'Last step in Muller iteration was ineffective'
!          exit
!       endif

      do j=1,3
         omega(j)=omega(j+1)
         fx(j)=fx(j+1)
      enddo

!       if(n.gt.40) then
!          write(logu, fmt= '(a,i0)') 'Error: Muller method did not converge for root: ', izero
!          exit
!       endif
      if(n > maxit)then
         vmit= [ vmit, izero ] ; exit
      end if

      n=n+1

   enddo
   
!solution of root finding procedure
! PRINT*, 'MULLER: W(4)=', OMEGA(4)
   solsk(izero)= omega(4)

end do

if(size(vmit) > 0)write(logu, fmt= '(a, i0, a, *(x, :, i0))')'Warning: Muller did not converge in ', &
                                                             maxit, ' steps for roots:', vmit
deallocate(vmit)
return

CONTAINS
   subroutine dflate(izero, zero, zeros, fzero, fzrdfl)
   integer, intent(in) :: izero
   complex, intent(in) :: zero
   complex, dimension(nroots), intent(inout) :: zeros
   complex, intent(out) :: fzero, fzrdfl
   integer :: j
   complex :: den
! PRINT*, 'DFLATE: ZERO=', ZERO
   fzero= disp_det_rg(zero, k, splcoeff1, splcoeff2)
! PRINT*, 'DFLATE: FZERO=', FZERO
   fzrdfl= fzero
   do j= 2, izero
! PRINT*, 'ENTERED LOOP'
      den= zero - zeros(j - 1)
      if(abs(den) == 0.0e0)then
         zeros(izero)= zero*1.001
         return
      else
         fzrdfl= fzrdfl/den
      end if
   end do
   return
   end subroutine dflate

end subroutine muller_rg
