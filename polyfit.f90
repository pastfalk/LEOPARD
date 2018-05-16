!> Fits a quadratic polynomial to previously derived solutions omega(k) using the Vandermonde Matrix and predicts frequency for the subsequent wavenumber to give reliable starting value for the following iteration
!! \param x array containing the wavenumbers of the previously derived solutions and the subsequent wavenumber
!! \param y array containing the complex frequencies of the previously derived solutions
!! \param omega_start predicted frequency for the subsequent wavenumber which is used as initial guess for following iteration
subroutine polyfit(x,y,omega_start)
  implicit none
  real, dimension(3,3) :: V, Vinv, adj
  real :: det

  real, dimension(4) :: x
  complex, dimension(3) :: y

  complex :: omega_start
  complex :: a1,a2,a3

  !note that x -> wavenumber k, y -> frequency omega

  !determine components of Vandermonde matrix
  V(1,1)=x(1)**2
  V(1,2)=x(1)
  V(1,3)=1.0
  V(2,1)=x(2)**2
  V(2,2)=x(2)
  V(2,3)=1.0
  V(3,1)=x(3)**2
  V(3,2)=x(3)
  V(3,3)=1.0

  !compute determinant of Vandermonde matrix
  det =  (V(1,1)*V(2,2)*V(3,3) + V(1,2)*V(2,3)*V(3,1) + V(1,3)*V(2,1)*V(3,2)) -&
       & (V(1,3)*V(2,2)*V(3,1) + V(1,1)*V(2,3)*V(3,2) + V(1,2)*V(2,1)*V(3,3))
  
  !compute adjugate of Vandermonde matrix
  adj(1,1) = V(2,2)*V(3,3)-V(2,3)*V(3,2)
  adj(1,2) = V(2,3)*V(3,1)-V(2,1)*V(3,3)
  adj(1,3) = V(2,1)*V(3,2)-V(2,2)*V(3,1)
  adj(2,1) = V(1,3)*V(3,2)-V(1,2)*V(3,3)
  adj(2,2) = V(1,1)*V(3,3)-V(1,3)*V(3,1)
  adj(2,3) = V(1,2)*V(3,1)-V(1,1)*V(3,2)
  adj(3,1) = V(1,2)*V(2,3)-V(1,3)*V(2,2)
  adj(3,2) = V(1,3)*V(2,1)-V(1,1)*V(2,3)
  adj(3,3) = V(1,1)*V(2,2)-V(1,2)*V(2,1)
  
  adj=transpose(adj)

  !determine inverse of Vandermonde matrix using relation V^-1 = 1/det * adj(V)
  Vinv = adj/det

  !want quadratic polynomial of form p = a1*x^2 + a2*x + a3 
  !have y = V(X) * a
  !want a = V^-1 * y

  !determine coefficients a
  a1=Vinv(1,1)*y(1)+Vinv(1,2)*y(2)+Vinv(1,3)*y(3)
  a2=Vinv(2,1)*y(1)+Vinv(2,2)*y(2)+Vinv(2,3)*y(3)
  a3=Vinv(3,1)*y(1)+Vinv(3,2)*y(2)+Vinv(3,3)*y(3)
 
  !determine polynomial and use it to get omega approximation for subsequent k
  omega_start=a1*x(4)**2+a2*x(4)+a3


end subroutine polyfit
