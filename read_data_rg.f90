!> Reads the input parameters from input.dat
!! \param omega_start first initial frequency guess for the initial wavenumber provided by input.dat
!! \param increment frequency increment suggested by the code user to determine the initial guesses for subsequent wavenumbers, before subroutine polyfit() takes over
!! \param kstart start value of the requested wavenumber interval
!! \param kend final value of the requested wavenumber interval
!! \param nk number of steps in the wavenumber interval provided by input.dat

! Modified in July/2022 by Rudi Gaelzer
!  Reads the input parameters from leopard***.ent, where * = 0 - 9
!  String set_num keeps the distribution grids fixed for different parameters, when it makes sense.

subroutine read_data_rg(inputf, set_num, nroots, omega_start, increment, kstart, kend, nk)
use param_mod
use leopardv3_mod, only: maxit
implicit none
character(len= *), intent(in) :: inputf
character(len= 3) :: set_num
real :: q_in, mu_in, dens_in
real :: drift_in, beta_para_in, beta_perp_in
integer :: mode_in, inputu, nroots
real :: kstart, kend
integer :: nk, n
complex :: omega_start, increment
real :: omega_r, omega_i, increment_r, increment_i

!define namelists and read input data using namelists

namelist /wavenumber/ &
   &  kstart, kend, nk

namelist /number_roots/ &
   &  nroots, maxit

namelist /initial_guess/ &
   &  omega_r, omega_i, increment_r, increment_i

! namelist /setup/ &
!    &  Nspecies, theta, delta
namelist /setup/ &
   &  set_num, Nspecies, theta, delta

namelist /accuracy/ &
   &  rf_error, eps_error

namelist /species/ &
      & mode_in, q_in, mu_in, dens_in,drift_in,&
      & beta_para_in, beta_perp_in

open(newunit= inputu, status='old', file= inputf)
read(inputu, wavenumber)
read(inputu, number_roots)
read(inputu, initial_guess)
read(inputu, setup)
read(inputu, accuracy)

allocate(mode(Nspecies),drift(Nspecies))
allocate(mu(Nspecies),q(Nspecies),dens(Nspecies))
allocate(beta_para(Nspecies),beta_perp(Nspecies),beta_ratio(Nspecies))

do n=1,Nspecies
   read(inputu, species)
   mode(n)=mode_in
   q(n)=q_in
   mu(n)=mu_in
   dens(n)=dens_in
   drift(n)=drift_in
   beta_para(n)=beta_para_in
   beta_perp(n)=beta_perp_in

   if(mode(n)/=1) then
      beta_ratio(n)=beta_perp(n)/beta_para(n)
   endif
   
enddo

close(inputu)

if(omega_r==0.0) omega_r=10.0**(-12) !algorithm cannot handle exact zero
if(omega_i==0.0) omega_i=10.0**(-12) !algorithm cannot handle exact zero
if(theta==0.0) theta=0.0001          !algorithm cannot handle exact zero


omega_start=omega_r+i*omega_i
increment=increment_r+i*increment_i
theta=theta*2*pi/360.0

end subroutine read_data_rg
