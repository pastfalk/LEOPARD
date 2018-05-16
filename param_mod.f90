!> Module containing all globally defined parameters
module param_mod
implicit none
complex, parameter :: i=(0.0,1.0)
real, parameter :: pi= 4*atan(1.0)
integer :: Nspecies
integer :: sign
real :: theta
integer, allocatable, dimension (:) :: mode
real, allocatable, dimension (:) :: mu
real, allocatable, dimension (:) :: q
real, allocatable, dimension (:) :: dens
real, allocatable, dimension (:) :: drift
real, allocatable, dimension (:) :: beta_para
real, allocatable, dimension (:):: beta_perp
real, allocatable, dimension (:) :: beta_ratio
real, allocatable, dimension (:,:,:) :: distribution
real, allocatable, dimension (:,:) :: vpara,vperp
integer, allocatable, dimension (:) :: npara, nperp
integer :: npara_max, nperp_max
integer :: narb
real :: delta 
real :: rf_error
real :: eps_error
end module param_mod
