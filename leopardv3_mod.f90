!> Module containing modified routines
module leopardv3_mod
implicit none
integer :: maxit= 40 ! Save # of Muller steps as the original version

CONTAINS
   include 'muller_rg.f90'
end module leopardv3_mod
