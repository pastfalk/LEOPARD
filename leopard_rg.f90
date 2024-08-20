!> Initializes the setup, scans through the requested wavenumber interval, computes corresponding frequencies, and prints dispersion relation to output file
! Modified in July/2022 by Rudi Gaelzer
program leopard_rg
use param_mod
use leopardv3_mod
implicit none
character(len= *), parameter :: leop= 'arq/leopardv3'
character(len= 3) :: file_num, set_num
character(len= 20) :: input_file, log_file, dat_file
complex :: omega_init, increment
integer :: nk,ik, iarb, comml, datu, logu, nroots, iw

real :: kstart, kend, dk
real, allocatable, dimension (:) :: krange !allocatable permite alocar apenas o espaco necessario, de acordo com a programacao
complex, allocatable, dimension (:,:) :: solution
complex, allocatable, dimension(:) :: omega_start, solsk

real, allocatable, dimension(:,:,:,:,:) :: splcoeff1, splcoeff2

real :: start, finish
real :: start2, finish2

! Get data file number and set up the file names.
if(command_argument_count() == 0)then
   print*, 'ERROR: You have to provide the number of the input data file.'
   print*, '   Command line: $>./leopard_rg ***  (where * = 0 - 9)'
   print*, '   Format of input file: leopard***.ent (where * = 0 - 9)'
   stop
end if
call get_command_argument(number= 1, value= file_num, length= comml)
if(comml /= 3)then
   print*, 'ERROR: The data file number must have 3 digits: *** (* = 0 - 9)'
   stop
end if
input_file= leop//file_num//'.ent'
log_file=   leop//file_num//'.log'
dat_file=   leop//file_num//'.dat'

! Proceed with the computation
open(newunit= datu, file= dat_file, status= 'unknown')
open(newunit= logu, file= log_file)

call cpu_time(start)

write(logu, *) 'Read input data'
call read_data_rg(input_file, set_num, nroots, omega_init, increment, kstart, kend, nk)
write(logu, *) '...done'

write(logu, *) 'Read velocity distributions from files'
call read_distr_rg(set_num, logu)
write(logu, *) '...done.'

! Initialize vectors of initial guesses and solutions
allocate(omega_start(nroots)) ; omega_start= omega_init
allocate(solsk(nroots))
allocate(solution(nk, nroots))

allocate(krange(nk))
dk=(kend-kstart)/(1.0*nk)
do ik=1,nk
   krange(ik)=kstart+(ik-1)*dk
enddo


!spline-interpolate the velocity distributions

allocate(splcoeff1(npara_max-1,nperp_max-1,4,3,narb))
allocate(splcoeff2(npara_max-1,nperp_max-1,4,3,narb))

do iarb=1,narb
   call get_splinecoeff(iarb,splcoeff1(:,:,:,:,iarb),splcoeff2(:,:,:,:,iarb))
enddo



!scan through wavenumber interval
do ik=1,nk

   write(logu, *) ' '
   write(logu, fmt= '(A7,I6,A10,F12.8)') '-------',ik,'------- k=', krange(ik)

   call cpu_time(start2)

   !use Muller method to iterate root of dispersion relation
   call muller_rg(logu, nroots, omega_start, krange(ik), solsk, splcoeff1, splcoeff2)
   solution(ik, :)= solsk

   call cpu_time(finish2)

!    write(logu, '(A9,E20.10,A9,E20.10)')  '   omega:', real(solution(ik)), '   gamma:',aimag(solution(ik))
   write(logu, *) 'time elapsed:', finish2-start2


   if ((ik .ge. 3).and.(ik .lt. nk))  then

      !if three subsequent solutions omega(k) are found, use quadratic polynomial fit 
      !to guess next starting frequency for Muller iteration
      do iw= 1, nroots
         call polyfit(krange(ik-2:ik+1), solution(ik-2:ik, iw), omega_start(iw))
      end do

   else

      !for the first two solution omega(k) guess next starting frequency for Muller iteration
      !by raising the computed omega by an increment which is provided by the user
      omega_start= solsk + increment

   end if

!    write(datu, fmt= '(F12.8,E20.10,E20.10)') krange(ik), real(solution(ik)), aimag(solution(ik))
   write(datu, fmt= '(F8.5, *(x, :, es12.5))') krange(ik), solsk

enddo

call cpu_time(finish)

write(logu, fmt= '(/, 2(a,g0))') 'Total time elapsed(s):', finish-start, '  (min): ', (finish-start)/60.0

! close(7)
close(logu) ; close(datu)

deallocate(krange, solution, omega_start, solsk)
deallocate(mu,q)
deallocate(beta_para,beta_perp,beta_ratio)
deallocate(splcoeff1,splcoeff2)
deallocate(mode, dens, drift)
deallocate(distribution,vpara,vperp,npara,nperp)

end program leopard_rg
