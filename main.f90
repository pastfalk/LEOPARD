!> Initializes the setup, scans through the requested wavenumber interval, computes corresponding frequencies, and prints dispersion relation to output file
program main
  use param_mod
  implicit none
  complex :: omega_start, increment
  integer :: nk,ik, iarb

  real :: kstart, kend, dk
  real, allocatable, dimension (:) :: krange
  complex, allocatable, dimension (:) :: solution

  real, allocatable, dimension(:,:,:,:,:) :: splcoeff1, splcoeff2

  real :: start, finish
  real :: start2, finish2

  open(unit=7,status='unknown',file='omega.dat')

  call cpu_time(start)


  write(*,*) 'Read input data'
  call read_data(omega_start, increment, kstart, kend, nk)
  write(*,*) '...done'

  write(*,*) 'Read velocity distributions from files'
  call read_distr
  write(*,*) '...done.'
  
  allocate(krange(nk),solution(nk))
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

     write(*,*) ' '
     write(*,'(A7,I6,A10,F12.8)') '-------',ik,'------- k=', krange(ik)

     call cpu_time(start2)

     !use Muller method to iterate root of dispersion relation
     call muller(omega_start,krange(ik),solution(ik),splcoeff1,splcoeff2)

     call cpu_time(finish2)

     write(*,'(A9,E20.10,A9,E20.10)')  '   omega:', real(solution(ik)), '   gamma:',aimag(solution(ik))
     write(*,*) 'time elapsed:', finish2-start2




     if ((ik .ge. 3).and.(ik < nk))  then

        !if three subsequent solutions omega(k) are found, use quadratic polynomial fit 
        !to guess next starting frequency for Muller iteration
        call polyfit(krange(ik-2:ik+1),solution(ik-2:ik),omega_start)

     else

        !for the first two solution omega(k) guess next starting frequency for Muller iteration
        !by raising the computed omega by an increment which is provided by the user
        omega_start=solution(ik)+increment

     end if


     write(7,'(F12.8,E20.10,E20.10)') krange(ik), real(solution(ik)), aimag(solution(ik))

  enddo

  call cpu_time(finish)
 
  write(*,*) 'Total time elapsed:', finish-start
 
  close(7)

  deallocate(krange,solution)
  deallocate(mu,q)
  deallocate(beta_para,beta_perp,beta_ratio)
  deallocate(splcoeff1,splcoeff2)
  deallocate(mode, dens, drift)
  deallocate(distribution,vpara,vperp,npara,nperp)

end program main
