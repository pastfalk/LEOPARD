!> Determines whether or not the parallel velocity integration requires arbitrary precision arithmetic, and how many significant digits have to be considered
!! \param dvpa resolution of the parallel velocity grid
!! \param zeta parameter showing up in the parallel velocity integration
!! \param ndp measure for the number of significant digits to be included for cases where arbitrary precision arithmetic is required (i.e. cases with ndp >0)
subroutine acc_Kvpa(dvpa,zeta,ndp)
  implicit none
  integer :: ndp
  real :: zeta
  real :: dvpa

  if(abs(zeta)<0.1) then
     ndp=0
  else if((abs(zeta).ge.0.1).and.(abs(zeta)<1.0)) then
     if(dvpa<0.001) then
        ndp=30
     else
        ndp=0
     endif

  else if((abs(zeta).ge.1.0).and.(abs(zeta)<10.0)) then
     if(dvpa<0.01) then
        ndp=30
     else
        ndp=0
     endif

  else if((abs(zeta).ge.10.0).and.(abs(zeta)<100.0)) then
     if(dvpa<0.001) then
        ndp=44
     else if((dvpa.ge.0.001).and.(dvpa<0.1)) then
        ndp=30
     else
        ndp=0
     endif

  else if((abs(zeta).ge.100.0).and.(abs(zeta)<1000.0)) then
     if(dvpa<0.01) then
        ndp=44
     else if((dvpa.ge.0.01).and.(dvpa<1.0)) then
        ndp=30
     else
        ndp=0
     endif

  else if((abs(zeta).ge.1000.0).and.(abs(zeta)<10000.0)) then
     if(dvpa<0.001) then
        ndp=58
     else if((dvpa.ge.0.001).and.(dvpa<0.1)) then
        ndp=44
     else if((dvpa.ge.0.1).and.(dvpa<10.0)) then
        ndp=30
     else
        ndp=0
     endif

  else if((abs(zeta).ge.10000.0).and.(abs(zeta)<100000.0)) then
     if(dvpa<0.01) then
        ndp=58
     else if((dvpa.ge.0.01).and.(dvpa<1.0)) then
        ndp=44
     else if((dvpa.ge.1.0).and.(dvpa<100.0)) then
        ndp=30
     else
        ndp=0
     endif

  else if((abs(zeta).ge.100000.0).and.(abs(zeta)<1000000.0)) then
     if(dvpa<0.1) then
        ndp=58
     else if((dvpa.ge.0.1).and.(dvpa<10.0)) then
        ndp=44
     else if((dvpa.ge.10.0).and.(dvpa<100.0)) then
        ndp=30
     else
        write(*,*) 'careful with Kvpa now, vpara is pretty large'
        ndp=0
     endif

  else if((abs(zeta).ge.1000000.0).and.(abs(zeta)<5000000.0)) then
     if(dvpa<0.001) then
        ndp=73
     else if((dvpa.ge.0.001).and.(dvpa<0.1)) then
        ndp=58
     else if((dvpa.ge.0.1).and.(dvpa<100.0)) then
        ndp=44
     else if((dvpa.ge.100.0).and.(dvpa<1000.0)) then
        ndp=30
     else
        write(*,*) 'careful with Kvpa now, vpara is pretty large'
        ndp=0
     endif

  else if((abs(zeta).ge.5000000.0).and.(abs(zeta)<10000000.0)) then
     if(dvpa<0.001) then
        ndp=73
     else if((dvpa.ge.0.001).and.(dvpa<1.0)) then
        ndp=58
     else if((dvpa.ge.1.0).and.(dvpa<100.0)) then
        ndp=44
     else if((dvpa.ge.100.0).and.(dvpa<1000.0)) then
        ndp=30
     else
        write(*,*) 'careful with Kvpa now, vpara is pretty large'
        ndp=0
     endif

  else
     write(*,*) 'Zeta too large!'
     write(*,*) 'Manually adjust ndp in acc_Kvpa().'
     stop
  endif

end subroutine acc_Kvpa
