!> Reads the velocity distribution data from the provided distribution files
subroutine read_distr
  use param_mod
  implicit none
  integer :: ios
  real :: start_pe
  real :: vpa, vpe
  real :: dist_value
  integer :: ipara, iperp, iarb
  integer :: j
  character(len=1024) :: filename

  !count number of included species with arbitrary velocity distribution

  narb=0
  do j=1,Nspecies
     if(mode(j)==1) then
        narb=narb+1
     endif
  enddo

  if(narb/=0) then

     allocate(npara(narb),nperp(narb))

     !determine dimension of the velocity grid for each velocity distribution

     do iarb=1,narb

        write(filename,'(A25,I1,A4)') 'distribution/distribution', iarb,'.dat'
        open(unit=17,status='old',file=filename)

        npara(iarb)=0
        nperp(iarb)=0

        read(17,*,iostat=ios) vpa, start_pe, dist_value
        rewind(17)

        do while(.true.)

           read(17,*,iostat=ios) vpa, vpe, dist_value

           if (ios/=0) exit

           if (vpe==start_pe) then
              npara(iarb)=npara(iarb)+1
              nperp(iarb)=0
           endif

           nperp(iarb)=nperp(iarb)+1

        enddo
        close(17)

        write(*,*) iarb, npara(iarb),nperp(iarb)

        if((iarb==1).or.(npara(iarb)>npara_max)) then 
           npara_max=npara(iarb)
        endif

        if((iarb==1).or.(nperp(iarb)>nperp_max))then 
           nperp_max=nperp(iarb)
        endif

     enddo


     !read distribution data from the files
     

     allocate(distribution(npara_max,nperp_max,narb))
     allocate(vpara(npara_max,narb),vperp(nperp_max,narb))

     do iarb=1,narb

        write(filename,'(A25,I1,A4)') 'distribution/distribution', iarb,'.dat'

        open(unit=17,status='old',file=filename)

        do ipara=1,npara(iarb)

           do iperp=1,nperp(iarb)

              read(17,*) vpa, vpe, dist_value
              distribution(ipara,iperp,iarb)=dist_value
              if(ipara==1) vperp(iperp,iarb)=vpe

           enddo

           vpara(ipara,iarb)=vpa

        enddo

        close(17)
     enddo

  endif


end subroutine read_distr
