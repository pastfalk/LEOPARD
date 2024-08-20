!> Reads the velocity distribution data from the provided distribution files

! Modified in July/2022 by Rudi Gaelzer
!  The grid file containing the distribution must now have the name: distribution***-n.dat,
!  where * = 0 - 9 and n = 1 - 9.

subroutine read_distr_rg(set_num, log_unit)
use param_mod
implicit none
character(len= 3), intent(in) :: set_num
integer, intent(in) :: log_unit
character(len= 1) :: ciarb
integer :: ios
real :: start_pe
real :: vpa, vpe
real :: dist_value
integer :: ipara, iperp, iarb
integer :: j, distu
character(len= 34) :: filename

!count number of included species with arbitrary velocity distribution

narb=0
do j=1,Nspecies
   if(mode(j).eq.1) then
      narb=narb+1
   endif
enddo

if(narb.ne.0) then

   allocate(npara(narb),nperp(narb))

   !determine dimension of the velocity grid for each velocity distribution

   do iarb=1,narb

      write(ciarb, fmt= '(i1)')iarb
      filename= 'distribution/distribution_'//set_num//'-'//ciarb//'.dat'
!       write(filename,'(A25,I1,A4)') 'distribution/distribution', iarb,'.dat'
      open(newunit= distu, status='old',file= filename, iostat= ios)
      if(ios /= 0)then
         write(*, '(2a)') 'Error during access to distribution file: ', filename
         stop
      end if

      npara(iarb)=0
      nperp(iarb)=0

      read(distu,*,iostat=ios) vpa, start_pe, dist_value
      if(ios /= 0)then
         write(*, '(2a)') 'Error during access to data in the distribution file: ', filename
         stop
      end if
      rewind(distu)

      do while(.true.)

         read(distu,*,iostat=ios) vpa, vpe, dist_value

         if(ios > 0)then
            write(*, '(2a)') 'Error during access to data in the distribution file: ', filename
            stop
         end if
!          if (ios.ne.0) exit
         if (ios < 0) exit

         if (vpe.eq.start_pe) then
            npara(iarb)=npara(iarb)+1
            nperp(iarb)=0
         endif

         nperp(iarb)=nperp(iarb)+1

      enddo
      close(distu)

      write(log_unit, *) iarb, npara(iarb),nperp(iarb)

      if((iarb.eq.1).or.(npara(iarb).gt.npara_max)) then 
         npara_max=npara(iarb)
      endif

      if((iarb.eq.1).or.(nperp(iarb).gt.nperp_max))then 
         nperp_max=nperp(iarb)
      endif

   enddo


   !read distribution data from the files
   

   allocate(distribution(npara_max,nperp_max,narb))
   allocate(vpara(npara_max,narb),vperp(nperp_max,narb))

   do iarb=1,narb

      write(ciarb, fmt= '(i1)')iarb
      filename= 'distribution/distribution'//set_num//'-'//ciarb//'.dat'
      open(newunit= distu, status='old',file= filename)
   
!       write(filename,'(A25,I1,A4)') 'distribution/distribution', iarb,'.dat'
!       open(unit=distu,status='old',file=filename)
      
      do ipara=1,npara(iarb)

         do iperp=1,nperp(iarb)

            read(distu,*) vpa, vpe, dist_value
            distribution(ipara,iperp,iarb)=dist_value
            if(ipara.eq.1) vperp(iperp,iarb)=vpe

         enddo

         vpara(ipara,iarb)=vpa

      enddo

      close(distu)
   enddo

endif

end subroutine read_distr_rg
