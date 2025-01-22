!> Determines how many terms have to be included for accurate evaluation of the generalized hypergeometric functions using a continued fraction formula, and whether arbitrary precision arithmetic has to be employed or not
!! \param z argument of the generalized hypergeometric function
!! \param ndp measure for the number of significant digits to be included for cases where arbitrary precision arithmetic is required (i.e. cases with ndp >0)
!! \param nfrac required number of terms to be included in continued fraction formula for an accurate evaluation of the generalized hypergeometric functions
subroutine acc_F(z,ndp,nfrac)
  implicit none
  integer :: ndp,nfrac
  real(kind=16) :: z

  if(abs(z)<0.04_16) then
     ndp=0
     nfrac=1
  else if((abs(z).ge.0.04_16).and.(abs(z)<0.1_16)) then
     ndp=0
     nfrac=2
  else if((abs(z).ge.0.1_16).and.(abs(z)<1.0_16)) then
     ndp=0
     nfrac=5
  else if((abs(z).ge.1.0_16).and.(abs(z)<5.0_16)) then
     ndp=0
     nfrac=9
  else if((abs(z).ge.5.0_16).and.(abs(z)<10.0_16)) then
     ndp=0
     nfrac=12
  else if((abs(z).ge.10.0_16).and.(abs(z)<30.0_16)) then
     ndp=0
     nfrac=19
  else if((abs(z).ge.30.0_16).and.(abs(z)<50.0_16)) then
     ndp=0
     nfrac=23
  else if((abs(z).ge.50.0_16).and.(abs(z)<100.0_16)) then
     ndp=0
     nfrac=31
  else if((abs(z).ge.100.0_16).and.(abs(z)<200.0_16)) then
     ndp=0
     nfrac=43
  else if((abs(z).ge.200.0_16).and.(abs(z)<400.0_16)) then
     ndp=0
     nfrac=57
  else if((abs(z).ge.400.0_16).and.(abs(z)<600.0_16)) then
     ndp=0
     nfrac=70
  else if((abs(z).ge.600.0_16).and.(abs(z)<800.0_16)) then
     ndp=30
     nfrac=82
  else if((abs(z).ge.800.0_16).and.(abs(z)<1000.0_16)) then
     ndp=30
     nfrac=91
  else if((abs(z).ge.1000.0_16).and.(abs(z)<1500.0_16)) then
     ndp=30
     nfrac=110
  else if((abs(z).ge.1500.0_16).and.(abs(z)<2000.0_16)) then
     ndp=44
     nfrac=126
  else if((abs(z).ge.2000.0_16).and.(abs(z)<3000.0_16)) then
     ndp=44
     nfrac=153
  else if((abs(z).ge.3000.0_16).and.(abs(z)<4000.0_16)) then
     ndp=58
     nfrac=176
  else if((abs(z).ge.4000.0_16).and.(abs(z)<6000.0_16)) then
     ndp=73
     nfrac=215
  else if((abs(z).ge.6000.0_16).and.(abs(z)<8000.0_16)) then
     ndp=73
     nfrac=247
  else if((abs(z).ge.8000.0_16).and.(abs(z)<10000.0_16)) then
     ndp=87
     nfrac=276
  else if((abs(z).ge.10000.0_16).and.(abs(z)<14000.0_16)) then
     ndp=102
     nfrac=326
  else if((abs(z).ge.14000.0_16).and.(abs(z)<20000.0_16)) then
     ndp=131
     nfrac=388
  else if((abs(z).ge.20000.0_16).and.(abs(z)<30000.0_16)) then
     ndp=159
     nfrac=474
  else if((abs(z).ge.30000.0_16).and.(abs(z)<40000.0_16)) then
     ndp=174
     nfrac=547
  else if((abs(z).ge.40000.0_16).and.(abs(z)<50000.0_16)) then
     ndp=203
     nfrac=611
  else if((abs(z).ge.50000.0_16).and.(abs(z)<70000.0_16)) then
     ndp=232
     nfrac=723
  else if((abs(z).ge.70000.0_16).and.(abs(z)<100000.0_16)) then
     ndp=275
     nfrac=863
  else if((abs(z).ge.100000.0_16).and.(abs(z)<150000.0_16)) then
     ndp=333
     nfrac=1056
  else if((abs(z).ge.150000.0_16).and.(abs(z)<200000.0_16)) then
     ndp=391
     nfrac=1219
  else if((abs(z).ge.200000.0_16).and.(abs(z)<300000.0_16)) then
     ndp=477
     nfrac=1492
  else if((abs(z).ge.300000.0_16).and.(abs(z)<400000.0_16)) then
     ndp=550
     nfrac=1722
  else if((abs(z).ge.400000.0_16).and.(abs(z)<500000.0_16)) then
     ndp=622
     nfrac=1925
  else if((abs(z).ge.500000.0_16).and.(abs(z)<700000.0_16)) then
     ndp=723
     nfrac=2277
  else
     write(*,*) 'Bessel argument very large!'
     write(*,*) 'Increase precision manually to ndp>477 and nfrac>1492'
     stop
  endif

end subroutine acc_F
