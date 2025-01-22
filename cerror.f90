!> Computes the error function for a complex argument 
!! \param z (complex) argument of the error function
!! \param cer value of the error function for the given argument
subroutine cerror( z_in, cer_sol )
  use param_mod
  !*****************************************************************************80
  !
  !! CERROR computes the error function for a complex argument.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    15 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  implicit none
  real(kind=16) :: a0
  complex(kind=16) :: c0
  complex(kind=16) :: cer
  complex :: cer_sol
  complex(kind=16) :: cl
  complex(kind=16) :: cr
  complex(kind=16) :: cs
  integer :: k2
  complex :: z_in
  complex(kind=16) :: z
  complex(kind=16) :: z1

  z=z_in

  a0 = abs ( z )
  c0 = exp ( - z * z )

  z1 = z

  if ( real (z) < 0.0_16 ) then
     z1 = - z
  end if

  if ( a0 <= 5.8_16 ) then    

     cs = z1
     cr = z1
     do k2 = 1, 120
        cr = cr * z1 * z1 / ( k2 + 0.5_16 )
        cs = cs + cr
        if ( abs ( cr / cs ) < 10.0_16**(-15) ) then
           exit
        end if
     end do

     cer = 2.0_16 * c0 * cs / sqrt ( 1.0_16*pi )

  else

     cl = 1.0_16 / z1              
     cr = cl
     do k2 = 1, 13
        cr = -cr * ( k2 - 0.5_16 ) / ( z1 * z1 )
        cl = cl + cr
        if ( abs ( cr / cl ) < 10.0_16**(-15) ) then
           exit
        end if
     end do

     cer = 1.0_16 - c0 * cl / sqrt ( 1.0_16*pi )

  end if

  if ( real ( z) < 0.0_16 ) then
     cer = -cer
  end if

  cer_sol=cer

end subroutine cerror
