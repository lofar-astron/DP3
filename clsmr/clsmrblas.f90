!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     File zlsmrBlasModule.f90
!
!     This file contains the following BLAS subroutines
!        scnrm2     sscal
!     required by cLSMR.
!
! 29 Jun 2013: Created complex BLAS module for 
!              complex LSMR.
! 20 Jan 2021: Single precision version
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*****************************************************************************
function scnrm2 ( n, x, incx )
!*****************************************************************************
!
! SCNRM2 returns the euclidean norm of a complex(4) vector.
!
!
!  Discussion:
!
!    SCNRM2 := sqrt ( sum ( conjg ( x(1:n) ) * x(1:n) ) )
!            = sqrt ( dot_product ( x(1:n), x(1:n) ) )
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, complex(4) X(*), the vector.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, real(4) DZNRM2, the norm of the vector.
!
  implicit none
!
  integer, intent(in)     :: incx
  integer                 :: ix
  integer, intent(in)     :: n
  real(4)                 :: norm
  real(4)                 :: scale
  real(4)                 :: scnrm2
  real(4)                 :: ssq
  real(4)                 :: temp
  real(4), parameter      :: one = 1.0
  real(4), parameter      :: zero = 0.0
  complex(4), intent(in)  :: x(*)

!
  if ( n < 1 .or. incx < 1 ) then

    norm  = zero

  else

    scale = zero
    ssq = one

    do ix = 1, 1 + ( n - 1 ) * incx, incx
      if ( real(x(ix), 4) /= zero ) then
        temp = abs ( real(x(ix), 4) )
        if ( scale < temp ) then
          ssq = one + ssq * ( scale / temp )**2
          scale = temp
        else
          ssq = ssq + ( temp / scale )**2
        end if
      end if

      if ( aimag ( x(ix) ) /= zero ) then
        temp = abs ( aimag ( x(ix) ) )
        if ( scale < temp ) then
          ssq = one + ssq * ( scale / temp )**2
          scale = temp
        else
          ssq = ssq + ( temp / scale )**2
        end if

      end if

    end do

    norm  = scale * sqrt ( ssq )

  end if

  scnrm2 = norm

  return
end function scnrm2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*******************************************************************
subroutine sscal(n, za, zx, incx)
!*****************************************************************************
!
! sscal returns the euclidean norm of a complex(4) vector.
!
!
!  Discussion:
!
!    SSCAL scales a vector by a constant.
!
!  Parameters:
!
!    Input, integer n, the number of entries in the vector.

!    Input, complex(4) za(*), the scaling parameter.

!    Input/Output, complex(4) zx(*), the vector.
!
!    Input, integer incx, the increment between successive entries of X.

  implicit none

  integer, intent(in)        :: incx
  integer                    :: i, nincx
  integer, intent(in)        :: n
  complex(4), intent(in)    :: za
  complex(4), intent(inout) :: zx(*)

  if (n .le. 0 .or. incx .le. 0) return
  if (incx .eq. 1) then
    do i = 1, n
      zx(i) = za * zx(i)
    end do
  else
    nincx = n * incx
    do i = 1, nincx, incx
      zx(i) = za * zx(i)
    end do
  end if
  return

end subroutine sscal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

