!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File lsmrblasInterface.f90
!
!    BLAS1 Interfaces:  dznrm2    zscal
!
! Maintained by Michael Saunders <saunders@stanford.edu>.
!
! 29 Jun 2013: LSMR version derived from LSQR equivalent.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module zlsmrblasInterface

  implicit none
  public   :: dznrm2, zscal

  interface                              ! Level 1 BLAS
     function dznrm2 (n,x,incx)
       use zlsmrDataModule, only : dp, ip
       integer(ip),  intent(in)    :: n,incx
       complex(dp), intent(in)     :: x(*)
       real(dp)                    :: dznrm2
     end function dznrm2

     subroutine zscal (n,za,zx,incx)
       use zlsmrDataModule, only : dp, ip
       integer(ip),  intent(in)       :: n,incx
       complex(dp), intent(in)        :: za
       complex(dp), intent(inout)     :: zx(*)
     end subroutine zscal
  end interface

end module zlsmrblasInterface
