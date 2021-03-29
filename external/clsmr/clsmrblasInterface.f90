!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File lsmrblasInterface.f90
!
!    BLAS1 Interfaces:  scnrm2    sscal
!
! Maintained by Michael Saunders <saunders@stanford.edu>.
!
! 29 Jun 2013: LSMR version derived from LSQR equivalent.
! 20 Jan 2021: single precision version
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module clsmrblasInterface

  implicit none
  public   :: scnrm2, sscal

  interface                              ! Level 1 BLAS
     function scnrm2 (n,x,incx)
       use clsmrDataModule, only : sp, ip
       integer(ip),  intent(in)    :: n,incx
       complex(sp), intent(in)     :: x(*)
       real(sp)                    :: scnrm2
     end function scnrm2

     subroutine sscal (n,sa,sx,incx)
       use clsmrDataModule, only : sp, ip
       integer(ip),  intent(in)       :: n,incx
       complex(sp), intent(in)        :: sa
       complex(sp), intent(inout)     :: sx(*)
     end subroutine sscal
  end interface

end module clsmrblasInterface
