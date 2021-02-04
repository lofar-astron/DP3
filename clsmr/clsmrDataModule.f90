!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File clsmrDataModule.f90
!
! Extends lsmrDataModule.f90 for use with complex numbers
! 29 Jun 2013: File created
! 20 Jan 2021: Single precision version
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module clsmrDataModule

  use    lsmrDataModule, only  :  dp, sp, ip, zero, one
  implicit none

  intrinsic                      ::      cmplx
  complex(sp), parameter, public :: zzero = cmplx(zero,zero,sp), zone = cmplx(one,zero,sp)

end module clsmrDataModule
