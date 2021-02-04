!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File lsmrDataModule.f90
!
! Defines real(dp) and a few constants for use in other modules.
!
! 24 Oct 2007: Allows floating-point precision dp to be defined
!              in exactly one place (here).  Note that we need
!                 use lsmrDataModule
!              at the beginning of modules AND inside interfaces.
!              zero and one are not currently used by LSMR,
!              but this shows how they should be declared
!              by a user routine that does need them.
! 16 Jul 2010: LSMR version derived from LSQR equivalent.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module lsmrDataModule

  implicit none

  intrinsic                   ::      selected_real_kind, selected_int_kind
  integer,  parameter, public :: dp = selected_real_kind(15)
  integer,  parameter, public :: sp    = selected_real_kind(6)
  integer, parameter, public  :: ip = selected_int_kind(9)       ! R: (-10^R, 10^R)
  real(sp), parameter, public :: zero = 0.0_sp, one = 1.0_sp, eps=epsilon(zero)


end module lsmrDataModule
