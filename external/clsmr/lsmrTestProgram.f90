!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File lsmrTestProgram.f90
!
!    lsmrTestProgram
!
! Main program for testing LSMR via subroutine lsmrtest in lsmrTestModule.
!
! Maintained by Michael Saunders <saunders@stanford.edu>.
!
! 16 Jul 2010: LSMR version derived from LSQR equivalent.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program lsmrTestProgram

  use   lsmrDataModule, only : dp, ip
  use   lsmrTestModule, only : lsmrtest
  implicit none

  !---------------------------------------------------------------------
  ! This program calls lsmrtest(...) to generate a series of test problems
  ! Ax = b or Ax ~= b and solve them with LSMR.
  ! The matrix A is m x n.  It is defined by routines in lsmrTestModule.
  !
  ! 23 Sep 2007: First version of lsmrTestProgram.f90.
  ! 24 Oct 2007: Use real(dp) instead of compiler option -r8.
  !---------------------------------------------------------------------

  ! Local variables
  integer(ip)  :: m,n,nbar,ndamp,nduplc,npower,nout
  integer(ip)  :: localSize  ! No. of vectors involved in local reorthogonalization (>= 0)
  real(dp)     :: damp


  nout   = 6
  open(nout,file='LSMR.txt',status='unknown')

  nbar   = 100  ! 1000
  nduplc =   4  ! 40


  m = 2*nbar        ! Over-determined systems
  n = nbar
  localSize = 0

  do ndamp = 2,6
     npower = ndamp
     damp   = 0.0
     if (ndamp > 2) damp   = 10.0**(-ndamp)
     call lsmrtest(m,n,nduplc,npower,damp,localSize,nout)
  end do

  localsize = 10    ! Repeat last test with local reorthogonalization
  call lsmrtest(m,n,nduplc,npower,damp,localSize,nout)


  m = nbar          ! Square systems
  n = nbar
  localSize = 0

  do ndamp = 2,6
     npower = ndamp
     damp   = 0.0
     if (ndamp > 2) damp   = 10.0**(-ndamp-6)
     call lsmrtest(m,n,nduplc,npower,damp,localSize,nout)
  end do

  localsize = 10    ! Repeat last test with local reorthogonalization
  call lsmrtest(m,n,nduplc,npower,damp,localSize,nout)


  m = nbar          ! Under-determined systems
  n = 2*nbar
  localSize = 0

  do ndamp = 2,6
     npower = ndamp
     damp   = 0.0
     if (ndamp > 2) damp   = 10.0**(-ndamp-6)
     call lsmrtest(m,n,nduplc,npower,damp,localSize,nout)
  end do

  localsize = 10    ! Repeat last test with local reorthogonalization
  call lsmrtest(m,n,nduplc,npower,damp,localSize,nout)

  close(nout)

end program lsmrTestProgram
