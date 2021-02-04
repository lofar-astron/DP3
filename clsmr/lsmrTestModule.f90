!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File lsmrTestModule.f90
!
!    Hprod   Aprod1   Aprod2    lstp     lsmrtest
!
! These routines define a class of least-squares test problems
! for testing algorithm LSMR (Paige and Saunders, ACM TOMS, 1982).
! Aprod1 and Aprod2 define matrix-vector products required by
! subroutine LSMR for a test matrix of the form  A = Y*D*Z,
! where Y and Z are Householder transformations and D is diagonal.
!
! This file illustrates how LSMR can call Aprod1 and Aprod2 with a
! short fixed parameter list, even if they need arbitrary other data.
! Here, they need arrays d, hx, hy created by subroutine lsmrtest
! in this module before LSMR is called (see below).
!
! Maintained by Michael Saunders <saunders@stanford.edu>.
!
! 21 Sep 2007: lsqrTestModule.f90 implemented.
! 16 Jul 2010: LSMR version derived from LSQR equivalent.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module lsmrTestModule

  use  lsmrDataModule,    only : dp, sp, ip, zero, one
  use  lsmrblasInterface, only : dnrm2, dscal
  use  lsmrModule,        only : LSMR
  use  lsmrCheckModule,   only : Acheck, xcheck

  implicit none
  public   :: lsmrtest
  private  :: Hprod, Aprod1, Aprod2, lstp

  ! DYNAMIC WORKSPACE DEFINED HERE.
  ! They are allocated in lsmrtest and used by Aprod1, Aprod2.

  real(dp), allocatable :: d(:), hy(:), hz(:) ! These define A = Y D Z.
  real(dp), allocatable :: wm(:), wn(:)       ! Work vectors.

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Hprod (n,z,x)

    integer(ip),  intent(in)    :: n
    real(dp), intent(in)    :: z(n)
    real(dp), intent(inout) :: x(n)

    !-------------------------------------------------------------------
    ! Hprod  applies a Householder transformation stored in z
    ! to return x = (I - 2*z*z')*x.
    !
    ! 23 Sep 2007: Fortran 90 version.
    !-------------------------------------------------------------------

    integer(ip)             :: i
    real(dp)                :: s

    s = zero
    do i = 1,n
       s = z(i)*x(i) + s
    end do

    s = s + s
    x = x - s*z
    
  end subroutine Hprod

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Aprod1(m,n,x,y)

    integer(ip),  intent(in)    :: m,n
    real(dp), intent(in)        :: x(n)
    real(dp), intent(inout)     :: y(m)

    !-------------------------------------------------------------------
    ! Aprod1 computes y = y + A*x without altering x,
    ! where A is a test matrix of the form  A = Y*D*Z,
    ! and the matrices D, Y, Z are represented by
    ! the allocatable vectors d, hy, hz in this module.
    !
    ! 23 Sep 2007: Fortran 90 version.
    !-------------------------------------------------------------------

    intrinsic           :: min
    integer(ip)         :: minmn

    minmn = min(m,n)
    wn    = x
    call Hprod (n,hz,wn)
    wm(1:minmn) = d(1:minmn)*wn(1:minmn)
    wm(n+1:m)   = zero
    call Hprod (m,hy,wm)
    y = y + wm

  end subroutine Aprod1

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Aprod2(m,n,x,y)

    integer(ip),  intent(in)    :: m,n
    real(dp), intent(inout)     :: x(n)
    real(dp), intent(in)        :: y(m)

    !-------------------------------------------------------------------
    ! Aprod2 computes x = x + A'*y without altering y,
    ! where A is a test matrix of the form  A = Y*D*Z,
    ! and the matrices D, Y, Z are represented by
    ! the allocatable vectors d, hy, hz in this module.
    !
    ! 23 Sep 2007: Fortran 90 version.
    !-------------------------------------------------------------------

    intrinsic               :: min
    integer(ip)             :: minmn

    minmn = min(m,n)
    wm    = y
    call Hprod (m,hy,wm)
    wn(1:minmn) = d(1:minmn)*wm(1:minmn)
    wn(m+1:n)   = zero
    call Hprod (n,hz,wn)
    x = x + wn

  end subroutine Aprod2

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine lstp  (m,n,nduplc,npower,damp,x,b,condA,normr)

    integer(ip),  intent(in)    :: m, n, nduplc, npower
    real(dp), intent(in)        :: damp
    real(dp), intent(inout)     :: x(n)
    real(dp), intent(out)       :: b(m)
    real(dp), intent(out)       :: condA, normr

    !-------------------------------------------------------------------
    ! lstp  generates a sparse least-squares test problem of the form
    !           (   A    )*x = ( b ) 
    !           ( damp*I )     ( 0 )
    ! for solution by LSMR, or a sparse underdetermined system
    !            Ax + damp*s = b
    ! for solution by LSMR or CRAIG.  The matrix A is m by n and is
    ! constructed in the form  A = Y*D*Z,  where D is an m by n
    ! diagonal matrix, and Y and Z are Householder transformations.
    !
    ! m and n may contain any positive values.
    ! If m >= n  or  damp = 0, the true solution is x as given.
    ! Otherwise, x is modified to contain the true solution.
    !
    ! 1982---1991: Various versions implemented.
    ! 06 Feb 1992: lstp generalized to allow any m and n.
    ! 07 Sep 2007: Line by line translation for Fortran 90 compilers
    !              by Eric Badel <badel@nancy.inra.fr>.
    ! 23 Sep 2007: Fortran 90 version with modules.
    !-------------------------------------------------------------------

    intrinsic           :: min, cos, sin, sqrt
    integer(ip)         :: i, j, minmn
    real(dp)            :: alfa, beta, dampsq, fourpi, t

    !-------------------------------------------------------------------
    ! Make two vectors of norm 1.0 for the Householder transformations.
    ! fourpi  need not be exact.
    !-------------------------------------------------------------------
    minmn  = min(m,n)
    dampsq = damp**2
    fourpi = 4.0 * 3.141592
    alfa   = fourpi / m
    beta   = fourpi / n

    do i = 1,m
       hy(i) = sin( alfa*i )
    end do

    do i = 1,n
       hz(i) = cos( beta*i )
    end do

    alfa   = dnrm2 ( m, hy, 1 )
    beta   = dnrm2 ( n, hz, 1 )
    call dscal ( m, (- one/alfa), hy, 1 )
    call dscal ( n, (- one/beta), hz, 1 )

    !-------------------------------------------------------------------
    ! Set the diagonal matrix D.  These are the singular values of A.
    !-------------------------------------------------------------------
    do i = 1,minmn
       j    = (i-1+nduplc) / nduplc
       t    =  j*nduplc
       t    =  t / minmn
       d(i) =  t**npower
    end do

    condA  = (d(minmn)**2 + dampsq) / (d(1)**2 + dampsq)
    condA  = sqrt( condA )

    !-------------------------------------------------------------------
    ! If m >=n, the input x will be the true solution.
    ! If m < n, reset x = different true solution.
    ! It must be of the form x = Z(w1) for some w1.
    !                             (0 )
    ! We get w1 from the top of w = Zx.
    !-------------------------------------------------------------------
    wn = x
    call Hprod (n,hz,wn)
    if (m < n) then
       wn(m+1:n) = zero
       call Hprod (n,hz,wn)
       x = wn
    end if

    ! Let r = Y rbar and xbar = Z x.
    ! Solve D r1bar = damp^2 x1bar, where r1bar is in wm(1:minmn).

    wm(1:minmn)   = dampsq * wn(1:minmn)/d(1:minmn)

    wm(minmn+1:m) = one      ! Set r2bar to be anything (empty if m <= n).
    call Hprod (m,hy,wm)     ! Form r = Y rbar  in wm.

    normr  = dnrm2 (m,wm,1)  ! Compute norm(r)
    b      = wm              ! and  b = r + Ax.
    call Aprod1(m,n,x,b)

  end subroutine lstp

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine lsmrtest(m,n,nduplc,npower,damp,localSize,nout)

    integer(ip),  intent(in) :: m, n, nduplc, npower, &
                                localSize,            & ! Local reorthogonalization
                                nout
    real(dp), intent(in)     :: damp

    !-------------------------------------------------------------------
    ! This is an example driver routine for running LSMR.
    ! It generates a test problem, solves it, and examines the results.
    !
    ! 27 Sep 2007: f90 version of lsqrtest.
    ! 17 Jul 2010: LSMR version derived from LSQR equivalent.
    !              localSize specifies that each bidiagonalization vector v
    !              should be reorthogonalized wrto the last "localSize" v's.
    !              localSize >= 0.
    !------------------------------------------------------------------------

    intrinsic       :: epsilon

    ! Local arrays and variables
    real(dp)        :: b(m), se(n), x(n), xtrue(n)
    integer(ip)     :: inform, istop, itn, itnlim, j, minmn, nprint
    real(dp)        :: atol, btol, conlim, normA, condA,             &
                       eps, normr, normAr, norme, etol, normw, normx

    ! Local constants
    character(len=*),  &
              parameter :: line = '----------------------------------'


    eps    = epsilon(eps)
    if (eps > 1e-9) then
       write(nout,*) ' '
       write(nout,*) 'WARNING: '
       write(nout,*) 'MACHINE PRECISION EPS =', eps, '  SEEMS TO BE INADEQUATE'
       write(nout,*) ' '
    end if

    !------------------------------------------------------------------------
    ! Generate the specified test problem
    ! and check that Aprod1, Aprod2 form y + Ax and x + A'y consistently.
    !------------------------------------------------------------------------
    do j = 1,n                                   ! Set the desired solution xtrue.
       xtrue(j) = 0.1*j                          ! For least-squares problems, this is it.
    end do                                       ! If m<n, lstp will alter it.

    minmn  = min(m,n)                            ! Allocate arrays.
    allocate( d(minmn), hy(m), hz(n) )           ! Vectors defining A = Y D Z.
    allocate( wm(m), wn(n) )                     ! Work vectors for Aprod1, Aprod2.

    call lstp  (m,n,nduplc,npower,damp, &        ! Generate test problem.
                xtrue,b,condA,normr)             ! If m<n, xtrue is altered.

    write(nout,1000) line,line,m,n,nduplc,npower,damp,condA,normr, &
                     line,line

    call Acheck(m,n,Aprod1,Aprod2,nout,inform)   ! Check Aprod1, Aprod2.
          
    if (inform > 0) then
       write(nout,'(a)') 'Check tol in subroutine Acheck in lsmrCheckModule'
       stop
    end if

    !------------------------------------------------------------------------
    ! Set input parameters for LSMR
    ! and solve the problem defined by Aprod1, Aprod2, b, damp.
    ! Next line would ask for standard errors only if they are well-defined.
    !------------------------------------------------------------------------
    atol   = eps**0.99                           ! Asks for high accuracy.
    btol   = atol
    conlim = 1000.0 * condA
    itnlim = 4*(m + n + 50)

    call LSMR  ( m, n, Aprod1, Aprod2, b, damp,               &
                 atol, btol, conlim, itnlim, localSize, nout, &
                 x, istop, itn, normA, condA, normr, normAr, normx )

    call xcheck( m, n, Aprod1, Aprod2, b, damp, x, &    ! Check x
                 normA, nout,                      &
                 inform )

    nprint = min(m,n,8)
    write(nout,2500)    (j, x(j), j=1,nprint)    ! Print some of the solution

    wn     = x - xtrue                           ! Print a clue about whether
    normw  = dnrm2 (n,wn,1)                      ! the solution looks OK.
    normx  = dnrm2 (n,xtrue,1)
    norme  = normw/(one + normx)
    etol   = 1.0e-3
    if (norme <= etol) then
       write(nout,3000) norme
    else
       write(nout,3100) norme
    end if

    deallocate( wn, wm, hz, hy, d )              ! Free arrays (in reverse)
    return

 1000 format(1p &
      // 1x, 2a &
      /  ' Least-Squares Test Problem      P(', 4i5, e12.2, ' )'     &
      /  ' Condition no. =', e12.4, '     Residual function =', e17.9&
      /  1x, 2a)
 2500 format(/' Solution  x:'          / 4(i6, g14.6))
 3000 format(1p &
      /  ' LSMR  appears to be successful.'                          &
      /  ' Relative error in  x  =', e10.2)
 3100 format(1p &
      /  ' LSMR  appears to have failed.  '                          &
      /  ' Relative error in  x  =', e10.2)

  end subroutine lsmrtest

end module lsmrTestModule
