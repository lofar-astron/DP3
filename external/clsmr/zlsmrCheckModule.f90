!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File zlsmrCheckModule.f90
!
!    zAcheck zxcheck
!
! zAcheck tests if a user's matrix-vector product routines for
! computing y + A*x and x + A'*y are working with the same A.
! zxcheck tests if a given x seems to be a solution of Ax = b or Ax ~= b.
!
! Extends lsmrCheckModule to complex numbers
!
! 29 Jun 2013: Derived from zlsqrCheckModule
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module zlsmrCheckModule

  use  zlsmrDataModule,       only : dp, ip, zzero, zero, one
  use  zlsmrblasInterface,    only : dznrm2, zscal

  implicit none
  private
  public      :: zAcheck, zxcheck
  intrinsic   :: abs, max, sqrt, dot_product, epsilon

contains

 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine zAcheck( m, n, Aprod1, Aprod2, nout, inform )

    integer(ip), intent(in)    :: m, n   ! No. of rows and cols of A
    integer(ip), intent(in)    :: nout   ! Output file number
    integer(ip), intent(out)   :: inform ! = 0 if Aprod1, Aprod2 seem ok
                                         ! = 1 otherwise
    interface
       subroutine Aprod1(m,n,x,y)                   ! y := y + A*x
         use zlsmrDataModule, only : dp, ip
         integer(ip), intent(in)       :: m,n
         complex(dp),    intent(in)    :: x(n)
         complex(dp),    intent(inout) :: y(m)
       end subroutine Aprod1

       subroutine Aprod2(m,n,x,y)                   ! x := x + A'*y
         use zlsmrDataModule, only : dp, ip
         integer(ip), intent(in)       :: m,n
         complex(dp),    intent(inout) :: x(n)
         complex(dp),    intent(in)    :: y(m)
       end subroutine Aprod2
    end interface

    !-------------------------------------------------------------------
    ! One-liner: Acheck checks Aprod1 and Aprod2 for LSMR.
    !
    ! Purpose:   Acheck tests the user subroutines Aprod1 and Aprod2
    !   called by LSMR.  For some m x n matrix A,
    !   Aprod1 computes y := y + A*x  from given x,y without altering x,
    !   Aprod2 computes x := x + A'*y from given x,y without altering y.
    !   Acheck tries to verify that A and A' refer to the same matrix.
    !
    ! Method:    We cook up some unlikely vectors x and y of unit length
    !   and test if  y'(y + Ax)  =  x'(x + A'y).
    !
    ! Parameter Constants:
    !   Param   Type   Description
    !   power   real   eps**power is the tolerance for judging if
    !                  y'(y + Ax) = x'(x + A'y) to sufficient accuracy.
    !                  power should be in the range (0.25, 0.9) say.
    !                  For example, power = 0.75 means we are happy
    !                  if 3/4 of the available digits agree.
    !                  power = 0.5 seems a reasonable requirement
    !                  (asking for half the digits to agree).
    !                    
    ! History:
    ! 04 Sep 1991  Initial design and code.
    !              Michael Saunders, Dept of Operations Research,
    !              Stanford University.
    ! 10 Feb 1992  Aprod added as parameter.
    !              tol defined via power.
    ! 10 Feb 1992: Acheck revised and xcheck implemented.
    ! 27 May 1993: Acheck and xcheck kept separate from test problems.
    ! 23 Sep 2007: Acheck implemented as part of this f90 module.
    ! 29 Jun 2013: zAcheck extends Acheck to complex numbers
    !-------------------------------------------------------------------

    ! Local arrays and variables
    complex(dp)                :: x(n), v(n), w(m), y(m)
    integer(ip)                :: i, j
    real(dp)                   :: alfa, beta, t, test1, test2, test3, tol, eps

    ! Local constants
    real(dp), parameter :: power = 0.5_dp
    eps = epsilon(eps)
    tol    = eps**power
    if (nout > 0) write(nout,1000)

    !===================================================================
    ! Cook up some unlikely vectors x and y of unit length.
    !===================================================================
    t = one
    do j=1,n
       t    = t + one
       x(j) = cmplx(sqrt(t),sqrt(t),dp)
    end do
 
    t = one
    do i=1,m
       t    = t + one
       y(i) = cmplx(one/sqrt(t),one/sqrt(t),dp)
    end do
 
    alfa = dznrm2 (n,x,1)
    beta = dznrm2 (m,y,1)
    call zscal (n, cmplx((one/alfa),zero,dp), x, 1)
    call zscal (m, cmplx((one/beta),zero,dp), y, 1)
      
    !===================================================================
    ! Test if y'(y + Ax) = x'(x + A'y).
    !===================================================================
    w(1:m) = y(1:m)               ! First set  w = y + Ax,  v = x + A'y.
    v(1:n) = x(1:n)
    call Aprod1(m,n,x,w)
    call Aprod2(m,n,v,y)
      
    alfa   = dot_product(y,w)    ! Now set    alfa = y'w,  beta = x'v.
    beta   = dot_product(x,v)
    test1  = abs(alfa - beta)            
    test2  = one + abs(alfa) + abs(beta)
    test3  = test1 / test2

    if (test3 <= tol) then        ! See if alfa and beta are essentially
       inform = 0                 ! the same.
       if (nout > 0) write(nout,1010) test3
    else
       inform = 1
       if (nout > 0) write(nout,1020) test3
    end if

    return

 1000 format(//' Enter Acheck.')
 1010 format(1p, &
        ' Aprod1, Aprod2 seem OK.  Relative error =', e10.1)
 1020 format(1p, &
        ' Aprod1, Aprod2 seem incorrect.  Relative error =', e10.1)
	          
  end subroutine zAcheck

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine zxcheck( m, n, Aprod1, Aprod2, b, damp, x, &
                     Anorm, tol, nout,                 &
                     inform )

    integer(ip),  intent(in)       :: m, n     ! No. of rows and cols of A
    integer(ip),  intent(in)       :: nout     ! Output file number
    integer(ip),  intent(out)      :: inform   ! = 0 if b = 0 and x = 0.
                                               ! = 1 2 or 3 if x seems to
                                               ! solve systems 1 2 or 3 below.
    real(dp), intent(in)           :: Anorm    ! An estimate of norm(A) or
                                               ! norm( A, delta*I ) if delta > 0.
                                               ! Provided by LSQR.
    real(dp), intent(in)           :: tol      ! tolerance for judging residuals.
                                               ! Typically the tol used for computing x.
    real(dp), intent(in)           :: damp     ! Defines problem 3 below.
    complex(dp), intent(in)        :: b(m)     ! The right-hand side of Ax ~= b.
    complex(dp), intent(in)        :: x(n)     ! The given solution estimate.

    interface
       subroutine Aprod1(m,n,x,y)                   ! y := y + A*x
         use zlsmrDataModule, only : dp, ip
         integer(ip),  intent(in)       :: m,n
         complex(dp), intent(in)        :: x(n)
         complex(dp), intent(inout)     :: y(m)
       end subroutine Aprod1

       subroutine Aprod2(m,n,x,y)                   ! x := x + A'*y
         use zlsmrDataModule, only : dp, ip
         integer(ip),  intent(in)       :: m,n
         complex(dp), intent(inout)     :: x(n)
         complex(dp), intent(in)        :: y(m)
       end subroutine Aprod2
    end interface

    !-------------------------------------------------------------------
    ! One-liner: xcheck tests if x solves a certain least-squares problem.
    !
    ! Purpose:   xcheck computes residuals and norms associated with
    ! the vector x and the least-squares problem solved by LSMR.
    ! It determines whether x seems to be a solution to any of three
    ! possible systems:  1. Ax = b
    !                    2. min norm(Ax - b)
    !                    3. min norm(Ax - b)^2 + damp^2 * norm(x)^2.
    !
    ! History:
    ! 07 Feb 1992  Initial design and code.
    !              Michael Saunders, Dept of Operations Research,
    !              Stanford University.
    ! 23 Sep 2007: xcheck implemented as part of this f90 module.
    ! 26 Oct 2012: tol is now an input parameter.
    ! 29 Jun 2013: zxcheck extends xcheck to complex numbers
    !-------------------------------------------------------------------
    intrinsic           :: epsilon
    ! Local variables and arrays
    complex(dp)         :: r(m), v(n)
    real(dp)            :: bnorm, dampsq, rho1, rho2, sigma1, sigma2, eps, &
                           test1, test2, test3, tol2, snorm, xnorm, xsnorm

    ! Local constants
    real(dp), parameter :: power = 0.5_dp
    eps    = epsilon(eps)
    dampsq = damp**2
    tol2   = max( tol, eps**power )

    r(1:m) = -b(1:m)        ! Compute the residual r = b - Ax
    call Aprod1(m,n,x,r)    ! via  r = -b + Ax,
    r(1:m) = -r(1:m)        !      r = -r. 

    v(1:n) = zzero           ! Compute v = A'r
    call Aprod2(m,n,v,r)    ! via  v = 0,  v = v + A'r. 

    bnorm  = dznrm2 (m,b,1)  ! Compute the norms of b, x, r, v.
    xnorm  = dznrm2 (n,x,1)
    rho1   = dznrm2 (m,r,1)
    sigma1 = dznrm2 (n,v,1)
      
    if (nout > 0) write(nout,2200) damp, xnorm, rho1, sigma1

    if (damp == zero) then
       rho2   = rho1
       sigma2 = sigma1
    else
       v(1:n) = v(1:n) - dampsq*x(1:n)  ! v = A'r - damp**2 x.
       rho2   = sqrt(rho1**2 + dampsq*xnorm**2)
       sigma2 = dznrm2 (n,v,1)
       snorm  = rho1/damp
       xsnorm = rho2/damp
       if (nout > 0) write(nout,2300) snorm, xsnorm, rho2, sigma2
    end if

    !-------------------------------------------------------------------
    ! See if x seems to solve Ax = b  or  min norm(Ax - b)
    ! or the damped least-squares system.
    !-------------------------------------------------------------------
    if (bnorm == zero  .and.  xnorm == zero) then
       inform = 0
       test1  = zero
       test2  = zero
       test3  = zero
    else
       inform = 4
       test1  = rho1 / (bnorm + Anorm*xnorm)
       test2  = zero
       if (rho1  > zero) test2  = sigma1 / (Anorm*rho1)
       test3  = test2
       if (rho2  > zero) test3  = sigma2 / (Anorm*rho2)
       
       if (test3 <= tol2) inform = 3
       if (test2 <= tol2) inform = 2
       if (test1 <= tol2) inform = 1
    end if

    if (nout > 0) write(nout,3000) inform, tol2, test1, test2, test3
    return

 2200 format(1p                                             &
      // ' Enter xcheck.     Does x solve Ax = b, etc?'     &
      /  '    damp            =', e10.3                     &
      /  '    norm(x)         =', e10.3                     &
      /  '    norm(r)         =', e15.8, ' = rho1'          &
      /  '    norm(A''r)       =',e10.3, '      = sigma1')
 2300 format(1p/  '    norm(s)         =', e10.3            &
      /  '    norm(x,s)       =', e10.3                     &
      /  '    norm(rbar)      =', e15.8, ' = rho2'          &
      /  '    norm(Abar''rbar) =',e10.3, '      = sigma2')
 3000 format(1p/  '    inform          =', i2               &
      /  '    tol             =', e10.3                     &
      /  '    test1           =', e10.3, ' (Ax = b)'        &
      /  '    test2           =', e10.3, ' (least-squares)' &
      /  '    test3           =', e10.3, ' (damped least-squares)')

  end subroutine zxcheck


end module zlsmrCheckModule
