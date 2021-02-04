!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File zlsmrTestModule.f90
!
!    Hprod   Aprod1   Aprod2    lstp     lsmrtest
!
! Complex version of lsmrTestModule
!
! 29 Jun 2013: Derived from zlsqr version
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module zlsmrTestModule

  use  zlsmrDataModule,     only  : dp, sp, ip, zzero, zero, one, zone
  use  zlsmrblasInterface,  only  : dznrm2, zscal
  use  zlsmrModule,         only  : zLSMR
  use  zlsmrCheckModule,    only  : zAcheck, zxcheck
  use  lsmrReadMtxModule,   only  : ReadMtxSize, ReadMtx, nnzmax

  implicit none
  private
  public      :: lsmrtest, zlsmrtestMtxCCH
  private     :: Hprod, Aprod1, Aprod2, lstp
  intrinsic   :: dot_product

  ! DYNAMIC WORKSPACE DEFINED HERE.
  ! They are allocated in lsmrtest and used by Aprod1, Aprod2.

  real(dp), allocatable    :: d(:)                             ! These define A = Y D Z.
  complex(dp), allocatable :: hy(:), hz(:), wm(:), wn(:)       ! Work vectors.

  integer(ip)              :: nnz      ! These ones are used by the
  integer(ip), allocatable :: indx(:)  ! Matrix Market Aprod routines
  integer(ip), allocatable :: jndx(:)  ! AprodMtxCDS, AprodMtxCPS, AprodMtxCRS
  real(dp)   , allocatable :: dval(:)  !
  integer(ip), allocatable :: ival(:)  !
  real(sp)   , allocatable :: rval(:)  !
  complex(dp), allocatable :: cval(:)  !

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Hprod (n,z,x)

    integer(ip),  intent(in)       :: n
    complex(dp), intent(in)        :: z(n)
    complex(dp), intent(inout)     :: x(n)

    !-------------------------------------------------------------------
    ! Hprod  applies a Householder transformation stored in z
    ! to return x = (I - 2*z*z')*x.
    !
    ! 23 Sep 2007: Fortran 90 version.
    !-------------------------------------------------------------------

    complex(dp) :: s

    s = dot_product( z,x )
    s = s + s
    x = x - s*z
  end subroutine Hprod

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine Aprod1(m,n,x,y)

    integer(ip),  intent(in)   :: m,n
    complex(dp), intent(in)    :: x(n)
    complex(dp), intent(inout) :: y(m)

    !-------------------------------------------------------------------
    ! Aprod1 computes y = y + A*x without altering x,
    ! where A is a test matrix of the form  A = Y*D*Z,
    ! and the matrices D, Y, Z are represented by
    ! the allocatable vectors d, hy, hz in this module.
    !
    ! 23 Sep 2007: Fortran 90 version.
    !-------------------------------------------------------------------

    intrinsic           :: min
    integer(ip)             :: minmn

    minmn = min(m,n)
    wn    = x
    call Hprod (n,hz,wn)
    wm(1:minmn) = d(1:minmn)*wn(1:minmn)
    wm(n+1:m)   = zzero
    call Hprod (m,hy,wm)
    y = y + wm

  end subroutine Aprod1

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Aprod2(m,n,x,y)

    integer(ip),  intent(in)   :: m,n
    complex(dp), intent(inout) :: x(n)
    complex(dp), intent(in)    :: y(m)

    !-------------------------------------------------------------------
    ! Aprod2 computes x = x + A'*y without altering y,
    ! where A is a test matrix of the form  A = Y*D*Z,
    ! and the matrices D, Y, Z are represented by
    ! the allocatable vectors d, hy, hz in this module.
    !
    ! 23 Sep 2007: Fortran 90 version.
    !-------------------------------------------------------------------

    intrinsic           :: min
    integer(ip)             :: minmn

    minmn = min(m,n)
    wm    = y
    call Hprod (m,hy,wm)
    wn(1:minmn) = d(1:minmn)*wm(1:minmn)
    wn(m+1:n)   = zzero
    call Hprod (n,hz,wn)
    x = x + wn

  end subroutine Aprod2

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine lstp  (m,n,nduplc,npower,damp,x,b,Acond,rnorm)

    integer(ip),  intent(in)   :: m, n, nduplc, npower
    real(dp), intent(in)       :: damp
    complex(dp), intent(inout) :: x(n)
    complex(dp), intent(out)   :: b(m)
    real(dp), intent(out)      :: Acond, rnorm

    !-------------------------------------------------------------------
    ! lstp  generates a sparse least-squares test problem of the form
    !           (   A    )*x = ( b ) 
    !           ( damp*I )     ( 0 )
    ! for by LSMR, or a sparse underdetermined system
    !            Ax + damp*s = b
    ! for by LSMR or CRAIG.  The matrix A is m by n and is
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
    ! 29 Jun 2013: Extended to complex A and b
    !-------------------------------------------------------------------

    intrinsic            :: min, cos, sin, sqrt
    integer(ip)          :: i, j, minmn
    real(dp)             :: alfa, beta, dampsq, fourpi
    real(dp)             :: t

    !-------------------------------------------------------------------
    ! Make two vectors of norm 1.0 for the Householder transformations.
    ! fourpi  need not be exact.
    !-------------------------------------------------------------------
    minmn  = min(m,n)
    dampsq = damp**2
    fourpi = 4.0_dp * 3.141592_dp
    alfa   = fourpi / m
    beta   = fourpi / n

!Make this test harder
    do i = 1,m
       hy(i) = cmplx(sin( alfa*i ), sin( alfa*i ), dp)
    end do

    do i = 1,n
       hz(i) = cmplx(cos( beta*i ), cos( beta*i ), dp)
    end do

    alfa   = dznrm2 ( m, hy, 1 )
    beta   = dznrm2 ( n, hz, 1 )
    call zscal ( m, cmplx((- one/alfa),0,dp), hy, 1 )
    call zscal ( n, cmplx((- one/beta),0,dp), hz, 1 )

    !-------------------------------------------------------------------
    ! Set the diagonal matrix D.  These are the singular values of A.
    !-------------------------------------------------------------------
    do i = 1,minmn
       j    = (i-1+nduplc) / nduplc
       t    =  j*nduplc
       t    =  t / minmn
       d(i) =  t**npower
    end do

    Acond  = (d(minmn)**2 + dampsq) / (d(1)**2 + dampsq)
    Acond  = sqrt( Acond )

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
       wn(m+1:n) = zzero
       call Hprod (n,hz,wn)
       x = wn
    end if

    ! Let r = Y rbar and xbar = Z x.
    ! Solve D r1bar = damp^2 x1bar, where r1bar is in wm(1:minmn).

    wm(1:minmn)   = dampsq * wn(1:minmn)/d(1:minmn)

    wm(minmn+1:m) = zone      ! Set r2bar to be anything (empty if m <= n).
    call Hprod (m,hy,wm)     ! Form r = Y rbar  in wm.

    rnorm  = dznrm2 (m,wm,1)  ! Compute norm(r)
    b      = wm              ! and  b = r + Ax.
    call Aprod1(m,n,x,b)

  end subroutine lstp

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine AprodMtxCCH (m,n,x,y)

    integer(ip), intent(in)  :: m, n
    complex(dp), intent(in)  :: x(n)
    complex(dp), intent(inout) :: y(n)

    !-------------------------------------------------------------------
    ! AprodMtxCCH computes y = A*x for some Hermitian matrix A.
    ! stored in Matrix Market CCH format (coordinate complex hermitian).
    ! Only subdiagonal and diagonal elements are in (indx, jndx, cval).
    !
    ! 16 Sep 2012: Diagonals treated as REAL here to guard against
    !              CCS format (coordinate complex symmetric), which
    !              might have complex diagonals.
    !-------------------------------------------------------------------

    intrinsic    :: conjg
    integer(ip)  :: i, j, k
    complex(dp)  :: d

    !y(1:n) = zzero

    do k = 1, nnz
       i = indx(k)
       j = jndx(k)
       d = cval(k)

       if (i > j) then          ! d = subdiagonal
          y(i) = y(i) + d*x(j)
          y(j) = y(j) + conjg(d)*x(i)
       else                     ! i = j, d = diagonal
          y(i) = y(i) + real(d)*x(i)
       end if

       ! if (k <= 10) write(*,*) '  ', i, ' ', j, ' ', d,  ' ', y(i)
    end do

  end subroutine AprodMtxCCH

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  subroutine AprodMtxCCH2 (m,n,x,y)

    integer(ip), intent(in)  :: m, n
    complex(dp), intent(inout)  :: x(n)
    complex(dp), intent(in) :: y(n)

    call AprodMtxCCH(m, n, y, x)

  end subroutine AprodMtxCCH2

    

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine zlsmrtestMtxCCH(input_file, consis, nout, tol)

    character(80),  intent(in)           :: input_file
    logical,        intent(in)           :: consis
    integer(ip),    intent(in)           :: nout
    real(dp),       intent(in), optional :: tol

    ! 30 Oct 2012: Use zxcheck to check computed x from zminresqlp.
    ! 02 Jan 2013: Print n in the "zminresqlp  appears to be ..." message
    !              to help identify the problem.
    ! 29 Jun 2013: Adapt zminresqlp version to zlsmr

    intrinsic      :: real, present, epsilon

    real(dp)       :: eps

    integer(ip)    :: input_unit, nrow, ncol
    character(14)  :: id
    character(10)  :: rep
    character( 6)  :: type
    character(7)   :: field
    character(19)  :: symm

    complex(dp), allocatable  :: b(:), x(:), r1(:), w(:)

    logical     :: checkA, disable, precon
    integer(ip) :: n, j, itn, itnlim, istop, inform, localSize
    real(dp)    :: shift, Anorm, Acond, Arnorm, rnorm, rtol, xnorm
    real(dp)    :: maxxnorm, TranCond, Acondlim
    real(dp)    :: realj, realn, relTol, test1, test2
    real(dp)    :: atol, btol, conlim, damp

    character(len=*), parameter :: headerStr =              &
       "(// '----------------------------------------'"  // &
       "  / ' Test of zLSMR on an MM CCH matrix '"  // &
       "  / '----------------------------------------')"

    write(nout, headerStr)
    write(  * , headerStr)

    call ReadMtxSize( input_file, input_unit, &
                      id, type, rep, field, symm, nrow, ncol, nnz )
    rewind( input_unit )

    ! Now we know the size of the problem.
    ! We should allocate only the arrays that will be used by the MM routines.
    ! For simplicity we allocate them all.

    nnzmax = nnz
    allocate( indx(nnz), jndx(nnz), ival(nnz) )
    allocate( rval(nnz) )   ! CRS
    allocate( dval(nnz) )   ! CDS
    allocate( cval(nnz) )   ! CCH

    call ReadMtx( input_file, input_unit, &
                  id, rep, field, symm, nrow, ncol, nnz, &
                  indx, jndx, ival, rval, dval, cval )

    n = nrow

    eps = epsilon(eps)
    relTol = 10_dp*eps
    if (present(tol)) then
       relTol = tol
    end if

    allocate( b(n) )
    allocate( x(n) )
    allocate( r1(n) )
    allocate( w(n) )

    realn = real(n,dp)
    do j = 1, n
       realj = real(j,dp)
       x(j)  = cmplx(realj,realj,dp) / realn
       b(j)  = zone
    end do

    write(nout,*) 'consis     = ', consis
    write(nout,*) 'input_file = ', trim(input_file)
    write(nout,*) 'n = ',  n, '  nnz = ', nnz

    checkA   = .true.          ! Set other parameters and solve.
    disable  = .false.
    precon   = .false.
    itnlim   = n*20
    rtol     = 1.0e-8_dp
    maxxnorm = 1.0e+8_dp
    TranCond = 1.0e+8_dp
    Acondlim = 1.0e+15_dp
    shift    = zero

    if (consis) then
      do j = 1, n
         b(j) = zzero
      end do
      call AprodMtxCCH (n,n,x,b)   ! b = A*x
      write(nout,*) ' '
      write(nout,*) 'norm(b) =', dznrm2(n,b,1)
      write(nout,*) 'Some of the x defining b'
      do j = 1, min(n,5)
         write(nout,*) j, x(j)
      end do
    end if

    write(nout,*) 'Some of b'
    do j = 1, min(n,5)
       write(nout,*) j, b(j)
    end do

!    if (debug) then
!       write(*,*)
!       write(*,*)  'n = ', n, ', shift = ', shift,  ', consis = ', consis, ', nout = ', nout
!       write(*,*)
!       write(*,*)  'checkA = ', checkA, 'itnlim = ', itnlim, ', nout = ', nout,                &
!                   ', maxxnorm = ', maxxnorm, ', TranCond = ', TranCond, ', Acondlim = ', Acondlim
!    end if

    atol      = 1e-14_dp                         ! Fixed accuracy.
    btol      = atol
    conlim    = 1e10_dp
    itnlim    = 4*(n + n + 50)
    damp      = 1e-6_dp
    localSize = 0

    do j = 1, n
       x(j) = zzero
    end do

    call zLSMR  (n, n, AprodMtxCCH, AprodMtxCCH2, b, damp,          &
                 atol, btol, conlim, itnlim, localSize, nout,       &
                 x, istop, itn, Anorm, Acond, rnorm, Arnorm, xnorm )

    call zxcheck(n, n, AprodMtxCCH, AprodMtxCCH2, b, damp, x, Anorm, &
                 tol, nout, inform)

    if (inform <= 3) then
       write(nout,3000) itn, rnorm, Arnorm
    else
       write(nout,3100) itn, rnorm, Arnorm
    end if

    write(nout,*) ' '
    write(nout,*) 'norm(x) =', dznrm2(n,x,1)
    write(nout,*) 'Some of the computed x'
    do j = 1, min(n,5)
       write(nout,*) x(j)
    end do

    deallocate( indx, jndx, ival, rval, dval, cval )
    deallocate( b, x, r1, w )

 3000 format(1p, " LSMR  appears to be successful.  Itns =", i7,  &
         "  norm(r) =", e9.2, "  norm(A'r) =", e9.2)
 3100 format(1p, " LSMR  appears to have failed.    Itns =", i7,  &
         "  norm(r) =", e9.2, "  norm(A'r) =", e9.2)

  end subroutine zlsmrtestMtxCCH

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine lsmrtest(m,n,nduplc,npower,damp,localSize,nout)

    integer(ip),  intent(in)    :: m, n, nduplc, npower, localSize, nout

    real(dp), intent(in)    :: damp

    !-------------------------------------------------------------------
    ! This is an example driver routine for running LSMR.
    ! It generates a test problem, solves it, and examines the results.
    !
    !------------------------------------------------------------------------

    intrinsic       :: epsilon
    ! Local arrays and variables
    complex(dp)     :: b(m), x(n), xtrue(n)
    integer(ip)     :: inform, istop, itn, itnlim, j, minmn, nprint
    real(dp)        :: atol, btol, conlim, Anorm, Acond,             &
                      rnorm, Arnorm, enorm, etol, wnorm, xnorm, eps

    ! Local constants

    character(len=*),  &
              parameter :: line = '----------------------------------'

    eps = epsilon(eps)

    if (eps > 1e-9_dp) then
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
       xtrue(j) = cmplx(0.1_dp*j, 0.1_dp*j,dp)   ! For least-squares problems, this is it.
    end do                                       ! If m<n, lstp will alter it.

    minmn  = min(m,n)                            ! Allocate arrays.
    allocate( d(minmn), hy(m), hz(n) )           ! Vectors defining A = Y D Z.
    allocate( wm(m), wn(n) )                     ! Work vectors for Aprod1, Aprod2.

    call lstp  (m,n,nduplc,npower,damp, &        ! Generate test problem.
                xtrue,b,Acond,rnorm)             ! If m<n, xtrue is altered.

    write(nout,1000) line,line,m,n,nduplc,npower,damp,Acond,rnorm, &
                     line,line

    call zAcheck(m,n,Aprod1,Aprod2,nout,inform)   ! Check Aprod1, Aprod2.
          
    if (inform > 0) then
       write(nout,'(a)') 'Check tol in subroutine zAcheck in zlsmrCheckModule'
       stop
    end if

    !------------------------------------------------------------------------
    ! Set input parameters for LSMR
    ! and solve the problem defined by Aprod1, Aprod2, b, damp.
    ! Next line would ask for standard errors only if they are well-defined.
    !------------------------------------------------------------------------

    atol   = 3.18e-16_dp                            ! Fixed accuracy.
    btol   = atol
    conlim = 1000.0_dp * Acond
    itnlim = 4*(m + n + 50)

    call zLSMR  (m, n, Aprod1, Aprod2, b, damp,     &
                 atol, btol, conlim, itnlim, localSize, nout,        &
                x, istop, itn, Anorm, Acond, rnorm, Arnorm, xnorm )

    call zxcheck(m, n, Aprod1, Aprod2, b, damp, x, &    ! Check x
                Anorm, atol, nout,                &
                inform)

    nprint = min(m,n,8)
    write(nout,2500)    (j, abs(x(j)), j=1,nprint)    ! Print some of the solution
    !write(nout,2500)    (j, abs(xtrue(j)), j=1,nprint)    ! Print some of the solution


    wn     = x - xtrue                           ! Print a clue about whether
    wnorm  = dznrm2 (n,wn,1)                      ! the solution looks OK.
    xnorm  = dznrm2 (n,xtrue,1)
    enorm  = wnorm/(one + xnorm)
    etol   = 1.0e-3_dp

    if (inform <= 3) then
       write(nout,3000) enorm
    else
       write(nout,3100) enorm
    end if

    deallocate( wn, wm, hz, hy, d )              ! Free arrays (in reverse)
    return

 1000 format(1p &
      // 1x, 2a &
      /  ' Least-Squares Test Problem      P(', 4i5, e12.2, ' )'     &
      /  ' Condition no. =', e12.4, '     Residual function =', e17.9&
      /  1x, 2a)
 2500 format(/' Solution  abs(x):'          / 4(i6, g14.6))
 3000 format(1p &
      /  ' LSMR  appears to be successful.'                          &
      /  ' Relative error in  x  =', e10.2)
 3100 format(1p &
      /  ' LSMR  appears to have failed.  '                          &
      /  ' Relative error in  x  =', e10.2)
  end subroutine lsmrtest
end module zlsmrTestModule
