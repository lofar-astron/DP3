!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File lsmrReadMtxModule.f90
!
! Calls mm_ioModule routines to input a sparse matrix stored in
! Matrix Market format.
!
! 13 Aug 2012: Testing MM_IO routines from John Burkardt
!              http://people.sc.fsu.edu/~jburkardt/f_src/mm_io/
!              Had to change some lines of mm_io.f90 to enable compilation.
!              Michael Saunders, SOL and ICME, Stanford University.
!
!              Intended for testing MINRESQLP.f90
!              (Sou-Cheng Choi, Chris Paige and Michael Saunders)
!
! 13 Aug 2012: Made mm_io.f90 into mm_ioModule.f90.
!              minresqlpReadMtxModule.f90 and mm_ioModule.f90 compile ok
! 28 Oct 2012: Debugged mm_ioModule.f90.  Can use it here now.
! 05 Jan 2013: ReadMtxSize reads the first few lines of data to find nnz.
!              ReadMtx     reads the whole file.
!              We can use this module for both real and complex MINRES-QLP,
!              at the expense of allocating space for 3 matrices instead of
!              of just one of dval, rval, cval.
! 29 Jun 2013: Adapted minresqlp version for lsmr
!--------------------------------------------------------------------------

module lsmrReadMtxModule
  use  lsmrDataModule, only : dp, sp, ip
  use  mm_ioModule, only : get_unit, mm_comment_read,          &
                           mm_file_read, mm_file_write,        &
                           mm_header_print, mm_header_read,    &
                           mm_size_print, mm_size_read_string, &
                           mm_values_print_some, timestamp

  implicit none

  private
  public              :: ReadMtxSize, ReadMtx
  integer(ip), public :: nnzmax

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine ReadMtxSize( input_file, input_unit, &
                          id, type, rep, field, symm, nrow, ncol, nnz )

    character(80), intent(in)    :: input_file
    integer(ip)  , intent(out)   :: input_unit
    character(14), intent(out)   :: id
    character(10), intent(out)   :: rep
    character( 6), intent(out)   :: type
    character( 7), intent(out)   :: field
    character(19), intent(out)   :: symm
    integer(ip)  , intent(out)   :: nrow, ncol, nnz

    ! ReadMtxSize mimics the beginning of mm_file_read
    ! to read past the first few lines of a Matrix Market *.mtx file
    ! to the line that says how big the matrix is.
    ! The calling routine should rewind input_unit
    ! before calling ReadMtx.
    !
    ! 05 Jan 2012: First version.

    integer(ip)     :: ios
    character(1024) :: chartmp

    call timestamp( )

    write(*, '(a)') ' '
    write(*, '(a)') 'lsmrReadMtxModule::ReadMtxSize():'
    write(*, '(a)') '  Reading ' // trim(input_file)
    write(*, '(a)') '  using the MM_IO library.'

    call get_unit( input_unit )

    open( unit = input_unit, file = input_file, status = 'old', iostat = ios )

    if (ios /= 0) then
       write(*, '(a)') ' '
       write(*, '(a)') 'ReadMtxSize - Fatal error!'
       write(*, '(a)') '  Could not open the input file.'
       write(*, *) ios
       stop
    end if

    ! Read mtx file header, then
    ! read through the comment lines.

    call mm_header_read ( input_unit, id, type, rep, field, symm )
    do
       call mm_comment_read ( input_unit, chartmp )
       if ( chartmp(1:1) /= '%' ) then
          exit
       end if
    end do

    ! Now read the size line
    call mm_size_read_string ( chartmp, rep, symm, nrow, ncol, nnz )

  end subroutine ReadMtxSize

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine ReadMtx( input_file, input_unit, &
                      id, rep, field, symm, nrow, ncol, nnz, &
                      indx, jndx, ival, rval, dval, cval )

    character(80), intent(in)    :: input_file
    integer(ip)  , intent(in)    :: input_unit
    integer(ip)  , intent(in)    :: nnz
    character(14), intent(out)   :: id
    character(10), intent(out)   :: rep
    character(7) , intent(out)   :: field
    character(19), intent(out)   :: symm
    integer(ip)  , intent(out)   :: nrow, ncol
    integer(ip)  , intent(out)   :: indx(nnz), jndx(nnzmax)
    integer(ip)  , intent(out)   :: ival(nnz)
    real(sp)     , intent(out)   :: rval(nnz)
    real(dp)     , intent(out)   :: dval(nnz)
    complex(dp)  , intent(out)   :: cval(nnz)

    ! ReadMtx reads a "coordinate real symmetric" matrix in MM format.
    ! The number of nonzeros (nnz) is assumed to be known,
    ! with help from ReadMtxSize.
    !
    ! 05 Jan 2013: ReadMtxSize is now called before ReadMtx.

    character( 6)  :: type
    integer(ip)    :: ihi
    integer(ip)    :: ilo
 

    call mm_file_read( input_unit, id, type, rep, field, symm, nrow, ncol, &
                       nnz, nnzmax, indx, jndx, ival, rval, dval, cval )

    close( unit = input_unit )

    call mm_header_print( input_file, id, type, rep, field, symm )

    call mm_size_print( input_file, rep, symm, nrow, ncol, nnz )

    ilo = 1
    ihi = 10

    call mm_values_print_some( rep, field, nnz, indx, jndx, ival, rval, &
                               dval, cval, ilo, ihi )

    !nout   = 6
    !open(nout,file='temp.txt',status='unknown')
    !call mm_file_write ( nout, id, type, rep, field, symm, nrow, &
    !  ncol, nnz, indx, jndx, ival, rval, dval, cval )

  end subroutine ReadMtx

end module lsmrReadMtxModule
