


program test_io_fixed
    !*  Unit tests for routines performing fixed-format text I/O.

    use, intrinsic :: iso_fortran_env

    use numfort_stats
    use numfort_common
    use numfort_common_testing
    use numfort_io

    use fcore_strings
    use fcore_testing

    implicit none

    integer, parameter :: PREC = real64

    call test_all ()


    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("IO unit tests for fixed-format files")

    call test_fixed (tests)

    call tests%print()

end subroutine


subroutine test_fixed (tests)
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    real (PREC), dimension(:), allocatable :: rwork1d
    real (PREC), dimension(:), allocatable :: dat1d_out, dat1d_in
    integer, dimension(:), allocatable :: idat1d_out, idat1d_in
    real (PREC), dimension(:,:), allocatable :: dat2d_out, dat2d_in
    real (PREC), dimension(:,:,:), allocatable :: dat3d_out, dat3d_in
    integer :: n, m, k
    character (100) :: path, msg
    character (20) :: fmt
    type (status_t) :: status

    tc => tests%add_test ("Unit tests for READ_FIXED / WRITE_FIXED")

    path = '/tmp/tests_numfort_io_fixed.txt'

    call set_seed (1234)

    ! === Test with REAL 1d-arrays ===
    m = 50
    allocate (dat1d_in(m), dat1d_out(m))

    call random_number (dat1d_out)

    fmt = '(*(5(f8.4),:,/))'
    call write_fixed (path, fmt, dat1d_out, msg=msg, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        "Writing 1d REAL data in flat form, reshaping via FMT")

    call read_fixed (path, fmt, dat1d_in, msg=msg, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (dat1d_out, dat1d_in, atol=1.0d-3), &
        "Reading 1d REAL data in flat form, reshaping via FMT")

    deallocate (dat1d_in, dat1d_out)

    ! === Test with INTEGER 1d-arrays ===
    m = 50
    allocate (idat1d_in(m), idat1d_out(m), rwork1d(m))

    call random_number (rwork1d)
    idat1d_out(:) = int(100.0d0 * rwork1d)

    fmt = '(*(5(i3),:,/))'
    call write_fixed (path, fmt, idat1d_out, msg=msg, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        "Writing 1d INTEGER data in flat form, reshaping via FMT")

    call read_fixed (path, fmt, idat1d_in, msg=msg, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all(idat1d_out == idat1d_in), &
        "Reading 1d INTEGER data in flat form, reshaping via FMT")

    deallocate (idat1d_out, idat1d_in, rwork1d)

    ! === Test with 2d-arrays ===
    n = 10
    m = 11
    allocate (dat2d_out(m,n), dat2d_in(m,n))

    call random_number (dat2d_out)

    ! Write in transposed form
    call write_fixed (path, '(*(f10.6))', dat2d_out, transform='transpose', &
        msg=msg, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        "Writing file in transposed order")

    call read_fixed (path, '(*(f10.6))', dat2d_in, transform='transpose', &
        msg=msg, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (dat2d_out, dat2d_in, atol=1.0d-5), &
        "Reading written file in transposed order")

    ! Write in non-transposed order such that one line in the file corresponds
    ! to one row in the array
    call write_fixed (path, '(*(f10.6))', dat2d_out, transform='none', &
        msg=msg, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        "Writing file in non-transposed order")

    call read_fixed (path, '(*(f10.6))', dat2d_in, transform='none', &
        msg=msg, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (dat2d_out, dat2d_in, atol=1.0d-5), &
        "Reading written file in non-transposed order")

    ! Write in flatten (column-major) order
    call write_fixed (path, '(*(f10.6))', dat2d_out, transform='flatten', &
        msg=msg, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        "Writing file in flatten form")

    call read_fixed (path, '(*(f10.6))', dat2d_in, transform='flatten', &
        msg=msg, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (dat2d_out, dat2d_in, atol=1.0d-5), &
        "Reading written file in flatten form")

    ! Test with flattened transform but apply format to reshape array
    fmt = '(*(2(f10.8),:,/))'
    call write_fixed (path, fmt, dat2d_out, transform='flatten', &
        msg=msg, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        "Writing file in flat form, reshaping via FMT")

    call read_fixed (path, fmt, dat2d_in, transform='flatten', &
        msg=msg, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (dat2d_out, dat2d_in, atol=1.0d-5), &
        "Reading written file in flat form, reshaping via FMT")

    ! === Test with 3d array ===
    m = 5
    n = 4
    k = 2
    allocate (dat3d_out(m,n,k), dat3d_in(m,n,k))

    call random_number (dat3d_out)

    fmt = '(*(4(f15.8),:,/))'
    call write_fixed (path, fmt, dat3d_out, transform='flatten', &
        msg=msg, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        "Writing 3d data in flat form, reshaping via FMT")

    call read_fixed (path, fmt, dat3d_in, transform='flatten', &
        msg=msg, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (dat3d_out, dat3d_in, atol=1.0d-7), &
        "Reading 3d data in flat form, reshaping via FMT")

end subroutine


end program

