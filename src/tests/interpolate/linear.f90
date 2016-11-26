program test_interp_linear

    use iso_fortran_env
    use corelib_testing
    use numfort_arrays
    use numfort_interpolate
    implicit none

    integer, parameter :: PREC = real64
    integer, parameter :: INTSIZE = int32

    call test_all ()

contains

subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("numfort_interpolate unit tests")

    call test_scalar (tests)
    ! call test_sub2ind (tests)

    ! print test statistics
    call tests%print ()

end subroutine


subroutine test_scalar (tests)

    class (test_suite) :: tests
    class (test_case), pointer :: tc

    real (PREC) :: x, fx
    real (PREC), dimension(10) :: xp, fp, fxx
    integer :: i

    tc => tests%add_test ("scalar arguments")

    ! identity function on [0, 1]
    call linspace (xp, 0.0d0, 1.0d0)
    fp = xp

    ! evaluate at all points in xp
    do i = 1, size(xp)
        call interp_linear (xp(i), xp, fp, fxx(i))
    end do

    call tc%assert_true (all(abs(fxx-fp)<1d-12), "Interpolate identity function")

end subroutine

end program
