program test_nf_integrate_basic
    !*   Unit tests for basic integration methods

    use, intrinsic :: iso_fortran_env
    use numfort_arrays
    use numfort_integrate

    use corelib_testing, only: test_suite, test_case
    use corelib_common, only: str, operator(//)

    use tests_nf_integrate_funcs

    implicit none

    integer, parameter :: PREC = real64

    integer, parameter :: NFUNCS = 4
 
    real (PREC), parameter :: BOUNDS(2,4) = reshape( &
        [0.0, 1.0, 1.0, 10.0, 0.0, 1.0, -1.0, 1.0], &
        shape=[2,4])
        !*  Integration boundaries for each function

   integer, parameter :: NPOINTS(*) = [4, 7, 10, 13, 100]


    call test_all()

contains

subroutine test_all ()
    type (test_suite) :: tests

    call tests%set_label ("numfort_interpolate_quadpack unit tests")

    call test_quad (tests)

    ! print test statistics
    call tests%print ()

end subroutine



subroutine test_quad (tests)
    !*  Unit tests for QUADPACK wrapper routines

    class (test_suite) :: tests
    class (test_case), pointer :: tc
    type (status_t) :: status
    type (str) :: msg

    real (PREC), parameter :: tol = 1d-7
        !   Error tolerance of Fortran vs. Scipy solution

    real (PREC), parameter :: epsabs = 1.49d-8, epsrel = 1.49d-8
        !   Default tol. in Scipy's wrapper

    ! Results for QUAD according to scipy's implementation
    real (PREC), parameter :: INT_SCIPY(size(NPOINTS), NFUNCS) = reshape( &
        [ real(PREC) :: &
        0.800001791067, 0.8, 0.8, 0.8, 0.8, &
        0.9, 0.9, 0.9, 0.9, 0.9, &
        1.718281828459, 1.718281828459, 1.718281828459, 1.718281828459, 1.718281828459, &
        0.551254597428, 0.551250154728, 0.551250001122, 0.55125, 0.55125 ], &
        shape=[size(NPOINTS), NFUNCS])

    real (PREC) :: res, diff
    procedure (fcn1), pointer :: ptr_fcn
    real (PREC) :: a, b
    integer :: i, k, nmax

    tc => tests%add_test ("QUAD routine")

    ! Test with degenerate interval
    status = NF_STATUS_UNDEFINED
    call quad (fcn1, 0.0_PREC, 0.0_PREC, res, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. (res == 0.0), &
        "Integrand: integrate over singleton set")

    ! Iterate over test functions; for each test function, test with
    ! all number of points.
    do i = 1, NFUNCS
        select case (i)
        case (1)
            ptr_fcn => fcn1
        case (2)
            ptr_fcn => fcn2
        case (3)
            ptr_fcn => fcn3
        case (4)
            ptr_fcn => fcn4
        case default
            ptr_fcn => null()
        end select

        do k = 1, size(NPOINTS)
            nmax = NPOINTS(k)
            a = BOUNDS(1, i)
            b = BOUNDS(2, i)

            ! Test with procedure argument
            res = 0.0
            status = NF_STATUS_UNDEFINED
            call quad (ptr_fcn, a, b, res, epsabs=epsabs, epsrel=epsrel, &
                    nmax=nmax, status=status)
            diff = abs(res - INT_SCIPY(k, i))
            msg = "Func. " // str(i, "i0") // "; res=" // str(res, "f0.8") &
                // "; diff=" // str(diff, "es9.2e3")
            call tc%assert_true (diff < tol, msg%to_char()) 

        end do
    end do


end subroutine

end program
