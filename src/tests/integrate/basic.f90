program test_nf_integrate_basic
    !*   Unit tests for basic integration methods

    use, intrinsic :: iso_fortran_env
    use numfort_arrays
    use numfort_integrate

    use fcore_testing, only: test_suite, test_case
    use fcore_common, only: str, operator(//)

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

    call tests%set_label ("numfort_interpolate_basic unit tests")

    call test_trapezoid (tests)

    ! print test statistics
    call tests%print ()

end subroutine



subroutine test_trapezoid (tests)
    !*  Unit tests for trapezoid rule

    class (test_suite) :: tests
    class (test_case), pointer :: tc
    type (status_t) :: status
    type (str) :: msg

    real (PREC), parameter :: tol = 1d-7
        !   Error tolerance of Fortran vs. Scipy solution

    ! Results for trapezoid rule according to scipy's implementation
    real (PREC), parameter :: INT_SCIPY(size(NPOINTS), NFUNCS) = reshape(&
        [ real(PREC) :: &
        0.721145896,  0.766453333,  0.779700174,  0.78579685 ,  0.798975959, &
        1.76372449 ,  1.192210267,  1.044767731,  0.985702522,  0.90137377, &
        1.73416246 ,  1.722257492,  1.720049244,  1.719276089,  1.718296438, &
        0.605555556,  0.558333333,  0.55617284 ,  0.554166667,  0.551254974 ], &
        shape=[size(NPOINTS), NFUNCS])

    real (PREC), dimension(:), allocatable :: x, fx
    real (PREC) :: res, diff
    procedure (fcn1), pointer :: ptr_fcn
    real (PREC) :: a, b
    integer :: i, j, k, n

    tc => tests%add_test ("Trapezoid rule")

    ! Tests with degenerate parameters
    allocate (x(0), fx(0))
    status = NF_STATUS_UNDEFINED
    call trapezoid (fx, res, x=x, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Zero-length input array")
    deallocate (x, fx)

    allocate (x(1), fx(0))
    x(1) = 0.0
    status = NF_STATUS_UNDEFINED
    call trapezoid (fx, res, x=x, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "fx, x of different length")
    deallocate (x, fx)

    allocate (fx(1), source=0.0_PREC)
    status = NF_STATUS_UNDEFINED
    call trapezoid (fx, res, dx=1.0_PREC, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Length-1 fx array")
    deallocate (fx)

    allocate (fx(2), source=0.0_PREC)
    status = NF_STATUS_UNDEFINED
    call trapezoid (fx, res, dx=0.0_PREC, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Invalid dx parameter")
    deallocate (fx)

    ! Test API that uses function argument and integration interval
    ! Test with too few evaluation NPOINTS, need n > 1
    status = NF_STATUS_UNDEFINED
    call trapezoid (fcn1, 0.0_PREC, 0.0_PREC, 1, res, status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "Integrand: n too small")
    
    ! Test with degenerate interval
    status = NF_STATUS_UNDEFINED
    call trapezoid (fcn1, 0.0_PREC, 0.0_PREC, 2, res, status)
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
            n = NPOINTS(k)
            a = BOUNDS(1, i)
            b = BOUNDS(2, i)
            
            ! Test array API with x-values 
            allocate (fx(n), x(n))
            call linspace (x, a, b)
            forall (j=1:n) fx(j) = ptr_fcn(x(j))
            status = NF_STATUS_UNDEFINED
            res = 0.0
            call trapezoid (fx, res, x=x, status=status)
            diff = abs(res - INT_SCIPY(k, i))

            msg = "Func. " // str(i, "i0") // "; array API; n=" // str(n, "i0") &
                // "; res=" // str(res, "f0.8")
            call tc%assert_true (diff < tol .and. status == NF_STATUS_OK, &
                msg%to_char()) 
            deallocate (fx, x)

            ! Test with procedure argument
            res = 0.0
            status = NF_STATUS_UNDEFINED
            call trapezoid (ptr_fcn, a, b, n, res, status)
            diff = abs(res - INT_SCIPY(k, i))
            msg = "Func. " // str(i, "i0") // "; fcn API; n=" // str(n, "i0") &
                // "; res=" // str(res, "f0.8")
            call tc%assert_true (diff < tol .and. status == NF_STATUS_OK, &
                msg%to_char()) 

        end do
    end do


end subroutine

end program
