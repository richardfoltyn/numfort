
program test_optimize_minpack
    !*  Unit tests for MINPACK wrapper.

    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env

    use fcore_common, FC_status_t => status_t
    use fcore_testing
    use numfort_arrays
    use numfort_common_testing
    use numfort_optimize, workspace => workspace_real64, &
        optim_result => optim_result_real64

    integer, parameter :: PREC = real64


    type, extends(args_data) :: args_linear
        !*  Stores matrix A and vector b of linear equation system Ax=b
        real (PREC), dimension(:,:), allocatable :: mat
        real (PREC), dimension(:), allocatable :: b
    end type

    interface dynamic_cast
        procedure cast_args_to_args_linear
    end interface

    call test_all ()

    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("MINPACK wrapper unit tests")

    call test_root_lm_linear (tests)
    call test_root_lm_poly (tests)

    call tests%print ()
end subroutine



subroutine test_root_lm_linear (tests)
    !*  Tests for MINPACK's LM root-finder with a linear equation system
    !   (no need for MINPACK routines in this case, of course).
    class (test_suite) :: tests

    class (test_case), pointer :: tc
    type (args_linear) :: args
    integer :: n, m
    type (workspace) :: work
    type (optim_result) :: res
    real (PREC), dimension(:), allocatable :: fx, x_desired, x

    tc => tests%add_test ("Tests for LM root-finder (linear objective)")

    ! Test 1: User-provided function evaluates function and its Jacobian
    n = 3
    m = 3
    allocate (args%mat(n,n), source=0.0_PREC)
    allocate (args%b(n), source=0.0_PREC)

    args%mat(1,1:n) = [1,2,3]
    args%mat(2,2:n) = [4,5]
    args%mat(3,3) = 6

    args%b(1:n) = [-10.0, 1.0, 5.0]

    allocate (x(n), fx(m), x_desired(m))

    x(:) = 0.0
    ! Solution from Scipy's wrapper of MINPACK's LMDER
    x_desired(:) = [-10.916666666666666d0, -0.79166666666666663d0, 0.83333333333333326d0]

    call root_lm (fcn_jac_linear, x, args, fx, work=work, res=res)

    ! General checks not specific not this test case: X and FX are identical
    ! with those stored in result object
    call tc%assert_true (all(res%x == x), "FCN and JAC given: X == RES%X")
    call tc%assert_true (all(res%fx == fx), "FCN and JAC given: f(X) == RES%FX")

    call tc%assert_true (res%status == NF_STATUS_OK, &
        "FCN and JAC given: exit status")
    call tc%assert_true (all_close (x, x_desired, atol=1.0d-8), &
        "FCN and JAC given: solution X close to Scipy")
    call tc%assert_true (maxval(abs(res%fx)) < 1.0d-6, &
        "FNC and JAC given: solution f(X) close to zero vector")

    ! Test 2: Use numeric differentiation
    x(:) = 0.0
    fx(:) = -100.0
    call root_lm (fcn_linear, x, args, fx, ndiff=.true., work=work, res=res)

    call tc%assert_true (res%status == NF_STATUS_OK, &
        "FCN given: exit status")
    call tc%assert_true (all_close (x, x_desired, atol=1.0d-8), &
        "FCN given: solution X close to Scipy")
    call tc%assert_true (maxval(abs(res%fx)) < 1.0d-6, &
        "FNC given: solution f(X) close to zero vector")

    deallocate (fx, x_desired, x)
    deallocate (args%mat, args%b)


end subroutine



subroutine test_root_lm_poly (tests)
    !*  Test MINPACK's LMDER with vector-valued polynomial function with
    !   domain R^4.
    class (test_suite) :: tests

    class (test_case), pointer :: tc

    type (workspace) :: work
    type (optim_result) :: res
    integer, parameter :: NDIM = 4
    real (PREC), dimension(NDIM) :: x, x_desired, fx

    tc => tests%add_test ("Tests for LM root-finder with polynomial objective")

    ! Result from Scipy's MINPACK wrapper
    x_desired = [ 0.310535466871667d0, -2.071117380200093d0,  &
        4.096662650196489d0,  2.365209283870028d0]

    ! Test with user-provided function and Jacobian evaluator
    x = 0.0
    call root_lm (fcn_jac_poly, x, fx, work=work, res=res)

    call tc%assert_true (res%status == NF_STATUS_OK, &
        "FCN and JAC given: exit status")
    call tc%assert_true (all_close (x, x_desired, atol=1.0d-8), &
        "FCN and JAC given: solution close to Scipy")
    call tc%assert_true (maxval(abs(res%fx)) < 1.0d-6, &
        "FCN and JAC given: solution f(X) close to zero vector")

    ! Test with numerical differentiation
    x = 0.0
    fx = -100.0
    ! Scipy result using numerical differentiation with default step size
    ! Note that this result is completely different from the one obtained
    ! with the analytically computed Jacobian.
    x_desired = [  9.999272040754798d0, 14.139850567306441d0, &
        -0.035508624284976d0, -0.020500913789493d0]
    call root_lm (fcn_poly, x, fx, ndiff=.true., work=work, res=res)

    call tc%assert_true (res%status == NF_STATUS_OK, &
        "FCN given: exit status")
    call tc%assert_true (all_close (x, x_desired, atol=1.0d-8), &
        "FCN given: solution close to Scipy")
    call tc%assert_true (maxval(abs(res%fx)) < 1.0d-6, &
        "FCN given: solution f(X) close to zero vector")

end subroutine


subroutine fcn_jac_linear (x, args, fx, fpx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), dimension(:), contiguous, optional :: fx
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: fpx

    type (args_linear), pointer :: ptr_args
    real (PREC), dimension(:), allocatable :: Ax
    integer :: m, n

    call dynamic_cast (args, ptr_args)

    m = size(ptr_args%mat, 1)
    n = size(ptr_args%mat, 2)

    allocate (Ax(m))

    Ax = matmul(ptr_args%mat, x)

    if (present(fx)) then
        ! Objective
        fx = Ax - ptr_args%b
    end if

    ! Jacobian
    if (present(fpx)) then
        fpx = ptr_args%mat
    end if
end subroutine



subroutine fcn_linear (x, args, fx)
    !*  Linear objective, no user-provided Jacobian
    real (PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), dimension(:), contiguous :: fx

    type (args_linear), pointer :: ptr_args
    real (PREC), dimension(:), allocatable :: Ax
    integer :: m

    call dynamic_cast (args, ptr_args)

    m = size(ptr_args%mat, 1)

    allocate (Ax(m))

    Ax = matmul(ptr_args%mat, x)

    ! Objective
    fx = Ax - ptr_args%b

end subroutine


subroutine cast_args_to_args_linear (tgt, ptr)
    class (args_data), intent(in), target :: tgt
    type (args_linear), intent(inout), pointer :: ptr

    nullify (ptr)

    select type (tgt)
    type is (args_linear)
        ptr => tgt
    end select

end subroutine


subroutine fcn_jac_poly (x, fx, fpx)
    !*  Returns vector-valued function and Jacobian for "polynomial" equation
    !   system.
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous, optional :: fx
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: fpx

    if (present(fx)) then
        fx(1) = 2.0 * x(1) ** 2.0 - x(2)**2.0 + x(3)
        fx(2) = x(1)*x(2)*x(3) - x(4) + 5.0
        fx(3) = x(3)**2.0 - 3.0 * x(4)**2.0
        fx(4) = x(1) + x(3) * x(4) - 10.0
    end if

    if (present(fpx)) then
        fpx(1,:) = [4.0 * x(1), -2.0*x(2), 1.0d0, 0.0d0]
        fpx(2,:) = [x(2)*x(3), x(1)*x(3), x(1)*x(2), -1.0d0]
        fpx(3,:) = [0.0d0, 0.0d0, 2.0*x(3), -6.0*x(4)]
        fpx(4,:) = [1.0d0, 0.0d0, x(4), x(3)]
    end if

end subroutine


subroutine fcn_poly (x, fx)
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx

    call fcn_jac_poly (x, fx)
end subroutine

end