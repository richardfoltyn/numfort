

program test_optimize_solver_map

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    use numfort_arrays
    use numfort_common
    use numfort_common_testing
    use numfort_optimize
    use numfort_stats, only: set_seed

    use fcore_strings
    use fcore_testing

    implicit none

    integer, parameter :: PREC = real64

    real (PREC) :: ninf
    real (PREC) :: inf

    ninf = ieee_value (0.0_PREC, IEEE_NEGATIVE_INF)
    inf = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)

    call test_all ()


    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ('optimize::solver_map unit tests')

    call test_linear (tests)
    call test_exp (tests)

    call tests%print ()

end subroutine



subroutine test_linear (tests)
    class (test_suite) :: tests

    class (test_case), pointer :: tc

        ! Note: y0 parameter is irrelevant for linear mappings
    real (PREC), parameter :: dx_all(*) = [1.0d-6, 0.1d0, 1.0d0, 2.34d0]
    real (PREC), parameter :: dy_all(*) = [-1.234d0, -0.123d0, 0.1d0, 1.0d0, 2.234d0]

    type (solver_map), dimension(:), allocatable :: maps
    type (solver_map) :: map
    real (PREC), dimension(:), allocatable :: x, y, x1, y1
    real (PREC), dimension(:), allocatable :: jac_diag, jac_inv_diag
    real (PREC), dimension(:,:), allocatable :: jac, jac_inv, eye, xx
    real (PREC) :: dx, dy
    integer :: i, j, k, nx, ny, n
    type (status_t) :: status
    type (str) :: msg

    tc => tests%add_test ('Linear mappings')

    ! Input checks with invalid arguments
    status = NF_STATUS_UNDEFINED
    call solver_map_init (map, lb=ninf, dy=0.0_PREC, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "INIT called with invalid dy=0")

    status = NF_STATUS_UNDEFINED
    call solver_map_init (map, lb=ninf, dx=0.0_PREC, dy=1.0_PREC, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "INIT called with invalid dx=0")

    nx = size(dx_all)
    ny = size(dy_all)
    n = nx * ny

    allocate (maps(n))

    do i = 1, size(dx_all)
        do j = 1, size(dy_all)
            k = (i-1) * ny + j

            dx = dx_all(i)
            dy = dy_all(j)

            status = NF_STATUS_UNDEFINED
            call solver_map_init (maps(k), lb=ninf, dy=dy, dx=dx, status=status)
            msg = 'INIT called with dx=' // str(dx, 'f6.3') // '; dy=' &
                // str(dy, 'f6.3')
            call tc%assert_true (status == NF_STATUS_OK, msg)
        end do
    end do

    ! Test mappings and its inverse
    allocate (x(n), y(n), x1(n), y1(n))

    call set_seed (1234)
    call random_number (x)
    x(:) = x * 10.0

    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps, x, y, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'EVAL called at values X')

    status = NF_STATUS_UNDEFINED
    call solver_map_eval_inverse (maps, y, x1, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'EVAL_INVERSE called at values f(X)')

    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (x1, x1, atol=1.0e-10_PREC, rtol=0.0_PREC), &
        'Checking X == f^-1(f(X))')

    status = NF_STATUS_UNDEFINED
    call random_number (y)
    call solver_map_eval_inverse (maps, y, x, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'EVAL_INVERSE called at values Y')

    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps, x, y1, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (y1, y, atol=1.0e-10_PREC, rtol=0.0_PREC), &
        'Checking Y == f(f^-1(Y))')

    ! Test diagonal Jacobians
    allocate (jac_diag(n), jac_inv_diag(n))
    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps, x, y, jac_diag, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'EVAL called with 1-d JAC argument')

    status = NF_STATUS_UNDEFINED
    call solver_map_eval_inverse (maps, y, x1, jac_inv_diag, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'EVAL_INVERSE called with 1-d JAC argument')

    call tc%assert_true (all_close (jac_inv_diag, 1.0_PREC/jac_diag), &
        'Checking inverse of Jacobian == Jacobian of inverse')

    deallocate (jac_diag, jac_inv_diag)

    ! Test Jacobians
    allocate (jac(n,n), jac_inv(n,n), eye(n,n), xx(n,n))
    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps, x, y, jac, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'EVAL called with 2-d JAC argument')

    status = NF_STATUS_UNDEFINED
    call solver_map_eval_inverse (maps, y, x1, jac_inv, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'EVAL_INVERSE called with 2-d JAC argument')

    xx = matmul (jac, jac_inv)
    call identity (eye)
    call tc%assert_true (all_close (xx, eye, atol=1.0e-10_PREC, rtol=0.0_PREC), &
        'Checking inverse of Jacobian == Jacobian of inverse')

    deallocate (jac, jac_inv)

end subroutine



subroutine test_exp (tests)
    class (test_suite) :: tests

    class (test_case), pointer :: tc

        ! Note: y0 parameter is irrelevant for linear mappings
    real (PREC), parameter :: dx_all(*) = [1.0d-6, 0.1d0, 1.34d0]
    real (PREC), parameter :: dy_all(*) = [-1.234d0, -0.123d0, 0.1d0]
    real (PREC), parameter :: y0_all(*) = [-1.234d0, 0.0d0, 0.123d0]
    real (PREC), parameter :: lb_all(*) = [-0.123d0, 0.0d0, 1.234d0]

    type (solver_map), dimension(:), allocatable :: maps
    type (solver_map) :: map
    real (PREC), dimension(:), allocatable :: x, y, x1, y1
    real (PREC), dimension(:), allocatable :: jac_diag, jac_inv_diag
    real (PREC), dimension(:,:), allocatable :: jac, jac_inv, eye, xx
    real (PREC) :: dx, dy, lb, y0
    integer :: i, j, k, l, n, nmax
    type (status_t) :: status
    type (str) :: msg

    tc => tests%add_test ('Linear mappings')

    ! Input checks with invalid arguments
    status = NF_STATUS_UNDEFINED
    call solver_map_init (map, lb=0.0_PREC, dy=0.0_PREC, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "INIT called with invalid dy=0")

    status = NF_STATUS_UNDEFINED
    call solver_map_init (map, lb=0.0_PREC, dx=0.0_PREC, dy=1.0_PREC, status=status)
    call tc%assert_true (status == NF_STATUS_INVALID_ARG, &
        "INIT called with invalid dx=0")

    nmax = size(dx_all) * size(dy_all) * size(y0_all) * size(lb_all)

    allocate (maps(nmax))

    n = 0
    do i = 1, size(dx_all)
        dx = dx_all(i)

        do j = 1, size(dy_all)
            dy = dy_all(j)

            do k = 1, size(y0_all)
                y0 = y0_all(k)

                do l = 1, size(lb_all)
                    lb = lb_all(l)

                    if (y0 <= lb .or. (y0 + dy) <= lb) cycle

                    n = n + 1

                    status = NF_STATUS_UNDEFINED
                    call solver_map_init (maps(n), lb=lb, y0=y0, dy=dy, dx=dx, &
                        status=status)
                    msg = 'INIT called with lb=' // str(lb, 'f6.3') &
                        // '; y0=' // str(y0, 'f6.3') //  '; dx=' &
                        // str(dx, 'f6.3') // '; dy=' // str(dy, 'f6.3')
                    call tc%assert_true (status == NF_STATUS_OK, msg)
                end do
            end do
        end do
    end do

    ! Test mappings and its inverse
    allocate (x(n), y(n), x1(n), y1(n))

    call set_seed (1234)
    call random_number (x)
    x(:) = x * 10.0

    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps(1:n), x, y, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'EVAL called at values X')

    status = NF_STATUS_UNDEFINED
    call solver_map_eval_inverse (maps(1:n), y, x1, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'EVAL_INVERSE called at values f(X)')

    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (x1, x1, atol=1.0e-10_PREC, rtol=0.0_PREC), &
        'Checking X == f^-1(f(X))')

    status = NF_STATUS_UNDEFINED
    call random_number (y)
    call solver_map_eval_inverse (maps(1:n), y, x, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'EVAL_INVERSE called at values Y')

    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps(1:n), x, y1, status=status)
    call tc%assert_true (status == NF_STATUS_OK .and. &
        all_close (y1, y, atol=1.0e-10_PREC, rtol=0.0_PREC), &
        'Checking Y == f(f^-1(Y))')

    ! Test diagonal Jacobians
    allocate (jac_diag(n), jac_inv_diag(n))
    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps(1:n), x, y, jac_diag, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'EVAL called with 1-d JAC argument')

    status = NF_STATUS_UNDEFINED
    call solver_map_eval_inverse (maps(1:n), y, x1, jac_inv_diag, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'EVAL_INVERSE called with 1-d JAC argument')

    call tc%assert_true (all_close (jac_inv_diag, 1.0_PREC/jac_diag), &
        'Checking inverse of Jacobian == Jacobian of inverse')

    deallocate (jac_diag, jac_inv_diag)

    ! Test Jacobians
    allocate (jac(n,n), jac_inv(n,n), eye(n,n), xx(n,n))
    status = NF_STATUS_UNDEFINED
    call solver_map_eval (maps(1:n), x, y, jac, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'EVAL called with 2-d JAC argument')

    status = NF_STATUS_UNDEFINED
    call solver_map_eval_inverse (maps(1:n), y, x1, jac_inv, status=status)
    call tc%assert_true (status == NF_STATUS_OK, &
        'EVAL_INVERSE called with 2-d JAC argument')

    xx = matmul (jac, jac_inv)
    call identity (eye)
    call tc%assert_true (all_close (xx, eye, atol=1.0e-10_PREC, rtol=0.0_PREC), &
        'Checking inverse of Jacobian == Jacobian of inverse')

    deallocate (jac, jac_inv)

end subroutine


end
