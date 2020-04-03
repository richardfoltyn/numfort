

subroutine fss_deriv (fcn, x, fpx, fx, eps, reps)
    !*  FSS_DERIV numerically differentiates a function f:R->R
    !   and returns its derivative
    procedure (fss) :: fcn
    real (PREC), intent(in) :: x
        !*  Point X at which derivative should be evaluated
    real (PREC), intent(out) :: fpx
        !*  Function derivative at point X
    real (PREC), intent(in), optional :: fx
        !*  If present, contains the function value evaluated at X so that
        !   it need not be computed inside the routine.
    real (PREC), intent(in), optional :: eps
        !*  If present, contains the step size used to compute the numerical
        !   derivatives, ie. the derivative of f(X)
        !   wrt. X is obtained as (f(X + eps) - f(X)) / eps.
        !   (default: square root of machine epsilon)
        !   Note: If both EPS and REPS are present, REPS is ignored.
    real (PREC), intent(in), optional :: reps
        !   If present, contains the relative step size used to compute the
        !   numerical derivatives, ie. the derivative of f(X) wrt. X
        !   is obtained as
        !       (f(X + h) - f(X) / h
        !   where h = reps * abs(X).
        !   Note: If both EPS and REPS are present, REPS is ignored.

    real (PREC) :: fx_h, lfx, h

    if (present(fx)) then
        lfx = fx
    else
        call fcn (x, lfx)
    end if

    h = get_step_size (x, eps, reps)

    call fcn (x + h, fx_h)
    fpx = (fx_h - lfx) / h
end subroutine



subroutine fss_deriv_args (fcn, x, args, fpx, fx, eps, reps)
    !*  FSS_DERIV_ARGS numerically differentiates a function f:R->R
    !   and returns its derivative
    procedure (fss_args) :: fcn
    real (PREC), intent(in) :: x
        !*  Point X at which derivative should be evaluated
    class (args_data), intent(inout) :: args
    real (PREC), intent(out) :: fpx
        !*  Function derivative at point X
    real (PREC), intent(in), optional :: fx
        !*  If present, contains the function value evaluated at X so that
        !   it need not be computed inside the routine.
    real (PREC), intent(in), optional :: eps
        !*  If present, contains the step size used to compute the numerical
        !   derivatives, ie. the derivative of f(X)
        !   wrt. X is obtained as (f(X + eps) - f(X)) / eps.
        !   (default: square root of machine epsilon)
        !   Note: If both EPS and REPS are present, REPS is ignored.
    real (PREC), intent(in), optional :: reps
        !   If present, contains the relative step size used to compute the
        !   numerical derivatives, ie. the derivative of f(X) wrt. X
        !   is obtained as
        !       (f(X + h) - f(X) / h
        !   where h = reps * abs(X).
        !   Note: If both EPS and REPS are present, REPS is ignored.

    real (PREC) :: fx_h, lfx, h

    if (present(fx)) then
        lfx = fx
    else
        call fcn (x, args, lfx)
    end if

    h = get_step_size (x, eps, reps)

    call fcn (x + h, args, fx_h)
    fpx = (fx_h - lfx) / h

end subroutine



! ------------------------------------------------------------------------------
! Numeric differentiation for functions that map vectors into scalars

subroutine fvs_deriv (fcn, x, fpx, fx, eps, reps)
    !*  FVS_DERIV numerically differentiates a function f:R^n->R^m
    !   and returns its m-by-n Jacobian.
    procedure (fvs_fcn) :: fcn
    real (PREC), intent(in), dimension(:), contiguous :: x
        !*  Point X at which gradient should be evaluated
    real (PREC), intent(out), dimension(:) :: fpx
        !*  Function gradient
    real (PREC), intent(in), optional :: fx
        !*  If present, contains the function value evaluated at X so that
        !   it need not be computed inside the routine.
    real (PREC), intent(in), optional :: eps
        !*  If present, contains the step size used to compute the numerical
        !   derivatives, ie. the derivative of f(X)
        !   wrt. X(j) is obtained as (f(X + e_j*eps) - f(X)) / eps.
        !   where e_j is the unit vector with e(j) = 1.
        !   (default: square root of machine epsilon)
        !   Note: If both EPS and REPS are present, REPS is ignored.
    real (PREC), intent(in), optional :: reps
        !   If present, contains the relative step size used to compute the
        !   numerical derivatives, ie. the derivative of f(X) wrt. X(j)
        !   is obtained as
        !       (f(X + e_j*h) - f(X) / h
        !   where h = reps * abs(X(j)) and e_j is the unit vector with e(j) = 1.
        !   Note: If both EPS and REPS are present, REPS is ignored.

    real (PREC), dimension(:), allocatable :: x_h
    real (PREC) :: fx_h, lfx, h
    integer :: j , n

    n = size(x)

    allocate (x_h(n))

    if (present(fx)) then
        lfx = fx
    else
        call fcn (x, lfx)
    end if

    do j = 1, n
        x_h(:) = x
        h = get_step_size (x(j), eps, reps)
        x_h(j) = x_h(j) + h

        call fcn (x_h, fx_h)

        fpx(j) = (fx_h - lfx) / h
    end do

    deallocate (x_h)

end subroutine



subroutine fvs_args_deriv (fcn, x, args, fpx, fx, eps, reps)
    !*  FVS_DERIV numerically differentiates a function f:R^n->R^m
    !   and returns its m-by-n Jacobian.
    procedure (fvs_fcn_args) :: fcn
    real (PREC), intent(in), dimension(:), contiguous :: x
        !*  Point X at which gradient should be evaluated
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), dimension(:) :: fpx
        !*  Function gradient
    real (PREC), intent(in), optional :: fx
        !*  If present, contains the function value evaluated at X so that
        !   it need not be computed inside the routine.
    real (PREC), intent(in), optional :: eps
        !*  If present, contains the step size used to compute the numerical
        !   derivatives, ie. the derivative of f(X)
        !   wrt. X(j) is obtained as (f(X + e_j*eps) - f(X)) / eps.
        !   where e_j is the unit vector with e(j) = 1.
        !   (default: square root of machine epsilon)
        !   Note: If both EPS and REPS are present, REPS is ignored.
    real (PREC), intent(in), optional :: reps
        !   If present, contains the relative step size used to compute the
        !   numerical derivatives, ie. the derivative of f(X) wrt. X(j)
        !   is obtained as
        !       (f(X + e_j*h) - f(X) / h
        !   where h = reps * abs(X(j)) and e_j is the unit vector with e(j) = 1.
        !   Note: If both EPS and REPS are present, REPS is ignored.

    real (PREC), dimension(:), allocatable :: x_h
    real (PREC) :: fx_h, lfx, h
    integer :: j , n

    n = size(x)

    allocate (x_h(n))

    if (present(fx)) then
        lfx = fx
    else
        call fcn (x, args, lfx)
    end if

    do j = 1, n
        x_h(:) = x
        h = get_step_size (x(j), eps, reps)
        x_h(j) = x_h(j) + h

        call fcn (x_h, args, fx_h)

        fpx(j) = (fx_h - lfx) / h
    end do

    deallocate (x_h)

end subroutine



! ------------------------------------------------------------------------------
! Numeric differentiation for functions that map vectors into vectors

subroutine fvv_deriv (fcn, x, fpx, fx, eps, reps)
    !*  JACOBIAN numerically differentiates a function f:R^n->R^m
    !   and returns its m-by-n Jacobian.
    procedure (fvv_fcn) :: fcn
    real (PREC), intent(in), dimension(:), contiguous :: x
        !*  Point X at this Jacobian should be evaluated
    real (PREC), intent(out), dimension(:,:) :: fpx
        !*  Contains m-by-n Jacobian on exit where N=SIZE(X) and M
        !   is the dimension of the function's range.
    real (PREC), intent(in), dimension(:), optional :: fx
        !*  If present, contains the function value evaluated at X so that
        !   it need not be computed inside the routine.
    real (PREC), intent(in), optional :: eps
        !*  If present, contains the step size used to compute the numerical
        !   derivatives, ie. the derivative of the i-th element of f(X)
        !   wrt. X(j) is obtained as (f_i(X + e_j*eps) - f_i(X)) / eps.
        !   where e_j is the unit vector with e(j) = 1.
        !   (default: square root of machine epsilon)
        !   Note: If both EPS and REPS are present, REPS is ignored.
    real (PREC), intent(in), optional :: reps
        !   If present, contains the relative step size used to compute the
        !   numerical derivatives, ie. the derivative of the i-th element
        !   of f(X) wrt. X(j) is obtained as
        !       (f_i(X + e_j*h) - f_i(X) / h
        !   where h = reps * abs(X(j)) and e_j is the unit vector with e(j) = 1.
        !   Note: If both EPS and REPS are present, REPS is ignored.

    real (PREC) :: h
    real (PREC), dimension(:), allocatable :: x_h, fx_h, lfx
    integer :: i, j, m, n

    ! Dimension of function range
    m = size(fpx,1)
    n = size(x)

    allocate (x_h(n), fx_h(m))

    if (present(fx)) then
        allocate (lfx(m), source=fx)
    else
        allocate (lfx(m))
        call fcn (x, lfx)
    end if

    do j = 1, n
        x_h(:) = x
        h = get_step_size (x(j), eps, reps)
        x_h(j) = x_h(j) + h

        call fcn (x_h, fx_h)

        do i = 1, m
            fpx(i,j) = (fx_h(i) - lfx(i)) / h
        end do
    end do

    deallocate (x_h, fx_h, lfx)

end subroutine


subroutine fvv_args_deriv (fcn, x, args, fpx, fx, eps, reps)
    !*  JACOBIAN numerically differentiates a function f:R^n->R^m
    !   and returns its m-by-n Jacobian.
    procedure (fvv_fcn_args) :: fcn
    real (PREC), intent(in), dimension(:), contiguous :: x
        !*  Point X at this Jacobian should be evaluated
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), dimension(:,:) :: fpx
        !*  Contains m-by-n Jacobian on exit where N=SIZE(X) and M
        !   is the dimension of the function's range.
    real (PREC), intent(in), dimension(:), optional :: fx
        !*  If present, contains the function value evaluated at X so that
        !   it need not be computed inside the routine.
    real (PREC), intent(in), optional :: eps
        !*  If present, contains the step size used to compute the numerical
        !   derivatives, ie. the derivative of the i-th element of f(X)
        !   wrt. X(j) is obtained as (f_i(X + e_j*eps) - f_i(X)) / eps.
        !   where e_j is the unit vector with e(j) = 1.
        !   (default: square root of machine epsilon)
        !   Note: If both EPS and REPS are present, REPS is ignored.
    real (PREC), intent(in), optional :: reps
        !   If present, contains the relative step size used to compute the
        !   numerical derivatives, ie. the derivative of the i-th element
        !   of f(X) wrt. X(j) is obtained as
        !       (f_i(X + e_j*h) - f_i(X) / h
        !   where h = reps * abs(X(j)) and e_j is the unit vector with e(j) = 1.
        !   Note: If both EPS and REPS are present, REPS is ignored.

    real (PREC) :: h
    real (PREC), dimension(:), allocatable :: x_h, fx_h, lfx
    integer :: i, j, m, n

    ! Dimension of function range
    m = size(fpx,1)
    n = size(x)

    allocate (x_h(n), fx_h(m))

    if (present(fx)) then
        allocate (lfx(m), source=fx)
    else
        allocate (lfx(m))
        call fcn (x, args, lfx)
    end if

    do j = 1, n
        x_h(:) = x
        h = get_step_size (x(j), eps, reps)
        x_h(j) = x_h(j) + h

        call fcn (x_h, args, fx_h)

        do i = 1, m
            fpx(i,j) = (fx_h(i) - lfx(i)) / h
        end do
    end do

    deallocate (x_h, fx_h, lfx)

end subroutine



pure function get_step_size (x, eps, reps) result(res)
    !*  GET_STEP_SIZE determines the forward-difference step size
    !   depending on (presence of) user-provided arguments.
    real (PREC), intent(in) :: x
    real (PREC), intent(in), optional :: eps
    real (PREC), intent(in), optional :: reps
    real (PREC) :: res

    real (PREC) :: h

    if (present(reps) .and. .not. present(eps)) then
        h = reps * abs(x)
        ! Take care of cases when x = 0.0
        if (h == 0.0_PREC) h = reps
    else if (present(eps)) then
        h = eps
    else
        ! Default step size
        h = sqrt(epsilon(0.0_PREC))
    end if

    res = h
end function
