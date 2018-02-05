

subroutine __APPEND(jacobian,__PREC) (fcn, x, fpx, fx, eps)
    !*  JACOBIAN numerically differentiates a function f:R^n->R^m
    !   and returns its m-by-n Jacobian.
    integer, parameter :: PREC = __PREC
    procedure (__APPEND(fvv,__PREC)) :: fcn
    real (PREC), intent(in), dimension(:) :: x
        !*  Point X at this Jacobian should be evaluated
    real (PREC), intent(out), dimension(:,:) :: fpx
        !*  Contains m-by-n Jacobian on exit where N=SIZE(X) and M
        !   is the dimension of the function's range.
    real (PREC), intent(in), dimension(:), optional :: fx
        !*  If present, contains the function value evaluated at X so that
        !   it need not be computed inside the routine.
    real (PREC), intent(in), optional :: eps
        !*  If present, contains the step size used to compute the numerical
        !   derivatives, ie. the derivative of the j-th element of f(X)
        !   wrt. to X(i) is obtained as (f_j(X(i) + eps) - f_j(X(i))) / eps.
        !   (default: square root of machine epsilon)

    real (PREC) :: leps
    real (PREC), dimension(:), allocatable :: x_eps, fx_eps, lfx
    integer :: i, j, m, n

    leps = sqrt(epsilon(0.0_PREC))
    if (present(eps)) leps = eps

    ! Dimension of function range
    m = size(fpx,1)
    n = size(x)

    allocate (x_eps(n), fx_eps(m))

    if (present(fx)) then
        allocate (lfx(m), source=fx)
    else
        allocate (lfx(m))
        call fcn (x, lfx)
    end if

    do j = 1, n
        x_eps(:) = x
        x_eps(j) = x_eps(j) + leps

        call fcn (x_eps, fx_eps)

        do i = 1, m
            fpx(i,j) = (fx_eps(i) - lfx(i)) / leps
        end do
    end do

    deallocate (x_eps, fx_eps, lfx)

end subroutine
