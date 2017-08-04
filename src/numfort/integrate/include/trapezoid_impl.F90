
!------------------------------------------------------------------------------- 
! INPUT CHECKING

pure subroutine __APPEND(check_inputs,__PREC) (a, b, n, status)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in) :: a, b
    integer, intent(in) :: n
    type (status_t), intent(out) :: status

    status = NF_STATUS_INVALID_ARG
    
    if (n <= 1) goto 100

    status = NF_STATUS_OK

100 continue
end subroutine


pure subroutine __APPEND(check_inputs_array,__PREC) (fx, x, dx, status)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:) :: fx
    real (PREC), intent(in), dimension(:), optional :: x
    real (PREC), intent(in), optional :: dx
    type (status_t), intent(out) :: status
    
    status = NF_STATUS_INVALID_ARG

    ! We need at least one interval, ie two points
    if (size(fx) < 2) goto 100

    ! At least one of x or dx must be present
    if (.not. present(x) .and. .not. present(dx)) goto 100

    if (present(x)) then
        if (size(fx) /= size(x)) goto 100
    end if

    if (present(dx)) then
        if (dx <= 0.0_PREC) goto 100
    end if

    status = NF_STATUS_OK

100 continue
end subroutine

!-------------------------------------------------------------------------------
! TRAPEZOID RULE

subroutine __APPEND(trapezoid,__PREC) (fcn, a, b, n, res, status)
    !*  Integrates a function using the composite trapezoid rule
    !   using equal-length intervals.

    integer, parameter :: PREC = __PREC
   
    abstract interface
        function f_integrand (x) result(fx)
            import __PREC
            integer, parameter :: PREC = __PREC
            real (PREC), intent(in) :: x
            real (PREC) :: fx
        end function
    end interface

    procedure (f_integrand) :: fcn
        !*  Function to be integrated
    real (PREC), intent(in) :: a
        !*  Lower bound of integration interval
    real (PREC), intent(in) :: b
        !*  Upper bound of integration interval (requires b >= a)  
    integer, intent(in) :: n
        !*  Number of (equidistant) points at which to evaluate 
        !   function. Composite trapezoid rule is thus applied to
        !   (n-1) sub-intervals of [a,b].
    real (PREC), intent(out) :: res
        !*  Integral of fcn
    type (status_t), intent(out), optional :: status
        !*  If present, contains status code on exit.

    type (status_t) :: lstatus
    real (PREC) :: dx
    integer :: i

    lstatus = NF_STATUS_OK
    call check_inputs (a, b, n, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100
    
    ! Handles integration over [a,b] with a > b automatically:
    ! In that case dx < 0, so it's equivalent to integrating
    ! over [b,a] scaled by (-1).
    dx = (b - a) / (n-1)
    res = 0.0_PREC

    if (a /= b) then
        res = fcn(a)
        do i = 1, n - 2
            res = res + 2.0 * fcn(a + dx * i)
        end do
        res = res + fcn(b)
        res = res * dx / 2.0
    end if
    
100 continue
    if (present(status)) status = lstatus

end subroutine


subroutine __APPEND(trapezoid_array,__PREC) (fx, res, x, dx, status)
    !*  Integrates a function defined by a set of points (x, fx) 
    !   using the composite trapezoid rule. Alternatively, if the
    !   intervals between consecutive elements of fx are identical,
    !   the integral is computing using the interval length dx.

    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:) :: fx
        !*  Array definint the function to be integrated at 
        !   discrete points.
    real (PREC), intent(out) :: res
        !*  Integral of f(x)
    real (PREC), intent(in), dimension(:), optional :: x
        !*  Array of values x at which the corresponding elements 
        !   in fx are evaluated. Must be omitted if dx is specified.
    real (PREC), intent(in), optional :: dx
        !*  Uniform spacing dx between values of x at which f(x)
        !   is evaluated. Ignored if x is present.
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: i, n

    lstatus = NF_STATUS_OK
    call check_inputs (fx, x, dx, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    n = size(fx)

    if (present(x)) then
        res = 0.0_PREC
        do i = 1, n-1
            res = res + (fx(i) + fx(i+1)) * (x(i+1)-x(i))
        end do
    else
        ! Uniform interval length, no need to 
        ! calculate dx for each interval.
        res = fx(1)
        do i = 2, n - 1
            res = res + 2.0 * fx(i)
        end do
        res = res + fx(n)
        res = res * dx
    end if

    res = res / 2.0
    
100 continue
    if (present(status)) status = lstatus

end subroutine



