


pure subroutine solver_map_init (self, lb, ub, y0, dy, dx, status)
    !*  MAPPING_INIT initializes the MAPPING object representing the function
    !   y = f(x) such that it maps the X on the real line into one of the
    !   intervals
    !       (1) [-inf, inf]
    !       (2) [y_lb, inf)
    !       (3) [-inf, y_ub] (not yet implemented)
    !       (4) [y_lb, y_ub]
    !
    !   Additionally, a scaling parameter can be computed that ensures
    !   a unit change in x results in a "reasonable" step size in y
    !   around x0 = f^{-1}(y0)
    !
    !   Linear map
    !   ==========
    !   For the mapping of type (1) we apply the transformation
    !       y = x/s
    !   and thus the scale parameter is determined as
    !       s = dx/dy
    !
    !   Exponential map
    !   ===============
    !   For mappings of type (2) we apply the transformation
    !       y = exp(x/s) + y_lb
    !   To pin down the scale parameter s, we use
    !       dy = exp((x0+dx)/s) - exp(x0/s)
    !          = exp(x0/s)*exp(dx/s) - exp(x0/s)
    !   where
    !       x0 = s * log(y0 - y_lb)
    !   and therefore
    !       dy = (y0 - y_lb) * exp(dx/s) - (y0 - y_lb)
    !       dx/s = log(y0 - y_lb + dy) - log(y0 - y_lb)
    !
    !   Logistic map
    !   ============
    !   For bounded mappings of type (4), we apply the transformation
    !       y = y_lb + F(x,s)(y_u - y_lb)
    !   where [y_lb, y_ub] are the bounds to be enforced for y and
    !   F(x,s) is the logistic CDF with the yet-to-be determined scale
    !   coefficient s,
    !       F(x,s) = 1/(1 + exp(-x/s))
    !
    !   We require that at y0
    !       dy = f(x0+dx) - f(x0)
    !   where dy is a "reasonable" step size for y at x0.
    !   Thus
    !       dy  = y_lb + F(x0+dx,s)(y_ub - y_lb) - y_lb - F(x0,s)(y_ub - y_lb)
    !           = (F(x0+dx,s)-F(x0,s))(y_ub - y_lb)
    !
    !   Additionally, the initial point y0 is given by
    !       (y0 - y_lb)/(y_ub - y_lb) = F(x0,s)
    !   so from the definition of the logistic CDF we get
    !       x0 = - s * log(1/z - 1)
    !   where z = (y0 - y_lb)/(y_ub - y_lb)
    !
    !   Combining these expressions, we find that s must satisfy
    !       dy/(y_ub - y_lb) = 1/(1 + (1/z-1)exp(-dx/s)) - z
    !   which implies that
    !       dx/s = - [log(y_ub - y0 - dy) - log(y0 - y_lb + dy)] +
    !               [log(y_ub - y0) - log(y0 - y_lb)]
    type (solver_map), intent(inout) :: self
    real (PREC), intent(in) :: lb
    real (PREC), intent(in), optional :: ub
    real (PREC), intent(in), optional :: y0
    real (PREC), intent(in), optional :: dx
    real (PREC), intent(in), optional :: dy
    type (status_t), intent(out), optional :: status

    real (PREC) :: inf, ninf
    integer :: transform
    real (PREC) :: ylb, yub, a1, a2, s, ldx
    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK

    inf = ieee_value (0.0_PREC, IEEE_POSITIVE_INF)
    ninf = ieee_value (0.0_PREC, IEEE_NEGATIVE_INF)

    ! Default range boundaries
    ylb = ninf
    yub = inf
    ! Default scale parameter
    s = 1.0_PREC
    ! Default transformation
    transform = TRANSFORM_LINEAR
    ldx = 1.0_PREC

    ylb = lb
    if (present(ub)) yub = ub
    if (present(dx)) ldx = dx

    ! Input validation
    call check_nonzero (0.0_PREC, dy, 'dy', status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_positive (0.0_PREC, dx, 'dx', status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (ieee_class (ylb) == IEEE_NEGATIVE_INF &
            .and. ieee_class (yub) == IEEE_POSITIVE_INF) then
        ! Maps from R -> R, possibly with rescaling
        transform = TRANSFORM_LINEAR
        if (present(dy)) then
            s = ldx / dy
        end if

    else if (ylb > ninf .and. ieee_class (yub) == IEEE_POSITIVE_INF) then
        ! Maps R -> [lb, inf); use exp-like transformation
        if (present(y0) .and. present(dy)) then
            ! Ensure that y0 in (y_lb, inf)
            if (y0 <= ylb .or. (y0 + dy) <= ylb) then
                lstatus = NF_STATUS_INVALID_ARG
                lstatus = lstatus + NF_STATUS_DOMAIN_ERROR
                goto 100
            end if

            s = ldx / (log(y0 + dy - ylb) - log(y0 - ylb))
        end if

        transform = TRANSFORM_EXP

    else if (ylb > ninf .and. yub < inf) then
        ! Maps R -> [lb,ub]; use logit CDF transformation
        if (present(y0) .and. present(dy)) then
            ! Ensure that y0 in (y_lb, y_ub)
            if (y0 <= ylb .or. yub <= y0 .or. (y0 + dy) <= ylb .or. &
                    (y0 + dy) >= yub) then
                lstatus = NF_STATUS_INVALID_ARG
                lstatus = lstatus + NF_STATUS_DOMAIN_ERROR
                goto 100
            end if

            a1 = log(y0 + dy - ylb) - log(yub - y0 - dy)
            a2 = log(yub - y0) - log(y0 - ylb)

            s = ldx/(a1 + a2)

            if (.not. ieee_is_finite (s)) then
                lstatus = NF_STATUS_INVALID_STATE
                goto 100
            end if
        end if

        transform = TRANSFORM_LOGISTIC
    end if

    self%transform = transform
    self%scale = s
    self%lb = ylb
    self%ub = yub

100 continue
    if (present(status)) status = lstatus

end subroutine



pure subroutine solver_map_eval_scalar (self, x, y, jac)
    type (solver_map), intent(in) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    real (PREC), intent(out), optional :: jac

    real (PREC) :: s, lb, ub, dydx, z

    s = real(self%scale, PREC)
    lb = real(self%lb, PREC)
    ub = real(self%ub, PREC)

    dydx = 0.0_PREC

    select case (self%transform)
    case (TRANSFORM_LINEAR)
        y = x / s
        dydx = 1.0_PREC / s
    case (TRANSFORM_EXP)
        y = exp(x/s) + lb
        dydx = exp(x/s) / s
    case (TRANSFORM_NEG_EXP)
        continue
    case (TRANSFORM_LOGISTIC)
        z = 1.0_PREC/(1.0_PREC + exp(-x/s))
        y = lb + z * (ub-lb)
        ! Prevent indeterminate expression for "large" X that lead
        ! to z = 0.0
        ! Otherwise leave dydx = 0 as initialized
        if (z > 0.0_PREC) then
            dydx = (ub-lb)*z**2.0_PREC * exp(-x/s) / s
        end if
    end select

    if (present(jac)) then
        jac = dydx
    end if
end subroutine



subroutine solver_map_eval_diag (map, x, y, jac, status)
    type (solver_map), intent(in), dimension(:) :: map
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(out), dimension(:) :: y
    real (PREC), intent(out), dimension(:), optional :: jac
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: i, n

    lstatus = NF_STATUS_OK

    n = size(map)

    if (size(y) /= n .or. size(x) /= n) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (present(jac)) then
        if (size(jac) /= n) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    if (present(jac)) then
        do i = 1, n
            call solver_map_eval (map(i), x(i), y(i), jac(i))
        end do
    else
        do i = 1, n
            call solver_map_eval (map(i), x(i), y(i))
        end do
    end if

100 continue
    if (present(status)) status = lstatus

end subroutine



subroutine solver_map_eval_matrix (map, x, y, jac, status)
    type (solver_map), intent(in), dimension(:) :: map
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(out), dimension(:) :: y
    real (PREC), intent(out), dimension(:,:) :: jac
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: i, n

    lstatus = NF_STATUS_OK

    n = size(map)

    if (size(y) /= n .or. size(x) /= n) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (size(jac,1) /= n .or. size(jac,1) /= n) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    jac = 0.0_PREC
    do i = 1, n
        call solver_map_eval (map(i), x(i), y(i), jac(i,i))
    end do

100 continue
    if (present(status)) status = lstatus

end subroutine



pure subroutine solver_map_eval_inverse_scalar (self, y, x, jac, status)
    type (solver_map), intent(in) :: self
    real (PREC), intent(in) :: y
    real (PREC), intent(out) :: x
    real (PREC), intent(out), optional :: jac
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    real (PREC) :: dxdy, s, a, lb, ub, z

    lstatus = NF_STATUS_OK
    s = real(self%scale,PREC)
    lb = real(self%lb, PREC)
    ub = real(self%ub, PREC)

    dxdy = 0.0_PREC

    select case (self%transform)
    case (TRANSFORM_LINEAR)
        x = y * s
        dxdy = s

    case (TRANSFORM_EXP)
        if (y < lb) then
            lstatus = NF_STATUS_INVALID_ARG
            lstatus = lstatus + NF_STATUS_DOMAIN_ERROR
            goto 100
        end if

        x = s * log(y - lb)
        dxdy = s / (y - lb)

    case (TRANSFORM_NEG_EXP)
        continue

    case (TRANSFORM_LOGISTIC)
        if (y < lb .or. y > ub) then
            lstatus = NF_STATUS_INVALID_ARG
            lstatus = lstatus + NF_STATUS_DOMAIN_ERROR
            goto 100
        end if

        z = (y-lb)/(ub-lb)
        a = 1.0_PREC/z - 1.0_PREC
        x = - s * log(a)
        ! Derivative : dx/dy = -s/a * da/dz * dz/dy
        ! where
        !   da/dz = -1/z^2
        !   dz/dy = 1/(ub-lb)
        dxdy = s / a / z**2.0_PREC / (ub-lb)
    end select

    if (present(jac)) then
        jac = dxdy
    end if

100 continue
    if (present(status)) status = lstatus

end subroutine



subroutine solver_map_eval_inverse_diag (map, y, x, jac, status)
    type (solver_map), intent(in), dimension(:) :: map
    real (PREC), intent(in), dimension(:) :: y
    real (PREC), intent(out), dimension(:) :: x
    real (PREC), intent(out), dimension(:), optional :: jac
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: i, n

    lstatus = NF_STATUS_OK

    n = size(map)

    if (size(y) /= n .or. size(x) /= n) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (present(jac)) then
        if (size(jac) /= n) then
            lstatus = NF_STATUS_INVALID_ARG
            goto 100
        end if
    end if

    if (present(jac)) then
        do i = 1, n
            call solver_map_eval_inverse (map(i), y(i), x(i), jac(i), lstatus)
            if (lstatus /= NF_STATUS_OK) goto 100
        end do
    else
        do i = 1, n
            call solver_map_eval_inverse (map(i), y(i), x(i), status=lstatus)
            if (lstatus /= NF_STATUS_OK) goto 100
        end do
    end if

100 continue
    if (present(status)) status = lstatus

end subroutine



subroutine solver_map_eval_inverse_matrix (map, y, x, jac, status)
    type (solver_map), intent(in), dimension(:) :: map
    real (PREC), intent(in), dimension(:) :: y
    real (PREC), intent(out), dimension(:) :: x
    real (PREC), intent(out), dimension(:,:) :: jac
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: i, n

    lstatus = NF_STATUS_OK

    n = size(map)

    if (size(y) /= n .or. size(x) /= n) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    if (size(jac,1) /= n .or. size(jac,1) /= n) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    jac = 0.0_PREC

    do i = 1, n
        call solver_map_eval_inverse (map(i), y(i), x(i), jac(i,i), lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end do


100 continue
    if (present(status)) status = lstatus

end subroutine


