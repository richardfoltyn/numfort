


pure subroutine __APPEND(solver_map_init,__PREC) (self, lb, ub, y0, dy, dx, status)
    integer, parameter :: PREC = __PREC
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
    !   where
    !       x0 = s * log(y0 - y_lb)
    !   and therefore
    !       dx/s = log(y0 - y_lb + dy) - log(y0 - y_lb)
    !
    !   Logistic map
    !   ============
    !   For bounded mappings of type (4), we apply the transformation
    !       y = y_lb + F(x,s)(y_u - y_lb)
    !   where [y_lb, y_ub] are the bounds to be enforced for y and
    !   F(x,s) is the logistic CDF with the yet-to-be determined scale
    !   coefficient s.
    !
    !   We require that at y0
    !       dy = f(x0+dx) - f(x0)
    !   where dy is in "reasonable" step size for y at x0.
    !   Thus
    !       dy  = y_lb + F(x0+dx,s)(y_ub - y_lb) - y_lb - F(x0,s)(y_ub - y_lb)
    !           = (F(x0+dx,s)-F(x0,s))(y_ub - y_lb)
    !
    !   Additionally, the initial point y0 is given by
    !       (y0 - y_lb)/(y_ub - y_lb) = F(x0,s)
    !   so from the definition of the logistic CDF we get
    !       x0 = s * log(1/z - 1)
    !   where z = (y0 - y_lb)/(y_ub - y_lb)
    !
    !   Combining these expressions, we find that s must satisfy
    !       dy/(y_ub - y_lb) = 1/(1 + (1/z-1)exp(dx/s)) - z
    !   which implies that
    !       dx/s = [log(y_ub - y0 - dy) - log(y0 - y_lb + dy)] /
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

            a1 = log(yub - y0 - dy) - log(y0 + dy - ylb)
            a2 = log(yub - y0) - log(y0 - ylb)

            s = ldx/(a1 - a2)

            if (.not. ieee_is_finite (s)) then
                lstatus = NF_STATUS_INVALID_STATE
                goto 100
            end if

            transform = TRANSFORM_LOGISTIC
        end if
    end if

    self%transform = transform
    self%scale = real(s,PREC)
    self%lb = real(ylb,PREC)
    self%ub = real(yub,PREC)

100 continue
    if (present(status)) status = lstatus

end subroutine



pure subroutine __APPEND(solver_map_eval_scalar,__PREC) (self, x, y, jac)
    integer, parameter :: PREC = __PREC
    type (solver_map), intent(in) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    real (PREC), intent(out), optional :: jac

    real (PREC) :: s, lb, ub, dydx

    s = real(self%scale, PREC)
    lb = real(self%lb, PREC)
    ub = real(self%ub, PREC)

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
        y = lb + (ub-lb)/(1.0_PREC + exp(x/s))
        dydx = (ub-lb)/(1.0_PREC + exp(x/s))**2.0_PREC * exp(x/s) / s
    end select

    if (present(jac)) then
        jac = dydx
    end if
end subroutine



subroutine __APPEND(solver_map_eval_diag,__PREC) (map, x, y, jac, status)
    integer, parameter :: PREC = __PREC
    type (solver_map), intent(in), dimension(:) :: map
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(out), dimension(:) :: y
    real (PREC), intent(out), dimension(:), optional :: jac
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    real (PREC), dimension(:), allocatable :: ly, ljac
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

    allocate (ly(n), ljac(n))

    do i = 1, n
        call solver_map_eval (map(i), x(i), ly(i), ljac(i))
    end do

    y = ly
    if (present(jac)) jac = ljac

100 continue
    if (present(status)) status = lstatus

end subroutine



subroutine __APPEND(solver_map_eval_matrix,__PREC) (map, x, y, jac, status)
    integer, parameter :: PREC = __PREC
    type (solver_map), intent(in), dimension(:) :: map
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(out), dimension(:) :: y
    real (PREC), intent(out), dimension(:,:) :: jac
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    real (PREC), dimension(:), allocatable :: ly, ljac
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

    allocate (ly(n), ljac(n))

    do i = 1, n
        call solver_map_eval (map(i), x(i), ly(i), ljac(i))
    end do

    y = ly
    ! Create diagonal Jacobian
    call diag (ljac, jac)

100 continue
    if (present(status)) status = lstatus

end subroutine



pure subroutine __APPEND(solver_map_eval_inverse_scalar,__PREC) (self, y, x, jac, status)
    integer, parameter :: PREC = __PREC
    type (solver_map), intent(in) :: self
    real (PREC), intent(in) :: y
    real (PREC), intent(out) :: x
    real (PREC), intent(out), optional :: jac
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    real (PREC) :: dxdy, s, a, lb, ub

    lstatus = NF_STATUS_OK
    s = real(self%scale,PREC)
    lb = real(self%lb, PREC)
    ub = real(self%ub, PREC)

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

        a = (ub-lb)/(y-lb) - 1.0_PREC
        x = s * log(a)
        dxdy = - s / a * (ub-lb) / (y-lb)**2.0_PREC
    end select

    if (present(jac)) then
        jac = dxdy
    end if

100 continue
    if (present(status)) status = lstatus

end subroutine



subroutine __APPEND(solver_map_eval_inverse_diag,__PREC) (map, y, x, jac, status)
    integer, parameter :: PREC = __PREC
    type (solver_map), intent(in), dimension(:) :: map
    real (PREC), intent(in), dimension(:) :: y
    real (PREC), intent(out), dimension(:) :: x
    real (PREC), intent(out), dimension(:), optional :: jac
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    real (PREC), dimension(:), allocatable :: lx, ljac
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

    allocate (lx(n), ljac(n))

    do i = 1, n
        call solver_map_eval_inverse (map(i), y(i), lx(i), ljac(i), lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end do

    x = lx
    if (present(jac)) jac = ljac

100 continue
    if (present(status)) status = lstatus

end subroutine



subroutine __APPEND(solver_map_eval_inverse_matrix,__PREC) (map, y, x, jac, status)
    integer, parameter :: PREC = __PREC
    type (solver_map), intent(in), dimension(:) :: map
    real (PREC), intent(in), dimension(:) :: y
    real (PREC), intent(out), dimension(:) :: x
    real (PREC), intent(out), dimension(:,:) :: jac
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    real (PREC), dimension(:), allocatable :: lx, ljac
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

    allocate (lx(n), ljac(n))

    do i = 1, n
        call solver_map_eval_inverse (map(i), y(i), lx(i), ljac(i), lstatus)
        if (lstatus /= NF_STATUS_OK) goto 100
    end do

    x = lx
    ! Create diagonal Jacobian
    call diag (ljac, jac)

100 continue
    if (present(status)) status = lstatus

end subroutine


