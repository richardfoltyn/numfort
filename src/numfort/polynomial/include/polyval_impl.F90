

pure subroutine __APPEND(polyval_check_input,__PREC) (x, y, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:) :: y
    type (status_t), intent(out) :: status

    status = NF_STATUS_INVALID_ARG

    if (size(y) /= size(x)) return

    status = NF_STATUS_OK
end subroutine


pure subroutine __APPEND(polyval,__PREC) (coefs, x, y, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(0:) :: coefs
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(out), dimension(:) :: y
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK

    call polyval_check_input (x, y, lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    call polyval_impl (coefs, x, y)

100 continue
    if (present(status)) status = lstatus
end subroutine



pure subroutine __APPEND(polyval_impl,__PREC) (coefs, x, y)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(0:) :: coefs
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(out), dimension(:) :: y

    integer :: deg, n, i, k
    real (PREC) :: yi, xp, xi

    deg = size(coefs) - 1
    if (deg < 0) then
        y = 0.0_PREC
        return
    end if

    n = size(x)

    do i = 1, n
        yi = coefs(0)
        xp = 1.0_PREC
        xi = x(i)
        do k = 1, deg
            xp = xp * xi
            yi = yi + coefs(k) * xp
        end do
        y(i) = yi
    end do
end subroutine


pure subroutine __APPEND(polyval_scalar,__PREC) (coefs, x, y, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(0:) :: coefs
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: y
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(1) :: x1d, y1d

    x1d(1) = x
    call polyval_impl (coefs, x1d, y1d)
    y = y1d(1)

    if (present(status)) status = NF_STATUS_OK
end subroutine
