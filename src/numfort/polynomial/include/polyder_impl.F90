

pure subroutine __APPEND(polyder_check_input,__PREC) (x, k, y, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: x
    integer,intent(in) :: k
    real (PREC), intent(in), dimension(:) :: y
    type (status_t), intent(out) :: status

    status = NF_STATUS_INVALID_ARG

    if (size(x) /= size(y)) return
    if (k < 0) return

    status = NF_STATUS_OK
end subroutine


pure subroutine __APPEND(polyder,__PREC) (coefs, x, k, y, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(0:) :: coefs
    real (PREC), intent(in), dimension(:) :: x
    integer,intent(in) :: k
    real (PREC), intent(out), dimension(:) :: y
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: deg, n, i, j
    real (PREC) :: xi, xp, yi
    real (PREC), dimension(:), allocatable :: fact
    
    lstatus = NF_STATUS_OK

    call polyder_check_input (x, k, y, lstatus)
    if (NF_STATUS_INVALID_ARG .in. lstatus) goto 100

    deg = size(coefs) - 1
    n = size(x)

    if (k > deg .or. deg < 0) then
        y = 0.0_PREC
        goto 100
    end if

    allocate (fact(k:deg))

    fact(k) = factorial (k)
    do j = k+1, deg
        fact(j) = fact(j-1) / (j-k) * j
    end do

    do i = 1, n
        yi = coefs(k) * fact(k)
        xp = 1.0_PREC
        xi = x(i)
        do j = k+1, deg
            xp = xp * xi
            yi = yi + coefs(j) * fact(j) * xp
        end do
        y(i) = yi
    end do

    deallocate (fact)

100 continue
    if (present(status)) status = lstatus

end subroutine


pure subroutine __APPEND(polyder_scalar,__PREC) (coefs, x, k, y, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(0:) :: coefs
    real (PREC), intent(in) :: x
    integer,intent(in) :: k
    real (PREC), intent(out) :: y
    type (status_t), intent(out), optional :: status

    real (PREC), dimension(1) :: x1d, y1d

    x1d(1) = x
    call polyder (coefs, x1d, k, y1d, status)
    y = y1d(1)
end subroutine

