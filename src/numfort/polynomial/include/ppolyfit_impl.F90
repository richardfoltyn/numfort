


pure subroutine __APPEND(bernstein_fit_deriv_check_input,__PREC) &
        (self, x, y, k, status)

    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:,:) :: y
    integer, intent(in) :: k
    type (status_t), intent(out) :: status

    integer :: n

    n = size(x)

    call check_nonneg (1, k, "k", status)
    if (status /= NF_STATUS_OK) goto 100

    ! Need at least one interval to fit a piecewise polynomial
    if (n < 2) goto 100

    if (size(x) /= size(y,2)) goto 100
    if (size(y,1) < ((k+1) / 2)) goto 100

    status = NF_STATUS_OK
    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine


pure subroutine __APPEND(bernstein_fit_deriv,__PREC) (self, x, y, k, &
        knots, coefs, status)
    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(inout) :: self
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(:,:) :: y
    integer, intent(in) :: k
    real (PREC), intent(out), dimension(:) :: knots
    real (PREC), intent(out), dimension(:), target, contiguous :: coefs
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: n

    lstatus = NF_STATUS_OK
    n = size(x)

    ! Perform input checking common to all routines processing Bernstein
    ! polynomials.
    call bernstein_check_input (self, knots, coefs, n, k, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    ! Perform input checking specifing to fitting polynomials to derivatives.
    call bernstein_fit_deriv_check_input (self, x, y, k, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call bernstein_fit_deriv_impl (self, x, y, k, knots, coefs, lstatus)

100 continue
    if (present(status)) status = lstatus
end subroutine


pure subroutine __APPEND(bernstein_fit_deriv_impl,__PREC) (self, x, y, k, &
        knots, coefs, status)
    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(inout) :: self
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in), dimension(0:,:) :: y
    integer, intent(in) :: k
    real (PREC), intent(out), dimension(:) :: knots
    real (PREC), intent(out), dimension(:), target, contiguous :: coefs
    type (status_t), intent(out) ::  status

    real (PREC), dimension(:,:), pointer, contiguous :: ptr_coefs
    integer :: n, i, ka, kb, q, j, nknots
    real (PREC) :: dx, cj, p
    real (PREC), dimension(:), allocatable :: ci

    status = NF_STATUS_OK

    n = size(x)
    nknots = n
    ptr_coefs(0:k,1:n-1) => coefs

    self%degree = k
    self%nknots = nknots

    ka = size(y, 1)
    kb = ka
    allocate (ci(0:k), source=0.0_PREC)

    do i = 1, n-1
        dx = x(i+1) - x(i)
        do q = 0, ka - 1
            p = poch (k+1-q, q)
            ci(q) = y(q,i) / p * dx ** q
            cj = 0.0
            do j = 0, q-1
                cj = cj - (-1)**(j+q) * comb(q, j) * ci(j)
            end do
            ci(q) = ci(q) + cj
        end do

        do q = 0, kb - 1
            p = poch (k+1-q, q)
            ci(k-q) = y(q,i+1) / p * (-1)**q * dx**q
            cj = 0.0
            do j = 0, q-1
                cj = cj - (-1)**(j+1) * comb (q, j+1) * ci(k+1-q+j)
            end do
            ci(k-q) = ci(k-q) + cj
        end do

        ptr_coefs(:,i) = ci
    end do

    deallocate (ci)

    ! For PCHIP-type
    knots = x

end subroutine
