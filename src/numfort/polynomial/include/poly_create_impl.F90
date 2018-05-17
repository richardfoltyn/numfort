

pure subroutine __APPEND(polyder_check_input,__PREC) (coefs, coefs_new, m, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: coefs
    real (PREC), intent(in), dimension(:) :: coefs_new
    integer, intent(in), optional :: m
    type (status_t), intent(out) :: status

    integer :: n
    ! Valid order of derivative
    if (m < 0) goto 100

    if (size(coefs) < 1) goto 100

    ! Check array size
    ! Note: Min. array size for new coefficients is 1, even if derivative
    ! is zero (ie the polynomial has all zero coefs, ie. one constant coef
    ! equal to zero)
    n = max(1, size(coefs) - m)
    if (size(coefs_new) < n) goto 100

    status = NF_STATUS_OK
    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine



pure subroutine __APPEND(polyder,__PREC) (coefs, coefs_new, m, status)
    !*  POLYDER differentiates a given polynomial of degree k and returns
    !   the coefficients of the resulting polynomial of degree (k-1).
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(0:) :: coefs
    real (PREC), intent(out), dimension(0:) :: coefs_new
    integer, intent(in), optional :: m
        !*  m-th derivative to compute.
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: lm, k, b_i, ik

    lstatus = NF_STATUS_OK

    call polyder_check_input (coefs, coefs_new, m, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lm = 1
    if (present(m)) lm = m

    k = size(coefs) - 1

    coefs_new = 0.0_PREC

    b_i = factorial (lm)
    do ik = 0, k-lm
        coefs_new(ik) = b_i * coefs(ik+lm)
        b_i = (ik + lm + 1) * b_i / (ik + 1)
    end do

100 continue

    if (present(status)) status = lstatus
end subroutine



pure subroutine __APPEND(polyint_check_input,__PREC) (coefs, coefs_new, x, y, &
        status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: coefs
    real (PREC), intent(in), dimension(:) :: coefs_new
    real (PREC), intent(in), optional :: x
    real (PREC), intent(in), optional :: y
    type (status_t), intent(out), optional :: status

    if (size(coefs_new) < (size(coefs) + 1)) goto 100

    if ((present(x) .and. .not. present(y)) .or. &
        (.not. present(x) .and. present(y))) goto 100

    status = NF_STATUS_OK
    return

100 continue

    status = NF_STATUS_INVALID_ARG

end subroutine



pure subroutine __APPEND(polyint,__PREC) (coefs, coefs_new, x, y, status)
    !*  POLYINT constructs a new polynomial as the integral of a given
    !   polynomial, defined by its coefficient array.
    !   The integration constant is optionally determined by providing
    !   an arbitrary point (x,y) from the graph of the antiderivative.
    integer, parameter ::  PREC = __PREC
    real (PREC), intent(in), dimension(0:) :: coefs
        !*  Coefficients of polynomial to be integrated
    real (PREC), intent(out), dimension(0:) :: coefs_new
        !*  Coefficients of polynomial representing the antiderivative
    real (PREC), intent(in), optional :: x
        !*  Optional point X used to determine the integration constant.
    real (PREC), intent(in), optional :: y
        !*  Optional point Y used to determine the integration constant,
        !   defined such that P(x) = Y where P is the antiderivative.
    type (status_t), intent(out), optional :: status
        !*  Optional status code.

    type (status_t) :: lstatus

    lstatus = NF_STATUS_OK

    call polyint_check_input (coefs, coefs_new, x, y, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call polyint_impl (coefs, coefs_new, x, y)

100 continue

    if (present(status)) status = lstatus

end subroutine


pure subroutine __APPEND(polyint_impl,__PREC) (coefs, coefs_new, x, y)
    !*  POLYINT_IMPL provides the actual integration functionality of
    !   POLYINT without the error-checking overhead.
    integer, parameter ::  PREC = __PREC
    real (PREC), intent(in), dimension(0:) :: coefs
        !*  Coefficients of polynomial to be integrated
    real (PREC), intent(out), dimension(0:) :: coefs_new
        !*  Coefficients of polynomial representing the antiderivative
    real (PREC), intent(in), optional :: x
        !*  Optional point X used to determine the integration constant.
    real (PREC), intent(in), optional :: y
        !*  Optional point Y used to determine the integration constant,
        !   defined such that P(x) = Y where P is the antiderivative.

    integer :: k, ik
    real (PREC) :: fx(1), x1(1)

    k = size(coefs) - 1

    coefs_new = 0.0
    do ik = 1, k+1
        coefs_new(ik) = coefs(ik-1) / ik
    end do

    if (present(x) .and. present(y)) then
        x1(1) = x
        call polyval_impl (coefs_new, x1, fx)
        coefs_new(0) = y - fx(1)
    end if

end subroutine



pure subroutine __APPEND(polyshift,__PREC) (coefs, x0, coefs_new, status)
    !*  POLYSHIFT computes the coefficients of a polynomial that results
    !   from shifting the domain of a given polynomial relative to X0.
    !
    !   That is, for a polynomial p(x) defined by the coefficient COEFS,
    !   this routine returns the coefficients of the polynomial p_s(s)
    !   where s = x - x0 and p(x) = p_s(x-x0) for all x.

    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(0:) :: coefs
        !*  Polynomial coefficients. The polynomial degree K is inferred as
        !   K = size(coefs) - 1.
    real (PREC), intent(in) :: x0
        !*  Value by which polynomial domain should be "shifted".
    real (PREC), intent(out), dimension(0:) :: coefs_new
        !*  Coefficients of "shifted" polynomial.
    type (status_t), intent(out), optional :: status
        !*  Optional exit code.

    type (status_t) :: lstatus
    integer :: k, i, j
    real (PREC), dimension(:), allocatable :: xp
    real (PREC) :: b_i

    lstatus = NF_STATUS_OK

    if (size(coefs) /= size(coefs_new)) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    k = size(coefs) - 1

    ! Procompute power series in x0
    allocate (xp(0:k))
    forall (i=0:k) xp(i) = x0 ** i

    ! Obtaining the coefficients of the shifted polynomial is an application
    ! of Pascal's triangle.
    ! We start with rewriting the original polynomial as
    !       p(x) = a_k (s+x0)^k + ... + a_1 (s+x0) + a_0
    ! The coefficient of the term s^i, i in {0,...,k} is a sum of expressions that
    ! show up after expanding the terms a_j (s+x0)^j in the original polynomial
    ! with j >= i.
    ! We can then derive that the coefficient for s^i has to be
    !       b_i = sum_{j>=i} a_j binom(j,i) x0^(j-i)
    ! Then the shifted polynomial can be written as
    !       p_s(s) = b_k * s^k + ... + b_1 s + b_0
    do i = 0, k
        b_i = 0.0
        do j = i, k
            b_i = b_i + coefs(j) * comb(j, i) * xp(j-i)
        end do
        coefs_new(i) = b_i
    end do

100 continue

    if (allocated(xp)) deallocate (xp)

    if (present(status)) status = lstatus

end subroutine

