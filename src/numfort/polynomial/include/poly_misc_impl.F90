


subroutine __APPEND(polyshift,__PREC) (coefs, x0, coefs_s, status)
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
    real (PREC), intent(out), dimension(0:) :: coefs_s
        !*  Coefficients of "shifted" polynomial.
    type (status_t), intent(out), optional :: status
        !*  Optional exit code.

    type (status_t) :: lstatus
    integer :: k, i, j
    real (PREC), dimension(:), allocatable :: xp
    real (PREC) :: b_i

    lstatus = NF_STATUS_OK

    if (size(coefs) /= size(coefs_s)) then
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
        coefs_s(i) = b_i
    end do

100 continue

    if (allocated(xp)) deallocate (xp)

    if (present(status)) status = lstatus

end subroutine

