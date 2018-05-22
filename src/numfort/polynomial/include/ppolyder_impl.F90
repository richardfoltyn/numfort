

pure subroutine __APPEND(bernstein_ppolyder,__PREC) (self, knots, coefs, &
        ppoly_out, coefs_out, m, status)
    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: knots
    real (PREC), intent(in), dimension(:), contiguous :: coefs
    integer, intent(in), optional :: m
    type (ppoly_bernstein), intent(out) :: ppoly_out
    real (PREC), intent(out), dimension(:), contiguous :: coefs_out
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: k1, lm, ncoefs

    lstatus = NF_STATUS_OK

    call bernstein_check_input (self, knots, coefs, status=lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call check_nonneg (1, m, "m", lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    lm = 1
    if (present(m)) lm = m

    ! TODO: Implement support for higher derivatives
    if (lm > 1) then
        lstatus = NF_STATUS_NOT_IMPLEMENTED
        goto 100
    end if

    k1 = self%degree - lm

    ncoefs = ppoly_get_ncoefs (self, k=k1)
    if (size(coefs_out) /= ncoefs) then
        lstatus = NF_STATUS_INVALID_ARG
        goto 100
    end if

    call ppolyder_impl (self, knots, coefs, lm, ppoly_out, coefs_out, lstatus)

100 continue
    if (present(status)) status = lstatus

end subroutine



pure subroutine __APPEND(bernstein_ppolyder_impl,__PREC) (self, knots, coefs, &
        m, ppoly_out, coefs_out, status)

    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: knots
    real (PREC), intent(in), dimension(0:), contiguous :: coefs
    integer, intent(in) :: m
    type (ppoly_bernstein), intent(out) :: ppoly_out
    real (PREC), intent(out), dimension(0:), contiguous :: coefs_out
    type (status_t), intent(out) :: status

    integer :: nknots, k1, i, j, k, ii0, ii1
    real (PREC) :: dx, dcoef

    status = NF_STATUS_OK

    nknots = self%nknots
    k = self%degree
    k1 = max(0, k - m)

    ppoly_out%degree = k1
    ppoly_out%nknots = nknots

    if (m == 0) then
        coefs_out(:) = coefs
        return
    end if

    ! Input polynomial was a constant function, thus derivative is 0
    if (k == 0) then
        coefs_out(:) = 0.0_PREC
        return
    end if

    ! Handle all *resulting* polynomial degrees >= 0.
    ii0 = 0
    ii1 = 0
    do i = 1, nknots - 1
        dx = knots(i+1) - knots(i)
        do j = 0, k1
            dcoef = coefs(ii0+j+1) - coefs(ii0+j)
            ! divide by DX to apply the chain rule as we are computing the
            ! derivative wrt. x, not s = (x-x_lb)/(x_ub-x_lb) on which the
            ! Bernstein polynomials are defined.
            coefs_out(ii1+j) = k * dcoef / dx
        end do

        ! Shift "pointer" to coefficient block corresponding to current segment
        ii0 = ii0 + (k+1)
        ii1 = ii1 + (k1+1)
    end do
end subroutine