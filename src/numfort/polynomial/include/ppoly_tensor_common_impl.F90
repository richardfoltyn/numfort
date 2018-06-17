


pure subroutine __APPEND(ppoly2d_check_input,__PREC) (self, n, k, knots, coefs, &
        status)
    integer, parameter :: PREC = __PREC
    type (ppoly2d), intent(in) :: self
    integer, intent(in), dimension(:) :: n
    integer, intent(in) :: k
    real (PREC), intent(in), dimension(:) :: knots
    real (PREC), intent(in), dimension(:) :: coefs
    type (status_t), intent(out) :: status

    integer :: ncoefs, nknots

    status = NF_STATUS_OK

    call check_nonneg (1, k, "k", status)
    if (status /= NF_STATUS_OK) goto 100

    nknots = ppoly_get_nknots (self, n, k)
    ncoefs = ppoly_get_ncoefs (self, n, k)

    if (size(knots) /= nknots) goto 100
    if (size(coefs) /= ncoefs) goto 100

    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine
