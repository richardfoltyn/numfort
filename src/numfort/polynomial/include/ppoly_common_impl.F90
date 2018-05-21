

pure subroutine __APPEND(bernstein_check_input,__PREC) (self, knots, coefs, &
        n, k, status)

    integer, parameter :: PREC = __PREC
    type (ppoly_bernstein), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: knots
        !*  Array of knots
    real (PREC), intent(in), dimension(:) :: coefs
        !*  Coefficient array
    integer, intent(in), optional :: n
        !*  Number of data points. Only relevant for routines that create
        !   piecewise polynomials.
    integer, intent(in), optional :: k
        !*  Polynomial degree. Only relevant for routines that create
        !   piecewise polynomials.
    type (status_t), intent(out) :: status

    integer :: nknots, ncoefs

    nknots = ppoly_get_nknots (self, n, k)
    ncoefs = ppoly_get_ncoefs (self, n, k)

    if (size(knots) < nknots) goto 100
    if (size(coefs) /= ncoefs) goto 100

    status = NF_STATUS_OK
    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine

