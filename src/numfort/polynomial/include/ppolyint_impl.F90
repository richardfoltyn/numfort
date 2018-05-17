

subroutine __APPEND(ppolyint_input_check,__PREC) (knots, coefs, k, &
        coefs_int, status)

    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(:) :: knots
    real (PREC), intent(in), dimension(:) :: coefs
    integer, intent(in) :: k
    real (PREC), intent(out), dimension(:) :: coefs_int
    type (status_t), intent(out), optional :: status

    integer :: n

    ! Polynomial degree
    if (k < 0) goto 100

    ! Require at least one segment, ie two knots
    if (size(knots) < 2) goto 100

    ! Check that coefficient array is large enough so that it actually does
    ! hold the required number of coefficients
    n = ppoly_size (size(knots), k)
    if (size(coefs) < n) goto 100

    ! Check that coefficient array of integrated polynomial is large enough
    n = ppoly_size (size(knots), k+1)
    if (size(coefs_int) < n) goto 100

    status = NF_STATUS_OK
    return

100 continue
    status = NF_STATUS_INVALID_ARG

end subroutine


subroutine __APPEND(ppolyint,__PREC) (knots, coefs, k, coefs_int, lbound, status)
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in), dimension(0:) :: knots
    real (PREC), intent(in), dimension(0:) :: coefs
    integer, intent(in) :: k
    real (PREC), intent(out), dimension(0:) :: coefs_int
    real (PREC), intent(in), optional :: lbound
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer :: n, i, ik, j0, j1
    logical :: has_bnd
    real (PREC) :: x0, y0, llbound

    lstatus = NF_STATUS_OK

    call ppolyint_input_check (knots, coefs, k, coefs_int, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    if (present(lbound)) llbound = lbound
    has_bnd = present(lbound)

    ! Iterate over segments, integrate each polynomial within and pin
    ! down the constant using the initial condition (lower bound) value.
    n = size(knots)
    do i = 0, n-2
        ! Segment offset "pointer" for array coefs, coefs_int
        j0 = (k+1) * i
        j1 = (k+2) * i
        ! Constant term is undetermined without some boundary condition,
        ! but we set it to zero. Any user-provided boundary condition
        ! will be applied later.
        coefs_int(j1) = 0.0
        do ik = 0, k
            coefs_int(j1+ik+1) = coefs(j0+ik) / real(ik+1, PREC)
        end do

        if (has_bnd) then
            ! Evaluate integrated polynomial at segment lower bound
            coefs_int(j1) = llbound
            ! Update lower bound for next segment, which is just the
            ! polynomial value at the segment endpoint
            x0 = knots(i+1) - knots(i)
            y0 = coefs_int(j1)
            do ik = 1, k+1
                y0 = y0 + coefs_int(j1+ik) * x0 ** ik
            end do
            ! Constant term
            llbound = y0
        end if
    end do

100 continue

    if (present(status)) status = lstatus

end subroutine