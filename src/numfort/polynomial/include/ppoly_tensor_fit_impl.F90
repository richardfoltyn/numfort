

pure subroutine __APPEND(ppoly2d_fit_check_input,__PREC) (self, x1, x2, y, &
        status)
    integer, parameter :: PREC = __PREC
    type (ppoly2d), intent(inout) :: self
    real (PREC), intent(in), dimension(:) :: x1
    real (PREC), intent(in), dimension(:) :: x2
    real (PREC), intent(in), dimension(:,:) :: y
    type (status_t), intent(out) :: status

    integer :: n1, n2

    n1 = size(x1)
    n2 = size(x2)

    ! Need at least one interval to fit piecewise polynomial
    if (n1 < 2 .or. n2 < 2) goto 100

    if (n1 /= size(y,1) .or. n2 /= size(y,2)) goto 100

    status = NF_STATUS_OK
    return

100 continue

    status = NF_STATUS_INVALID_ARG

end subroutine



pure subroutine __APPEND(ppoly2d_fit,__PREC) (self, x1, x2, y, k, &
        knots, coefs, status)

    integer, parameter :: PREC = __PREC
    type (ppoly2d), intent(inout) :: self
    real (PREC), intent(in), dimension(:) :: x1
    real (PREC), intent(in), dimension(:) :: x2
    real (PREC), intent(in), dimension(:,:) :: y
    integer, intent(in) :: k
    real (PREC), intent(out), dimension(:) :: knots
    real (PREC), intent(out), dimension(:), contiguous :: coefs
    type (status_t), intent(out), optional :: status

    type (status_t) :: lstatus
    integer, dimension(2) :: n

    lstatus = NF_STATUS_OK

    n(1) = size(x1)
    n(2) = size(x2)

    ! Perform input checking common to all routines processing bivariate
    ! piecewise polynomials.
    call ppoly2d_check_input (self, n, k, knots, coefs, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    call ppoly2d_fit_check_input (self, x1, x2, y, lstatus)
    if (lstatus /= NF_STATUS_OK) goto 100

    select case (k)
    case (1)
        call ppoly2d_fit_bilinear (self, x1, x2, y, k, knots, coefs, lstatus)
    case default
        ! TODO: implement bivariate polynomials for orders higher than bilinear
        lstatus = NF_STATUS_NOT_IMPLEMENTED
    end select


100 continue
    if (present(status)) status = lstatus

end subroutine




pure subroutine __APPEND(ppoly2d_fit_bilinear,__PREC) (self, x1, x2, y, k, &
        knots, coefs, status)
    integer, parameter :: PREC = __PREC
    type (ppoly2d), intent(inout) :: self
    real (PREC), intent(in), dimension(:) :: x1
    real (PREC), intent(in), dimension(:) :: x2
    real (PREC), intent(in), dimension(:,:) :: y
    integer, intent(in) :: k
    real (PREC), intent(out), dimension(:) :: knots
    real (PREC), intent(out), dimension(:), target, contiguous :: coefs
    type (status_t), intent(out) :: status

    real (PREC), dimension(:,:,:), pointer, contiguous :: ptr_coefs

    integer :: n1, n2, i1, i2
    real (PREC) :: a0, a1, a2, a3

    status = NF_STATUS_OK

    n1 = size(x1)
    n2 = size(x2)
    ptr_coefs(0:3,1:n1-1,1:n2-1) => coefs

    knots(1:n1) = x1
    knots(n1+1:n1+n2) = x2

    self%degree = k
    self%nknots(1) = n1
    self%nknots(2) = n2

    do i2 = 1, n2 - 1
        do i1 = 1, n1 - 1
            a0 = y(i1,i2)
            a1 = y(i1+1,i2) - a0
            a2 = y(i1,i2+1) - a0
            a3 = y(i1+1,i2+1) - a1 - a2 + 3 * a0

            ptr_coefs(0,i1,i2) = a0
            ptr_coefs(1,i1,i2) = a1
            ptr_coefs(2,i1,i2) = a2
            ptr_coefs(3,i1,i2) = a3
        end do
    end do


end subroutine