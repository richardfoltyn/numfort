

pure subroutine __APPEND(update_scalar_scalar,__PREC) (res, x, fx, status, nit, nfev, msg)
    integer, parameter :: PREC = __PREC
    type (__APPEND(optim_result,__PREC)), intent(in out) :: res
    real (PREC) :: x, fx
    integer :: nit, nfev
    type (status_t), intent(in), optional :: status
    character (len=*) :: msg

    intent (in) :: x, fx, nit, nfev, msg
    optional :: nit, nfev, msg

    real (PREC), dimension(1) :: x1, fx1

    x1(1) = x
    fx1(1) = fx

    call result_update (res, x1, fx1, status, nit, nfev, msg)
end subroutine

pure subroutine __APPEND(update_vec_scalar,__PREC) (res, x, fx, status, nit, nfev, msg)
    integer, parameter :: PREC = __PREC
    type (__APPEND(optim_result,__PREC)), intent(in out) :: res
    real (PREC) :: x(:), fx
    type (status_t), intent(in), optional :: status
    integer :: nit, nfev
    character (len=*) :: msg

    intent (in) :: x, fx, nit, nfev, msg
    optional :: nit, nfev, msg

    real (PREC), dimension(1) :: fx1

    fx1(1) = fx

    call result_update (res, x, fx1, status, nit, nfev, msg)
end subroutine

pure subroutine __APPEND(update,__PREC) (res, x, fx, status, nit, nfev, msg)
    integer, parameter :: PREC = __PREC
    type (__APPEND(optim_result,__PREC)), intent(in out) :: res
    real (PREC), dimension(:) :: x, fx
    type (status_t), intent(in), optional :: status
    integer :: nit, nfev
    character (len=*) :: msg

    intent (in) :: x, fx, nit, nfev, msg
    optional :: x, fx, nit, nfev, msg

    call result_reset (res)

    call alloc_assign (x, res%x)
    call alloc_assign (fx, res%fx)

    if (present(status)) then
        res%success = status_contains (status, NF_STATUS_OK)
        res%status = status
    end if

    if (present(msg)) res%msg = msg

    if (present(nit)) res%nit = nit
    if (present(nfev)) res%nfev = nfev

end subroutine

pure subroutine __APPEND(reset,__PREC) (res)
    integer, parameter :: PREC = __PREC
    type (__APPEND(optim_result,__PREC)), intent(in out) :: res

    res%nit = UNINITIALIZED_COUNTER
    res%nfev = UNINITIALIZED_COUNTER
    res%msg = ""
    res%x = 0.0_PREC
    res%fx = 0.0_PREC
    call status_set (res%status, NF_STATUS_UNDEFINED)
    res%success = .false.
end subroutine
