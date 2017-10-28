

pure subroutine __APPEND(update_ss,__PREC) (res, x, fx, status, nit, nfev, msg)
    integer, parameter :: PREC = __PREC
    type (__APPEND(optim_result,__PREC)), intent(in out) :: res
    real (PREC), intent(in) :: x, fx
    type (status_t), intent(in), optional :: status
    integer, intent(in), optional :: nit
    integer, intent(in), optional :: nfev
    character (len=*), intent(in), optional :: msg

    real (PREC), dimension(1) :: x1, fx1

    x1(1) = x
    fx1(1) = fx

    call result_update (res, x1, fx1, status, nit, nfev, msg)
end subroutine

pure subroutine __APPEND(update_vs,__PREC) (res, x, fx, status, nit, nfev, msg)
    integer, parameter :: PREC = __PREC
    type (__APPEND(optim_result,__PREC)), intent(in out) :: res
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in) :: fx
    type (status_t), intent(in), optional :: status
    integer, intent(in), optional :: nit
    integer, intent(in), optional :: nfev
    character (len=*), intent(in), optional :: msg

    real (PREC), dimension(1) :: fx1

    fx1(1) = fx

    call result_update (res, x, fx1, status, nit, nfev, msg)
end subroutine

pure subroutine __APPEND(update,__PREC) (res, x, fx, status, nit, nfev, msg)
    integer, parameter :: PREC = __PREC
    type (__APPEND(optim_result,__PREC)), intent(in out) :: res
    real (PREC), intent(in), dimension(:), optional :: x
    real (PREC), intent(in), dimension(:), optional :: fx
    type (status_t), intent(in), optional :: status
    integer, intent(in), optional :: nit
    integer, intent(in), optional :: nfev
    character (len=*), intent(in), optional :: msg

    if (present(x)) call copy_alloc (x, res%x)
    if (present(fx)) call copy_alloc (fx, res%fx)

    if (present(status)) then
        res%status = status
    end if

    res%success = (NF_STATUS_OK .in. res%status)

    if (present(msg)) res%msg = msg

    if (present(nit)) res%nit = nit
    if (present(nfev)) res%nfev = nfev

end subroutine


pure subroutine __APPEND(update_int_status,__PREC) (res, x, fx, status, nit, &
        nfev, msg, istatus)
    !*  UPDATE_INT_STATUS updates OPTIM_RESULT object attributes with given
    !   values, accepting an integer status value instead of STATUS_T.
    integer, parameter :: PREC = __PREC
    type (__APPEND(optim_result,__PREC)), intent(in out) :: res
    real (PREC), intent(in), dimension(:), optional :: x
    real (PREC), intent(in), dimension(:), optional :: fx
    integer (NF_ENUM_KIND), intent(in) :: status
    integer, intent(in), optional :: nit
    integer, intent(in), optional :: nfev
    character (len=*), intent(in), optional :: msg
    integer, intent(in), optional :: istatus
        !*  Optional "original" integer status returned by the underlying
        !   implementation.

    type (status_t) :: lstatus

    lstatus = status
    if (present(istatus)) lstatus%code_orig = istatus

    call result_update (res, x, fx, lstatus, nit, nfev, msg)
end subroutine


pure subroutine __APPEND(reset,__PREC) (res)
    integer, parameter :: PREC = __PREC
    type (__APPEND(optim_result,__PREC)), intent(in out) :: res

    res%nit = UNINITIALIZED_COUNTER
    res%nfev = UNINITIALIZED_COUNTER
    res%msg = ""
    if (allocated(res%x)) res%x = 0.0_PREC
    if (allocated(res%fx)) res%fx = 0.0_PREC
    res%status = NF_STATUS_UNDEFINED
    res%success = .false.
end subroutine


pure subroutine __APPEND(assert_alloc_ptr,__PREC) (res, ptr)
    type (__APPEND(optim_result,__PREC)), intent(in out), target, optional :: res
    type (__APPEND(optim_result,__PREC)), pointer, intent(in out) :: ptr

    if (present(res)) then
        ptr => res
    else
        allocate (ptr)
    end if
end subroutine


pure subroutine __APPEND(assert_dealloc_ptr,__PREC) (res, ptr)
    type (__APPEND(optim_result,__PREC)), intent(in out), target, optional :: res
    type (__APPEND(optim_result,__PREC)), pointer, intent(in out) :: ptr

    if (associated(ptr)) then
        if (.not. present(res)) then
            deallocate (ptr)
            nullify (ptr)
        else if (.not. associated(ptr,res)) then
            deallocate (ptr)
            nullify (ptr)
        end if
    end if
end subroutine
