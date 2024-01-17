

pure subroutine update_ss (res, x, fx, status, nit, nfev, msg)
    type (optim_result), intent(inout) :: res
    real (PREC), intent(in) :: x, fx
    type (status_t), intent(in), optional :: status
    integer, intent(in), optional :: nit
    integer, intent(in), optional :: nfev
    character (len=*), intent(in), optional :: msg

    real (PREC), dimension(1) :: x1, fx1

    x1(1) = x
    fx1(1) = fx

    call result_update (res, x1, fx1, status=status, nit=nit, nfev=nfev, msg=msg)
end subroutine



pure subroutine update_vs (res, x, fx, status, nit, nfev, msg)
    type (optim_result), intent(inout) :: res
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(in) :: fx
    type (status_t), intent(in), optional :: status
    integer, intent(in), optional :: nit
    integer, intent(in), optional :: nfev
    character (len=*), intent(in), optional :: msg

    real (PREC), dimension(1) :: fx1

    fx1(1) = fx

    call result_update (res, x, fx1, status=status, nit=nit, nfev=nfev, msg=msg)
end subroutine



pure subroutine update (res, x, fx, jac_inv, status, nit, nfev, msg)
    type (optim_result), intent(inout) :: res
    real (PREC), intent(in), dimension(:), optional :: x
    real (PREC), intent(in), dimension(:), optional :: fx
    real (PREC), intent(in), dimension(:,:), optional :: jac_inv
    type (status_t), intent(in), optional :: status
    integer, intent(in), optional :: nit
    integer, intent(in), optional :: nfev
    character (len=*), intent(in), optional :: msg

    if (present(x)) call copy_alloc (x, res%x)
    if (present(fx)) call copy_alloc (fx, res%fx)
    if (present(jac_inv)) call copy_alloc (jac_inv, res%jac_inv)

    if (present(status)) then
        res%status = status
    end if

    res%success = (NF_STATUS_OK .in. res%status)

    if (present(msg)) res%msg = msg

    if (present(nit)) res%nit = nit
    if (present(nfev)) res%nfev = nfev

end subroutine



pure subroutine reset (res)
    type (optim_result), intent(inout) :: res

    res%nit = UNINITIALIZED_COUNTER
    res%nfev = UNINITIALIZED_COUNTER
    res%msg = ""
    if (allocated(res%x)) res%x = 0.0_PREC
    if (allocated(res%fx)) res%fx = 0.0_PREC
    if (allocated(res%jac_inv)) deallocate (res%jac_inv)
    res%status = NF_STATUS_UNDEFINED
    res%success = .false.
end subroutine



pure subroutine assert_alloc_ptr (res, ptr)
    type (optim_result), intent(inout), target, optional :: res
    type (optim_result), pointer, intent(inout) :: ptr

    if (present(res)) then
        ptr => res
    else
        allocate (ptr)
    end if
end subroutine



pure subroutine assert_dealloc_ptr (res, ptr)
    type (optim_result), intent(inout), target, optional :: res
    type (optim_result), pointer, intent(inout) :: ptr

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



subroutine result_finalize (self)
    !*  RESULT_FINALIZE deallocates any allocated (allocatable) attributes.
    type (optim_result), intent(inout) :: self

    if (allocated(self%x)) deallocate (self%x)
    if (allocated(self%fx)) deallocate (self%fx)

end subroutine
