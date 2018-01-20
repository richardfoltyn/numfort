
subroutine __APPEND(lm_data_update,__PREC) (self, model, coefs, nobs, nvars, &
        ncomp, add_const, trans_rhs, var_expl)
    
    integer, parameter :: PREC  = __PREC
    type (__APPEND(lm_data,__PREC)), intent(in out) :: self
    integer, intent(in), optional :: model
    real (PREC), intent(in), dimension(:), optional :: coefs
    integer, intent(in), optional :: nobs
    integer, intent(in), optional :: nvars
    integer, intent(in), optional :: ncomp
    logical, intent(in), optional :: add_const
    logical, intent(in), optional :: trans_rhs
    real (PREC), intent(in), optional :: var_expl
    
    if (present(coefs)) then
        call copy_alloc (coefs, self%coefs)
    end if
    
    if (present(model)) self%model = model
    if (present(nobs)) self%nobs = nobs
    if (present(nvars)) self%nvars = nvars
    if (present(ncomp)) self%ncomp = ncomp
    if (present(add_const)) self%add_const = add_const
    if (present(trans_rhs)) self%trans_rhs = trans_rhs
    if (present(var_expl)) self%var_expl = var_expl
end subroutine


subroutine __APPEND(lm_data_finalize,__PREC) (self)
    type (__APPEND(lm_data,__PREC)), intent(in out) :: self
    
    if (allocated(self%coefs)) deallocate (self%coefs)
end subroutine


subroutine __APPEND(lm_data_assert_alloc_ptr,__PREC) (self, ptr)
    type (__APPEND(lm_data,__PREC)), intent(in out), target, optional :: self
    type (__APPEND(lm_data,__PREC)), intent(in out), pointer :: ptr
    
    if (present(self)) then
        ptr => self
    else
        allocate (ptr)
    end if
end subroutine


subroutine __APPEND(lm_data_assert_dealloc_ptr,__PREC) (self, ptr)
    type (__APPEND(lm_data,__PREC)), intent(in), target, optional :: self
    type (__APPEND(lm_data,__PREC)), intent(in out), pointer :: ptr

    if (associated(ptr)) then
        if (.not. present(self)) then
            call finalize (ptr)
            deallocate (ptr)
            nullify (ptr)
        else if (.not. associated (ptr, self)) then
            call finalize (ptr)
            deallocate (ptr)
            nullify (ptr)
        end if
    end if
end subroutine
