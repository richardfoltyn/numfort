
subroutine __APPEND(fss_args,__PREC) (x, args, fx)
    !*  Abstract interface for a scalar-valued function f:R->R
    !   that takes additional arguments.
    import __PREC
    import args_data
    real (__PREC), intent(in)  :: x
    class (args_data), intent(inout) :: args
    real (__PREC), intent(out) :: fx
end subroutine

subroutine __APPEND(fss,__PREC) (x, fx)
    !*  Abstract interface for a scalar-valued function f:R->R.
    import __PREC
    real (__PREC), intent(in)  :: x
    real (__PREC), intent(out) :: fx
end subroutine

subroutine __APPEND(fss_fcn_jac,__PREC) (x, fx, fpx)
    !*  Abstract interface for a scalar-valued function f:R->R
    !   which also returns the first derivative.
    import __PREC
    real (__PREC), intent(in) :: x
    real (__PREC), intent(out) :: fx
    real (__PREC), intent(out) :: fpx
end subroutine

subroutine __APPEND(fss_fcn_jac_args,__PREC) (x, args, fx, fpx)
    !*  Abstract interface for a scalar-valued function f:R->R
    !   that takes additional arguments and also returns the first derivative.
    import __PREC
    import args_data
    real (__PREC), intent(in) :: x
    class (args_data), intent(inout) :: args
    real (__PREC), intent(out) :: fx
    real (__PREC), intent(out) :: fpx
end subroutine

subroutine __APPEND(fss_fcn_jac_opt,__PREC) (x, fx, fpx)
    !*  Abstract interface for a scalar-valued function f:R->R
    !   which (optionally) returns returns the function value and/or
    !   first derivative.
    import __PREC
    real (__PREC), intent(in) :: x
    real (__PREC), intent(out), optional :: fx
    real (__PREC), intent(out), optional :: fpx
end subroutine

subroutine __APPEND(fss_fcn_jac_opt_args,__PREC) (x, args, fx, fpx)
    !*  Abstract interface for a scalar-valued function f:R->R
    !   which (optionally) returns returns the function value and/or
    !   first derivative and takes additional arguments.
    import __PREC
    import args_data
    real (__PREC), intent(in) :: x
    class (args_data), intent(inout) :: args
    real (__PREC), intent(out), optional :: fx
    real (__PREC), intent(out), optional :: fpx
end subroutine


! ------------------------------------------------------------------------------
! Interfaces for function mapping vectors into scalars

subroutine __APPEND(fvs_fcn,__PREC) (x, fx)
    import __PREC
    real (__PREC), intent(in), dimension(:), contiguous :: x
    real (__PREC), intent(out) :: fx
end subroutine

subroutine __APPEND(fvs_jac,__PREC) (x, fpx)
    import __PREC
    real (__PREC), intent(in), dimension(:), contiguous :: x
    real (__PREC), intent(out), dimension(:), contiguous :: fpx
end subroutine

subroutine __APPEND(fvs_fcn_jac,__PREC) (x, fx, fpx)
    import __PREC
    real (__PREC), intent(in), dimension(:), contiguous :: x
    real (__PREC), intent(out) :: fx
    real (__PREC), intent(out), dimension(:), contiguous :: fpx
end subroutine

subroutine __APPEND(fvs_fcn_jac_opt,__PREC) (x, fx, fpx)
    import __PREC
    real (__PREC), intent(in), dimension(:), contiguous :: x
    real (__PREC), intent(out), optional :: fx
    real (__PREC), intent(out), dimension(:), contiguous, optional :: fpx
end subroutine

subroutine __APPEND(fvs_fcn_args,__PREC) (x, args, fx)
    import __PREC
    import args_data
    real (__PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (__PREC), intent(out) :: fx
end subroutine

subroutine __APPEND(fvs_jac_args,__PREC) (x, args, fpx)
    import __PREC
    import args_data
    real (__PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (__PREC), intent(out), dimension(:), contiguous, optional :: fpx
end subroutine

subroutine __APPEND(fvs_fcn_jac_args,__PREC) (x, args, fx, fpx)
    import __PREC
    import args_data
    real (__PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (__PREC), intent(out) :: fx
    real (__PREC), intent(out), dimension(:), contiguous :: fpx
end subroutine

subroutine __APPEND(fvs_fcn_jac_opt_args,__PREC) (x, args, fx, fpx)
    import __PREC
    import args_data
    real (__PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (__PREC), intent(out), optional :: fx
    real (__PREC), intent(out), dimension(:), contiguous, optional :: fpx
end subroutine


! ------------------------------------------------------------------------------
! Interfaces for functions mapping vectors into vectors

subroutine __APPEND(fvv_fcn,__PREC) (x, fx)
    import __PREC
    real (__PREC), intent(in), dimension(:), contiguous :: x
    real (__PREC), intent(out), dimension(:), contiguous :: fx
end subroutine

subroutine __APPEND(fvv_jac,__PREC) (x, fpx)
    import __PREC
    real (__PREC), intent(in), dimension(:), contiguous :: x
    real (__PREC), intent(out), dimension(:,:), contiguous :: fpx
end subroutine

subroutine __APPEND(fvv_fcn_jac,__PREC) (x, fx, fpx)
    import __PREC
    real (__PREC), intent(in), dimension(:), contiguous :: x
    real (__PREC), intent(out), dimension(:), contiguous :: fx
    real (__PREC), intent(out), dimension(:,:), contiguous :: fpx
end subroutine

subroutine __APPEND(fvv_fcn_jac_opt,__PREC) (x, fx, fpx)
    import __PREC
    real (__PREC), intent(in), dimension(:), contiguous :: x
    real (__PREC), intent(out), dimension(:), contiguous, optional :: fx
    real (__PREC), intent(out), dimension(:,:), contiguous, optional :: fpx
end subroutine

subroutine __APPEND(fvv_fcn_args,__PREC) (x, args, fx)
    import __PREC
    import args_data
    real (__PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (__PREC), intent(out), dimension(:), contiguous :: fx
end subroutine

subroutine __APPEND(fvv_jac_args,__PREC) (x, args, fpx)
    import __PREC
    import args_data
    real (__PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (__PREC), intent(out), dimension(:,:), contiguous, optional :: fpx
end subroutine

subroutine __APPEND(fvv_fcn_jac_args,__PREC) (x, args, fx, fpx)
    import __PREC
    import args_data
    real (__PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (__PREC), intent(out), dimension(:), contiguous :: fx
    real (__PREC), intent(out), dimension(:,:), contiguous :: fpx
end subroutine

subroutine __APPEND(fvv_fcn_jac_opt_args,__PREC) (x, args, fx, fpx)
    import __PREC
    import args_data
    real (__PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (__PREC), intent(out), dimension(:), contiguous, optional :: fx
    real (__PREC), intent(out), dimension(:,:), contiguous, optional :: fpx
end subroutine

