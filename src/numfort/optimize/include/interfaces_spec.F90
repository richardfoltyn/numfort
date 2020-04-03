
subroutine fss_args (x, args, fx)
    !*  Abstract interface for a scalar-valued function f:R->R
    !   that takes additional arguments.
    import PREC
    import args_data
    real (PREC), intent(in)  :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out) :: fx
end subroutine

subroutine fss (x, fx)
    !*  Abstract interface for a scalar-valued function f:R->R.
    import PREC
    real (PREC), intent(in)  :: x
    real (PREC), intent(out) :: fx
end subroutine

subroutine fss_fcn_jac (x, fx, fpx)
    !*  Abstract interface for a scalar-valued function f:R->R
    !   which also returns the first derivative.
    import PREC
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx
    real (PREC), intent(out) :: fpx
end subroutine

subroutine fss_fcn_jac_args (x, args, fx, fpx)
    !*  Abstract interface for a scalar-valued function f:R->R
    !   that takes additional arguments and also returns the first derivative.
    import PREC
    import args_data
    real (PREC), intent(in) :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out) :: fx
    real (PREC), intent(out) :: fpx
end subroutine

subroutine fss_fcn_jac_opt (x, fx, fpx)
    !*  Abstract interface for a scalar-valued function f:R->R
    !   which (optionally) returns returns the function value and/or
    !   first derivative.
    import PREC
    real (PREC), intent(in) :: x
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), optional :: fpx
end subroutine

subroutine fss_fcn_jac_opt_args (x, args, fx, fpx)
    !*  Abstract interface for a scalar-valued function f:R->R
    !   which (optionally) returns returns the function value and/or
    !   first derivative and takes additional arguments.
    import PREC
    import args_data
    real (PREC), intent(in) :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), optional :: fpx
end subroutine


! ------------------------------------------------------------------------------
! Interfaces for function mapping vectors into scalars

subroutine fvs_fcn (x, fx)
    import PREC
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out) :: fx
end subroutine

subroutine fvs_jac (x, fpx)
    import PREC
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fpx
end subroutine

subroutine fvs_fcn_jac (x, fx, fpx)
    import PREC
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out) :: fx
    real (PREC), intent(out), dimension(:), contiguous :: fpx
end subroutine

subroutine fvs_fcn_jac_opt (x, fx, fpx)
    import PREC
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), dimension(:), contiguous, optional :: fpx
end subroutine

subroutine fvs_fcn_args (x, args, fx)
    import PREC
    import args_data
    real (PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out) :: fx
end subroutine

subroutine fvs_jac_args (x, args, fpx)
    import PREC
    import args_data
    real (PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), dimension(:), contiguous, optional :: fpx
end subroutine

subroutine fvs_fcn_jac_args (x, args, fx, fpx)
    import PREC
    import args_data
    real (PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out) :: fx
    real (PREC), intent(out), dimension(:), contiguous :: fpx
end subroutine

subroutine fvs_fcn_jac_opt_args (x, args, fx, fpx)
    import PREC
    import args_data
    real (PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), dimension(:), contiguous, optional :: fpx
end subroutine


! ------------------------------------------------------------------------------
! Interfaces for functions mapping vectors into vectors

subroutine fvv_fcn (x, fx)
    import PREC
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx
end subroutine

subroutine fvv_jac (x, fpx)
    import PREC
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:,:), contiguous :: fpx
end subroutine

subroutine fvv_fcn_jac (x, fx, fpx)
    import PREC
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx
    real (PREC), intent(out), dimension(:,:), contiguous :: fpx
end subroutine

subroutine fvv_fcn_jac_opt (x, fx, fpx)
    import PREC
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous, optional :: fx
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: fpx
end subroutine

subroutine fvv_fcn_args (x, args, fx)
    import PREC
    import args_data
    real (PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), dimension(:), contiguous :: fx
end subroutine

subroutine fvv_jac_args (x, args, fpx)
    import PREC
    import args_data
    real (PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: fpx
end subroutine

subroutine fvv_fcn_jac_args (x, args, fx, fpx)
    import PREC
    import args_data
    real (PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), dimension(:), contiguous :: fx
    real (PREC), intent(out), dimension(:,:), contiguous :: fpx
end subroutine

subroutine fvv_fcn_jac_opt_args (x, args, fx, fpx)
    import PREC
    import args_data
    real (PREC), intent(in), dimension(:), contiguous :: x
    class (args_data), intent(inout) :: args
    real (PREC), intent(out), dimension(:), contiguous, optional :: fx
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: fpx
end subroutine

