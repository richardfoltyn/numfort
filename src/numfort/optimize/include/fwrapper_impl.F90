
!-------------------------------------------------------------------------------
! Routines for FWRAPPER_SS

subroutine __APPEND(fss_init,__PREC) (self, fcn, jac, fcn_jac, fcn_jac_opt, &
        fcn_args, jac_args, fcn_jac_args, fcn_jac_opt_args, args, eps)
    !*  FSS_INIT intializes a wrapper found a scalar-valued function f:R->R.
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_ss,__PREC)), intent(inout) :: self
    procedure (__APPEND(fss,__PREC)), optional :: fcn
        !*  Pointer to function that returns the function value
    procedure (__APPEND(fss,__PREC)), optional :: jac
        !*  Pointer to function that returns the first derivative
    procedure (__APPEND(fss_fcn_jac,__PREC)), optional :: fcn_jac
        !*  Pointer to function that returns the function value and the
        !   first derivative in a single call.
    procedure (__APPEND(fss_fcn_jac_opt,__PREC)), optional :: fcn_jac_opt
        !*  Pointer to function that returns the function value and/or the
        !   first derivative in a single call.
    procedure (__APPEND(fss_args,__PREC)), optional :: fcn_args
        !*  Pointer to function that returns the function value and accepts
        !   additional arguments.
    procedure (__APPEND(fss_args,__PREC)), optional :: jac_args
        !*  Pointer to function that returns the first derivative and accepts
        !   additional arguments.
    procedure (__APPEND(fss_fcn_jac_args,__PREC)), optional :: fcn_jac_args
        !*  Pointer to function that returns the function value and the
        !   first derivative in a single call and accepts additional arguments.
    procedure (__APPEND(fss_fcn_jac_opt_args,__PREC)), optional :: fcn_jac_opt_args
        !*  Pointer to function that returns the function value and/or the
        !   first derivative in a single call and accepts additional arguments.
    real (PREC), intent(in), dimension(:), optional, target :: args
        !*  Additional arguments that need to be passed to function.
    real (PREC), intent(in), optional :: eps
        !*  If present, sets the step size used when computing numerical
        !   derivative if no other way to obtain the derivative is avaiable.

    nullify (self%ptr_args)
    nullify (self%fcn, self%jac, self%fcn_jac, self%fcn_jac_opt)
    nullify (self%fcn_args, self%jac_args, self%fcn_jac_args, self%fcn_jac_opt_args)
    self%nfev = 0
    self%eps = sqrt(epsilon(0.0_PREC))

    if (present(args)) then
        if (present(fcn_args)) self%fcn_args => fcn_args
        if (present(jac_args)) self%jac_args => jac_args
        if (present(fcn_jac_args)) self%fcn_jac_args => fcn_jac_args
        if (present(fcn_jac_opt_args)) self%fcn_jac_opt_args => fcn_jac_opt_args
        self%ptr_args => args
    else
        if (present(fcn)) self%fcn => fcn
        if (present(jac)) self%jac => jac
        if (present(fcn_jac)) self%fcn_jac => fcn_jac
        if (present(fcn_jac_opt)) self%fcn_jac_opt => fcn_jac_opt
    end if

    if (present(eps)) then
        self%eps = eps
    end if

end subroutine

recursive subroutine __APPEND(fss_dispatch_fcn,__PREC) (self, x, fx)
    !*  DISPATCH_FCN returns the value of the wrapped function at a given point.
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_ss,__PREC)), intent(inout) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx

    real (PREC) :: fpx
        !*  First derivative; not used.

    if (associated(self%ptr_args)) then
        if (associated(self%fcn_args)) then
            call self%fcn_args (x, self%ptr_args, fx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt_args)) then
            call self%fcn_jac_opt_args (x, self%ptr_args, fx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_args)) then
            call self%fcn_jac_args (x, self%ptr_args, fx, fpx)
            self%nfev = self%nfev + 1
        end if
    else
        if (associated(self%fcn)) then
            call self%fcn (x, fx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt)) then
            call self%fcn_jac_opt (x, fx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac)) then
            call self%fcn_jac (x, fx, fpx)
            self%nfev = self%nfev + 1
        end if
    end if

end subroutine


recursive subroutine __APPEND(fss_dispatch_jac,__PREC) (self, x, fpx)
    !*  DISPATH_JAC returns the first derivative of the wrapped function
    !   at a given point.
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_ss,__PREC)), intent(inout) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fpx

    real (PREC) :: fx
        !   Function value, not used.

    if (associated(self%ptr_args)) then
        if (associated(self%jac_args)) then
            call self%jac_args (x, self%ptr_args, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt_args)) then
            call self%fcn_jac_opt_args (x, self%ptr_args, fpx=fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_args)) then
            call self%fcn_jac_args (x, self%ptr_args, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_args)) then
            call num_diff (self%fcn_args, x, self%ptr_args, fpx, eps=self%eps)
            self%nfev = self%nfev + 2
        end if
    else
        if (associated(self%jac)) then
            call self%jac (x, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt)) then
            call self%fcn_jac_opt (x, fpx=fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac)) then
            call self%fcn_jac (x, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn)) then
            call num_diff (self%fcn, x, fpx, eps=self%eps)
            self%nfev = self%nfev + 2
        end if
    end if

end subroutine


recursive subroutine __APPEND(fss_dispatch_fcn_jac,__PREC) (self, x, fx, fpx)
    !*  DISPATCH_FCN_JAC returns the function value and first derivative
    !   of the wrapped function at a given point.
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_ss,__PREC)), intent(inout) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx
    real (PREC), intent(out) :: fpx

    ! Check whether attribute to jointly compute function value and derivative
    ! is present
    if (associated(self%ptr_args)) then
        if (associated(self%fcn_jac_args)) then
            call self%fcn_jac_args (x, self%ptr_args, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt_args)) then
            call self%fcn_jac_opt_args (x, self%ptr_args, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_args)) then
            call self%fcn_args (x, self%ptr_args, fx)
            self%nfev = self%nfev + 1

            if (associated(self%jac_args)) then
                call self%jac_args (x, self%ptr_args, fpx)
                self%nfev = self%nfev + 1
            else
                call num_diff (self%fcn_args, x, self%ptr_args, fpx, fx, self%eps)
                self%nfev = self%nfev + 1
            end if
        end if
    else
        if (associated(self%fcn_jac)) then
            call self%fcn_jac (x, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt)) then
            call self%fcn_jac_opt (x, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn)) then
            call self%fcn (x, fx)
            self%nfev = self%nfev + 1

            if (associated(self%jac)) then
                call self%jac (x, fpx)
                self%nfev = self%nfev + 1
            else
                call num_diff (self%fcn, x, fpx, fx, eps=self%eps)
                self%nfev = self%nfev + 1
            end if
        end if
    end if
end subroutine


pure function __APPEND(fss_is_associated,__PREC) (self) result(res)
    !*  IS_ASSOCIATED returns TRUE if the wrapper object is associated 
    !   with any user-provided routine.
    type (__APPEND(fwrapper_ss,__PREC)), intent(in) :: self
    logical :: res 
    
    res = associated(self%fcn) .or. associated(self%fcn_args) &
        .or. associated(self%fcn_jac) .or. associated(self%fcn_jac_args) &
        .or. associated(self%fcn_jac_opt) .or. associated(self%fcn_jac_opt_args)
end function



!-------------------------------------------------------------------------------
! Routines for FWRAPPER_VS

subroutine __APPEND(fvs_init,__PREC) (self, fcn, jac, fcn_jac, fcn_jac_opt, &
        fcn_args, jac_args, fcn_jac_args, fcn_jac_opt_args, args, eps)
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_vs,__PREC)), intent(inout) :: self
    procedure (__APPEND(fvs_fcn,__PREC)), optional :: fcn
        !*  Pointer to function that returns the function value
    procedure (__APPEND(fvs_jac,__PREC)), optional :: jac
        !*  Pointer to function that returns the first derivative
    procedure (__APPEND(fvs_fcn_jac,__PREC)), optional :: fcn_jac
        !*  Pointer to function that returns the function value and the
        !   first derivative in a single call.
    procedure (__APPEND(fvs_fcn_jac_opt,__PREC)), optional :: fcn_jac_opt
        !*  Pointer to function that returns the function value and/or the
        !   first derivative in a single call.
    procedure (__APPEND(fvs_fcn_args,__PREC)), optional :: fcn_args
        !*  Pointer to function that returns the function value and accepts
        !   additional arguments.
    procedure (__APPEND(fvs_jac_args,__PREC)), optional :: jac_args
        !*  Pointer to function that returns the first derivative and accepts
        !   additional arguments.
    procedure (__APPEND(fvs_fcn_jac_args,__PREC)), optional :: fcn_jac_args
        !*  Pointer to function that returns the function value and the
        !   first derivative in a single call and accepts additional arguments.
    procedure (__APPEND(fvs_fcn_jac_opt_args,__PREC)), optional :: fcn_jac_opt_args
        !*  Pointer to function that returns the function value and/or the
        !   first derivative in a single call and accepts additional arguments.
    real (PREC), intent(in), dimension(:), optional, target :: args
        !*  Additional arguments that need to be passed to function.
    real (PREC), intent(in), optional :: eps
        !*  If present, sets the step size used when computing numerical
        !   derivative if no other way to obtain the derivative is avaiable.
        
    nullify (self%ptr_args)
    nullify (self%fcn, self%jac, self%fcn_jac, self%fcn_jac_opt)
    nullify (self%fcn_args, self%jac_args, self%fcn_jac_args, self%fcn_jac_opt_args)
    self%nfev = 0
    self%eps = sqrt(epsilon(0.0_PREC))

    if (present(args)) then
        if (present(fcn_args)) self%fcn_args => fcn_args
        if (present(jac_args)) self%jac_args => jac_args
        if (present(fcn_jac_args)) self%fcn_jac_args => fcn_jac_args
        if (present(fcn_jac_opt_args)) self%fcn_jac_opt_args => fcn_jac_opt_args
        self%ptr_args => args
    else
        if (present(fcn)) self%fcn => fcn
        if (present(jac)) self%jac => jac
        if (present(fcn_jac)) self%fcn_jac => fcn_jac
        if (present(fcn_jac_opt)) self%fcn_jac_opt => fcn_jac_opt
    end if
    
    if (present(eps)) then
        self%eps = eps 
    end if
end subroutine



recursive subroutine __APPEND(fvs_dispatch_fcn,__PREC) (self, x, fx)
    !*  DISPATCH_FCN returns the value of the wrapped function at a given point.
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_vs,__PREC)), intent(inout) :: self
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out) :: fx

    real (PREC), dimension(:), allocatable :: fpx
        !*  Local variable to store gradient, not returned to caller.

    if (associated(self%ptr_args)) then
        if (associated(self%fcn_args)) then
            call self%fcn_args (x, self%ptr_args, fx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt_args)) then
            call self%fcn_jac_opt_args (x, self%ptr_args, fx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_args)) then
            allocate (fpx(size(x)))
            call self%fcn_jac_args (x, self%ptr_args, fx, fpx)
            deallocate (fpx)
            self%nfev = self%nfev + 1
        end if
    else
        if (associated(self%fcn)) then
            call self%fcn (x, fx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt)) then
            call self%fcn_jac_opt (x, fx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac)) then
            allocate (fpx(size(x)))
            call self%fcn_jac (x, fx, fpx)
            self%nfev = self%nfev + 1
            deallocate (fpx)
        end if
    end if

end subroutine


recursive subroutine __APPEND(fvs_dispatch_jac,__PREC) (self, x, fpx)
    !*  DISPATH_JAC returns the first derivative of the wrapped function
    !   at a given point.
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_vs,__PREC)), intent(inout) :: self
    real (PREC), intent(in), dimension(:), contiguous  :: x
    real (PREC), intent(out), dimension(:), contiguous :: fpx

    real (PREC) :: fx
        !   Function value; not returned to caller

    if (associated(self%ptr_args)) then
        if (associated(self%jac_args)) then
            call self%jac_args (x, self%ptr_args, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt_args)) then
            call self%fcn_jac_opt_args (x, self%ptr_args, fpx=fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_args)) then
            call self%fcn_jac_args (x, self%ptr_args, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_args)) then
            call num_diff (self%fcn_args, x, self%ptr_args, fpx, eps=self%eps)
            self%nfev = self%nfev + size(x) + 1
        end if
    else
        if (associated(self%jac)) then
            call self%jac (x, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt)) then
            call self%fcn_jac_opt (x, fpx=fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac)) then
            call self%fcn_jac (x, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn)) then
            call num_diff (self%fcn, x, fpx, eps=self%eps)
            self%nfev = self%nfev + size(x) + 1
        end if
    end if

end subroutine


recursive subroutine __APPEND(fvs_dispatch_fcn_jac,__PREC) (self, x, fx, fpx)
    !*  DISPATCH_FCN_JAC returns the function value and first derivative
    !   of the wrapped function at a given point.
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_vs,__PREC)), intent(inout) :: self
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out) :: fx
    real (PREC), intent(out), dimension(:), contiguous :: fpx

    ! Check whether attribute to jointly compute function value and derivative
    ! is present
    if (associated(self%ptr_args)) then
        if (associated(self%fcn_jac_args)) then
            call self%fcn_jac_args (x, self%ptr_args, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt_args)) then
            call self%fcn_jac_opt_args (x, self%ptr_args, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_args)) then
            call self%fcn_args (x, self%ptr_args, fx)
            self%nfev = self%nfev + 1
            
            if (associated(self%jac_args)) then
                call self%jac_args (x, self%ptr_args, fpx)
                self%nfev = self%nfev + 1
            else
                call num_diff (self%fcn_args, x, self%ptr_args, fpx, fx, self%eps)
                self%nfev = self%nfev + size(x)
            end if
        end if
    else
        if (associated(self%fcn_jac)) then
            call self%fcn_jac (x, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt)) then
            call self%fcn_jac_opt (x, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn)) then
            call self%fcn (x, fx)
            self%nfev = self%nfev + 1
            
            if (associated(self%jac)) then
                call self%jac (x, fpx)
                self%nfev = self%nfev + 1
            else
                call num_diff (self%fcn, x, fpx, fx, eps=self%eps)
                self%nfev = self%nfev + size(x)
            end if
        end if
    end if
end subroutine

pure function __APPEND(fvs_is_associated,__PREC) (self) result(res)
    !*  IS_ASSOCIATED returns TRUE if the wrapper object is associated 
    !   with any user-provided routine.
    type (__APPEND(fwrapper_vs,__PREC)), intent(in) :: self
    logical :: res 
    
    res = associated(self%fcn) .or. associated(self%fcn_args) &
        .or. associated(self%fcn_jac) .or. associated(self%fcn_jac_args) &
        .or. associated(self%fcn_jac_opt) .or. associated(self%fcn_jac_opt_args)
end function

!-------------------------------------------------------------------------------
! Routines for FWRAPPER_VV

subroutine __APPEND(fvv_init,__PREC) (self, fcn, jac, fcn_jac, fcn_jac_opt, &
        fcn_args, jac_args, fcn_jac_args, fcn_jac_opt_args, args, eps)
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_vv,__PREC)), intent(inout) :: self
    procedure (__APPEND(fvv_fcn,__PREC)), optional :: fcn
        !*  Pointer to function that returns the function value
    procedure (__APPEND(fvv_jac,__PREC)), optional :: jac
        !*  Pointer to function that returns the first derivative
    procedure (__APPEND(fvv_fcn_jac,__PREC)), optional :: fcn_jac
        !*  Pointer to function that returns the function value and the
        !   first derivative in a single call.
    procedure (__APPEND(fvv_fcn_jac_opt,__PREC)), optional :: fcn_jac_opt
        !*  Pointer to function that returns the function value and/or the
        !   first derivative in a single call.
    procedure (__APPEND(fvv_fcn_args,__PREC)), optional :: fcn_args
        !*  Pointer to function that returns the function value and accepts
        !   additional arguments.
    procedure (__APPEND(fvv_jac_args,__PREC)), optional :: jac_args
        !*  Pointer to function that returns the first derivative and accepts
        !   additional arguments.
    procedure (__APPEND(fvv_fcn_jac_args,__PREC)), optional :: fcn_jac_args
        !*  Pointer to function that returns the function value and the
        !   first derivative in a single call and accepts additional arguments.
    procedure (__APPEND(fvv_fcn_jac_opt_args,__PREC)), optional :: fcn_jac_opt_args
        !*  Pointer to function that returns the function value and/or the
        !   first derivative in a single call and accepts additional arguments.
    real (PREC), intent(in), dimension(:), optional, target :: args
        !*  Additional arguments that need to be passed to function.
    real (PREC), intent(in), optional :: eps
        !*  If present, sets the step size used when computing numerical
        !   derivative if no other way to obtain the derivative is avaiable.
        
    nullify (self%ptr_args)
    nullify (self%fcn, self%jac, self%fcn_jac, self%fcn_jac_opt)
    nullify (self%fcn_args, self%jac_args, self%fcn_jac_args, self%fcn_jac_opt_args)
    self%nfev = 0
    self%eps = sqrt(epsilon(0.0_PREC))

    if (present(args)) then
        if (present(fcn_args)) self%fcn_args => fcn_args
        if (present(jac_args)) self%jac_args => jac_args
        if (present(fcn_jac_args)) self%fcn_jac_args => fcn_jac_args
        if (present(fcn_jac_opt_args)) self%fcn_jac_opt_args => fcn_jac_opt_args
        self%ptr_args => args
    else
        if (present(fcn)) self%fcn => fcn
        if (present(jac)) self%jac => jac
        if (present(fcn_jac)) self%fcn_jac => fcn_jac
        if (present(fcn_jac_opt)) self%fcn_jac_opt => fcn_jac_opt
    end if
    
    if (present(eps)) then
        self%eps = eps 
    end if
end subroutine



recursive subroutine __APPEND(fvv_dispatch_fcn,__PREC) (self, x, fx)
    !*  DISPATCH_FCN returns the value of the wrapped function at a given point.
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_vv,__PREC)), intent(inout) :: self
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx

    real (PREC), dimension(:,:), allocatable :: fpx
        !*  Local array used to store Jacobian; not used by caller.

    if (associated(self%ptr_args)) then
        if (associated(self%fcn_args)) then
            call self%fcn_args (x, self%ptr_args, fx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt_args)) then
            call self%fcn_jac_opt_args (x, self%ptr_args, fx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_args)) then
            allocate (fpx(size(fx),size(x)))
            call self%fcn_jac_args (x, self%ptr_args, fx, fpx)
            self%nfev = self%nfev + 1
            deallocate (fpx)
        end if
    else
        if (associated(self%fcn)) then
            call self%fcn (x, fx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt)) then
            call self%fcn_jac_opt (x, fx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac)) then
            allocate (fpx(size(fx),size(x)))
            call self%fcn_jac (x, fx, fpx)
            self%nfev = self%nfev + 1
            deallocate (fpx)
        end if
    end if

end subroutine


recursive subroutine __APPEND(fvv_dispatch_jac,__PREC) (self, x, fpx)
    !*  DISPATH_JAC returns the first derivative of the wrapped function
    !   at a given point.
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_vv,__PREC)), intent(inout) :: self
    real (PREC), intent(in), dimension(:), contiguous  :: x
    real (PREC), intent(out), dimension(:,:), contiguous :: fpx

    real (PREC), dimension(:), allocatable :: fx
        !   Local array to store function value; not returned to caller.

    if (associated(self%ptr_args)) then
        if (associated(self%jac_args)) then
            call self%jac_args (x, self%ptr_args, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt_args)) then
            call self%fcn_jac_opt_args (x, self%ptr_args, fpx=fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_args)) then
            allocate (fx(size(fpx,1)))
            call self%fcn_jac_args (x, self%ptr_args, fx, fpx)
            self%nfev = self%nfev + 1
            deallocate (fx)
        else if (associated(self%fcn_args)) then
            call num_diff (self%fcn_args, x, self%ptr_args, fpx, eps=self%eps)
            self%nfev = self%nfev + size(x) + 1
        end if
    else
        if (associated(self%jac)) then
            call self%jac (x, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt)) then
            call self%fcn_jac_opt (x, fpx=fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac)) then
            allocate (fx(size(fpx,1)))
            call self%fcn_jac (x, fx, fpx)
            self%nfev = self%nfev + 1
            deallocate (fx)
        else if (associated(self%fcn)) then
            call num_diff (self%fcn, x, fpx, eps=self%eps)
            self%nfev = self%nfev + size(x) + 1
        end if
    end if

end subroutine


recursive subroutine __APPEND(fvv_dispatch_fcn_jac,__PREC) (self, x, fx, fpx)
    !*  DISPATCH_FCN_JAC returns the function value and first derivative
    !   of the wrapped function at a given point.
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_vv,__PREC)), intent(inout) :: self
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous :: fx
    real (PREC), intent(out), dimension(:,:), contiguous :: fpx

    ! Check whether attribute to jointly compute function value and derivative
    ! is present
    if (associated(self%ptr_args)) then
        if (associated(self%fcn_jac_args)) then
            call self%fcn_jac_args (x, self%ptr_args, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt_args)) then
            call self%fcn_jac_opt_args (x, self%ptr_args, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_args)) then
            call self%fcn_args (x, self%ptr_args, fx)
            self%nfev = self%nfev + 1
            
            if (associated(self%jac_args)) then
                call self%jac_args (x, self%ptr_args, fpx)
                self%nfev = self%nfev + 1
            else
                call num_diff (self%fcn_args, x, self%ptr_args, fpx, fx, self%eps)
                self%nfev = self%nfev + size(x)
            end if
        end if
    else
        if (associated(self%fcn_jac)) then
            call self%fcn_jac (x, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt)) then
            call self%fcn_jac_opt (x, fx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn)) then
            call self%fcn (x, fx)
            self%nfev = self%nfev + 1
            
            if (associated(self%jac)) then
                call self%jac (x, fpx)
                self%nfev = self%nfev + 1
            else
                call num_diff (self%fcn, x, fpx, fx, eps=self%eps)
                self%nfev = self%nfev + size(x)
            end if
        end if
    end if
end subroutine


pure function __APPEND(fvv_is_associated,__PREC) (self) result(res)
    !*  IS_ASSOCIATED returns TRUE if the wrapper object is associated 
    !   with any user-provided routine.
    type (__APPEND(fwrapper_vv,__PREC)), intent(in) :: self
    logical :: res 
    
    res = associated(self%fcn) .or. associated(self%fcn_args) &
        .or. associated(self%fcn_jac) .or. associated(self%fcn_jac_args) &
        .or. associated(self%fcn_jac_opt) .or. associated(self%fcn_jac_opt_args)
end function
