
!-------------------------------------------------------------------------------
! Routines for FWRAPPER_SS

subroutine fss_init (self, fcn, jac, fcn_jac, fcn_jac_opt, &
        fcn_args, jac_args, fcn_jac_args, fcn_jac_opt_args, args, eps, reps)
    !*  FSS_INIT intializes a wrapper found a scalar-valued function f:R->R.
    type (fwrapper_ss), intent(inout) :: self
    procedure (fss), optional :: fcn
        !*  Pointer to function that returns the function value
    procedure (fss), optional :: jac
        !*  Pointer to function that returns the first derivative
    procedure (fss_fcn_jac), optional :: fcn_jac
        !*  Pointer to function that returns the function value and the
        !   first derivative in a single call.
    procedure (fss_fcn_jac_opt), optional :: fcn_jac_opt
        !*  Pointer to function that returns the function value and/or the
        !   first derivative in a single call.
    procedure (fss_args), optional :: fcn_args
        !*  Pointer to function that returns the function value and accepts
        !   additional arguments.
    procedure (fss_args), optional :: jac_args
        !*  Pointer to function that returns the first derivative and accepts
        !   additional arguments.
    procedure (fss_fcn_jac_args), optional :: fcn_jac_args
        !*  Pointer to function that returns the function value and the
        !   first derivative in a single call and accepts additional arguments.
    procedure (fss_fcn_jac_opt_args), optional :: fcn_jac_opt_args
        !*  Pointer to function that returns the function value and/or the
        !   first derivative in a single call and accepts additional arguments.
    class (args_data), intent(in), optional, target :: args
        !*  Additional arguments that need to be passed to function.
    real (PREC), intent(in), optional :: eps
        !*  If present, sets the step size used when computing numerical
        !   derivatives if no other way to obtain the derivative is avaiable.
        !   Note: REPS is ignored if EPS is present.
    real (PREC), intent(in), optional :: reps
        !*  If present, sets the relative step size used when computing numerical
        !   derivatives by forward differencing.
        !   Note: REPS is ignored if EPS is present.

    nullify (self%ptr_args)
    nullify (self%fcn, self%jac, self%fcn_jac, self%fcn_jac_opt)
    nullify (self%fcn_args, self%jac_args, self%fcn_jac_args, self%fcn_jac_opt_args)
    self%nfev = 0
    self%eps = sqrt(epsilon(0.0_PREC))
    self%reps = sqrt(epsilon(0.0_PREC))
    self%rel_diff = .false.

    if (present(args)) then
        if (present(fcn_args)) self%fcn_args => fcn_args
        if (present(jac_args)) self%jac_args => jac_args
        if (present(fcn_jac_args)) self%fcn_jac_args => fcn_jac_args
        if (present(fcn_jac_opt_args)) self%fcn_jac_opt_args => fcn_jac_opt_args
        self%ptr_args => args

        ! Set flag whether numerical differentiation will be performed when
        ! Jacobian needs to be computed.
        self%num_diff = .not. (present(jac_args) .or. present(fcn_jac_args) &
            .or. present(fcn_jac_opt_args))
    else
        if (present(fcn)) self%fcn => fcn
        if (present(jac)) self%jac => jac
        if (present(fcn_jac)) self%fcn_jac => fcn_jac
        if (present(fcn_jac_opt)) self%fcn_jac_opt => fcn_jac_opt

        ! Set flag whether numerical differentiation will be performed when
        ! Jacobian needs to be computed.
        self%num_diff = .not. (present(jac) .or. present(fcn_jac) .or. &
            present(fcn_jac_opt))
    end if

    if (present(eps)) then
        self%eps = eps
    end if

    if (present(reps) .and. .not. present(eps)) then
        self%rel_diff = .true.
        self%reps = reps
    end if

end subroutine



recursive subroutine fss_dispatch_fcn (self, x, fx)
    !*  DISPATCH_FCN returns the value of the wrapped function at a given point.
    type (fwrapper_ss), intent(inout) :: self
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



recursive subroutine fss_dispatch_jac (self, x, fpx, fx)
    !*  DISPATH_JAC returns the first derivative of the wrapped function
    !   at a given point.
    type (fwrapper_ss), intent(inout) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fpx
    real (PREC), intent(in), optional :: fx
        !*  Function value at point X. If present, numerical differentiation
        !   takes one less function evaluation to compute. Ignored if
        !   derivative is computed via a user-provided function.

    real (PREC) :: lfx
        !   Function value, used to store return value from joint function/
        !   derivative routine.

    if (associated(self%ptr_args)) then
        if (associated(self%jac_args)) then
            call self%jac_args (x, self%ptr_args, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt_args)) then
            call self%fcn_jac_opt_args (x, self%ptr_args, fpx=fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_args)) then
            call self%fcn_jac_args (x, self%ptr_args, lfx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_args)) then
            if (.not. self%rel_diff) then
                call num_diff (self%fcn_args, x, self%ptr_args, fpx, fx, self%eps)
            else
                call num_diff (self%fcn_args, x, self%ptr_args, fpx, fx, &
                    reps=self%reps)
            end if
            ! Number of function evaluations requierd to compute f(X+eps)
            self%nfev = self%nfev + 1
            ! Add function evaluated required to obtain f(X)
            if (.not. present(fx)) self%nfev = self%nfev + 1
        end if
    else
        if (associated(self%jac)) then
            call self%jac (x, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt)) then
            call self%fcn_jac_opt (x, fpx=fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac)) then
            call self%fcn_jac (x, lfx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn)) then
            if (.not. self%rel_diff) then
                call num_diff (self%fcn, x, fpx, fx, self%eps)
            else
                call num_diff (self%fcn, x, fpx, fx, reps=self%reps)
            end if
            ! Number of function evaluations requierd to compute f(X+eps)
            self%nfev = self%nfev + 1
            ! Add function evaluated required to obtain f(X)
            if (.not. present(fx)) self%nfev = self%nfev + 1
        end if
    end if

end subroutine



recursive subroutine fss_dispatch_fcn_jac (self, x, fx, fpx)
    !*  DISPATCH_FCN_JAC returns the function value and first derivative
    !   of the wrapped function at a given point.
    type (fwrapper_ss), intent(inout) :: self
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
                if (.not. self%rel_diff) then
                    call num_diff (self%fcn_args, x, self%ptr_args, fpx, fx, self%eps)
                else
                    call num_diff (self%fcn_args, x, self%ptr_args, fpx, fx, &
                        reps=self%reps)
                end if
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
                if (.not. self%rel_diff) then
                    call num_diff (self%fcn, x, fpx, fx, eps=self%eps)
                else
                    call num_diff (self%fcn, x, fpx, fx, reps=self%reps)
                end if
                self%nfev = self%nfev + 1
            end if
        end if
    end if
end subroutine



pure function fss_is_associated (self) result(res)
    !*  IS_ASSOCIATED returns TRUE if the wrapper object is associated 
    !   with any user-provided routine.
    type (fwrapper_ss), intent(in) :: self
    logical :: res 
    
    res = associated(self%fcn) .or. associated(self%fcn_args) &
        .or. associated(self%fcn_jac) .or. associated(self%fcn_jac_args) &
        .or. associated(self%fcn_jac_opt) .or. associated(self%fcn_jac_opt_args)
end function



!-------------------------------------------------------------------------------
! Routines for FWRAPPER_VS

subroutine fvs_init (self, fcn, jac, fcn_jac, &
        fcn_jac_opt, fcn_args, jac_args, fcn_jac_args, fcn_jac_opt_args, &
        args, eps, reps)
    type (fwrapper_vs), intent(inout) :: self
    procedure (fvs_fcn), optional :: fcn
        !*  Pointer to function that returns the function value
    procedure (fvs_jac), optional :: jac
        !*  Pointer to function that returns the first derivative
    procedure (fvs_fcn_jac), optional :: fcn_jac
        !*  Pointer to function that returns the function value and the
        !   first derivative in a single call.
    procedure (fvs_fcn_jac_opt), optional :: fcn_jac_opt
        !*  Pointer to function that returns the function value and/or the
        !   first derivative in a single call.
    procedure (fvs_fcn_args), optional :: fcn_args
        !*  Pointer to function that returns the function value and accepts
        !   additional arguments.
    procedure (fvs_jac_args), optional :: jac_args
        !*  Pointer to function that returns the first derivative and accepts
        !   additional arguments.
    procedure (fvs_fcn_jac_args), optional :: fcn_jac_args
        !*  Pointer to function that returns the function value and the
        !   first derivative in a single call and accepts additional arguments.
    procedure (fvs_fcn_jac_opt_args), optional :: fcn_jac_opt_args
        !*  Pointer to function that returns the function value and/or the
        !   first derivative in a single call and accepts additional arguments.
    class (args_data), intent(in), optional, target :: args
        !*  Additional arguments that need to be passed to function.
    real (PREC), intent(in), optional :: eps
        !*  If present, sets the step size used when computing numerical
        !   derivative if no other way to obtain the derivative is avaiable.
        !   Note: REPS is ignored if EPS is present.
    real (PREC), intent(in), optional :: reps
        !*  If present, sets the relative step size used when computing numerical
        !   derivatives by forward differencing.
        !   Note: REPS is ignored if EPS is present.

    nullify (self%ptr_args)
    nullify (self%fcn, self%jac, self%fcn_jac, self%fcn_jac_opt)
    nullify (self%fcn_args, self%jac_args, self%fcn_jac_args, self%fcn_jac_opt_args)
    self%nfev = 0
    self%eps = sqrt(epsilon(0.0_PREC))
    self%reps = sqrt(epsilon(0.0_PREC))
    self%rel_diff = .false.

    if (present(args)) then
        if (present(fcn_args)) self%fcn_args => fcn_args
        if (present(jac_args)) self%jac_args => jac_args
        if (present(fcn_jac_args)) self%fcn_jac_args => fcn_jac_args
        if (present(fcn_jac_opt_args)) self%fcn_jac_opt_args => fcn_jac_opt_args
        self%ptr_args => args

        ! Set flag whether numerical differentiation will be performed when
        ! Jacobian needs to be computed.
        self%num_diff = .not. (present(jac_args) .or. present(fcn_jac_args) &
            .or. present(fcn_jac_opt_args))
    else
        if (present(fcn)) self%fcn => fcn
        if (present(jac)) self%jac => jac
        if (present(fcn_jac)) self%fcn_jac => fcn_jac
        if (present(fcn_jac_opt)) self%fcn_jac_opt => fcn_jac_opt

        ! Set flag whether numerical differentiation will be performed when
        ! Jacobian needs to be computed.
        self%num_diff = .not. (present(jac) .or. present(fcn_jac) .or. &
            present(fcn_jac_opt))
    end if
    
    if (present(eps)) then
        self%eps = eps 
    end if

    if (present(reps) .and. .not. present(eps)) then
        self%rel_diff = .true.
        self%reps = reps
    end if
end subroutine



recursive subroutine fvs_dispatch_fcn (self, x, fx)
    !*  DISPATCH_FCN returns the value of the wrapped function at a given point.
    type (fwrapper_vs), intent(inout) :: self
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



recursive subroutine fvs_dispatch_jac (self, x, fpx, fx)
    !*  DISPATH_JAC returns the first derivative of the wrapped function
    !   at a given point.
    type (fwrapper_vs), intent(inout) :: self
    real (PREC), intent(in), dimension(:), contiguous  :: x
    real (PREC), intent(out), dimension(:), contiguous :: fpx
    real (PREC), intent(in), optional :: fx
        !*  Function value at point X. If present, numerical differentiation
        !   takes one less function evaluation to compute. Ignored if
        !   derivative is computed via a user-provided function.

    real (PREC) :: lfx
        !   Function value, used to store return value from joint function/
        !   derivative routine.

    if (associated(self%ptr_args)) then
        if (associated(self%jac_args)) then
            call self%jac_args (x, self%ptr_args, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt_args)) then
            call self%fcn_jac_opt_args (x, self%ptr_args, fpx=fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_args)) then
            call self%fcn_jac_args (x, self%ptr_args, lfx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_args)) then
            if (.not. self%rel_diff) then
                call num_diff (self%fcn_args, x, self%ptr_args, fpx, fx, self%eps)
            else
                call num_diff (self%fcn_args, x, self%ptr_args, fpx, fx, &
                    reps=self%reps)
            end if
            ! Number of function evaluations requierd to compute f(X+eps)
            self%nfev = self%nfev + size(x)
            ! Add function evaluated required to obtain f(X)
            if (.not. present(fx)) self%nfev = self%nfev + 1
        end if
    else
        if (associated(self%jac)) then
            call self%jac (x, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt)) then
            call self%fcn_jac_opt (x, fpx=fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac)) then
            call self%fcn_jac (x, lfx, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn)) then
            if (.not. self%rel_diff) then
                call num_diff (self%fcn, x, fpx, fx, self%eps)
            else
                call num_diff (self%fcn, x, fpx, fx, reps=self%reps)
            end if
            ! Number of function evaluations requierd to compute f(X+eps)
            self%nfev = self%nfev + size(x)
            ! Add function evaluated required to obtain f(X)
            if (.not. present(fx)) self%nfev = self%nfev + 1
        end if
    end if

end subroutine



recursive subroutine fvs_dispatch_fcn_jac (self, x, fx, fpx)
    !*  DISPATCH_FCN_JAC returns the function value and first derivative
    !   of the wrapped function at a given point.
    type (fwrapper_vs), intent(inout) :: self
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
                if (.not. self%rel_diff) then
                    call num_diff (self%fcn_args, x, self%ptr_args, fpx, fx, self%eps)
                else
                    call num_diff (self%fcn_args, x, self%ptr_args, fpx, fx, &
                        reps=self%reps)
                end if
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
                if (.not. self%rel_diff) then
                    call num_diff (self%fcn, x, fpx, fx, eps=self%eps)
                else
                    call num_diff (self%fcn, x, fpx, fx, reps=self%reps)
                end if
                self%nfev = self%nfev + size(x)
            end if
        end if
    end if
end subroutine



pure function fvs_is_associated (self) result(res)
    !*  IS_ASSOCIATED returns TRUE if the wrapper object is associated 
    !   with any user-provided routine.
    type (fwrapper_vs), intent(in) :: self
    logical :: res 
    
    res = associated(self%fcn) .or. associated(self%fcn_args) &
        .or. associated(self%fcn_jac) .or. associated(self%fcn_jac_args) &
        .or. associated(self%fcn_jac_opt) .or. associated(self%fcn_jac_opt_args)
end function

!-------------------------------------------------------------------------------
! Routines for FWRAPPER_VV

subroutine fvv_init (self, fcn, jac, fcn_jac, fcn_jac_opt, &
        fcn_args, jac_args, fcn_jac_args, fcn_jac_opt_args, args, eps, reps)
    type (fwrapper_vv), intent(inout) :: self
    procedure (fvv_fcn), optional :: fcn
        !*  Pointer to function that returns the function value
    procedure (fvv_jac), optional :: jac
        !*  Pointer to function that returns the first derivative
    procedure (fvv_fcn_jac), optional :: fcn_jac
        !*  Pointer to function that returns the function value and the
        !   first derivative in a single call.
    procedure (fvv_fcn_jac_opt), optional :: fcn_jac_opt
        !*  Pointer to function that returns the function value and/or the
        !   first derivative in a single call.
    procedure (fvv_fcn_args), optional :: fcn_args
        !*  Pointer to function that returns the function value and accepts
        !   additional arguments.
    procedure (fvv_jac_args), optional :: jac_args
        !*  Pointer to function that returns the first derivative and accepts
        !   additional arguments.
    procedure (fvv_fcn_jac_args), optional :: fcn_jac_args
        !*  Pointer to function that returns the function value and the
        !   first derivative in a single call and accepts additional arguments.
    procedure (fvv_fcn_jac_opt_args), optional :: fcn_jac_opt_args
        !*  Pointer to function that returns the function value and/or the
        !   first derivative in a single call and accepts additional arguments.
    class (args_data), intent(in), optional, target :: args
        !*  Additional arguments that need to be passed to function.
    real (PREC), intent(in), optional :: eps
        !*  If present, sets the step size used when computing numerical
        !   derivative if no other way to obtain the derivative is avaiable.
        !   Note: REPS is ignored if EPS is present.
    real (PREC), intent(in), optional :: reps
        !*  If present, sets the relative step size used when computing numerical
        !   derivatives by forward differencing.
        !   Note: REPS is ignored if EPS is present.

    nullify (self%ptr_args)
    nullify (self%fcn, self%jac, self%fcn_jac, self%fcn_jac_opt)
    nullify (self%fcn_args, self%jac_args, self%fcn_jac_args, self%fcn_jac_opt_args)
    self%nfev = 0
    self%eps = sqrt(epsilon(0.0_PREC))
    self%reps = sqrt(epsilon(0.0_PREC))
    self%rel_diff = .false.

    if (present(args)) then
        if (present(fcn_args)) self%fcn_args => fcn_args
        if (present(jac_args)) self%jac_args => jac_args
        if (present(fcn_jac_args)) self%fcn_jac_args => fcn_jac_args
        if (present(fcn_jac_opt_args)) self%fcn_jac_opt_args => fcn_jac_opt_args
        self%ptr_args => args

        ! Set flag whether numerical differentiation will be performed when
        ! Jacobian needs to be computed.
        self%num_diff = .not. (present(jac_args) .or. present(fcn_jac_args) &
            .or. present(fcn_jac_opt_args))
    else
        if (present(fcn)) self%fcn => fcn
        if (present(jac)) self%jac => jac
        if (present(fcn_jac)) self%fcn_jac => fcn_jac
        if (present(fcn_jac_opt)) self%fcn_jac_opt => fcn_jac_opt

        ! Set flag whether numerical differentiation will be performed when
        ! Jacobian needs to be computed.
        self%num_diff = .not. (present(jac) .or. present(fcn_jac) .or. &
            present(fcn_jac_opt))
    end if
    
    if (present(eps)) then
        self%eps = eps 
    end if

    if (present(reps) .and. .not. present(eps)) then
        self%rel_diff = .true.
        self%reps = reps
    end if
end subroutine



recursive subroutine fvv_dispatch_fcn (self, x, fx)
    !*  DISPATCH_FCN returns the value of the wrapped function at a given point.
    type (fwrapper_vv), intent(inout) :: self
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



recursive subroutine fvv_dispatch_jac (self, x, fpx, fx)
    !*  DISPATH_JAC returns the first derivative of the wrapped function
    !   at a given point.
    type (fwrapper_vv), intent(inout) :: self
    real (PREC), intent(in), dimension(:), contiguous  :: x
    real (PREC), intent(out), dimension(:,:), contiguous :: fpx
    real (PREC), intent(in), dimension(:), optional, contiguous :: fx
        !*  Function value at point X. If present, numerical differentiation
        !   takes one less function evaluation to compute. Ignored if
        !   derivative is computed via a user-provided function.

    real (PREC), dimension(:), allocatable :: lfx
        !   Function value, used to store return value from joint function/
        !   derivative routine.

    if (associated(self%ptr_args)) then
        if (associated(self%jac_args)) then
            call self%jac_args (x, self%ptr_args, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt_args)) then
            call self%fcn_jac_opt_args (x, self%ptr_args, fpx=fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_args)) then
            allocate (lfx(size(fpx,1)))
            call self%fcn_jac_args (x, self%ptr_args, lfx, fpx)
            self%nfev = self%nfev + 1
            deallocate (lfx)
        else if (associated(self%fcn_args)) then
            if (.not. self%rel_diff) then
                call num_diff (self%fcn_args, x, self%ptr_args, fpx, fx, self%eps)
            else
                call num_diff (self%fcn_args, x, self%ptr_args, fpx, fx, &
                    reps=self%reps)
            end if
            ! Number of function evaluations requierd to compute f(X+eps)
            self%nfev = self%nfev + size(x)
            ! Add function evaluated required to obtain f(X)
            if (.not. present(fx)) self%nfev = self%nfev + 1
        end if
    else
        if (associated(self%jac)) then
            call self%jac (x, fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac_opt)) then
            call self%fcn_jac_opt (x, fpx=fpx)
            self%nfev = self%nfev + 1
        else if (associated(self%fcn_jac)) then
            allocate (lfx(size(fpx,1)))
            call self%fcn_jac (x, lfx, fpx)
            self%nfev = self%nfev + 1
            deallocate (lfx)
        else if (associated(self%fcn)) then
            if (.not. self%rel_diff) then
                call num_diff (self%fcn, x, fpx, fx, self%eps)
            else
                call num_diff (self%fcn, x, fpx, fx, reps=self%reps)
            end if
            ! Number of function evaluations requierd to compute f(X+eps)
            self%nfev = self%nfev + size(x)
            ! Add function evaluated required to obtain f(X)
            if (.not. present(fx)) self%nfev = self%nfev + 1
        end if
    end if

end subroutine



recursive subroutine fvv_dispatch_fcn_jac (self, x, fx, fpx)
    !*  DISPATCH_FCN_JAC returns the function value and first derivative
    !   of the wrapped function at a given point.
    type (fwrapper_vv), intent(inout) :: self
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
                if (.not. self%rel_diff) then
                    call num_diff (self%fcn_args, x, self%ptr_args, fpx, fx, self%eps)
                else
                    call num_diff (self%fcn_args, x, self%ptr_args, fpx, fx, &
                        reps=self%reps)
                end if
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
                if (.not. self%rel_diff) then
                    call num_diff (self%fcn, x, fpx, fx, eps=self%eps)
                else
                    call num_diff (self%fcn, x, fpx, fx, reps=self%reps)
                end if
                self%nfev = self%nfev + size(x)
            end if
        end if
    end if
end subroutine



pure function fvv_is_associated (self) result(res)
    !*  IS_ASSOCIATED returns TRUE if the wrapper object is associated 
    !   with any user-provided routine.
    type (fwrapper_vv), intent(in) :: self
    logical :: res 
    
    res = associated(self%fcn) .or. associated(self%fcn_args) &
        .or. associated(self%fcn_jac) .or. associated(self%fcn_jac_args) &
        .or. associated(self%fcn_jac_opt) .or. associated(self%fcn_jac_opt_args)
end function
