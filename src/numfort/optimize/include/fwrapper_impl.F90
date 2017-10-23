
!-------------------------------------------------------------------------------
! Routines for FWRAPPER_SS

subroutine __APPEND(fwrapper_ss_wrap,__PREC) (self, fcn, fcn_args)
    type (__APPEND(fwrapper_ss,__PREC)), intent(in out) :: self
    procedure (__APPEND(fss,__PREC)), optional :: fcn
    procedure (__APPEND(fss_args,__PREC)), optional :: fcn_args

    if (present(fcn_args)) then
        self%fcn_args => fcn_args
    else if (present(fcn)) then
        self%fcn => fcn
    end if
end subroutine

subroutine __APPEND(fwrapper_ss_dispatch,__PREC) (self, x, fx, args)
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_ss,__PREC)), intent(in out) :: self
    real (PREC), intent(in) :: x
    real (PREC), intent(out) :: fx
    real (PREC), intent(in), dimension(:), optional :: args

    if (associated(self%fcn_args)) then
        call self%fcn_args (x, fx, args)
        self%nfev = self%nfev + 1
    else if (associated(self%fcn)) then
        call self%fcn (x, fx)
        self%nfev = self%nfev + 1
    end if

end subroutine

pure function __APPEND(fwrapper_ss_is_present,__PREC) (self) result(res)
    type (__APPEND(fwrapper_ss,__PREC)), intent(in) :: self
    logical :: res

    res = associated(self%fcn_args) .or. associated(self%fcn)
end function

!-------------------------------------------------------------------------------
! Routines for FWRAPPER_VS_JAC

subroutine __APPEND(fwrapper_vs_jac_wrap,__PREC) (self, fcn, fcn_args)
    type (__APPEND(fwrapper_vs_jac,__PREC)), intent(in out) :: self
    procedure (__APPEND(fvs_jac,__PREC)), optional :: fcn
    procedure (__APPEND(fvs_jac_args,__PREC)), optional :: fcn_args

    if (present(fcn_args)) then
        self%fcn_args => fcn_args
    else if (present(fcn)) then
        self%fcn => fcn
    end if
end subroutine

subroutine __APPEND(fwrapper_vs_jac_dispatch,__PREC) (self, x, fx, fpx, args)
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_vs_jac,__PREC)), intent(in out) :: self
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), dimension(:), contiguous, optional :: fpx
    real (PREC), intent(in), dimension(:), optional :: args

    if (associated(self%fcn_args)) then
        call self%fcn_args (x, fx, fpx, args)
        self%nfev = self%nfev + 1
    else if (associated(self%fcn)) then
        call self%fcn (x, fx, fpx)
        self%nfev = self%nfev + 1
    end if

end subroutine


pure function __APPEND(fwrapper_vs_jac_is_present,__PREC) (self) result(res)
    type (__APPEND(fwrapper_vs_jac,__PREC)), intent(in) :: self
    logical :: res

    res = associated(self%fcn_args) .or. associated(self%fcn)
end function

!-------------------------------------------------------------------------------
! Routines for FWRAPPER_VV_JAC

subroutine __APPEND(fwrapper_vv_jac_wrap,__PREC) (self, fcn, fcn_args)
    type (__APPEND(fwrapper_vv_jac,__PREC)), intent(in out) :: self
    procedure (__APPEND(fvv_jac,__PREC)), optional :: fcn
    procedure (__APPEND(fvv_jac_args,__PREC)), optional :: fcn_args

    if (present(fcn_args)) then
        self%fcn_args => fcn_args
    else if (present(fcn)) then
        self%fcn => fcn
    end if
end subroutine


subroutine __APPEND(fwrapper_vv_jac_dispatch,__PREC) (self, x, fx, fpx, args)
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_vv_jac,__PREC)), intent(in out) :: self
    real (PREC), intent(in), dimension(:), contiguous :: x
    real (PREC), intent(out), dimension(:), contiguous, optional :: fx
    real (PREC), intent(out), dimension(:,:), contiguous, optional :: fpx
    real (PREC), intent(in), dimension(:), optional :: args

    if (associated(self%fcn_args)) then
        call self%fcn_args (x, fx, fpx, args)
        self%nfev = self%nfev + 1
    else if (associated(self%fcn)) then
        call self%fcn (x, fx, fpx)
        self%nfev = self%nfev + 1
    end if

end subroutine


pure function __APPEND(fwrapper_vv_jac_is_present,__PREC) (self) result(res)
    type (__APPEND(fwrapper_vv_jac,__PREC)), intent(in) :: self
    logical :: res

    res = associated(self%fcn_args) .or. associated(self%fcn)
end function

