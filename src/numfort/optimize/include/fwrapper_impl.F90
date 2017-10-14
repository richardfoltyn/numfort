
subroutine __APPEND(fwrapper_vec_scalar_wrap,__PREC) (self, fcn, fcn_args)
    type (__APPEND(fwrapper_vec_scalar,__PREC)), intent(in out) :: self
    procedure (__APPEND(fvec_scalar,__PREC)), optional :: fcn
    procedure (__APPEND(fvec_scalar_args,__PREC)), optional :: fcn_args

    if (present(fcn)) then
        self%fcn => fcn
    else if (present(fcn_args)) then
        self%fcn_args => fcn_args
    end if
end subroutine


subroutine __APPEND(fwrapper_vec_scalar_dispatch,__PREC) (self, x, fx, fpx, args)
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_vec_scalar,__PREC)), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(out), optional :: fx
    real (PREC), intent(out), dimension(:), optional :: fpx
    real (PREC), intent(in), dimension(:), optional :: args

    if (associated(self%fcn_args)) then
        call self%fcn_args (x, fx, fpx, args)
    else if (associated(self%fcn)) then
        call self%fcn (x, fx, fpx)
    end if

end subroutine


pure function __APPEND(fwrapper_vec_scalar_is_present,__PREC) (self) result(res)
    type (__APPEND(fwrapper_vec_scalar,__PREC)), intent(in) :: self
    logical :: res

    res = associated(self%fcn_args) .or. associated(self%fcn)
end function


subroutine __APPEND(fwrapper_vec_vec_wrap,__PREC) (self, fcn, fcn_args)
    type (__APPEND(fwrapper_vec_vec,__PREC)), intent(in out) :: self
    procedure (__APPEND(fvec_vec,__PREC)), optional :: fcn
    procedure (__APPEND(fvec_vec_args,__PREC)), optional :: fcn_args

    if (present(fcn)) then
        self%fcn => fcn
    else if (present(fcn_args)) then
        self%fcn_args => fcn_args
    end if
end subroutine


subroutine __APPEND(fwrapper_vec_vec_dispatch,__PREC) (self, x, fx, fpx, args)
    integer, parameter :: PREC = __PREC
    type (__APPEND(fwrapper_vec_vec,__PREC)), intent(in) :: self
    real (PREC), intent(in), dimension(:) :: x
    real (PREC), intent(out), dimension(:), optional :: fx
    real (PREC), intent(out), dimension(:,:), optional :: fpx
    real (PREC), intent(in), dimension(:), optional :: args

    if (associated(self%fcn_args)) then
        call self%fcn_args (x, fx, fpx, args)
    else if (associated(self%fcn)) then
        call self%fcn (x, fx, fpx)
    end if

end subroutine


pure function __APPEND(fwrapper_vec_vec_is_present,__PREC) (self) result(res)
    type (__APPEND(fwrapper_vec_vec,__PREC)), intent(in) :: self
    logical :: res

    res = associated(self%fcn_args) .or. associated(self%fcn)
end function

