abstract interface

    subroutine __APPEND(fcn_der2_args,__PREC) (x, args, fx, fpx, fppx)
        import __PREC
        integer, parameter :: PREC = __PREC
        real (PREC), intent(in) :: x
        real (PREC), intent(in out), dimension(:) :: args
        real (PREC), intent(out) :: fx, fpx, fppx
    end subroutine

    subroutine __APPEND(fcn_der2,__PREC) (x, fx, fpx, fppx)
        import __PREC
        integer, parameter :: PREC = __PREC
        real (PREC), intent(in) :: x
        real (PREC), intent(out) :: fx, fpx, fppx
    end subroutine

end interface

interface root_newton
    procedure __APPEND(root_newton_args,__PREC)
end interface

interface root_newton
    procedure __APPEND(root_newton_jac_args,__PREC)
end interface

interface root_newton
    procedure __APPEND(root_newton_fcn_jac_args,__PREC)
end interface

interface root_newton
    procedure __APPEND(root_newton,__PREC)
end interface

interface root_newton
    procedure __APPEND(root_newton_jac,__PREC)
end interface

interface root_newton
    procedure __APPEND(root_newton_fcn_jac,__PREC)
end interface

interface root_halley
    procedure __APPEND(root_halley,__PREC)
end interface

interface root_halley
    procedure __APPEND(root_halley_args,__PREC)
end interface

interface root_newton_impl
    procedure __APPEND(root_newton_impl,__PREC)
end interface

interface root_halley_impl
    procedure __APPEND(root_halley_impl,__PREC)
end interface

interface check_inputs
    procedure __APPEND(check_inputs,__PREC)
end interface

interface set_defaults
    procedure __APPEND(set_defaults,__PREC)
end interface

interface root_newton_bisect
    procedure __APPEND(newton_bisect,__PREC)
end interface

interface root_newton_bisect
    procedure __APPEND(newton_bisect_jac,__PREC)
end interface

interface root_newton_bisect
    procedure __APPEND(newton_bisect_fcn_jac,__PREC)
end interface

interface root_newton_bisect
    procedure __APPEND(newton_bisect_args,__PREC)
end interface

interface root_newton_bisect
    procedure __APPEND(newton_bisect_jac_args,__PREC)
end interface

interface root_newton_bisect
    procedure __APPEND(newton_bisect_fcn_jac_args,__PREC)
end interface

interface newton_bisect_impl
    procedure __APPEND(newton_bisect_impl,__PREC)
end interface

interface check_bracket
    procedure __APPEND(check_bracket,__PREC)
end interface

