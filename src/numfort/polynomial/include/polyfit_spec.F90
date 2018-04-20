interface polyfit
    procedure &
        __APPEND(polyfit,__PREC), &
        __APPEND(polyfit_1d,__PREC)
end interface

interface polyfit_deriv
    procedure &
        __APPEND(polyfit_deriv_scalar,__PREC)
end interface

interface polyfit_exact
    procedure &
        __APPEND(polyfit_exact,__PREC)
end interface

interface polyfit_check_input
    procedure &
        __APPEND(polyfit_check_input,__PREC)
end interface

interface polyfit_deriv_check_input
    procedure &
        __APPEND(polyfit_deriv_check_input,__PREC)
end interface
