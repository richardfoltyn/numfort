
interface polyval
    procedure &
        __APPEND(polyval_scalar,__PREC), &
        __APPEND(polyval,__PREC)
end interface


interface polyval_impl
    procedure &
        __APPEND(polyval_impl,__PREC)
end interface


interface polyval_check_input
    procedure &
        __APPEND(polyval_check_input,__PREC)
end interface
