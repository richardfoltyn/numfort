

interface interp_linear
    procedure __APPEND(interp_linear_scalar,__PREC), &
        __APPEND(interp_linear_1d,__PREC)
end interface

interface interp_linear_eval
    procedure __APPEND(interp_linear_eval_scalar,__PREC), &
        __APPEND(interp_linear_eval_1d,__PREC)
end interface

interface interp_linear_eval_impl
    procedure __APPEND(interp_linear_eval_impl_scalar,__PREC), &
        __APPEND(interp_linear_eval_impl_1d,__PREC)
end interface

interface interp_linear_check_input
    procedure __APPEND(interp_linear_check_input,__PREC)
end interface

interface interp_linear_eval_check_input
    procedure __APPEND(interp_linear_eval_check_input,__PREC)
end interface

interface interp_bilinear_impl
    procedure __APPEND(interp_bilinear_impl,__PREC)
end interface

interface interp_bilinear_check_input
    procedure __APPEND(interp_bilinear_check_input,__PREC)
end interface

interface interp_bilinear
    procedure __APPEND(interp_bilinear_scalar,__PREC)
end interface

interface interp_bilinear
    procedure __APPEND(interp_bilinear_1d,__PREC)
end interface
