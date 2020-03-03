

interface power_ppoly_check_input
    procedure __APPEND(power_ppoly_check_input,__PREC)
end interface

interface ppolyval
    procedure __APPEND(power_ppolyval_scalar,__PREC)
end interface

interface ppolyval
    procedure __APPEND(power_ppolyval_1d,__PREC)
end interface

interface ppolyval_eval
    procedure __APPEND(power_ppolyval_eval_scalar,__PREC)
end interface

interface ppolyval_eval
    procedure __APPEND(power_ppolyval_eval_1d,__PREC)
end interface

interface ppolyval_eval_impl
    procedure __APPEND(power_ppolyval_eval_impl_scalar,__PREC)
end interface

interface ppolyval_eval_impl
    procedure __APPEND(power_ppolyval_eval_impl_1d,__PREC)
end interface

interface ppolyval_eval_impl_deg2
    procedure __APPEND(power_ppolyval_eval_impl_deg2_1d,__PREC)
end interface

interface ppolyval_eval_impl_deg3
    procedure __APPEND(power_ppolyval_eval_impl_deg3_1d,__PREC)
end interface

interface ppolyval_eval_impl_degk
    procedure __APPEND(power_ppolyval_eval_impl_degk_1d,__PREC)
end interface

interface power_polyval
    procedure __APPEND(power_polyval,__PREC)
end interface

interface ppolyder
    procedure __APPEND(power_ppolyder,__PREC)
end interface

