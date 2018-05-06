

interface interp_pchip_fit
    procedure __APPEND(interp_pchip_fit,__PREC)
end interface

interface interp_pchip_eval
    procedure __APPEND(interp_pchip_eval,__PREC)
end interface

interface pchip_fit_input_check
    procedure __APPEND(pchip_fit_input_check,__PREC)
end interface

interface pchip_slope_end
    procedure __APPEND(pchip_slope_end,__PREC)
end interface

interface pchip_eval_input_check
    procedure __APPEND(pchip_eval_input_check,__PREC)
end interface
