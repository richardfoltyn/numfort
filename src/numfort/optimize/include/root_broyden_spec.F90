

interface root_broyden_check_input
    procedure __APPEND(root_broyden_check_input,__PREC)

end interface

interface dumb_line_search
    procedure __APPEND(dumb_line_search,__PREC)
end interface

interface root_broyden
    procedure __APPEND(root_broyden,__PREC)
end interface

interface root_broyden
    procedure __APPEND(root_broyden_jac,__PREC)
end interface

interface root_broyden
    procedure __APPEND(root_broyden_fcn_jac_opt,__PREC)
end interface

interface root_broyden
    procedure __APPEND(root_broyden_args,__PREC)
end interface

interface root_broyden
    procedure __APPEND(root_broyden_jac_args,__PREC)
end interface

interface root_broyden
    procedure __APPEND(root_broyden_fcn_jac_opt_args,__PREC)
end interface

interface root_broyden_impl
    procedure __APPEND(root_broyden_impl,__PREC)
end interface
