

interface ppolyval
    procedure __APPEND(ppolyval,__PREC)
end interface

interface ppolyval
    procedure __APPEND(ppolyval_scalar,__PREC)
end interface

interface ppolyval_check_input
    procedure __APPEND(ppolyval_check_input,__PREC)
end interface

interface ppolyval
    procedure __APPEND(bernstein_ppolyval,__PREC)
end interface

interface ppolyval
    procedure __APPEND(bernstein_ppolyval_scalar,__PREC)
end interface

interface ppolyval_impl
    procedure __APPEND(bernstein_ppolyval_impl,__PREC)
end interface

interface ppolyval_check_input
    procedure __APPEND(bernstein_ppolyval_check_input,__PREC)
end interface

interface bernstein_polyval
    procedure __APPEND(bernstein_polyval,__PREC)
end interface
