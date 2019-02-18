

interface ppolyval
    procedure __APPEND(ppolyval,__PREC)
end interface

interface ppolyval
    procedure __APPEND(ppolyval_scalar,__PREC)
end interface

interface ppolyval_check_input
    procedure __APPEND(ppolyval_check_input,__PREC)
end interface
