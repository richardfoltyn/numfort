

interface ppolyval
    procedure __APPEND(ppoly2d_val,__PREC)
end interface

interface ppolyval
    procedure __APPEND(ppoly2d_val_scalar,__PREC)
end interface

interface ppoly2d_val_check_input
    procedure __APPEND(ppoly2d_val_check_input,__PREC)
end interface

interface ppoly2d_val_bilinear
    procedure __APPEND(ppoly2d_val_bilinear,__PREC)
end interface
