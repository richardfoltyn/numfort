

interface root_bisect
    procedure __APPEND(root_bisect,__PREC)
end interface

interface root_bisect
    procedure __APPEND(root_bisect_args,__PREC)
end interface

interface root_bisect_impl
    procedure __APPEND(root_bisect_impl,__PREC)
end interface

interface bisect_check_inputs
    procedure __APPEND(bisect_check_inputs,__PREC)
end interface
