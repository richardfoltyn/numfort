

interface polybasis_complete
    procedure __APPEND(polybasis,__PREC)
end interface


interface polybasis_complete
    procedure __APPEND(polybasis_scalar,__PREC)
end interface

interface polybasis_complete_impl
    procedure __APPEND(polybasis_impl,__PREC)
end interface

interface polybasis_check_input
    procedure __APPEND(polybasis_check_input_array,__PREC)
end interface

interface polybasis_check_input
    procedure __APPEND(polybasis_check_input_scalar,__PREC)
end interface

interface polybasis_jac_complete
    procedure __APPEND(polybasis_jac,__PREC)
end interface

interface polybasis_jac_complete_impl
    procedure __APPEND(polybasis_jac_impl,__PREC)
end interface
