
!*  Contains the type and interface definitions for function wrapper
!   that are precision specific.
!   To be included in module specification part.

type :: __APPEND(fwrapper_vec_scalar,__PREC)
    private
    procedure (__APPEND(fvec_scalar_args,__PREC)), pointer, nopass :: &
        fcn_args => null()
    procedure (__APPEND(fvec_scalar,__PREC)), pointer, nopass :: &
        fcn => null()
end type

interface dispatch
    module procedure __APPEND(fwrapper_vec_scalar_dispatch,__PREC)
end interface

interface is_present
    module procedure __APPEND(fwrapper_vec_scalar_is_present,__PREC)
end interface

interface wrap_procedure
    module procedure __APPEND(fwrapper_vec_scalar_wrap,__PREC)
end interface



type :: __APPEND(fwrapper_vec_vec,__PREC)
    private
    procedure (__APPEND(fvec_vec_args,__PREC)), pointer, nopass :: &
        fcn_args => null()
    procedure (__APPEND(fvec_vec,__PREC)), pointer, nopass :: &
        fcn => null()
end type

interface dispatch
    module procedure __APPEND(fwrapper_vec_vec_dispatch,__PREC)
end interface

interface is_present
    module procedure __APPEND(fwrapper_vec_vec_is_present,__PREC)
end interface

interface wrap_procedure
    module procedure __APPEND(fwrapper_vec_vec_wrap,__PREC)
end interface

