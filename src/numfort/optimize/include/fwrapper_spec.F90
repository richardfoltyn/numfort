
!*  Contains the type and interface definitions for function wrapper
!   that are precision specific.
!   To be included in module specification part.
!
!   Naming convention:
!   fwrapper_vs_jac_args defines a wrapper around ta function that
!       1. Maps a vector into a scalar
!       2. Takes an (optional) Jacobian as dummy argument
!       3. Takes (optional) ARGS array as dummy argument

type :: __APPEND(fwrapper_vs_jac,__PREC)
    private
    procedure (__APPEND(fvs_jac_args,__PREC)), pointer, nopass :: &
        fcn_args => null()
    procedure (__APPEND(fvs_jac,__PREC)), pointer, nopass :: &
        fcn => null()
end type

interface dispatch
    module procedure __APPEND(fwrapper_vs_jac_dispatch,__PREC)
end interface

interface is_present
    module procedure __APPEND(fwrapper_vs_jac_is_present,__PREC)
end interface

interface wrap_procedure
    module procedure __APPEND(fwrapper_vs_jac_wrap,__PREC)
end interface



type :: __APPEND(fwrapper_vv_jac,__PREC)
    private
    procedure (__APPEND(fvv_jac_args,__PREC)), pointer, nopass :: &
        fcn_args => null()
    procedure (__APPEND(fvv_jac,__PREC)), pointer, nopass :: &
        fcn => null()
end type

interface dispatch
    module procedure __APPEND(fwrapper_vv_jac_dispatch,__PREC)
end interface

interface is_present
    module procedure __APPEND(fwrapper_vv_jac_is_present,__PREC)
end interface

interface wrap_procedure
    module procedure __APPEND(fwrapper_vv_jac_wrap,__PREC)
end interface

