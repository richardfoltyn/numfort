
!*  Contains the type and interface definitions for function wrapper
!   that are precision specific.
!   To be included in module specification part.
!
!   Naming convention:
!   fwrapper_vs_jac_args defines a wrapper around ta function that
!       1. Maps a vector into a scalar
!       2. Takes an (optional) Jacobian as dummy argument
!       3. Takes (optional) ARGS array as dummy argument

type, public :: __APPEND(fwrapper_ss,__PREC)
    !*  FWRAPPER_SS implements a function wrapper for a scalar-valued function
    !   f:R->R that can optionally take additional arguments or compute
    !   its first derivative.
    private
    procedure (__APPEND(fss_args,__PREC)), pointer, nopass :: &
        fcn_args => null()
    procedure (__APPEND(fss,__PREC)), pointer, nopass :: &
        fcn => null()
    procedure (__APPEND(fss_jac,__PREC)), pointer, nopass :: &
        fcn_jac => null()
    procedure (__APPEND(fss_jac_args,__PREC)), pointer, nopass :: &
        fcn_jac_args => null()
    procedure (__APPEND(fss,__PREC)), pointer, nopass :: &
        jac => null()
    procedure (__APPEND(fss_args,__PREC)), pointer, nopass :: &
        jac_args => null()
    integer, public :: nfev = 0
        !*  Function evaluation counter
    real (__PREC), dimension(:), pointer :: ptr_args => null()
        !*  Pointer to (optional) additional argument array
    real (__PREC) :: eps
        !*  Step size to use for numerical differentiation
end type

interface dispatch_fcn_jac
    procedure __APPEND(fss_dispatch_fcn_jac,__PREC)
end interface

interface dispatch_jac
    procedure __APPEND(fss_dispatch_jac,__PREC)
end interface

interface dispatch
    procedure __APPEND(fss_dispatch_fcn,__PREC)
end interface

interface wrap_procedure
    procedure __APPEND(fss_init,__PREC)
end interface

interface is_associated
    procedure __APPEND(fss_is_associated,__PREC)
end interface

! ------------------------------------------------------------------------------
! Wrapper for functions mapping vectors into scalars


type, public :: __APPEND(fwrapper_vs,__PREC)
    private
    procedure (__APPEND(fvs_fcn,__PREC)), pointer, nopass :: fcn => null()
    procedure (__APPEND(fvs_jac,__PREC)), pointer, nopass :: jac => null()
    procedure (__APPEND(fvs_fcn_jac,__PREC)), pointer, nopass :: fcn_jac => null()
    procedure (__APPEND(fvs_fcn_args,__PREC)), pointer, nopass :: fcn_args => null()
    procedure (__APPEND(fvs_jac_args,__PREC)), pointer, nopass :: jac_args => null()
    procedure (__APPEND(fvs_fcn_jac_args,__PREC)), pointer, nopass :: fcn_jac_args => null()
    integer, public :: nfev = 0
        !*  Function evaluation counter
    real (__PREC), dimension(:), pointer :: ptr_args => null()
        !*  Pointer to (optional) additional argument array
    real (__PREC) :: eps
        !*  Step size to use for numerical differentiation
end type

interface wrap_procedure
    procedure __APPEND(fvs_init,__PREC)
end interface

interface dispatch
    module procedure __APPEND(fvs_dispatch_fcn,__PREC)
end interface

interface dispatch_jac
    procedure __APPEND(fvs_dispatch_jac,__PREC)
end interface

interface dispatch_fcn_jac
    procedure __APPEND(fvs_dispatch_fcn_jac,__PREC)
end interface

interface is_associated
    procedure __APPEND(fvs_is_associated,__PREC)
end interface

! ------------------------------------------------------------------------------
! Wrapper for function mapping vectors into vectors

type, public :: __APPEND(fwrapper_vv,__PREC)
    private
    procedure (__APPEND(fvv_fcn,__PREC)), pointer, nopass :: fcn => null()
    procedure (__APPEND(fvv_jac,__PREC)), pointer, nopass :: jac => null()
    procedure (__APPEND(fvv_fcn_jac,__PREC)), pointer, nopass :: fcn_jac => null()
    procedure (__APPEND(fvv_fcn_args,__PREC)), pointer, nopass :: fcn_args => null()
    procedure (__APPEND(fvv_jac_args,__PREC)), pointer, nopass :: jac_args => null()
    procedure (__APPEND(fvv_fcn_jac_args,__PREC)), pointer, nopass :: fcn_jac_args => null()
    integer, public :: nfev = 0
        !*  Function evaluation counter
    real (__PREC), dimension(:), pointer :: ptr_args => null()
        !*  Pointer to (optional) additional argument array
    real (__PREC) :: eps
        !*  Step size to use for numerical differentiation
end type

interface wrap_procedure
    procedure __APPEND(fvv_init,__PREC)
end interface

interface dispatch
    module procedure __APPEND(fvv_dispatch_fcn,__PREC)
end interface

interface dispatch_jac 
    procedure __APPEND(fvv_dispatch_jac,__PREC)
end interface

interface dispatch_fcn_jac
    procedure __APPEND(fvv_dispatch_fcn_jac,__PREC)
end interface

interface is_associated
    procedure __APPEND(fvv_is_associated,__PREC)
end interface

