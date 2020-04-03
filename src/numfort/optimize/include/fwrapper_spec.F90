
!*  Contains the type and interface definitions for function wrapper
!   that are precision specific.
!   To be included in module specification part.
!

type, public :: fwrapper_ss
    !*  FWRAPPER_SS implements a function wrapper for a scalar-valued function
    !   f:R->R that can optionally take additional arguments or compute
    !   its first derivative.
    private
    procedure (fss), pointer, nopass :: fcn => null()
    procedure (fss), pointer, nopass :: jac => null()
    procedure (fss_fcn_jac), pointer, nopass :: fcn_jac => null()
    procedure (fss_fcn_jac_opt), pointer, nopass :: fcn_jac_opt => null()
    procedure (fss_args), pointer, nopass :: fcn_args => null()
    procedure (fss_args), pointer, nopass :: jac_args => null()
    procedure (fss_fcn_jac_args), pointer, nopass :: fcn_jac_args => null()
    procedure (fss_fcn_jac_opt_args), pointer, nopass :: fcn_jac_opt_args => null()
    integer, public :: nfev = 0
        !*  Function evaluation counter
    class (args_data), pointer :: ptr_args => null()
        !*  Pointer to (optional) additional arguments
    real (PREC) :: eps = sqrt(epsilon(0.0))
        !*  Step size to use for numerical differentiation
    real (PREC) :: reps = sqrt(epsilon(0.0))
        !*  Relative step size for numerical differentiation.
    logical :: rel_diff = .false.
        !*  Flag determining whether absolute or relative step size will be
        !   used for computing forward differences.
    logical, public :: num_diff = .true.
        !*  Flag indicating whether derivatives are obtained by numerical
        !   differentiation.
end type

interface dispatch_fcn_jac
    procedure fss_dispatch_fcn_jac
end interface

interface dispatch_jac
    procedure fss_dispatch_jac
end interface

interface dispatch
    procedure fss_dispatch_fcn
end interface

interface wrap_procedure
    procedure fss_init
end interface

interface is_associated
    procedure fss_is_associated
end interface

! ------------------------------------------------------------------------------
! Wrapper for functions mapping vectors into scalars


type, public :: fwrapper_vs
    private
    procedure (fvs_fcn), pointer, nopass :: fcn => null()
    procedure (fvs_jac), pointer, nopass :: jac => null()
    procedure (fvs_fcn_jac), pointer, nopass :: fcn_jac => null()
    procedure (fvs_fcn_jac_opt), pointer, nopass :: fcn_jac_opt => null()
    procedure (fvs_fcn_args), pointer, nopass :: fcn_args => null()
    procedure (fvs_jac_args), pointer, nopass :: jac_args => null()
    procedure (fvs_fcn_jac_args), pointer, nopass :: fcn_jac_args => null()
    procedure (fvs_fcn_jac_opt_args), pointer, nopass :: fcn_jac_opt_args => null()
    integer, public :: nfev = 0
        !*  Function evaluation counter
    class (args_data), pointer :: ptr_args => null()
        !*  Pointer to (optional) additional arguments
    real (PREC) :: eps = sqrt(epsilon(0.0))
        !*  Step size to use for numerical differentiation
    real (PREC) :: reps = sqrt(epsilon(0.0))
        !*  Relative step size for numerical differentiation.
    logical :: rel_diff = .false.
        !*  Flag determining whether absolute or relative step size will be
        !   used for computing forward differences.
    logical, public :: num_diff = .true.
    !*  Flag indicating whether derivatives are obtained by numerical
    !   differentiation.
end type

interface wrap_procedure
    procedure fvs_init
end interface

interface dispatch
    module procedure fvs_dispatch_fcn
end interface

interface dispatch_jac
    procedure fvs_dispatch_jac
end interface

interface dispatch_fcn_jac
    procedure fvs_dispatch_fcn_jac
end interface

interface is_associated
    procedure fvs_is_associated
end interface

! ------------------------------------------------------------------------------
! Wrapper for function mapping vectors into vectors

type, public :: fwrapper_vv
    private
    procedure (fvv_fcn), pointer, nopass :: fcn => null()
    procedure (fvv_jac), pointer, nopass :: jac => null()
    procedure (fvv_fcn_jac), pointer, nopass :: fcn_jac => null()
    procedure (fvv_fcn_jac_opt), pointer, nopass :: fcn_jac_opt => null()
    procedure (fvv_fcn_args), pointer, nopass :: fcn_args => null()
    procedure (fvv_jac_args), pointer, nopass :: jac_args => null()
    procedure (fvv_fcn_jac_args), pointer, nopass :: fcn_jac_args => null()
    procedure (fvv_fcn_jac_opt_args), pointer, nopass :: fcn_jac_opt_args => null()
    integer, public :: nfev = 0
        !*  Function evaluation counter
    class (args_data), pointer :: ptr_args => null()
        !*  Pointer to (optional) additional arguments
    real (PREC) :: eps = sqrt(epsilon(0.0))
        !*  Step size to use for numerical differentiation
    real (PREC) :: reps = sqrt(epsilon(0.0))
        !*  Relative step size for numerical differentiation.
    logical :: rel_diff = .false.
        !*  Flag determining whether absolute or relative step size will be
        !   used for computing forward differences.
    logical, public :: num_diff = .true.
    !*  Flag indicating whether derivatives are obtained by numerical
    !   differentiation.
end type

interface wrap_procedure
    procedure fvv_init
end interface

interface dispatch
    module procedure fvv_dispatch_fcn
end interface

interface dispatch_jac 
    procedure fvv_dispatch_jac
end interface

interface dispatch_fcn_jac
    procedure fvv_dispatch_fcn_jac
end interface

interface is_associated
    procedure fvv_is_associated
end interface

