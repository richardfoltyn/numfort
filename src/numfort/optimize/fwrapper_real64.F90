

module numfort_optimize_fwrapper_real64

    use, intrinsic :: iso_fortran_env

    use numfort_optimize_interfaces_real64
    use numfort_optimize_diff_real64

    implicit none

    private

    public :: wrap_procedure
    public :: is_associated
    public :: dispatch
    public :: dispatch_jac
    public :: dispatch_fcn_jac

    integer, parameter :: PREC = real64

#include "fwrapper_spec.F90"

    contains

#include "fwrapper_impl.F90"

end module
