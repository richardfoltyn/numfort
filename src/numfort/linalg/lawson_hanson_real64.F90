
module numfort_linalg_lh95_real64


    use, intrinsic :: iso_fortran_env, only: real64

    use numfort_common_enums
    use numfort_common_status

    use numfort_linalg_lh95_common

    use blas_interfaces, only: NRM2
    use lapack_interfaces, only: LARTGP

    implicit none

    private

    public :: h12
    public :: hfti, hfti_query
    public :: nnls, nnls_query
    public :: ldp, ldp_query

    integer, parameter :: PREC = real64

    interface h12
        procedure h12_1d, h12_2d
    end interface

    interface h12_apply
        procedure h12_apply_1d, h12_apply_2d
    end interface

    interface h12_compute
        procedure h12_compute
    end interface

    interface hfti
        procedure hfti_1d, hfti_2d
    end interface

    interface nnls
        procedure nnls
    end interface

    interface ldp
        procedure ldp
    end interface


    contains

#include "lawson_hanson_impl.f90"


end
