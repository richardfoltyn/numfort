#include <numfort.h>

module numfort_stats_dmvnorm

    use, intrinsic :: iso_fortran_env

    use numfort_common
    use numfort_common_input_checks

    use random, only: random_normal_real32, random_normal_real64

    use lapack_interfaces, only: LAPACK_GESDD => GESDD
    use blas_interfaces, only: BLAS_GEMM => GEMM

    use numfort_core

    implicit none
    private

    public :: rvs
    public :: dist_set_params
    public :: dmvnorm_real64

#include <numfort_real64.h>
#include "cont/dmvnorm_spec.F90"

contains

#include <numfort_real64.h>
#include "cont/dmvnorm_impl.F90"

end module
