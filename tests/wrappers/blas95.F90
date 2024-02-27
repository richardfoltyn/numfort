

#include <numfort.h>

program test_wrappers_blas95

    use, intrinsic :: iso_fortran_env

    use blas95
    use numfort_common
    use numfort_common_testing
    use numfort_stats, only: set_seed

    use fcore_testing, only: test_suite, test_case
    use fcore_strings

    implicit none


    call test_all ()


    contains



subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("Unit test for Fortran 95 BLAS wrappers")

    call test_gemv_real32 (tests)
    call test_gemv_real64 (tests)

    call test_gemm_real32 (tests)
    call test_gemm_real64 (tests)

    call tests%print ()

end subroutine


#include <numfort_real32.h>
#include "blas95_gemv.F90"
#include "blas95_gemm.F90"


#include <numfort_real64.h>
#include "blas95_gemv.F90"
#include "blas95_gemm.F90"

end program
