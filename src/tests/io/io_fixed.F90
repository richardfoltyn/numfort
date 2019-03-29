
#include <numfort.h>

program test_io_fixed
    !*  Unit tests for routines performing fixed-format text I/O.

    use, intrinsic :: iso_fortran_env

    use numfort_stats
    use numfort_common
    use numfort_common_testing
    use numfort_io

    use fcore_strings
    use fcore_testing

    implicit none

    call test_all ()


    contains


subroutine test_all ()

    type (test_suite) :: tests

    call tests%set_label ("IO unit tests for fixed-format files")

    call test_fixed_real32 (tests)
    call test_fixed_real64 (tests)

    call test_fixed_alloc_real32 (tests)
    call test_fixed_alloc_real64 (tests)

    call tests%print()

end subroutine

#include <numfort_real32.h>
#include "io_fixed_impl.F90"

#include <numfort_real64.h>
#include "io_fixed_impl.F90"



end program

