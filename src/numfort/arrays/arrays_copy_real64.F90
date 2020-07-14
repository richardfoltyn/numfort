


module numfort_arrays_copy_real64

    use, intrinsic :: iso_fortran_env

    use numfort_common_status
    use numfort_common_input_checks

    implicit none

    private

    public :: pack_indexed
    public :: unpack_indexed

    interface pack_indexed
        procedure pack_indexed_1d, pack_indexed_2d, pack_indexed_dim_2d
    end interface

    interface unpack_indexed
        procedure unpack_indexed_1d, unpack_indexed_2d, unpack_indexed_dim_2d
    end interface

    integer, parameter :: PREC = real64

    contains

#include "include/copy_impl.f90"


end module
