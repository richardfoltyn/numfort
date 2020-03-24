

module numfort_interpolate_search_real64

    use, intrinsic :: iso_fortran_env

    use numfort_common_kinds, only: NF_ENUM_KIND
    use numfort_common_status
    use numfort_interpolate_common

    implicit none

    private

    public :: bsearch
    public :: bsearch_cached

    public :: interp_find
    public :: interp_find_impl

    public :: interp_find_decr
    public :: interp_find_decr_impl
    public :: interp_find_incr_impl

    integer, parameter :: PREC = real64

#include "search_spec.F90"


    contains


#include "search_impl.F90"


end module
