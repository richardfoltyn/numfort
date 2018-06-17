
#include <numfort.h>

module numfort_interpolate_search

    use, intrinsic :: iso_fortran_env

    implicit none

    private

    public :: bsearch
    public :: interp_find

#include <numfort_real32.h>
#include "search_spec.F90"

#include <numfort_real64.h>
#include "search_spec.F90"

   contains

#include <numfort_real32.h>
#include "search_impl.F90"

#include <numfort_real64.h>
#include "search_impl.F90"


end module