
#include <numfort.h>

module numfort_interpolate_search

    use, intrinsic :: iso_fortran_env

    implicit none

    private

    public :: bsearch
    public :: interp_find

    public :: search_cache
    public :: bsearch_cached
    public :: interp_find_cached

    type :: search_cache
        private
        integer :: i = 1
    end type

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