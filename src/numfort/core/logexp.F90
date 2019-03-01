

#include <numfort.h>


module numfort_core_logexp

    use, intrinsic :: iso_fortran_env

    use numfort_core_libc, only: log1p

    implicit none
    private


    public :: logaddexp

#include <numfort_real32.h>
#include "logexp_spec.F90"

#include <numfort_real64.h>
#include "logexp_spec.F90"

    contains


#include <numfort_real32.h>
#include "logexp_impl.F90"


#include <numfort_real64.h>
#include "logexp_impl.F90"


end module
