
type, public :: __APPEND(optim_result,__PREC)
    integer :: nfev = UNINITIALIZED_COUNTER
        !*  Number of function evaluations performed
    integer :: nit = UNINITIALIZED_COUNTER
        !*  Number of iterations performed
    type (status_t) :: status
        !*  Detailed exit status
    logical :: success = .false.
        !*  Exit status flag
    character (100) :: msg
    real (__PREC), dimension(:), allocatable :: x
    real (__PREC), dimension(:), allocatable :: fx
end type

interface result_update
    procedure __APPEND(update_ss,__PREC)
end interface

interface result_update
    procedure __APPEND(update_vs,__PREC)
end interface

interface result_update
    procedure __APPEND(update,__PREC)
end interface

interface result_update
    procedure __APPEND(update_int_status,__PREC)
end interface

interface result_reset
    procedure __APPEND(reset,__PREC)
end interface

interface assert_alloc_ptr
    procedure __APPEND(assert_alloc_ptr,__PREC)
end interface

interface assert_dealloc_ptr
    procedure __APPEND(assert_dealloc_ptr,__PREC)
end interface

interface result_finalize
    procedure __APPEND(result_finalize,__PREC)
end interface
