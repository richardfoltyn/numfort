

type, public, extends(args_data) :: __APPEND(args_default,__PREC)
    real (__PREC), dimension(:), allocatable :: rdata
    integer, dimension(:), allocatable :: idata
end type


interface dynamic_cast
   procedure __APPEND(cast_to_args_default,__PREC)
end interface

interface cond_alloc
    procedure __APPEND(cond_alloc_args_default,__PREC)
end interface
