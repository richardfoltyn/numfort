

type, public :: __APPEND(workspace,__PREC)
    integer (int64), private :: roffset = 0
    integer (int64), private :: ioffset = 0
    integer (int64), private :: coffset = 0
    integer (int64), private :: loffset = 0

    real (__PREC), dimension(:), allocatable :: rwrk
    integer, dimension(:), allocatable :: iwrk
    logical, dimension(:), allocatable :: lwrk
    character (:), allocatable :: cwrk
end type

interface assert_alloc
    procedure __APPEND(ws_assert_alloc_int32,__PREC)
end interface

interface assert_alloc
    procedure __APPEND(ws_assert_alloc_int64,__PREC)
end interface

interface assert_alloc_ptr
    procedure __APPEND(ws_assert_alloc_ptr,__PREC)
end interface

interface assert_dealloc_ptr
    procedure __APPEND(ws_assert_dealloc_ptr,__PREC)
end interface

interface workspace_finalize
    procedure __APPEND(ws_finalize,__PREC)
end interface

interface workspace_get_ptr
    procedure __APPEND(ws_get_rptr_1d_int32,__PREC)
end interface

interface workspace_get_ptr
    procedure __APPEND(ws_get_rptr_2d_int32,__PREC)
end interface

interface workspace_get_ptr
    procedure __APPEND(ws_get_iptr_1d_int32,__PREC)
end interface

interface workspace_reset
    procedure __APPEND(ws_reset,__PREC)
end interface
