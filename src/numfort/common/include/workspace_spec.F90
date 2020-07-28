

type, public :: workspace
    integer (int64), private :: roffset = 0
    integer (int64), private :: ioffset = 0
    integer (int64), private :: coffset = 0
    integer (int64), private :: loffset = 0

    real (PREC), dimension(:), allocatable :: rwrk
    integer, dimension(:), allocatable :: iwrk
    logical, dimension(:), allocatable :: lwrk
    character (:), allocatable :: cwrk
end type

interface assert_alloc
    procedure ws_assert_alloc_int32
end interface

interface assert_alloc
    procedure ws_assert_alloc_int64
end interface

interface assert_alloc_ptr
    procedure ws_assert_alloc_ptr
end interface

interface assert_dealloc_ptr
    procedure ws_assert_dealloc_ptr
end interface

interface workspace_get_ptr
    procedure ws_get_rptr_1d_int32
end interface

interface workspace_get_ptr
    procedure ws_get_rptr_2d_int32
end interface

interface workspace_get_ptr
    procedure ws_get_iptr_1d_int32
end interface

interface workspace_reset
    procedure ws_reset
end interface
