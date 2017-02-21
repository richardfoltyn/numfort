module numfort_common_enums

    use numfort_common_kinds, only: NF_ENUM_KIND
    implicit none

    integer (NF_ENUM_KIND), parameter :: NF_STATUS_UNDEFINED = 0

    integer (NF_ENUM_KIND), parameter :: NF_STATUS_OK = ishft(1, 0)
    integer (NF_ENUM_KIND), parameter :: NF_STATUS_INVALID_ARG = ishft(1, 10)
    integer (NF_ENUM_KIND), parameter :: NF_STATUS_UNKNOWN = ishft(1, 11)
    integer (NF_ENUM_KIND), parameter :: NF_STATUS_UNSUPPORTED_OP = ishft(1, 12)
    integer (NF_ENUM_KIND), parameter :: NF_STATUS_INVALID_STATE = ishft(1, 13)
    integer (NF_ENUM_KIND), parameter :: NF_STATUS_BOUNDS_ERROR = ishft(1, 14)
    integer (NF_ENUM_KIND), parameter :: NF_STATUS_MAX_ITER = ishft(1, 15)
    integer (NF_ENUM_KIND), parameter :: NF_STATUS_MAX_EVAL = ishft(1, 16)
    integer (NF_ENUM_KIND), parameter :: NF_STATUS_NOT_CONVERGED = ishft(1, 17)
    integer (NF_ENUM_KIND), parameter :: NF_STATUS_STORAGE_ERROR = ishft(1, 18)
    integer (NF_ENUM_KIND), parameter :: NF_STATUS_OTHER = ishft(1, 19)
    integer (NF_ENUM_KIND), parameter :: NF_STATUS_APPROX = ishft(1, 20)

    integer (NF_ENUM_KIND), parameter :: NF_PRINT_NONE = 0
    integer (NF_ENUM_KIND), parameter :: NF_PRINT_MINIMAL = 10
    integer (NF_ENUM_KIND), parameter :: NF_PRINT_VERBOSE = 20
    integer (NF_ENUM_KIND), parameter :: NF_PRINT_ALL = 30


end module
