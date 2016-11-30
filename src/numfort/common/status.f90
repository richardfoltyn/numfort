module numfort_common_status

    use iso_fortran_env, only: int32
    implicit none

    integer, parameter :: STATUS_KIND = int32

    integer (STATUS_KIND), parameter :: STATUS_OK = 0
    integer (STATUS_KIND), parameter :: STATUS_INVALID_INPUT = 2 ** 10
    integer (STATUS_KIND), parameter :: STATUS_UNKNOWN = 2 ** 11
    integer (STATUS_KIND), parameter :: STATUS_UNSUPPORTED_OPERATION = 2 * 12
    integer (STATUS_KIND), parameter :: STATUS_INVALID_STATE = 2 ** 13


end module
