

module numfort_interpolate_result

    use numfort_common, only: ENUM_KIND
    implicit none

    integer (ENUM_KIND), parameter :: INTERP_STATUS_SUCCESS = 0
    integer (ENUM_KIND), parameter :: INTERP_STATUS_INVALID_INPUT = 2 ** 29
    integer (ENUM_KIND), parameter :: INTERP_STATUS_UNKNOWN = 2 ** 30

end module
