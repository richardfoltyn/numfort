

module numfort_linalg_lh95_common

    use numfort_common_enums

    implicit none

    ! Error codes originally defined in the implementation by
    ! Lawson / Hanson (1995)

    integer (NF_ENUM_KIND), parameter :: STATUS_OK = 1
    integer (NF_ENUM_KIND), parameter :: STATUS_INVALID_DIMS = 2
    integer (NF_ENUM_KIND), parameter :: STATUS_MAX_ITER = 3
    integer (NF_ENUM_KIND), parameter :: STATUS_INCOMPAT_CONSTR = 4

end module
