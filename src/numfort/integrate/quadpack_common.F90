
module numfort_integrate_quadpack_common

    use numfort_common_enums

    ! Custom error codes for QUADPACK wrappers
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_MAX_INTERVALS = ishft(1, 1)
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_INTEGRAND_ERROR = ishft(1, 2)

end module
