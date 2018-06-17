module numfort_interpolate_common

    use iso_fortran_env
    use numfort_common, only: NF_ENUM_KIND

    implicit none
    private

    integer (NF_ENUM_KIND), public, parameter :: NF_INTERP_EVAL_EXTRAPOLATE = 0
    integer (NF_ENUM_KIND), public, parameter :: NF_INTERP_EVAL_ZERO = ishft(1, 0)
    integer (NF_ENUM_KIND), public, parameter :: NF_INTERP_EVAL_ERROR = ishft(1, 1)
    integer (NF_ENUM_KIND), public, parameter :: NF_INTERP_EVAL_BOUNDARY = ishft(1, 2)
    integer (NF_ENUM_KIND), public, parameter :: NF_INTERP_EVAL_CONST = ishft(1, 3)



end module
