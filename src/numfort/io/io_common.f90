

module numfort_io_common

    use numfort_common

    implicit none

    integer (NF_ENUM_KIND), public, parameter :: NF_IO_TRANSFORM_NONE = 1
    integer (NF_ENUM_KIND), public, parameter :: NF_IO_TRANSFORM_TRANSPOSE = 2
    integer (NF_ENUM_KIND), public, parameter :: NF_IO_TRANSFORM_FLATTEN = 4


end module