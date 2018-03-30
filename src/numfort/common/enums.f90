module numfort_common_enums

    use numfort_common_kinds, only: NF_ENUM_KIND
    implicit none

    integer (NF_ENUM_KIND), parameter :: NF_PRINT_NONE = 0
    integer (NF_ENUM_KIND), parameter :: NF_PRINT_MINIMAL = 10
    integer (NF_ENUM_KIND), parameter :: NF_PRINT_VERBOSE = 20
    integer (NF_ENUM_KIND), parameter :: NF_PRINT_ALL = not(0)

    integer (NF_ENUM_KIND), parameter :: NF_LINESEARCH_BACKTRACK = 0
    integer (NF_ENUM_KIND), parameter :: NF_LINESEARCH_EXACT = 1

end module
