


module numfort_optimize_root_broyden_common

    use numfort_common_kinds
    use numfort_common_enums

    integer (NF_ENUM_KIND), parameter :: PRINT_LSEARCH = 1
    integer (NF_ENUM_KIND), parameter :: PRINT_STEP = 2
    integer (NF_ENUM_KIND), parameter :: PRINT_JAC = 4

    integer (NF_ENUM_KIND), parameter :: PRINT_NONE = NF_PRINT_NONE
    integer (NF_ENUM_KIND), parameter :: PRINT_ALL = NF_PRINT_ALL

    ! Public fully qualified parameter names
    integer (NF_ENUM_KIND), public, parameter :: &
        NF_ROOT_BROYDEN_PRINT_LSEARCH = PRINT_LSEARCH
    integer (NF_ENUM_KIND), public, parameter :: &
        NF_ROOT_BROYDEN_PRINT_STEP = PRINT_STEP
    integer (NF_ENUM_KIND), public, parameter :: &
        NF_ROOT_BROYDEN_PRINT_JAC = PRINT_JAC
    integer (NF_ENUM_KIND), public, parameter :: &
        NF_ROOT_BROYDEN_PRINT_ALL = PRINT_ALL
    integer (NF_ENUM_KIND), public, parameter :: &
        NF_ROOT_BROYEN_PRINT_NONE = PRINT_NONE

    integer, parameter :: LINESEARCH_MAX_STEPS = 4
        !*  Max. number of evaluations during line search (ie. at most
        !   LINESEARCH_MAX_STEPS-1 backtracking steps will be performed)

end module
