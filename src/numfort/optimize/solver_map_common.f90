

module numfort_optimize_solver_map_common

    use numfort_common_kinds

    integer (NF_ENUM_KIND), parameter :: TRANSFORM_LINEAR = 1
    integer (NF_ENUM_KIND), parameter :: TRANSFORM_EXP = 2
    integer (NF_ENUM_KIND), parameter :: TRANSFORM_NEG_EXP = 3
    integer (NF_ENUM_KIND), parameter :: TRANSFORM_LOGISTIC = 4

end module
