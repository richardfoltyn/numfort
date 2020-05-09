


module numfort_optimize_slsqp_ng_common

    use numfort_common_enums

    implicit none

    ! Internal status codes used only by SLSQP-NG
    ! These are the same codes as used in the original implementation
    integer (NF_ENUM_KIND), parameter :: STATUS_OK = 1
    integer (NF_ENUM_KIND), parameter :: STATUS_INVALID_DIMS = 2
    integer (NF_ENUM_KIND), parameter :: STATUS_LSQ_MAX_ITER = 3
    integer (NF_ENUM_KIND), parameter :: STATUS_INCOMPAT_CONSTR = 4
    integer (NF_ENUM_KIND), parameter :: STATUS_LSI_SINGULAR = 5
    integer (NF_ENUM_KIND), parameter :: STATUS_LSEI_SINGULAR_CONSTR = 6
    integer (NF_ENUM_KIND), parameter :: STATUS_HFTI_RANK_DEFECT = 7


    character (*), parameter :: MSG_OBJECTIVE_NONFINITE &
        = 'SLSQP: Objective function returned non-finite values'

    character (*), parameter :: MSG_CONSTRAINTS_NONFINITE &
        = 'SLSQP: Constraints function returned non-finite values'




end module
