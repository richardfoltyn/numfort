

module numfort_optimize

    use numfort_common_enums
    use numfort_common_kinds
    use numfort_common_status
    use numfort_common, only : workspace_real32, workspace_real64

    use numfort_optimize_interfaces_common, only: args_data
    use numfort_optimize_interfaces_real32, only: dynamic_cast, cond_alloc, &
        args_default_real32 => args_default
    use numfort_optimize_interfaces_real64, only: dynamic_cast, cond_alloc, &
        args_default_real64 => args_default

    use numfort_optimize_result_real32, optim_result_real32 => optim_result
    use numfort_optimize_result_real64, optim_result_real64 => optim_result

    use numfort_optimize_fminbound_real32
    use numfort_optimize_fminbound_real64

    use numfort_optimize_lbfgsb_real64

    use numfort_optimize_simplex_real32, only: minimize_simplex
    use numfort_optimize_simplex_real64, only: minimize_simplex

    use numfort_optimize_dfls_real32
    use numfort_optimize_dfls_real64

    use numfort_optimize_bisect_real32
    use numfort_optimize_bisect_real64

    use numfort_optimize_brent_real32
    use numfort_optimize_brent_real64

    use numfort_optimize_root_broyden_real32
    use numfort_optimize_root_broyden_real64

    use numfort_optimize_newton_real32
    use numfort_optimize_newton_real64

    use numfort_optimize_minpack_real64

    use numfort_optimize_slsqp_real64

    use numfort_optimize_slsqp_ng_real64

    use numfort_optimize_solver_map

    implicit none

end module
