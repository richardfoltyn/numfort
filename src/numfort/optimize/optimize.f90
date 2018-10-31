

module numfort_optimize

    use numfort_common_enums
    use numfort_common_kinds
    use numfort_common_status
    use numfort_common, only : workspace_real32, workspace_real64, workspace_finalize

    use numfort_optimize_result, only: optim_result_real32, optim_result_real64, &
        result_finalize, result_update

    use numfort_optimize_interfaces, only: args_data, args_default_real32, &
        args_default_real64, dynamic_cast, cond_alloc

    use numfort_optimize_diff

    use numfort_optimize_fminbound
    use numfort_optimize_lbfgsb, only: minimize_lbfgsb
    use numfort_optimize_simplex, only: minimize_simplex
    use numfort_optimize_dfls, only: minimize_dfls

    use numfort_optimize_bisect
    use numfort_optimize_brent, only: root_brentq
    use numfort_optimize_root_broyden
    use numfort_optimize_newton, only: root_newton, root_halley, root_newton_bisect

    use numfort_optimize_minpack

    use numfort_optimize_slsqp

    implicit none

end module
