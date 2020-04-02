


module numfort_optimize_solver_map

    ! Import status so status codes can be passed without user having
    ! to USE numfort_common_status.
    use numfort_common_status

    use numfort_optimize_solver_map_common
    use numfort_optimize_solver_map_real32, solver_map_real32 => solver_map
    use numfort_optimize_solver_map_real64, solver_map_real64 => solver_map

end module
