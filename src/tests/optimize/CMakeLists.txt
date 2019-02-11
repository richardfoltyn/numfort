
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(${COMPONENT_NAME}_FILES
    minimize_dfls.f90
    minimize_lbfgsb.f90
    minimize_simplex.f90
    minimize_slsqp.f90
    minpack.f90
    root_broyden.f90
    root_newton.f90
    solver_map.f90
    PARENT_SCOPE
)
