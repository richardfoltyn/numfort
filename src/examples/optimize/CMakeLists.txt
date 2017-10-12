get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(${COMPONENT_NAME}_FILES
    root_brent.f90
    minimize_lbfgsb.f90
    minimize_simplex.f90
    minimize_slsqp.f90
    minpack_lm.f90
    minpack_hybr.f90
    root_newton.f90
    PARENT_SCOPE
)
