
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(${COMPONENT_NAME}_FILES
    minimize_slsqp.f90
    root_newton.f90
    PARENT_SCOPE
)
