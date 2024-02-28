# CMake configuration for numfort examples

get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(${COMPONENT_NAME}_FILES
    concon.f90
    curfit.f90
    minimize_lbfgsb.f90
    minimize_slsqp.f90
	PARENT_SCOPE
)