
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(SOURCE_FILES
    optimize.f90
    diff_helpers_real32.F90
    diff_helpers_real64.F90
    interfaces_common.f90
    interfaces_real32.F90
    interfaces_real64.F90
    fwrapper_real32.F90
    fwrapper_real64.F90
    optim_result_real32.F90
    optim_result_real64.F90
    minimize_bounded_real32.F90
    minimize_bounded_real64.F90
    minimize_dfls_real32.F90
    minimize_dfls_real64.F90
    minimize_lbfgsb_real64.F90
    minimize_simplex_real32.F90
    minimize_simplex_real64.F90
    minimize_slsqp_real64.f90
    minimize_slsqp_ng_common.f90
    minimize_slsqp_ng_real64.F90
    minpack_real64.F90
    root_bisect_real32.F90
    root_bisect_real64.F90
    root_brent_real32.F90
    root_brent_real64.F90
    root_broyden_common.f90
    root_broyden_real32.F90
    root_broyden_real64.F90
    root_newton_real32.F90
    root_newton_real64.F90
    solver_map.f90
    solver_map_common.f90
    solver_map_real32.F90
    solver_map_real64.F90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common arrays)

# external libraries
set(EXTERNAL interfaces simplex_csiro lbfgsb_bmnz minpack newuoa2 slsqp)

foreach (_ext IN LISTS EXTERNAL)
    add_dependencies(${COMPONENT_NAME} ${_ext})
    target_include_directories(${COMPONENT_NAME} PRIVATE
        "${EXTERNAL_MODULE_DIRECTORY}/${_ext}"
    )
endforeach()
