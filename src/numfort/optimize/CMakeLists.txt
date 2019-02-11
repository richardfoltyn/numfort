
set(COMPONENT_NAME optimize)

set(SOURCE_FILES
    optimize.f90
    diff_helpers.F90
    interfaces.F90
    fwrapper.F90
    optim_result.F90
    minimize_bounded.F90
    minimize_dfls.F90
    minimize_lbfgsb.f90
    minimize_simplex.F90
    minimize_slsqp.f90
    minpack.F90
    root_bisect.F90
    root_brent.F90
    root_broyden.F90
    root_newton.F90
    solver_map.F90
)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common arrays)
# external libraries
add_dependencies(${COMPONENT_NAME}
    interfaces
    simplex_csiro
    lbfgsb_bmnz
    minpack
    newuoa2
    slsqp
)
