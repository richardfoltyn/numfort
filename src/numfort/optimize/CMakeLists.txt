
set(COMPONENT_NAME optimize)

set(SOURCE_FILES
    optimize.f90
    interfaces.f90
    fwrapper.f90
    optim_result.F90
    minimize_lbfgsb.f90
    minimize_simplex.f90
    minimize_slsqp.f90
    brent.F90
    minpack.f90
    root_newton.F90
)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common)
# external libraries
add_dependencies(${COMPONENT_NAME}
    interfaces
    simplex_csiro
    lbfgsb_bmnz
    minpack
    slsqp
)
