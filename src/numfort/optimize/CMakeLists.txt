
set(COMPONENT_NAME optimize)

set(SOURCE_FILES
    optimize.f90
    common.f90
    optim_result.f90
    minimize_lbfgsb.f90
    minimize_simplex.f90
)

include_directories(${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS})
add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common)
