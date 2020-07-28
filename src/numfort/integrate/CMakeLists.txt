# src/numfort/numfort/integrate/CMakeLists.txt

set(COMPONENT_NAME integrate)

set(SOURCE_FILES
    integrate.f90
    integrate_basic.F90
    quadpack.f90
    quadpack_common.F90
    quadpack_real64.F90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common quadpack)
