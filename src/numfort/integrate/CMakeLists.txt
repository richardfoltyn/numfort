# src/numfort/numfort/integrate/CMakeLists.txt

set(COMPONENT_NAME integrate)

set(SOURCE_FILES
    integrate.f90
    integrate_basic.F90
    quadpack.F90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common quadpack)
