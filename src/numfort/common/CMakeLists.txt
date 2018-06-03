# src/numfort/numfort/interp/CMakeLists.txt

set(COMPONENT_NAME common)

set(SOURCE_FILES
    alloc.f90
    cond_alloc.F90
    copy_alloc.F90
    common.f90
    enums.f90
    kinds.f90
    status.f90
    workspace.F90
    input_checks.f90
    shape.F90
    strings.F90
    swap.F90
    testing.F90
)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
