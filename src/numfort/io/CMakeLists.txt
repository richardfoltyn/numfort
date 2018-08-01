# src/numfort/numfort/interp/CMakeLists.txt

set(COMPONENT_NAME io)

set(SOURCE_FILES
    io.f90
    io_common.f90
    io_fixed.F90
)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common)
