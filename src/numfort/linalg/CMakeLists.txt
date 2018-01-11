# src/numfort/numfort/arrays/CMakeLists.txt

set(COMPONENT_NAME linalg)

set(SOURCE_FILES
    linalg.f90
	linalg_core.F90
)

include_directories(
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} interfaces common)
