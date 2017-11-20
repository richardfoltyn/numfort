# src/numfort/numfort/arrays/CMakeLists.txt

set(COMPONENT_NAME arrays)

set(SOURCE_FILES
    arrays.f90
    arrays_sort.f90
    arrays_index.f90
    arrays_copy.f90
    arrays_create.f90
    arrays_norms.F90
    arrays_setops.F90
    arrays_shape.f90
)

include_directories(
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} core common orderpack)
