# src/numfort/numfort/arrays/CMakeLists.txt

set(COMPONENT_NAME arrays)

set(SOURCE_FILES
    arrays.f90
    arrays_sort.f90
    arrays_index.F90
    arrays_create.F90
    arrays_manipulate.f90
    arrays_norms.F90
    arrays_setops.F90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} core common orderpack)
