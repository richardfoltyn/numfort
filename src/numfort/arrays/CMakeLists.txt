# src/numfort/numfort/arrays/CMakeLists.txt

set(COMPONENT_NAME arrays)

set(SOURCE_FILES
    arrays.f90
    arrays_sort.F90
    arrays_index.F90
    arrays_create.F90
    arrays_grid.F90
    arrays_grid_real32.F90
    arrays_grid_real64.F90
    arrays_copy.f90
    arrays_copy_logical.F90
    arrays_copy_real32.F90
    arrays_copy_real64.F90
    arrays_manipulate.f90
    arrays_norms.F90
    arrays_setops.F90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} core common orderpack)
