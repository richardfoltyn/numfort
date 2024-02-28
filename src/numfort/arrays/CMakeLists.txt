
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

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
add_dependencies(${COMPONENT_NAME} core common)

# External dependencies
target_link_libraries(${COMPONENT_NAME} PRIVATE orderpack)
