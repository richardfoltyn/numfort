
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(SOURCE_FILES
    alloc.f90
    cond_alloc.F90
    copy_alloc.F90
    copy_masked.F90
    common.f90
    enums.f90
    kinds.f90
    status.f90
    input_checks.F90
    print_helpers.f90
    shape.F90
    strings.F90
    swap.F90
    testing.f90
    testing_real32.F90
    testing_real64.F90
    transpose.f90
    workspace.f90
    workspace_real32.F90
    workspace_real64.F90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
