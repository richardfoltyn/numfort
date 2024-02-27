
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(${COMPONENT_NAME}_FILES
    copy.f90
    delete.f90
    indexing.f90
    insert.f90
    kron.f90
    setops.f90
    sort.F90
    vander.f90
    PARENT_SCOPE
)
