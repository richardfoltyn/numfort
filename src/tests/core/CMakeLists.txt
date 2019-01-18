
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(${COMPONENT_NAME}_FILES
    combinatorics.f90
    cephes_stats.f90
    PARENT_SCOPE
)
