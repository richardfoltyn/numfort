
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(${COMPONENT_NAME}_FILES
    cond_alloc.f90
    copy_alloc.f90
    has_shape.f90
    PARENT_SCOPE
)
