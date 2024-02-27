
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(${COMPONENT_NAME}_FILES
    basic.f90
    quadpack.f90
    PARENT_SCOPE
)

set(${COMPONENT_NAME}_SHARED_FILES
    funcs.f90
    PARENT_SCOPE
)
