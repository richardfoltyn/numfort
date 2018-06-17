
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(${COMPONENT_NAME}_FILES
    concon.f90
    curfit.f90
    linear.f90
    pchip.f90
    search.f90
    PARENT_SCOPE
)
