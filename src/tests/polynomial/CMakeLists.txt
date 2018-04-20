
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(${COMPONENT_NAME}_FILES
    polyfit.f90
    polyroots.f90
    polyval.f90
    poly_complete.f90
    PARENT_SCOPE
)
