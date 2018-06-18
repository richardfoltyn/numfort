
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(${COMPONENT_NAME}_FILES
    polyfit.f90
    polyroots.f90
    polyval.f90
    poly_complete.f90
    poly_create.f90
    ppoly.f90
    ppoly_tensor.f90
    PARENT_SCOPE
)
