
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(SOURCE_FILES
    hermite_e.F90
    polynomial.f90
    poly_create.F90
    polyfit.F90
    polyroots.F90
    polyval.F90
    poly_complete.f90
    polyval_complete.F90
    ppoly.F90
    ppoly_tensor.F90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} arrays common core linalg)
