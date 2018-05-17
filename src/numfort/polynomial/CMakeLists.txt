
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(SOURCE_FILES
    polynomial.f90
    poly_create.F90
    polyfit.F90
    polyroots.F90
    polyval.F90
    poly_complete.f90
    polyval_complete.F90
    ppoly.F90
    ppolyval.F90
    ppolyint.F90
)

include_directories(
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} arrays common core linalg)
