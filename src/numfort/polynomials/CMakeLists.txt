
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(SOURCE_FILES
    polynomial.f90
    polyfit.F90
    polyroot.F90
    polyval.F90
)

include_directories(
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common core linalg)
