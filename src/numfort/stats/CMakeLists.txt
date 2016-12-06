
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(SOURCE_FILES
    combinatorics.F90
    descriptive.F90
    dcont.F90
    disc_dist.F90
    dnorm.F90
    duniform.F90
    drandint.F90
    stats.f90
)

include_directories(
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common core)
# external libraries
add_dependencies(${COMPONENT_NAME} cdflib90 random)
