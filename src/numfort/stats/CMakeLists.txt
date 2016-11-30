
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(SOURCE_FILES
    cont_dist.f90
    disc_dist.F90
    dist_uniform.f90
    dist_randint.f90
)

include_directories(${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS})
add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common)
# external libraries
add_dependencies(${COMPONENT_NAME} cdflib90 random)