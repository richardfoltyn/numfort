

set(COMPONENT_NAME "wrappers_blas95")


set(SOURCE_FILES
    blas95.F90
)

include_directories(
    ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}
)

message(STATUS "Building BLAS95 wrappers")

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} interfaces)

# append to list of object library targets
set_property(GLOBAL APPEND
    PROPERTY ${CMAKE_PROJECT_NAME}_COMPONENT_TARGETS ${COMPONENT_NAME}
)
