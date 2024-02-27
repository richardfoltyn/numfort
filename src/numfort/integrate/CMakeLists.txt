
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(SOURCE_FILES
    integrate.f90
    integrate_basic.F90
    quadpack.f90
    quadpack_common.F90
    quadpack_real64.F90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common quadpack)

target_include_directories(${COMPONENT_NAME} PRIVATE
    "${EXTERNAL_MODULE_DIRECTORY}/quadpack"
)
