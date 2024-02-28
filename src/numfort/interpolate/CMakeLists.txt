
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(SOURCE_FILES
	interpolate.f90
    common.F90
    interp_linear.F90
	interp_pchip.F90
	fitpack.f90
	search.F90
	search_real32.F90
	search_real64.F90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common core)

# External dependencies
target_link_libraries(${COMPONENT_NAME} PRIVATE fitpack)
