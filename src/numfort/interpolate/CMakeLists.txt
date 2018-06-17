# src/numfort/numfort/interp/CMakeLists.txt

set(COMPONENT_NAME interpolate)

set(SOURCE_FILES
	interpolate.f90
	common.f90
    interp_linear.F90
	interp_pchip.F90
	fitpack.f90
	search.F90
)

include_directories("${CMAKE_CURRENT_LIST_DIR}/include")
add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common core fitpack)
