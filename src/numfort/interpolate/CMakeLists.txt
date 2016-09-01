# src/numfort/numfort/interp/CMakeLists.txt

set(COMPONENT_NAME interpolate)

set(SOURCE_FILES
	interpolate.f90
	common.f90
	interp_result.f90
	interp_linear.f90
	fitpack_interfaces.f90
	fitpack08.f90
	fitpack/concon.f90
)

include_directories(${${PROJECT_NAME}_MODDIR})
add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common)
