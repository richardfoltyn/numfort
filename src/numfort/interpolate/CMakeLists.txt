# src/numfort/numfort/interpolate/CMakeLists.txt

set(COMPONENT_NAME interpolate)

set(SOURCE_FILES
	interpolate.f90
    common.F90
    interp_linear.F90
	interp_pchip.F90
	fitpack.f90
	search.F90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common core fitpack)
