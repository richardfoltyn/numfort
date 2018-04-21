# src/numfort/numfort/interp/CMakeLists.txt

set(COMPONENT_NAME interpolate)

set(SOURCE_FILES
	interpolate.f90
	common.f90
	interp_result.f90
    interp_linear.F90
	fitpack.f90
)

#include_directories(${${PROJECT_NAME}_MODDIR})
include_directories("${CMAKE_CURRENT_LIST_DIR}/include")
add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common core fitpack)
