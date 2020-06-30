# src/numfort/numfort/arrays/CMakeLists.txt

set(COMPONENT_NAME linalg)

set(SOURCE_FILES
    linalg.f90
	linalg_core.f90
	linalg_core_real32.F90
	linalg_core_real64.F90
	lawson_hanson_common.f90
	lawson_hanson_real64.F90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} interfaces common)
