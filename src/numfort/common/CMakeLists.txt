# src/numfort/numfort/interp/CMakeLists.txt

set(COMPONENT_NAME common)

set(SOURCE_FILES
	common.f90
	enums.f90
	kinds.f90
	status.f90
	workspace.F90
    input_checks.f90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
