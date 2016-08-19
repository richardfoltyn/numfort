# src/numfort/numfort/interp/CMakeLists.txt

set(COMPONENT_NAME common)

set(SOURCE_FILES
	common.f90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
set_target_properties(${COMPONENT_NAME}
	PROPERTIES COMPILE_FLAGS "${${PROJECT_NAME}_FFLAGS}"
)
