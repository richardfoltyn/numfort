# src/numfort/numfort/interp/CMakeLists.txt

set(COMPONENT_NAME common)

set(SOURCE_FILES
	common.F90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
