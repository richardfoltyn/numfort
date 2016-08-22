# src/numfort/numfort/interp/CMakeLists.txt

set(COMPONENT_NAME core)

set(SOURCE_FILES
	core.f90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
