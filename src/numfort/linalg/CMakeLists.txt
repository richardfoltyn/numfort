# src/numfort/numfort/arrays/CMakeLists.txt

set(COMPONENT_NAME linalg)

set(SOURCE_FILES
	linalg.f90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
