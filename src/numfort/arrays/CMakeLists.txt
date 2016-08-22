# src/numfort/numfort/arrays/CMakeLists.txt

set(COMPONENT_NAME arrays)

set(SOURCE_FILES
	arrays.f90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
