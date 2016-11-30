# src/numfort/numfort/arrays/CMakeLists.txt

set(COMPONENT_NAME arrays)

set(SOURCE_FILES
	arrays.f90
	arrays_sort.f90
	arrays_index.f90
	arrays_create.f90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common orderpack)
