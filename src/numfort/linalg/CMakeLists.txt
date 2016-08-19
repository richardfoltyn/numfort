# src/numfort/numfort/arrays/CMakeLists.txt

set(COMPONENT_NAME linalg)

set(SOURCE_FILES
	linalg.f90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})

# set_source_files_properties(array_creation.f90 PROPERTIES
#	Fortran_MODULE_DIRECTORY "${NUMFORT_MODDIR}")
