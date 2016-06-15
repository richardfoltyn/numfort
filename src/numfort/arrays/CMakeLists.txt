# src/numfort/numfort/arrays/CMakeLists.txt

set(SOURCE_FILES
	${SOURCE_FILES}
	${CMAKE_CURRENT_SOURCE_DIR}/array_creation.f90
	${CMAKE_CURRENT_SOURCE_DIR}/array_indexing.f90
	PARENT_SCOPE)
	
# set_source_files_properties(array_creation.f90 PROPERTIES 
#	Fortran_MODULE_DIRECTORY "${NUMFORT_MODDIR}")