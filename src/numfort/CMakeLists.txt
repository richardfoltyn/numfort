# src/numfort/CMakeLists.txt
# CMake configuration for numfort library

set(SOURCE_FILES 
	${SOURCE_FILES}
	${CMAKE_CURRENT_SOURCE_DIR}/numfort.f90 
	${CMAKE_CURRENT_SOURCE_DIR}/numfort_linalg.f90)

add_library(numfort_static STATIC ${SOURCE_FILES})
add_library(numfort SHARED ${SOURCE_FILES})

if (NOT WIN32)
	target_link_libraries(numfort ${NUMFORT_LAPACK_LIBRARIES})
endif (NOT WIN32)

# set mod-file output directory to subproject-specific folder
set_target_properties(numfort PROPERTIES Fortran_MODULE_DIRECTORY "${NUMFORT_MODDIR}")
set_target_properties(numfort_static PROPERTIES Fortran_MODULE_DIRECTORY "${NUMFORT_MODDIR}")
