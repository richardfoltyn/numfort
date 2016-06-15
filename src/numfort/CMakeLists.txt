# src/numfort/CMakeLists.txt
# CMake configuration for numfort library

unset(SOURCE_FILES)

# include relevant subdirectories
add_subdirectory(arrays)

set(SOURCE_FILES 
	${SOURCE_FILES}
	numfort.f90 
	numfort_linalg.f90)

# always build additional static library to be used when linking tests and benchmarks
add_library(numfort_static STATIC ${SOURCE_FILES})
add_library(numfort ${SOURCE_FILES})

if (NOT USE_MKL)
	target_link_libraries(numfort ${NUMFORT_LAPACK_LIBRARIES})
endif (NOT USE_MKL)

# set mod-file output directory
set_target_properties(numfort PROPERTIES Fortran_MODULE_DIRECTORY "${NUMFORT_MODDIR}")
# set_target_properties(numfort_static PROPERTIES Fortran_MODULE_DIRECTORY "${NUMFORT_MODDIR}")
