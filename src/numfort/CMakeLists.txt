# src/numfort/CMakeLists.txt
# CMake configuration for numfort library

unset(SOURCE_FILES)

# include relevant subdirectories
add_subdirectory(arrays)
add_subdirectory(interp)

set(SOURCE_FILES
	${SOURCE_FILES}
	numfort.f90
	numfort_linalg.f90)

# always build additional static library to be used when linking tests and benchmarks

if (BUILD_TESTS OR BUILD_BENCHMARKS)
	add_library(${LIBRARY_NAME}_static STATIC ${SOURCE_FILES})
	set_target_properties(${LIBRARY_NAME}_static PROPERTIES OUTPUT_NAME ${LIBRARY_NAME})
endif()

add_library(${LIBRARY_NAME} ${SOURCE_FILES})
target_link_libraries(${LIBRARY_NAME} ${NUMFORT_LAPACK_LIBRARIES})

# set mod-file output directory
set_target_properties(${LIBRARY_NAME} PROPERTIES
    Fortran_MODULE_DIRECTORY "${NUMFORT_MODDIR}"
)

install(TARGETS ${LIBRARY_NAME}
    RUNTIME DESTINATION "${INSTALL_LIBRARY_PREFIX}"
    LIBRARY DESTINATION "${INSTALL_LIBRARY_PREFIX}"
    ARCHIVE DESTINATION "${INSTALL_LIBRARY_PREFIX}"
)

install(DIRECTORY ${NUMFORT_MODDIR}/ DESTINATION "${INSTALL_INCLUDE_PREFIX}")
