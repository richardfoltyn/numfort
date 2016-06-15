# src/tests/CMakeLists.txt
# Configuration file for numfort unit tests

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

set(SOURCE_FILES 
	main.f90
	tests_numfort.f90
	tests_cumsum.f90)

add_executable(numfort_tests ${SOURCE_FILES})
target_link_libraries(numfort_tests numfort_static)

if (NOT USE_MKL)
	target_link_libraries(numfort_tests ${NUMFORT_LAPACK_LIBRARIES})
endif (NOT USE_MKL)
