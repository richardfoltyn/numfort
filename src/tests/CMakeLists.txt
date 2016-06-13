# src/tests/CMakeLists.txt
# Configuration file for numfort unit tests

set(SOURCE_FILES 
	main.f90
	tests_numfort.f90
	tests_cumsum.f90)

add_executable(numfort_tests ${SOURCE_FILES})
target_link_libraries(numfort_tests numfort_static)

if (NOT WIN32)
	target_link_libraries(numfort_tests ${NUMFORT_LAPACK_LIBRARIES})
endif (NOT WIN32)