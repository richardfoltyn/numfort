# src/benchmarks/CMakeLists.txt
# CMake configuration for numfort benchmarks

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/benchmarks)

set(CUMSUM_SOURCE_FILES 
	cumsum.f90)
	
add_executable(cumsum ${CUMSUM_SOURCE_FILES})
target_link_libraries(cumsum numfort_static)

if (NOT USE_MKL)
	target_link_libraries(cumsum ${NUMFORT_LAPACK_LIBRARIES})
endif (NOT USE_MKL)
