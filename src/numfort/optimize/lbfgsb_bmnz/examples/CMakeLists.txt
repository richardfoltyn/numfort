# src/fitpack/fitpack/CMakeLists.txt
# CMake configuration for numfort library

set(SOURCE_FILES
	driver1.f
	driver2.f
	driver3.f
	driver1.f90
	driver2.f90
	driver3.f90
)

include_directories("${${PROJECT_NAME}_MODDIR}")

# set output directory to test directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
	"${${PROJECT_NAME}_EXAMPLES_DIRECTORY}/optimize/lbfgsb-bmnz"
)

# Add all tests as separate executables
foreach (_file IN LISTS SOURCE_FILES)
	get_filename_component(exe_name "${_file}" NAME_WE)
	get_filename_component(_ext "${_file}" EXT)

	# add different suffixes based on Fortran 77/90 to avoid name collissions
	if (_ext STREQUAL ".f90")
		set(exe_name ${exe_name}_F90)
	else (_ext STREQUAL ".f")
		set(exe_name ${exe_name}_F77)
	endif ()

	add_executable("${exe_name}" "${_file}")
	target_link_libraries("${exe_name}"
		${LIBRARY_NAME}_static ${${PROJECT_NAME}_BLA_LIBRARIES}
	)
endforeach ()

# copy sample data from samples/ subdirectory to test directory
# file(GLOB FITPACK_SAMPLE_DATA "../samples/*.txt")
# file(COPY ${FITPACK_SAMPLE_DATA} DESTINATION "${${PROJECT_NAME}_TESTS_DIRECTORY}")
