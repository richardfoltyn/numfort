# src/fitpack/fitpack/CMakeLists.txt
# CMake configuration for numfort library

set(SOURCE_FILES
	t_minim.f90
)

# set output directory to test directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${${PROJECT_NAME}_EXAMPLES_DIRECTORY}")

# include_directories("${_MODDIR}")

# Add all tests as separate executables
foreach (_file IN LISTS SOURCE_FILES)
	get_filename_component(exe_name "${_file}" NAME_WE)
	set(_target csiro_${exe_name})

	add_executable(${_target} "${_file}")
	target_link_libraries(${_target} ${${PROJECT_NAME}_LIBRARIES})
	set_target_properties(${_target} PROPERTIES OUTPUT_NAME ${exe_name})

endforeach (_file)
