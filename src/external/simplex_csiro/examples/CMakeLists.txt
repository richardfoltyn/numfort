# src/external/simplex_csiro/examples/CMakeLists.txt

set(SOURCE_FILES
	t_minim.f90
)

# set output directory to examples directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${EXAMPLES_DIRECTORY}")

include_directories("${MODULE_DIRECTORY}")

# Add all tests as separate executables
foreach (_file IN LISTS SOURCE_FILES)
	get_filename_component(exe_name "${_file}" NAME_WE)
	set(_target ${PROJECT_NAME}_${exe_name})

	add_executable(${_target} "${_file}"
		$<TARGET_OBJECTS:${LIBRARY_NAME}>
	)

	set_target_properties(${_target} PROPERTIES OUTPUT_NAME ${exe_name})

endforeach (_file)
