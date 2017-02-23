# CMake configuration for numfort examples

# List of subdirectories to search for example files
list(APPEND _MODULES interpolate optimize)

include_directories("${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}")

find_package(OpenMP)

# Add all tests as separate executables
foreach (_dir IN LISTS _MODULES)
    # subdirectory stores source files (without directory prefix) in
    # the following variable
    unset(${_dir}_FILES)
    add_subdirectory(${_dir})

	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
		"${${PROJECT_NAME}_EXAMPLES_DIRECTORY}/${_dir}"
	)

	foreach (_file IN LISTS ${_dir}_FILES)
		get_filename_component(exe_name "${_file}" NAME_WE)
		# set target name that does not possibly conflict with other stuff
		set(_target example_${_dir}_${exe_name})
		add_executable(${_target} "${_dir}/${_file}")
		target_link_libraries(${_target} ${${CMAKE_PROJECT_NAME}_LIBRARIES})
		set_target_properties(${_target} PROPERTIES
			OUTPUT_NAME ${exe_name}
		)
	endforeach ()

    # Process example files that require OpenMP
    if (OPENMP_FOUND)
        foreach (_file IN LISTS ${_dir}_FILES_OMP)
    		get_filename_component(exe_name "${_file}" NAME_WE)
    		# set target name that does not possibly conflict with other stuff
    		set(_target example_${_dir}_${exe_name})
    		add_executable(${_target} "${_dir}/${_file}")
    		target_link_libraries(${_target} ${${CMAKE_PROJECT_NAME}_LIBRARIES})
    		set_target_properties(${_target} PROPERTIES
    			OUTPUT_NAME ${exe_name}
    		)
            target_compile_options(${_target} PRIVATE ${OpenMP_Fortran_FLAGS})
    	endforeach ()
    endif()
endforeach ()
