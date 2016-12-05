# src/benchmarks/CMakeLists.txt
# CMake configuration for numfort benchmarks

include_directories(${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS})

list(APPEND COMPONENTS stats)
set(LIBS ${${CMAKE_PROJECT_NAME}_LIBRARIES})

# keep only child items that are directories
foreach (_dir IN LISTS COMPONENTS)
	# subdirectory stores source files (without directory prefix) in
	# the following variable
	unset(${_dir}_FILES)
	add_subdirectory(${_dir})

	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
		"${${PROJECT_NAME}_BENCHMARKS_DIRECTORY}/${_dir}"
	)

	foreach (test_source IN LISTS ${_dir}_FILES)
        get_filename_component(exe_name "${test_source}" NAME_WE)
        set(target_name bm_${_dir}_${exe_name})
    	add_executable(${target_name} "${_dir}/${test_source}")
    	target_link_libraries(${target_name} ${LIBS})
        set_target_properties(${target_name} PROPERTIES
            OUTPUT_NAME ${exe_name}
        )
    endforeach()

endforeach()
