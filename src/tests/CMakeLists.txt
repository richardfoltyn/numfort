# src/tests/CMakeLists.txt
# Configuration file for numfort unit tests

find_package(FCorelib REQUIRED)
set(LIBS ${${CMAKE_PROJECT_NAME}_LIBRARIES} ${FCorelib_LIBRARIES})
include_directories(${FCorelib_INCLUDE_DIRS} ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS})

list(APPEND COMPONENTS arrays interpolate stats)

# keep only child items that are directories
foreach (_dir IN LISTS COMPONENTS)
	# subdirectory stores source files (without directory prefix) in
	# the following variable
	unset(${_dir}_FILES)
	add_subdirectory(${_dir})

	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
		"${${PROJECT_NAME}_TESTS_DIRECTORY}/${_dir}"
	)

	foreach (test_source IN LISTS ${_dir}_FILES)
        get_filename_component(exe_name "${test_source}" NAME_WE)
        set(target_name ${_dir}_${exe_name})
    	add_executable(${target_name} "${_dir}/${test_source}")
    	target_link_libraries(${target_name} ${LIBS})
        set_target_properties(${target_name} PROPERTIES
            OUTPUT_NAME ${exe_name}
        )
    	add_test(NAME ${exe_name} COMMAND ${exe_name}
            WORKING_DIRECTORY "${${PROJECT_NAME}_TESTS_DIRECTORY}/${_dir}"
        )
    endforeach()

endforeach()

add_subdirectory(openmp)
