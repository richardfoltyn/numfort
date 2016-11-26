# src/tests/CMakeLists.txt
# Configuration file for numfort unit tests

find_package(FCorelib REQUIRED)

set(TEST_DIRS arrays interpolate)

set(LIBS ${${CMAKE_PROJECT_NAME}_LIBRARIES} ${FCorelib_LIBRARIES})

include_directories(${FCorelib_INCLUDE_DIRS} ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS})

foreach (dir IN LISTS TEST_DIRS)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
		"${${CMAKE_PROJECT_NAME}_TESTS_DIRECTORY}/${dir}"
	)

    unset(SOURCE_FILES)
    file(GLOB SOURCE_FILES RELATIVE "${CMAKE_CURRENT_LIST_DIR}" "${dir}/*.f90")

    foreach (test_source IN LISTS SOURCE_FILES)
        get_filename_component(stub "${test_source}" NAME_WE)
        set(exe_name ${dir}_${stub})
    	add_executable(${exe_name} "${test_source}")
    	target_link_libraries(${exe_name} ${LIBS})
    	add_test(NAME ${exe_name} COMMAND ${exe_name}
            WORKING_DIRECTORY "${${PROJECT_NAME}_TESTS_DIRECTORY}/${dir}")
    endforeach()
endforeach()
