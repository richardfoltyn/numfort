# src/tests/CMakeLists.txt
# Configuration file for numfort unit tests

set(TEST_DIRS arrays)
set(TEST_SOURCE_FILES_ARRAYS arrays/indexing.f90)

set(LIBS ${${CMAKE_PROJECT_NAME}_LIBRARIES} ${FCore_LIBRARIES})

include_directories(${FCore_INCLUDE_DIRS} ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS})

foreach (dir IN LISTS TEST_DIRS)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
		"${${CMAKE_PROJECT_NAME}_TESTS_DIRECTORY}/${dir}"
	)

    string(TOUPPER "${dir}" DIR_UPPER)
    set(TEST_SOURCE_FILES ${TEST_SOURCE_FILES_${DIR_UPPER}})

    foreach (test_source IN LISTS TEST_SOURCE_FILES)
        get_filename_component(stub "${test_source}" NAME_WE)
        set(exe_name ${dir}_${stub})
    	add_executable(${exe_name} "${test_source}")
    	target_link_libraries(${exe_name} ${LIBS})
    	add_test(NAME ${exe_name} COMMAND ${exe_name}
            WORKING_DIRECTORY "${${PROJECT_NAME}_TESTS_DIRECTORY}/${dir}")
    endforeach()
endforeach()
