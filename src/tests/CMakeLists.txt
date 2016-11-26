# src/tests/CMakeLists.txt
# Configuration file for numfort unit tests

find_package(FCorelib REQUIRED)

# get subdirectories
file(GLOB _children RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/*
)

# keep only child items that are directories
foreach (_c IN LISTS _children)
	if (IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${_c}")
		list(APPEND _COMPONENTS "${_c}")
	endif ()
endforeach()

set(LIBS ${${CMAKE_PROJECT_NAME}_LIBRARIES} ${FCorelib_LIBRARIES})

include_directories(${FCorelib_INCLUDE_DIRS} ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS})

foreach (dir IN LISTS _COMPONENTS)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
		"${${CMAKE_PROJECT_NAME}_TESTS_DIRECTORY}/${dir}"
	)

    unset(SOURCE_FILES)
    file(GLOB SOURCE_FILES RELATIVE "${CMAKE_CURRENT_LIST_DIR}" "${dir}/*.f90")

    foreach (test_source IN LISTS SOURCE_FILES)
        get_filename_component(exe_name "${test_source}" NAME_WE)
        set(target_name ${dir}_${exe_name})
    	add_executable(${target_name} "${test_source}")
    	target_link_libraries(${target_name} ${LIBS})
        set_target_properties(${target_name} PROPERTIES
            OUTPUT_NAME ${exe_name}
        )
    	add_test(NAME ${exe_name} COMMAND ${exe_name}
            WORKING_DIRECTORY "${${PROJECT_NAME}_TESTS_DIRECTORY}/${dir}"
        )
    endforeach()
endforeach()
