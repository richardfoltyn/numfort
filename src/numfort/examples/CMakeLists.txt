# src/fitpack/fitpack/CMakeLists.txt
# CMake configuration for numfort library

# get subdirectories
file(GLOB _children RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/*
)

# keep only child items that are directories
foreach (_c IN LISTS _children)
	if (IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${_c}")
		list(APPEND _MODULES "${_c}")
	endif ()
endforeach()

include_directories("${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}")

# Add all tests as separate executables
foreach (_dir IN LISTS _MODULES)

	file(GLOB _FILES "${_dir}/*.f90")

	set (_OUTPUT_DIR "${${PROJECT_NAME}_EXAMPLES_DIRECTORY}/${_dir}")

	foreach (_file IN LISTS _FILES)
		get_filename_component(exe_name "${_file}" NAME_WE)
		# set target name that does not possibly conflict with other stuff
		set(_target example_${_dir}_${exe_name})
		add_executable(${_target} "${_file}")
		set_target_properties(${_target} PROPERTIES
			RUNTIME_OUTPUT_DIRECTORY "${_OUTPUT_DIR}"
			OUTPUT_NAME ${exe_name}
		)
		target_link_libraries(${_target} ${${CMAKE_PROJECT_NAME}_LIBRARIES})
	endforeach ()
endforeach ()

# copy sample data from samples/ subdirectory to test directory
# file(GLOB FITPACK_SAMPLE_DATA "../samples/*.txt")
# file(COPY ${FITPACK_SAMPLE_DATA} DESTINATION "${${PROJECT_NAME}_TESTS_DIRECTORY}")
