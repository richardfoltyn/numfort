# src/external/simplex_csiro/examples/CMakeLists.txt

set(SOURCE_FILES
	t_minim.f90
)

# set output directory to examples directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${EXAMPLES_DIRECTORY}")

include_directories("${MODULE_DIRECTORY}")

if (CMAKE_Fortran_FLAGS_ORIG)
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_ORIG}")
endif()

if (CMAKE_Fortran_FLAGS_DEBUG_ORIG)
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG_ORIG}")
endif()

# Add all tests as separate executables
foreach (_file IN LISTS SOURCE_FILES)
	get_filename_component(exe_name "${_file}" NAME_WE)
	set(_target ${PROJECT_NAME}_${exe_name})

	if (BUILD_STANDALONE)
		add_executable(${_target} "${_file}")
		target_link_libraries(${_target} ${LIBRARY_NAME})
	else ()
		add_executable(${_target} "${_file}"
			$<TARGET_OBJECTS:${LIBRARY_NAME}>
		)
	endif ()

	set_target_properties(${_target} PROPERTIES OUTPUT_NAME ${exe_name})

endforeach (_file)
