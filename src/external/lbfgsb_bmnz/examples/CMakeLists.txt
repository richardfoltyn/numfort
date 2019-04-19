# Example drivers
# Note: These are the original example files shipped with the L-BFGS-B
# code. They do not use modules and instead call the global function SETULB,
# so we compile them in the original way instead of modifying the code
# to use the module.

set(SOURCE_FILES
	driver1.f
	driver2.f
	driver3.f
	driver1.f90
	driver2.f90
	driver3.f90
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
	get_filename_component(_ext "${_file}" EXT)

	# add different suffixes based on Fortran 77/90 to avoid name collissions
	if (_ext STREQUAL ".f90")
		set(exe_name ${exe_name}_F90)
	else (_ext STREQUAL ".f")
		set(exe_name ${exe_name}_F77)
	endif ()

	# Additional source files required for each driver.
	# We use whatever is shipped with the original source code instead of
	# linking to BLAS, etc.
	set(SOURCE_FILES
		${_file}
		${CMAKE_CURRENT_SOURCE_DIR}/../lbfgsb.f
		${CMAKE_CURRENT_SOURCE_DIR}/../linpack.f
        ${CMAKE_CURRENT_SOURCE_DIR}/../timer.f
		${CMAKE_CURRENT_SOURCE_DIR}/../blas.f
	)

	set(_target ${PROJEC_NAME}_${exe_name})
    add_executable(${_target} ${SOURCE_FILES})

	target_link_libraries(${_target} ${${PROJECT_NAME}_LIBRARIES})
	set_target_properties(${_target} PROPERTIES OUTPUT_NAME ${exe_name})

endforeach (_file)
