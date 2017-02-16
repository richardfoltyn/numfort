
set(SOURCE_FILES
	driver1.f
	driver2.f
	driver3.f
	driver1.f90
	driver2.f90
	driver3.f90
)

# set output directory to test directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${${PROJECT_NAME}_EXAMPLES_DIRECTORY}")

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

	set(_target lbfgsb_bmnz_${exe_name})
	add_executable(${_target} "${_file}"
		$<TARGET_OBJECTS:${${PROJECT_NAME}_LIBRARY_NAME}>
	)
	target_link_libraries(${_target} ${${PROJECT_NAME}_LIBRARIES})
	set_target_properties(${_target} PROPERTIES OUTPUT_NAME ${exe_name})

endforeach (_file)
