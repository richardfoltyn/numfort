# CMake configuration for numfort examples

include_directories("${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}")

set(SOURCE_FILES
    concon.f90
    curfit.f90
    minimize_lbfgsb.f90
)

foreach (_file IN LISTS SOURCE_FILES)
	get_filename_component(exe_name "${_file}" NAME_WE)
	# set target name that does not possibly conflict with other stuff
	set(_target example_omp_${exe_name})
	add_executable(${_target} "${_file}")
	target_link_libraries(${_target} ${${CMAKE_PROJECT_NAME}_LIBRARIES})
	set_target_properties(${_target} PROPERTIES
		OUTPUT_NAME ${exe_name}
	)
endforeach ()
