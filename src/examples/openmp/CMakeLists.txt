# CMake configuration for numfort examples

get_filename_component(_dir "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
include_directories("${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}")

set(SOURCE_FILES
    concon.f90
    curfit.f90
    minimize_lbfgsb.f90
    minimize_slsqp.f90
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    "${${PROJECT_NAME}_EXAMPLES_DIRECTORY}/${_dir}"
)

foreach (_file IN LISTS SOURCE_FILES)
	get_filename_component(exe_name "${_file}" NAME_WE)
	# set target name that does not possibly conflict with other stuff
	set(_target "example_${_dir}_${exe_name}")
	add_executable(${_target} "${_file}")
	target_link_libraries(${_target} ${${CMAKE_PROJECT_NAME}_LIBRARIES})
	set_target_properties(${_target} PROPERTIES
		OUTPUT_NAME ${exe_name}
	)
endforeach ()
