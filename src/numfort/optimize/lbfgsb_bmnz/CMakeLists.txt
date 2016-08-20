
set(F77_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/lbfgsb.f
    ${CMAKE_CURRENT_SOURCE_DIR}/linpack.f
    ${CMAKE_CURRENT_SOURCE_DIR}/timer.f
)

set(SOURCE_FILES
    ${SOURCE_FILES}
    ${F77_SOURCE_FILES}
    PARENT_SCOPE
)

# set less restrictive compile options for external code
set_source_files_properties(${F77_SOURCE_FILES}
    PROPERTIES COMPILE_FLAGS "${${PROJECT_NAME}_FFLAGS_EXTERNAL}"
)

if (BUILD_EXAMPLES)
    set(EXAMPLE_SOURCE_FILES
    	driver1.f
    	driver2.f
    	driver3.f
    	driver1.f90
    	driver2.f90
    	driver3.f90
    )

    # set output directory to example directory
    set(_OUTPUT_DIR "${${PROJECT_NAME}_EXAMPLES_DIRECTORY}/optimize/lbfgsb-bmnz")

    # Add all tests as separate executables
    foreach (_file IN LISTS EXAMPLE_SOURCE_FILES)
    	get_filename_component(exe_name "${_file}" NAME_WE)
    	get_filename_component(_ext "${_file}" EXT)

    	# add different suffixes based on Fortran 77/90 to avoid name collissions
    	if (_ext STREQUAL ".f90")
    		set(exe_name ${exe_name}_F90)
    	else (_ext STREQUAL ".f")
    		set(exe_name ${exe_name}_F77)
    	endif ()

        set(_target lbfgsb_bmnz_${exe_name})
    	add_executable(${_target} "${_file}")
    	target_link_libraries(${_target}
    		${LIBRARY_NAME}_static ${${PROJECT_NAME}_BLA_LIBRARIES}
    	)
        set_target_properties(${_target} PROPERTIES
            OUTPUT_NAME ${exe_name}
            RUNTIME_OUTPUT_DIRECTORY "${_OUTPUT_DIR}"
        )
    endforeach ()
endif()
