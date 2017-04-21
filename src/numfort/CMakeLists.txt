# src/numfort/CMakeLists.txt
# CMake configuration for numfort library

project(numfort C Fortran)

################################################################################
# Target definitions

# add individual numfort components as object libraries
set(COMPONENTS core common arrays interpolate linalg optimize stats)
# set(COMPONENTS core arrays interpolate linalg)

foreach (_comp IN LISTS COMPONENTS)
	add_subdirectory(${_comp})
	# set mod-file output directory
	set_target_properties(${_comp} PROPERTIES
	    Fortran_MODULE_DIRECTORY "${${PROJECT_NAME}_MODDIR}"
	)

    # we need to make sure that -fPIC etc. is passed to compiler when building
    # object libraries that will be used to create a shared lib
    if (BUILD_SHARED_LIBS)
        set_target_properties(${_comp} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif ()

	# append to list of object library targets
	set_property(GLOBAL APPEND
		PROPERTY ${CMAKE_PROJECT_NAME}_COMPONENT_TARGETS ${_comp}
	)
endforeach()
