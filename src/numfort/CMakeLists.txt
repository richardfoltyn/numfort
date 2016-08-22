# src/numfort/CMakeLists.txt
# CMake configuration for numfort library

project(numfort C Fortran)

################################################################################
# Compiler FLAGS

set(_FFLAGS)
set(_FFLAGS_DEBUG)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    # explicitly add -cpp to run C preprocessor if required
    set(GNU_WARN "-Wall -Wextra -Wimplicit-interface -Wimplicit-procedure \
      -Warray-temporaries -Wrealloc-lhs -pedantic")
    set(_FFLAGS "${GNU_WARN} -std=f2008")
    set(_FFLAGS_DEBUG "-fcheck=all -Og")
    # enable optimization that does not break debugging
    # set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Og")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (WIN32)
        set (_FFLAGS "/warn:all /standard-semantics /stand:f08")
        set(_FFLAGS_DEBUG "/check:all")
    else ()
        set (_FFLAGS "-warn all -standard-semantics -std08")
        set(_FFLAGS_DEBUG "-check all")
    endif (WIN32)
endif ()

message(STATUS "Default compiler flags: ${CMAKE_Fortran_FLAGS}")
message(STATUS "Project-specific compiler flags: ${_FFLAGS}")
# Make sure build type is upper case, otherwise the wrong flags will be shown
string(TOUPPER "${CMAKE_BUILD_TYPE}" _CONFIG)
message(STATUS "Compiler flags for ${CMAKE_BUILD_TYPE}: ${CMAKE_Fortran_FLAGS_${_CONFIG}}")

# set FFLAGS for all code in numfort/* directories
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_FFLAGS}")
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${_FFLAGS_DEBUG}")

################################################################################
# Target definitions

# add individual numfort components as object libraries
set(COMPONENTS core common arrays interpolate linalg optimize)
# set(COMPONENTS core arrays interpolate linalg)

include_directories("${${PROJECT_NAME}_MODDIR}")
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

if (BUILD_TESTS)
    add_subdirectory(tests)
endif()

if (BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()
