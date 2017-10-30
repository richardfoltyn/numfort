# src/numfort/CMakeLists.txt
# CMake configuration for numfort library

project(numfort C Fortran)

################################################################################
# Additional compiler flags for core numfort library (can be more restrictive
# than flags for external code, examples, unit tests, etc.

unset(_FFLAGS)
unset(_FFLAGS_DEBUG)
unset(_FFLAGS_RELEASE)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    # explicitly add -cpp to run C preprocessor if required
    set(GNU_WARN "-Wextra -Warray-temporaries -Wrealloc-lhs -pedantic -Wsurprising")
    set(_FFLAGS "${GNU_WARN} -std=f2008")
    set(_FFLAGS_DEBUG "-save-temps")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (WIN32)
        set(_FFLAGS "/warn:all /standard-semantics /stand:f08")
        set(_FFLAGS_DEBUG "/Qsave-temps")
        #        set(_FLAGS_RELEASE "/Qipo")

        # NOTE: When setting -standard-semantics, resets any previously
        # given -assume std_mod_proc_name, so we need to append it again!
        if (NOT IFORT_STD_MOD_PROC_NAME)
            set(_FFLAGS "${_FFLAGS} /assume:nostd_mod_proc_name")
        else ()
            set(_FFLAGS "${_FFLAGS} /assume:std_mod_proc_name")
        endif ()
    else ()
        set(_FFLAGS "-warn all -standard-semantics -std08")
        set(_FFLAGS_DEBUG "-save-temps")
        #set(_FFLAGS_RELEASE "-ipo")

        if (NOT IFORT_STD_MOD_PROC_NAME)
            set(_FFLAGS "${_FFLAGS} -assume nostd_mod_proc_name")
        else ()
            set(_FFLAGS "${_FFLAGS} -assume std_mod_proc_name")
        endif ()
    endif (WIN32)
endif ()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_FFLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${_FFLAGS_DEBUG}")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${_FFLAGS_RELEASE}")

#set_property(DIRECTORY PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)

string(TOUPPER "${CMAKE_BUILD_TYPE}" _CONFIG)
message(STATUS "Flags for NUMFORT library: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${_CONFIG}}")

################################################################################
# Target definitions

# add individual numfort components as object libraries
set(COMPONENTS
    core
    common
    arrays
    interpolate
    integrate
    linalg
    optimize
    stats
)

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
