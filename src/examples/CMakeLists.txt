# CMake configuration for numfort examples

################################################################################
# Additional compiler flags for examples

unset(_FFLAGS)
unset(_FFLAGS_DEBUG)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    # explicitly add -cpp to run C preprocessor if required
    set(_FFLAGS "-std=f2008")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (WIN32)
        set(_FFLAGS "/standard-semantics /stand:f08")
        if (NOT IFORT_STD_MOD_PROC_NAME)
            set(_FFLAGS "${_FFLAGS} /assume:nostd_mod_proc_name")
        else ()
            set(_FFLAGS "${_FFLAGS} /assume:std_mod_proc_name")
        endif ()
    else ()
        set (_FFLAGS "-standard-semantics -std08")
        if (NOT IFORT_STD_MOD_PROC_NAME)
            set(_FFLAGS "${_FFLAGS} -assume nostd_mod_proc_name")
        else ()
            set(_FFLAGS "${_FFLAGS} -assume std_mod_proc_name")
        endif ()
    endif (WIN32)
endif ()

# add standard semantics to compiler flags even for F77 files, in case
# these reference any submodules
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_FFLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${_FFLAGS_DEBUG}")

string(TOUPPER "${CMAKE_BUILD_TYPE}" _CONFIG)
message(STATUS "Flags for EXAMPLES: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${_CONFIG}}")

################################################################################



# List of subdirectories to search for example files
list(APPEND _MODULES interpolate optimize stats)

include_directories("${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}")

# Add all tests as separate executables
foreach (_dir IN LISTS _MODULES)
    # subdirectory stores source files (without directory prefix) in
    # the following variable
    unset(${_dir}_FILES)
    add_subdirectory(${_dir})

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
        "${${PROJECT_NAME}_EXAMPLES_DIRECTORY}/${_dir}"
    )

    foreach (_file IN LISTS ${_dir}_FILES)
        get_filename_component(exe_name "${_file}" NAME_WE)
        # set target name that does not possibly conflict with other stuff
        set(_target example_${_dir}_${exe_name})
        add_executable(${_target} "${_dir}/${_file}")
        target_link_libraries(${_target} ${${CMAKE_PROJECT_NAME}_LIBRARIES})
        
        # Add additional include directories if applicable
        set(_incdir "${CMAKE_CURRENT_SOURCE_DIR}/${_dir}/include")
        if (IS_DIRECTORY "${_incdir}")
            target_include_directories(${_target} PUBLIC "${_dir}/include")
        endif()
        
        set_target_properties(${_target} PROPERTIES
            OUTPUT_NAME ${exe_name}
        )
    endforeach ()

endforeach ()

find_package(OpenMP)

if (OPENMP_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    add_subdirectory(openmp)
endif()
