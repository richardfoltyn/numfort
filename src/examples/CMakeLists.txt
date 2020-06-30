# CMake configuration for numfort examples
find_package(OpenMP)

if (OPENMP_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

################################################################################
# Additional compiler flags for examples

unset(_FFLAGS)
unset(_FFLAGS_DEBUG)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    # explicitly add -cpp to run C preprocessor if required
    set(GNU_WARN "-Wall -pedantic -Wuninitialized -Wno-unused")

    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 6)
        # breaks compilation with gfortran 5.x, possibly with 6.x as well.
        set(GNU_WARN "${GNU_WARN} -Wno-unused-dummy-argument")
    endif()
    # Disable "surprising" warnings for older GFORTRAN version as this
    # produces endless warnings about array finalizers
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 7)
        set(GNU_WARN "${GNU_WARN} -Wno-surprising")
    endif()
    set(_FFLAGS "${GNU_WARN} --std=f2008")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (WIN32)
        set(_FFLAGS "/warn:all /standard-semantics /stand:f08")
        # NOTE: When setting -standard-semantics, resets any previously
        # given -assume std_mod_proc_name, so we need to append it again!
        if (NOT IFORT_STD_MOD_PROC_NAME)
            set(_FFLAGS "${_FFLAGS} /assume:nostd_mod_proc_name")
        else ()
            set(_FFLAGS "${_FFLAGS} /assume:std_mod_proc_name")
        endif ()
    else ()
        set(_FFLAGS "-warn all -standard-semantics -std08")
        # NOTE: When setting -standard-semantics, resets any previously
        # given -assume std_mod_proc_name, so we need to append it again!
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

message(STATUS "Building examples")

string(TOUPPER "${CMAKE_BUILD_TYPE}" _CONFIG)
message(STATUS "Flags for EXAMPLES: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${_CONFIG}}")

################################################################################



# List of subdirectories to search for example files
list(APPEND _MODULES interpolate optimize stats)


# Add all tests as separate executables
foreach (_comp IN LISTS _MODULES)
    # subdirectory stores source files (without directory prefix) in
    # the following variable
    unset(${_comp}_FILES)
    add_subdirectory(${_comp})

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
        "${NUMFORT_EXAMPLES_DIRECTORY}/${_comp}"
    )

    foreach (_file IN LISTS ${_comp}_FILES)
        get_filename_component(exe_name "${_file}" NAME_WE)
        # set target name that does not possibly conflict with other stuff
        set(_target example_${_comp}_${exe_name})
        add_executable(${_target} "${_comp}/${_file}")

        target_link_libraries(${_target} PRIVATE
            ${NUMFORT_LIBRARY_NAME}
            ${BLA_LIBRARIES}
        )

        if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/${_comp}/include")
            target_include_directories(${_target} PRIVATE "${_comp}/include")
        endif()

        target_include_directories(${_target} PRIVATE
            "${NUMFORT_HEADER_INCLUDE_DIR}"
            "${NUMFORT_MODULE_DIRECTORY}"
            ${BLA_INCLUDE_DIRS}
        )

        set_target_properties(${_target} PROPERTIES
            OUTPUT_NAME ${exe_name}
        )
    endforeach ()

endforeach ()


if (OPENMP_FOUND)
    add_subdirectory(openmp)
endif()
