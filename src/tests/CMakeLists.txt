# src/tests/CMakeLists.txt
# Configuration file for numfort unit tests

find_package(fcore REQUIRED)

################################################################################
# Additional compiler flags for unit tests

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
message(STATUS "Flags for TESTS: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${_CONFIG}}")

################################################################################
# MKL include directories
# Unit test code is permitted to rely on MKL Fortran95 interfaces even though
# we avoid this dependency in the main library.

if (USE_MKL)
    include_directories(${MKL_BLAS95_INCLUDE_DIR} ${MKL_LAPACK95_INCLUDE_DIR})
endif ()


################################################################################


list(APPEND COMPONENTS
    arrays
    common
    core
    integrate
    interpolate
    io
    optimize
    polynomial
    stats
)

if (NOT USE_MKL)
    list(APPEND COMPONENTS wrappers)
endif()

# keep only child items that are directories
foreach (_comp IN LISTS COMPONENTS)
    # subdirectory stores source files (without directory prefix) in
    # the following variable
    unset(${_comp}_FILES)
    add_subdirectory(${_comp})

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
        "${${PROJECT_NAME}_TESTS_DIRECTORY}/${_comp}"
    )

    # Collect additional source files for this component and prepend
    # component directory
    if (${_comp}_SHARED_FILES)
        foreach (f IN LISTS ${_comp}_SHARED_FILES)
            list(APPEND ${_comp}_aux_source "${_comp}/${f}")
        endforeach()
    endif()

    # Create component-specific directory to store module files to
    # avoid potential naming collisions.
    set(CMAKE_Fortran_MODULE_DIRECTORY "${${PROJECT_NAME}_MODDIR}/tests/${_comp}")
    include_directories(BEFORE ${CMAKE_Fortran_MODULE_DIRECTORY})

    foreach (test_source IN LISTS ${_comp}_FILES)
        get_filename_component(exe_name "${test_source}" NAME_WE)
        set(target_name "test_${_comp}_${exe_name}")
        add_executable(${target_name}
            "${_comp}/${test_source}"
            ${${_comp}_aux_source}
        )

        if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/${_comp}/include")
            target_include_directories(${target_name} PRIVATE "${_comp}/include")
        endif()

        target_link_libraries(${target_name} ${${CMAKE_PROJECT_NAME}_LIBRARIES})
        target_link_libraries(${target_name} fcore::fcore)
        set_target_properties(${target_name} PROPERTIES
            OUTPUT_NAME ${exe_name}
        )

        add_test(NAME ${exe_name} COMMAND ${exe_name}
            WORKING_DIRECTORY "${${PROJECT_NAME}_TESTS_DIRECTORY}/${_comp}"
        )
    endforeach()

endforeach()

add_subdirectory(openmp)
