# src/tests/CMakeLists.txt
# Configuration file for numfort unit tests

find_package(fcore REQUIRED)

################################################################################
# Additional compiler flags for unit tests

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

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_FFLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${_FFLAGS_DEBUG}")

message(STATUS "Building unit tests")

string(TOUPPER "${CMAKE_BUILD_TYPE}" _CONFIG)
message(STATUS "Flags for TESTS: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${_CONFIG}}")


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

if (BUILD_BLAS95)
    list(APPEND COMPONENTS wrappers)
endif()

# keep only child items that are directories
foreach (_comp IN LISTS COMPONENTS)
    # subdirectory stores source files (without directory prefix) in
    # the following variable
    unset(${_comp}_FILES)
    add_subdirectory(${_comp})

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
        "${NUMFORT_TESTS_DIRECTORY}/${_comp}"
    )

    # Collect additional source files for this component and prepend
    # component directory
    if (${_comp}_SHARED_FILES)
        foreach (f IN LISTS ${_comp}_SHARED_FILES)
            list(APPEND ${_comp}_aux_source "${_comp}/${f}")
        endforeach()
    endif()

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

        target_link_libraries(${target_name}
            ${NUMFORT_LIBRARY_NAME}
            fcore::fcore
            ${BLA_LIBRARIES}
        )
        target_include_directories(${target_name} PRIVATE
            "${NUMFORT_HEADER_INCLUDE_DIR}"
            "${NUMFORT_MODULE_DIRECTORY}"
            "${BLA_INCLUDE_DIRS}"
        )

        set_target_properties(${target_name} PROPERTIES
            OUTPUT_NAME ${exe_name}
        )

        add_test(NAME ${exe_name} COMMAND ${exe_name}
            WORKING_DIRECTORY "${NUMFORT_TESTS_DIRECTORY}/${_comp}"
        )
    endforeach()

endforeach()

add_subdirectory(openmp)
