# Configuration file for numfort unit tests

find_package(fcore REQUIRED)

if (USE_OPENMP)
    find_package(OpenMP)

    if (OPENMP_FOUND)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    endif()
endif()

################################################################################
# Additional compiler flags for all unit tests

unset(_FFLAGS)
unset(_FFLAGS_DEBUG)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    set(GNU_WARN "-Wall -pedantic -Wuninitialized -Wno-unused")
    set(_FFLAGS "${GNU_WARN} --std=f2018 -fimplicit-none")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel(LLVM)?")
    if (WIN32)
        set(_FFLAGS "/warn:all,nounused,nousage /stand:f18")
    else ()
        set(_FFLAGS "-warn all,nounused,nousage -std18")
    endif (WIN32)
endif ()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_FFLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${_FFLAGS_DEBUG}")

message(STATUS "Building unit tests")


################################################################################

list(APPEND COMPONENTS
    arrays
    common
    core
    integrate
    interpolate
    io
    linalg
    optimize
    polynomial
    stats
)

if (BUILD_BLAS95)
    list(APPEND COMPONENTS wrappers)
endif()

foreach (_comp IN LISTS COMPONENTS)
    # subdirectory stores source files (without directory prefix) in
    # the following variable
    unset(${_comp}_FILES)
    add_subdirectory(${_comp})

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

        if (IS_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/${_comp}/include")
            target_include_directories(${target_name} PRIVATE "${_comp}/include")
        endif()

        target_link_libraries(${target_name} PRIVATE
            ${NUMFORT_LIBRARY_NAME}_static
            fcore::fcore
            ${BLA_LIBRARIES}
            ${OpenMP_Fortran_LIBRARIES}
        )

        target_include_directories(${target_name} PRIVATE
            ${BLA_INCLUDE_DIRS}
            ${OpenMP_Fortran_INCLUDE_DIRS}
        )

        set_target_properties(${target_name} PROPERTIES
            OUTPUT_NAME ${exe_name}
            RUNTIME_OUTPUT_DIRECTORY "${NUMFORT_TESTS_DIRECTORY}/${_comp}"
        )

        add_test(NAME ${exe_name} COMMAND ${exe_name}
            WORKING_DIRECTORY "${NUMFORT_TESTS_DIRECTORY}/${_comp}"
        )
    endforeach()

endforeach()

add_subdirectory(openmp)
