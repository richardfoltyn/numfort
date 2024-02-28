# CMake configuration for numfort examples

# Note: Even examples not in openmp/ might need to be linked accordingly
# because the routines they call might use OpenMP if USE_OPENMP=ON
if (USE_OPENMP)
    find_package(OpenMP)

    if (OPENMP_FOUND)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    endif()
endif()

################################################################################
# Additional compiler flags for all examples

unset(_FFLAGS)
unset(_FFLAGS_Debug)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    set(GNU_WARN "-Wall -pedantic -Wuninitialized -Wno-unused")
    set(_FFLAGS "${GNU_WARN} --std=f2018")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel(LLVM)?")
    if (WIN32)
        set(_FFLAGS "/warn:all,nounused,nousage /stand:f18")
    else ()
        set(_FFLAGS "-warn all,nounused,nousage -std18")
    endif (WIN32)
endif ()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_FFLAGS}")
set(CMAKE_Fortran_FLAGS_Debug "${CMAKE_Fortran_FLAGS_DEBUG} ${_FFLAGS_Debug}")

message(STATUS "Building examples")


################################################################################


# List of subdirectories to search for example files
list(APPEND COMPONENTS interpolate openmp optimize stats)


foreach (_comp IN LISTS COMPONENTS)
    # subdirectory stores source files (without directory prefix) in
    # the following variable
    unset(${_comp}_FILES)
    add_subdirectory(${_comp})

    foreach (_file IN LISTS ${_comp}_FILES)
        get_filename_component(exe_name "${_file}" NAME_WE)
        # set target name that does not possibly conflict with other stuff
        set(target_name example_${_comp}_${exe_name})
        add_executable(${target_name} "${_comp}/${_file}")

        target_link_libraries(${target_name} PRIVATE
            ${NUMFORT_LIBRARY_NAME}_static
            ${BLA_LIBRARIES}
            ${OpenMP_Fortran_LIBRARIES}
        )

        if (IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${_comp}/include")
            target_include_directories(${target_name} PRIVATE "${_comp}/include")
        endif()

        target_include_directories(${target_name} PRIVATE
            ${BLA_INCLUDE_DIRS}
            ${OpenMP_Fortran_INCLUDE_DIRS}
        )

        set_target_properties(${target_name} PROPERTIES
            OUTPUT_NAME ${exe_name}
            RUNTIME_OUTPUT_DIRECTORY "${NUMFORT_EXAMPLES_DIRECTORY}/${_comp}"
        )
    endforeach ()

endforeach ()
