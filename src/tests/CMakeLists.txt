# src/tests/CMakeLists.txt
# Configuration file for numfort unit tests

find_package(fcore REQUIRED)

list(APPEND COMPONENTS
    arrays
    integrate
    interpolate
    stats
)

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
        set(target_name ${_comp}_${exe_name})
        add_executable(${target_name}
            "${_comp}/${test_source}"
            ${${_comp}_aux_source}
        )
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
