# src/tests/CMakeLists.txt

get_filename_component(_dir "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

find_package(fcore REQUIRED)
find_package(OpenMP)

set(SOURCE_FILES
    linear_omp.f90
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    "${${PROJECT_NAME}_TESTS_DIRECTORY}/${_dir}"
)

if (OPENMP_FOUND)
    foreach (test_source IN LISTS SOURCE_FILES)
        get_filename_component(exe_name "${test_source}" NAME_WE)
        set(target_name ${_dir}_${exe_name})
        add_executable(${target_name} "${test_source}")
        target_link_libraries(${target_name} ${${CMAKE_PROJECT_NAME}_LIBRARIES})
        target_link_libraries(${target_name} fcore::fcore)
        set_target_properties(${target_name} PROPERTIES
            OUTPUT_NAME ${exe_name}
        )
        add_test(NAME ${exe_name} COMMAND ${exe_name}
            WORKING_DIRECTORY "${${PROJECT_NAME}_TESTS_DIRECTORY}/${_dir}"
        )
    endforeach()
endif()
