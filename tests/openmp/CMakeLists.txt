# src/tests/CMakeLists.txt

get_filename_component(_dir "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

find_package(fcore REQUIRED)
find_package(OpenMP)

set(SOURCE_FILES
    linear_omp.f90
    minpack.f90
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    "${NUMFORT_TESTS_DIRECTORY}/${_dir}"
)

if (OPENMP_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    foreach (test_source IN LISTS SOURCE_FILES)
        get_filename_component(exe_name "${test_source}" NAME_WE)
        set(target_name test_${_dir}_${exe_name})
        add_executable(${target_name} "${test_source}")

        target_link_libraries(${target_name} ${NUMFORT_LIBRARY_NAME}
            fcore::fcore
            ${BLA_LIBRARIES}
        )

        target_include_directories(${target_name} PRIVATE
            ${CMAKE_Fortran_MODULE_DIRECTORY}
            ${BLA_INCLUDE_DIRS}
        )

        set_target_properties(${target_name} PROPERTIES
            OUTPUT_NAME ${exe_name}
        )

        add_test(NAME ${exe_name} COMMAND ${exe_name}
            WORKING_DIRECTORY "${${PROJECT_NAME}_TESTS_DIRECTORY}/${_dir}"
        )
    endforeach()
endif()
