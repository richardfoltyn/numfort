
set(F77_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/lbfgsb.f
    ${CMAKE_CURRENT_SOURCE_DIR}/linpack.f
    ${CMAKE_CURRENT_SOURCE_DIR}/timer.f
)

set(SOURCE_FILES
    ${SOURCE_FILES}
    ${F77_SOURCE_FILES}
    PARENT_SCOPE
)

# set less restrictive compile options for external code
set_source_files_properties(${F77_SOURCE_FILES}
    PROPERTIES COMPILE_FLAGS "${${PROJECT_NAME}_FFLAGS_EXTERNAL}"
)

if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
