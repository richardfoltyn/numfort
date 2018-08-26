
set(SOURCE_FILES
    dfovec.f
    dfoxs.f
	main.f
)

# set output directory to examples directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${${PROJECT_NAME}_EXAMPLES_DIRECTORY}")

set(_target example_external_${PROJECT_NAME})
add_executable(${_target} ${SOURCE_FILES}
    $<TARGET_OBJECTS:${${PROJECT_NAME}_LIBRARY_NAME}>
)
target_link_libraries(${_target} ${${PROJECT_NAME}_LIBRARIES})

