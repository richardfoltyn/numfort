
set(SOURCE_FILES
    dfovec.f
    dfoxs.f
	main.f
)

# set output directory to examples directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${EXAMPLES_DIRECTORY}")

include_directories("${MODULE_DIRECTORY}")

set(_target example_external_${PROJECT_NAME})

add_executable(${_target} ${SOURCE_FILES}
	$<TARGET_OBJECTS:${LIBRARY_NAME}>
)
