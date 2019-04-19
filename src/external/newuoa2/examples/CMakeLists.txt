
set(SOURCE_FILES
    dfovec.f
    dfoxs.f
	main.f
)

# set output directory to examples directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${EXAMPLES_DIRECTORY}")

include_directories("${MODULE_DIRECTORY}")

if (CMAKE_Fortran_FLAGS_ORIG)
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_ORIG}")
endif()

if (CMAKE_Fortran_FLAGS_DEBUG_ORIG)
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG_ORIG}")
endif()

set(_target example_external_${PROJECT_NAME})

if (BUILD_STANDALONE)
	add_executable(${_target} ${SOURCE_FILES})
	target_link_libraries(${_target} ${LIBRARY_NAME})
else ()
	add_executable(${_target} ${SOURCE_FILES}
    	$<TARGET_OBJECTS:${LIBRARY_NAME}>
	)
endif ()
