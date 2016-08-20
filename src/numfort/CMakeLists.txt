# src/numfort/CMakeLists.txt
# CMake configuration for numfort library

set(COMPONENTS core common arrays interpolate linalg optimize)

foreach (_comp IN LISTS COMPONENTS)
	add_subdirectory(${_comp})
	# set mod-file output directory
	set_target_properties(${_comp} PROPERTIES
	    Fortran_MODULE_DIRECTORY "${${PROJECT_NAME}_MODDIR}"
	)
endforeach()

set(COMPONENT_OBJECTS)
foreach (_comp IN LISTS COMPONENTS)
	list(APPEND COMPONENT_OBJECTS $<TARGET_OBJECTS:${_comp}>)
endforeach()

set(_SHARED STATIC)
if (BUILD_SHARED_LIBS)
	set(_SHARED SHARED)
endif ()
add_library(${LIBRARY_NAME} ${_SHARED} ${COMPONENT_OBJECTS})
target_link_libraries(${LIBRARY_NAME} ${${PROJECT_NAME}_BLA_LIBRARIES})
set_target_properties(${LIBRARY_NAME} PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
	LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
	ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

# for benchmarks, examples, etc
add_library(${LIBRARY_NAME}_static STATIC ${COMPONENT_OBJECTS})

install(TARGETS ${LIBRARY_NAME}
    RUNTIME DESTINATION "${INSTALL_LIBRARY_PREFIX}"
    LIBRARY DESTINATION "${INSTALL_LIBRARY_PREFIX}"
    ARCHIVE DESTINATION "${INSTALL_LIBRARY_PREFIX}"
)

install(DIRECTORY ${${PROJECT_NAME}_MODDIR}/
	DESTINATION "${INSTALL_INCLUDE_PREFIX}"
)
