# src/numfort/CMakeLists.txt
# CMake configuration for numfort library

set(COMPONENTS core arrays interpolate linalg)

foreach (_comp IN LISTS COMPONENTS)
	add_subdirectory(${_comp})
	# set mod-file output directory
	set_target_properties(${_comp} PROPERTIES
	    Fortran_MODULE_DIRECTORY "${NUMFORT_MODDIR}"
	)
endforeach()

set(COMPONENT_OBJECTS)
foreach (_comp IN LISTS COMPONENTS)
	list(APPEND COMPONENT_OBJECTS $<TARGET_OBJECTS:${_comp}>)
endforeach()

add_library(${LIBRARY_NAME} ${COMPONENT_OBJECTS})
target_link_libraries(${LIBRARY_NAME} ${NUMFORT_LAPACK_LIBRARIES})

# for benchmarks, examples, etc
add_library(${LIBRARY_NAME}_static STATIC ${COMPONENT_OBJECTS})
target_link_libraries(${LIBRARY_NAME} ${NUMFORT_LAPACK_LIBRARIES})

install(TARGETS ${LIBRARY_NAME}
    RUNTIME DESTINATION "${INSTALL_LIBRARY_PREFIX}"
    LIBRARY DESTINATION "${INSTALL_LIBRARY_PREFIX}"
    ARCHIVE DESTINATION "${INSTALL_LIBRARY_PREFIX}"
)

install(DIRECTORY ${NUMFORT_MODDIR}/ DESTINATION "${INSTALL_INCLUDE_PREFIX}")
