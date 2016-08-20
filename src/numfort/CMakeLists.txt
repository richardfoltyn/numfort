# src/numfort/CMakeLists.txt
# CMake configuration for numfort library

# number of version components
set(NUMFORT_VERSION_COUNT 3)
set(NUMFORT_VERSION_MAJOR 0)
set(NUMFORT_VERSION_MINOR 1)
set(NUMFORT_VERSION_PATCH 0)
set(NUMFORT_VERSION
	"${NUMFORT_VERSION_MAJOR}.${NUMFORT_VERSION_MINOR}.${NUMFORT_VERSION_PATCH}"
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.f90.in
	${CMAKE_BINARY_DIR}/version.f90
)

set(SOURCE_FILES
	${CMAKE_BINARY_DIR}/version.f90
	numfort.f90
)

# add individual numfort components as object libraries
set(COMPONENTS core common arrays interpolate linalg optimize)
# set(COMPONENTS core arrays interpolate linalg)

foreach (_comp IN LISTS COMPONENTS)
	add_subdirectory(${_comp})
	# set mod-file output directory
	set_target_properties(${_comp} PROPERTIES
	    Fortran_MODULE_DIRECTORY "${${PROJECT_NAME}_MODDIR}"
	)

    # we need to make sure that -fPIC etc. is passed to compiler when building
    # object libraries that will be used to create a shared lib
    if (BUILD_SHARED_LIBS)
        set_target_properties(${_comp} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif ()
endforeach()

set(COMPONENT_OBJECTS)
foreach (_comp IN LISTS COMPONENTS)
	list(APPEND COMPONENT_OBJECTS $<TARGET_OBJECTS:${_comp}>)
endforeach()

include_directories("${${PROJECT_NAME}_MODDIR}")

add_library(${LIBRARY_NAME} ${SOURCE_FILES} ${COMPONENT_OBJECTS})
target_link_libraries(${LIBRARY_NAME} ${${PROJECT_NAME}_BLA_LIBRARIES})
set_target_properties(${LIBRARY_NAME} PROPERTIES
	Fortran_MODULE_DIRECTORY "${${PROJECT_NAME}_MODDIR}"
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
	LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
	ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

# for benchmarks, examples, etc
add_library(${LIBRARY_NAME}_static STATIC ${SOURCE_FILES} ${COMPONENT_OBJECTS})

install(TARGETS ${LIBRARY_NAME}
    RUNTIME DESTINATION "${INSTALL_LIBRARY_PREFIX}"
    LIBRARY DESTINATION "${INSTALL_LIBRARY_PREFIX}"
    ARCHIVE DESTINATION "${INSTALL_LIBRARY_PREFIX}"
)

install(DIRECTORY ${${PROJECT_NAME}_MODDIR}/
	DESTINATION "${INSTALL_INCLUDE_PREFIX}"
)
