# fitpack/CMakeLists.txt
# CMake configuration for FITPACK library

cmake_minimum_required(VERSION 3.0)

project(simplex_csiro Fortran)
set(PROJECT_NAME_PRETTY "SIMPLEX (CSIRO implementation)")
set(${PROJECT_NAME}_LIBRARY_NAME ${PROJECT_NAME})

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

################################################################################
# User-definable options
if (WIN32)
    option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
else ()
    option(BUILD_SHARED_LIBS "Build shared libraries" ON)
endif ()
option(BUILD_EXAMPLES "Build examples/demo code" OFF)

################################################################################
# Numfort subproject setup

if (${CMAKE_PROJECT_NAME} STREQUAL "numfort")
	set(BUILD_OBJECT_LIBS ON)
	set(${PROJECT_NAME}_EXAMPLES_DIRECTORY
		"${${CMAKE_PROJECT_NAME}_EXAMPLES_DIRECTORY}/external/${PROJECT_NAME}"
	)
    # link against final numfort library
	set(${PROJECT_NAME}_LIBRARIES ${${CMAKE_PROJECT_NAME}_LIBRARIES})
    # add to list of object libraries that will be used to create ultimate
    # numfort library
    set_property(GLOBAL APPEND
        PROPERTY ${CMAKE_PROJECT_NAME}_COMPONENT_TARGETS
        ${${PROJECT_NAME}_LIBRARY_NAME}
    )
    message(STATUS "Building ${PROJECT_NAME_PRETTY} as sub-project")
else ()
	set(BUILD_OBJECT_LIBS OFF)
	set(${PROJECT_NAME}_EXAMPLES_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
	set(${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_LIBRARY_NAME})
	message(STATUS "Building ${PROJECT_NAME_PRETTY} as standalone project")
endif ()

################################################################################
# Set Fortran compiler FLAGS

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # enable optimization that does not break debugging
    set(CMAKE_Fortran_FFLAGS_DEBUG "-Og")
endif ()

################################################################################
# Target definitions

set(SOURCE_FILES
    simplex_real64.f90
    simplex_real32.f90
    simplex.f90
)

if (BUILD_OBJECT_LIBS)
	set(_OBJECT OBJECT)
endif ()

add_library(${${PROJECT_NAME}_LIBRARY_NAME} ${_OBJECT} ${SOURCE_FILES})
if (BUILD_SHARED_LIBS)
    set_target_properties(${${PROJECT_NAME}_LIBRARY_NAME}
        PROPERTIES POSITION_INDEPENDENT_CODE ON
    )
endif ()

if (BUILD_EXAMPLES)
	add_subdirectory(examples)
endif ()
