# src/external/newuoa2/CMakeLists.txt
# CMake configuration for modified NEWUOA library

cmake_minimum_required(VERSION 3.0)

project(newuoa2 Fortran)
set(PROJECT_NAME_PRETTY NEWUOA2)
set(LIBRARY_NAME ${PROJECT_NAME})

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
# Default standalone build

set(BUILD_STANDALONE ON)
set(EXAMPLES_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include")

################################################################################
# Numfort subproject setup

if (${CMAKE_PROJECT_NAME} STREQUAL "numfort")
    set(BUILD_STANDALONE OFF)
	set(EXAMPLES_DIRECTORY
		"${${CMAKE_PROJECT_NAME}_EXAMPLES_DIRECTORY}/external/${PROJECT_NAME}"
	)
    # add to list of object libraries that will be used to create ultimate
    # numfort library
    set_property(GLOBAL APPEND
        PROPERTY ${CMAKE_PROJECT_NAME}_COMPONENT_TARGETS
        ${LIBRARY_NAME}
    )

    # Store module files in global module directory
    set(MODULE_DIRECTORY ${${CMAKE_PROJECT_NAME}_MODDIR}/external)

    message(STATUS "Building ${PROJECT_NAME_PRETTY} as sub-project")
else ()
	message(STATUS "Building ${PROJECT_NAME_PRETTY} as standalone project")
endif ()

################################################################################
# Project-specific compiler flags

# The Fortran 77 code uses implicit variable declarations, turn off IFORT
# warnings
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (WIN32)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /warn:nodeclarations")
    else ()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn nodeclarations")
    endif ()
endif ()

################################################################################
# Target definitions

set(SOURCE_FILES
    newuoa2_real32.F
    newuoa2_real64.F
)

set(_name ${LIBRARY_NAME})

# build shared/static library from object library only if standalone build.
# Otherwise, parent project will use object library directly.
if (BUILD_STANDALONE)
    add_library(${_name} ${SOURCE_FILES})
    if (BUILD_SHARED_LIBS)
        set_target_properties(${_name} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif ()
else ()
    add_library(${_name} OBJECT ${SOURCE_FILES})
endif ()

set_target_properties(${_name} PROPERTIES
    Fortran_MODULE_DIRECTORY ${MODULE_DIRECTORY}
)

target_include_directories(${_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()
