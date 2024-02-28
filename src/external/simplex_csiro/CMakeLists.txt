# src/external/simplex_csiro/CMakeLists.txt
# CMake configuration for simplex (CSIRO) library

cmake_minimum_required(VERSION 3.12)

project(simplex_csiro Fortran)
set(PROJECT_NAME_PRETTY "SIMPLEX (CSIRO implementation)")
set(LIBRARY_NAME ${PROJECT_NAME})

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif ()

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
		"${NUMFORT_EXAMPLES_DIRECTORY}/external/${PROJECT_NAME}"
	)
    # add to list of object libraries that will be used to create ultimate
    # numfort library
    set_property(GLOBAL APPEND
        PROPERTY ${CMAKE_PROJECT_NAME}_COMPONENT_TARGETS
        ${LIBRARY_NAME}
    )

    # Store module files in global module directory
    set(MODULE_DIRECTORY "${EXTERNAL_MODULE_DIRECTORY}/${LIBRARY_NAME}")

    message(STATUS "Building ${PROJECT_NAME_PRETTY} as sub-project")
else ()
	message(STATUS "Building ${PROJECT_NAME_PRETTY} as standalone project")
endif ()

################################################################################
# Target definitions

set(SOURCE_FILES
    simplex_real32.F90
    simplex_real64.F90
    simplex.f90
)

# build shared/static library from object library only if standalone build.
# Otherwise, parent project will use object library directly.
if (BUILD_STANDALONE)
    add_library(${LIBRARY_NAME} ${SOURCE_FILES})
    if (BUILD_SHARED_LIBS)
        set_target_properties(${LIBRARY_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif ()
else ()
    add_library(${LIBRARY_NAME} OBJECT ${SOURCE_FILES})
endif ()

set_target_properties(${LIBRARY_NAME} PROPERTIES
    Fortran_MODULE_DIRECTORY "${MODULE_DIRECTORY}"
)

target_include_directories(${LIBRARY_NAME} PUBLIC
    "${MODULE_DIRECTORY}"
)

target_compile_options(${LIBRARY_NAME} PRIVATE
    ${EXTERNAL_FFLAGS}
    $<$<CONFIG:Debug>:${EXTERNAL_FFLAGS_DEBUG}>
)

if (USE_LTO)
    set_target_properties(${LIBRARY_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()
